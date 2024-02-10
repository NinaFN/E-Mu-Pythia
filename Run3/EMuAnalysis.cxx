// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "CommonConstants/MathConstants.h"
#include "TDatabasePDG.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/MixingHandler.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/MixingLibrary.h"

#include <TH1F.h>
#include <TH3F.h>
#include <TList.h>
#include <THashList.h>
#include <TString.h>
#include <algorithm>

#include <iostream>
#include <vector>
#include <string>
#include <cstdio>
#include <tuple>
#include <stdlib.h>
#include "cmath"

using namespace o2;
using namespace o2::framework;
using namespace constants::math;
using namespace o2::framework::expressions;
using namespace o2::aod;
//using namespace o2::soa;

namespace o2::aod
{
  namespace extrainfo
  {
    DECLARE_SOA_COLUMN(IsEventSelected, isEventSelected, int);
    DECLARE_SOA_COLUMN(MixingHashTest, mixingHash, int);
    DECLARE_SOA_COLUMN(IsBarrelSelected, isBarrelSelected, int);
    DECLARE_SOA_COLUMN(IsMuonSelected, isMuonSelected, int);
    DECLARE_SOA_COLUMN(IsBarrelSelectedPrefilter, isBarrelSelectedPrefilter, int);
    DECLARE_SOA_COLUMN(IsPrefiltered, isPrefiltered, int);

  }

  DECLARE_SOA_TABLE(EventCuts, "AOD", "EVENTCUTS", extrainfo::IsEventSelected);
  DECLARE_SOA_TABLE(MixingHashes, "AOD", "MIXINGHASHES", extrainfo::MixingHashTest);
  //DECLARE_SOA_TABLE(BarrelTrackCuts, "AOD", "BARRELTRACKCUTS", extrainfo::IsBarrelSelected, extrainfo::IsBarrelSelectedPrefilter);
  DECLARE_SOA_TABLE(BarrelSelected, "AOD", "BARRELSELECTED", extrainfo::IsBarrelSelected);
  DECLARE_SOA_TABLE(MuonTrackCuts, "AOD", "MUONTRACKCUTS", extrainfo::IsMuonSelected);
  DECLARE_SOA_TABLE(Prefilter, "AOD", "PREFILTER", extrainfo::IsPrefiltered);
} 

using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended>;
using MyEventsSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts>;
using MyEventsHashSelected = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::EventCuts, aod::MixingHashes>;


using MyBarrelTracks = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID>;
using MyBarrelTracksSelected = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID, aod::BarrelSelected>;
using MyMuonTracks = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra>;
using MyMuonTracksSelected = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::MuonTrackCuts>;

constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelPID;
constexpr static uint32_t gkMuonFillMap = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra;

//////Trying to follow structure of other DQ tasks (tableReader, dileptonMuMu) - I'm frankensteining it, most of the code comes from there with minor adaptations

//Function definition to compute the difference in phi between two particles
  //Taken from the O2 Tutorial session run on 13/10/22
Double_t ComputeDeltaPhi(Double_t phi1, Double_t phi2) 
{
  //To be completely sure, use inner products
  Double_t x1, y1, x2, y2;
  x1 = TMath::Cos( phi1 );
  y1 = TMath::Sin( phi1 );
  x2 = TMath::Cos( phi2 );
  y2 = TMath::Sin( phi2 );
  
  Double_t lInnerProd = x1*x2 + y1*y2;
  Double_t lVectorProd = x1*y2 - x2*y1;
  Double_t lReturnVal = 0;

  if( lVectorProd > 1e-8 ) {lReturnVal = TMath::ACos(lInnerProd);}
  if( lVectorProd < -1e-8 ) {lReturnVal = -TMath::ACos(lInnerProd);}
  if( lReturnVal < -TMath::Pi()/2. ) {lReturnVal += 2.*TMath::Pi();}

  return lReturnVal;
}

void DefineHistograms(HistogramManager* histMan, TString histClasses, Configurable<std::string> configVar);

struct EventSelection
{
  //Produces
  Produces<aod::EventCuts> eventSel;
  Produces<aod::MixingHashes> hash;
  //Configurables
  Configurable<string> fConfigMixingVariables{"cfgMixingVars", "", "Mixing configs separated by a comma, default no mixing"};
  Configurable<string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<std::string> fConfigAddEventHistogram{"cfgAddEventHistogram", "", "Comma separated list of histograms"};

  //Other
  OutputObj<THashList> fOutputList{"output"};
  HistogramManager* fHistMan = nullptr;
  MixingHandler* fMixHandler = nullptr;
  AnalysisCompositeCut* fEventCut;
  float* fValues;


  void init(o2::framework::InitContext&)
  {
    fValues = new float[VarManager::kNVars];
    VarManager::SetDefaultVarNames();

    //Cuts
    DefineCuts(); //seperated this as in dileptonMuMu, not sure if there's a huge coding benefit but it feels more organized since I'm also including the additional mixing cuts

    //Histo Setup
    if (fConfigQA) //requirement from tableReader but not in dileptonMuMu, decided to use it 
    {
      fHistMan = new HistogramManager("analysisHistos", "", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);
      DefineHistograms(fHistMan, "Event_BeforeCuts;Event_AfterCuts;EventMixing_BeforeCuts;EventMixing_AfterCuts;", fConfigAddEventHistogram); // define all histograms
      VarManager::SetUseVars(fHistMan->GetUsedVars());                                           // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }

    //Mixing Setup
    TString mixVarsString = fConfigMixingVariables.value;
    std::unique_ptr<TObjArray> objArray(mixVarsString.Tokenize(","));
    if (objArray->GetEntries() > 0) {
      fMixHandler = new MixingHandler("mixingHandler", "mixing handler");
      fMixHandler->Init();
      for (int iVar = 0; iVar < objArray->GetEntries(); ++iVar) {
        dqmixing::SetUpMixing(fMixHandler, objArray->At(iVar)->GetName());
      }
    }

  }


  void DefineCuts()
  {
    fEventCut = new AnalysisCompositeCut(true);
    TString eventCutStr = fConfigEventCuts.value; //extract cut info from configurable
    fEventCut->AddCut(dqcuts::GetAnalysisCut(eventCutStr.Data()));
    VarManager::SetUseVars(AnalysisCut::fgUsedVars);
  }

  void process(MyEvents::iterator const& event)
  {
    // Reset the fValues array
    VarManager::ResetValues(0, VarManager::kNEventWiseVariables, fValues);


    VarManager::FillEvent<gkEventFillMap>(event);

    if (fConfigQA) {
      fHistMan->FillHistClass("Event_BeforeCuts", fValues); // automatically fill all the histograms in the class Event
      fHistMan->FillHistClass("EventMixing_BeforeCuts", fValues);
    }

    if (fEventCut->IsSelected(fValues)) {
      if (fConfigQA) {
        fHistMan->FillHistClass("Event_AfterCuts", fValues);
      }
      eventSel(1);
    } else {eventSel(0);}

    if (fMixHandler != nullptr) {
      int hh = fMixHandler->FindEventCategory(fValues);
      hash(hh);
    }
  }

};


struct BarrelSelection
{
  //Produces
  Produces<aod::BarrelSelected> trackSel;

  //Configurables
  //Configurable<string> fConfigTrackCuts{"cfgTrackCuts", "", "Comma separated list of barrel track cuts"};
  Configurable<int> fConfigPrefilterCutId{"cfgPrefilterCutId", 32, "Id of the Prefilter track cut (starting at 0)"}; // In order to create another column prefilter (should be temporary before improving cut selection in configurables, then displaced to AnalysisPrefilterSelection)
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<string> fConfigAddBarrelHistogram{"cfgAddBarrelHistogram", "", "Comma separated list of histograms"};
  Configurable<float> fPtCut{"cfgPtCut", 5.0, "Minimum track pt"};
  Configurable<bool> fPID{"cfgPID", false, "If true, apply electron PID cuts"};
  Configurable<float> fTPCCut{"cfgTPCCut", 2, "Cut on electron Nsigma from TPC"};
  Configurable<float> fTOFCut{"cfgTOFCut", 2, "Cut on electron Nsigma from TOF"};
  Configurable<bool> fPIDRej{"cfgPIDRej", false, "If true, apply electron PID rejection cuts"};
  //Configurable<float> fTPCCut{"cfgTOFCut", 2, "Cut on electron Nsigma from TOF"};


  
  //Other
  OutputObj<THashList> fOutputList{"outputTracks"};
  HistogramManager* fHistMan;
  float* fValues;

  AnalysisCut* fTrackCut = new AnalysisCut("fTrackCut","");

  void init(o2::framework::InitContext&)
  {
    fTrackCut->AddCut(VarManager::kPt, 0.0, fPtCut, true);
    fTrackCut->AddCut(VarManager::kEta, -0.9, 0.9);

    if(fPID == true)
    {
      fTrackCut->AddCut(VarManager::kTPCnSigmaEl, -fTPCCut, fTPCCut);
      fTrackCut->AddCut(VarManager::kTOFnSigmaEl, -fTOFCut, fTOFCut);
    }

    /*if(fPIDRej == true)
    {

    }*/

    //cut->AddCut(VarManager::kTPCchi2, 0.0, 4.0);
    //cut->AddCut(VarManager::kTPCncls, 70, 161.);
    //cut->AddCut(VarManager::kTPCsignal, 70., 100.);

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
    fValues = new float[VarManager::kNVars];

    //Add CCDB stuff?

  }

  void process(MyEvents::iterator const& event, MyBarrelTracks const& tracks)
  {
    VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);

    // fill event information which might be needed in histograms/cuts that combine track and event properties
    VarManager::FillEvent<gkEventFillMap>(event,fValues);

    //potentially CCDB stuff here

    //Barrel Selection
    trackSel.reserve(tracks.size());
    uint32_t filterMapBarrel = 0;
    bool prefilterSelected = false;
    int iCut = 0;

    for (auto& track : tracks) {
      filterMapBarrel = 0;
      prefilterSelected = false;
      VarManager::FillTrack<gkTrackFillMap>(track,fValues);
      if (fConfigQA) { 
        fHistMan->FillHistClass("TrackBarrel_BeforeCuts", fValues);
      }
      iCut = 0;
      //for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, iCut++) {
        if (fTrackCut->IsSelected(fValues)) {
          //if (iCut != fConfigPrefilterCutId) {
          //  filterMapBarrel |= (uint32_t(1) << iCut);
          //}
          //if (iCut == fConfigPrefilterCutId) {
            prefilterSelected = true;
          //}
          //if (fConfigQA) { 
          //  fHistMan->FillHistClass(Form("TrackBarrel_%s", (*cut)->GetName()), VarManager::fgValues);
          //}
        }
      //}

      trackSel(static_cast<int>(prefilterSelected));//static_cast<int>(filterMapBarrel), static_cast<int>(prefilterSelected)
    } // end loop over tracks

  }

};


struct MuonSelection
{
  //Produces
  Produces<aod::MuonTrackCuts> muonSel;


  //Configurables
  //Configurable<string> fConfigMuonCuts{"cfgMuonCuts", "muonQualityCuts", "Comma separated list of muon cuts"};
  Configurable<int> fConfigPrefilterCutId{"cfgPrefilterCutId", 32, "Id of the Prefilter track cut (starting at 0)"}; // In order to create another column prefilter (should be temporary before improving cut selection in configurables, then displaced to AnalysisPrefilterSelection)
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<string> fConfigAddMuonHistogram{"cfgAddMuonHistogram", "", "Comma separated list of histograms"};
  Configurable<float> fPtCut{"cfgPtCut", 2.0, "Minimum muon pt"};


  //Other
  OutputObj<THashList> fOutputList{"outputMuons"};
  HistogramManager* fHistMan;
  float* fValues;
  //std::vector<AnalysisCompositeCut> fMuonCuts;

  AnalysisCut* fMuonCut = new AnalysisCut("fMuonCut","");
  //std::vector<AnalysisCut*> fMuonCuts;

  void init(o2::framework::InitContext&)
  {

    fMuonCut->AddCut(VarManager::kPt, 0.0, fPtCut, true); 
    fMuonCut->AddCut(VarManager::kEta, -4.0, -2.5);
    fMuonCut->AddCut(VarManager::kMuonChi2MatchMCHMID, 0.0, 20); // matching MCH-MID
    fMuonCut->AddCut(VarManager::kMuonChi2MatchMCHMFT, 0.0, 45); // matching MFT-MCH
    
    /*fMuonCut->AddCut(VarManager::kEta, -4.0, -2.5);
    fMuonCut->AddCut(VarManager::kMuonRAtAbsorberEnd, 17.6, 89.5);
    fMuonCut->AddCut(VarManager::kMuonPDca, 0.0, 594.0, false, VarManager::kMuonRAtAbsorberEnd, 17.6, 26.5);
    fMuonCut->AddCut(VarManager::kMuonPDca, 0.0, 324.0, false, VarManager::kMuonRAtAbsorberEnd, 26.5, 89.5);
    fMuonCut->AddCut(VarManager::kMuonChi2, 0.0, 1e6);
    fMuonCut->AddCut(VarManager::kMuonChi2MatchMCHMID, 0.0, 1e6); // matching MCH-MID
    fMuonCut->AddCut(VarManager::kMuonChi2MatchMCHMFT, 0.0, 1e6); // matching MFT-MCH*/
    //fMuonCuts.push_back(fMuonPtCut);

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
    fValues = new float[VarManager::kNVars];

  /*
    //Histogram Setup
    if (fConfigQA) {
      VarManager::SetDefaultVarNames();
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

      // set one histogram directory for each defined track cut
      TString histMuonDirNames = "TrackMuon_BeforeCuts;";
      for (auto& cut : fMuonCuts) {
        histMuonDirNames += Form("TrackMuon_%s;", cut->GetName()); 
      }

      DefineHistograms(fHistMan, histMuonDirNames.Data(),fConfigAddMuonHistogram);
      VarManager::SetUseVars(fHistMan->GetUsedVars());                           // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }
  */
    //Add CCDB stuff?

  }

  void process(MyEvents::iterator const& event, MyMuonTracks const& muons)
  {
    VarManager::ResetValues(0, VarManager::kNMuonTrackVariables);
    
    // fill event information which might be needed in histograms/cuts that combine track and event properties
    VarManager::FillEvent<gkEventFillMap>(event);

    //potentially CCDB stuff here
  
    //Muon Selection
    muonSel.reserve(muons.size());
    //uint32_t filterMapMuon = 0;
    bool prefilterSelected = false;
    //int iCut = 0;

    for (auto& muon : muons) {
      //filterMapMuon = 0;
      //iCut = 0;
      prefilterSelected = false;
      VarManager::FillTrack<gkMuonFillMap>(muon);

      //for (auto cut = fMuonCuts.begin(); cut != fMuonCuts.end(); cut++, iCut++) {
        if ((fMuonCut)->IsSelected(VarManager::fgValues)) {
          //if (iCut != fConfigPrefilterCutId) {
          //  filterMapMuon |= (uint32_t(1) << iCut);
          //}
          prefilterSelected = true;
          //if (fConfigQA) { 
          //  fHistMan->FillHistClass(Form("TrackMuon_%s", (*cut)->GetName()), VarManager::fgValues);
          //}
        }
      //}

      muonSel(static_cast<int>(prefilterSelected));
    } // end loop over muons
    

  }


  //void process()
};

/*struct EventMixing
{
  Configurable<string> fTrackPt{"cfgTrackPt", "", "barrel track cuts"};
  Configurable<string> fMuonPt{"cfgMuonPt", "", "muon cuts"};
  Configurable<int> fConfigMixingDepth{"cfgMixingDepth", 100, "Number of Events stored for event mixing"};

  //Filters
  //Filter filterEventSelected = aod::extrainfo::isEventSelected == 1;
  Filter filterTrackSelected = aod::extrainfo::isBarrelSelected > 0;
  Filter filterMuonTrackSelected = aod::extrainfo::isMuonSelected > 0;

  OutputObj<THashList> fOutputList{"output"};
  HistogramManager* fHistMan;
  float* fValues;

  void init(o2::framework::InitContext&)
  {

  }

  void process(soa::Filtered<MyEventsHashSelected>& events, soa::Filtered<MyBarrelTracksSelected> const& tracks, soa::Filtered<MyMuonTracksSelected> const& muons)
  {
    events.bindExternalIndices(&muons);
    auto tracksTuple = std::make_tuple(tracks);
    auto muonsTuple = std::make_tuple(muons);
    GroupSlicer slicerTracks(events, tracksTuple);
    GroupSlicer slicerMuons(events, muonsTuple);
    

    //////////////////////////////////////

    uint32_t twoTrackFilter = 0;
    for (auto& track : tracks) 
    {
      for (auto& muon : muons) 
      {
        twoTrackFilter = uint32_t(track.isBarrelSelected()) & uint32_t(muon.isMuonSelected()) & fTwoTrackFilterMask;
        if (!twoTrackFilter) { 
          continue;
        }
        
        VarManager::FillPairME<kElectronMuon>(track, muon); 

        if (track.sign() * muon.sign() < 0) {
        //LOGF(info, "Pair Type: Opp");
        PhiCorrelations.get<TH1>(HIST("PhiCorrelation_OS"))->Fill(dPhi);
        }

        else{
          //LOGF(info, "Pair Type: Like");
          PhiCorrelations.get<TH1>(HIST("PhiCorrelation_LS"))->Fill(dPhi);
        }

      }
    }


  }

}*/


struct AnalysisLight
{
  Configurable<int> fBins{"cfgBins", 70, "Number of bins used for histogram"};

  HistogramManager* fHistMan;
  std::vector<std::vector<TString>> fTrackMuonHistNames;
  int count = 0;
  int countFilled = 0;
  int TrigAssocRatio = 0;

  //Filters to implement cuts chosen in event/barrel/muon selections
  //Filter filterEventSelected = aod::extrainfo::isEventSelected == 1;
  Filter filterBarrelTrackSelected = aod::extrainfo::isBarrelSelected > 0;
  Filter filterMuonTrackSelected = aod::extrainfo::isMuonSelected > 0;


  float* fValues = new float[VarManager::kNVars];

  HistogramRegistry EventStats{
    "EventStats",
    {
      {"VtxZ", "Z Vertex position in cm", {HistType::kTH1F, {{fBins, -10, 10,"Vtx Z (cm)"}}}},
      {"TrigAssocRatio", "Ratio of recorded trigger to associated particles per event", {HistType::kTH1F, {{fBins, 0, 25,"N_{trig (#mu)}/N_{assoc (barrel)}"}}}}
    }
    };

  HistogramRegistry UncorrelatedTracks{
    "UncorrelatedTracks",
    {

      {"BarrelPt", "Barrel Tracks Transverse Momentum", {HistType::kTH1F, {{fBins, 0, 50,"p_{T}^{e} (GeV)"}}}},
      {"BarrelEta", "Barrel Tracks Pseudorapidity", {HistType::kTH1F, {{fBins, -5, 5,"#eta^{e}"}}}},
      {"MuonPt", "Muon Track Transverse Momentums", {HistType::kTH1F, {{fBins, 0, 30,"p_{T}^{#mu} (GeV)"}}}},
      {"MuonEta", "Muon Tracks Pseudorapidity", {HistType::kTH1F, {{fBins, -5, 5,"#eta^{#mu}"}}}}
    }
    };

  HistogramRegistry PhiCorrelations{
    "PhiCorrelations",
    {
      {"PhiCorrelation_All", "#Delta #phi (All Pairs)", {HistType::kTH1F, {{fBins, 0, 2*M_PI,"#Delta #phi"}}}},
      {"PhiCorrelation_LS", "#Delta #phi (Like-Sign Pairs)", {HistType::kTH1F, {{fBins, 0, 2*M_PI,"#Delta #phi"}}}},
      {"PhiCorrelation_PP", "#Delta #phi (e^{+} - #mu^{+})", {HistType::kTH1F, {{fBins, 0, 2*M_PI,"#Delta #phi"}}}},
      {"PhiCorrelation_MM", "#Delta #phi (e^{-} - #mu^{-})", {HistType::kTH1F, {{fBins, 0, 2*M_PI,"#Delta #phi"}}}},
      {"PhiCorrelation_OS", "#Delta #phi (Opposite-Sign Pairs)", {HistType::kTH1F, {{fBins, 0, 2*M_PI,"#Delta #phi"}}}},
      {"PhiCorrelation_PM", "#Delta #phi (e^{+} - #mu^{-})", {HistType::kTH1F, {{fBins, 0, 2*M_PI,"#Delta #phi"}}}},
      {"PhiCorrelation_MP", "#Delta #phi (e^{-} - #mu^{+})", {HistType::kTH1F, {{fBins, 0, 2*M_PI,"#Delta #phi"}}}}
    }
    };

  
  void process(MyEvents::iterator const& event, soa::Filtered<MyBarrelTracksSelected> const& tracks, soa::Filtered<MyMuonTracksSelected> const& muons)
  {
    //LOGF(info, "Event Number: %d", count);
    count++;

    if (tracks.size()==0 || muons.size()==0) {
      return;
    }

    /*LOGF(info, " ");
    LOGF(info, "Event Number With Contributions: %d", countFilled);
    LOGF(info, "Number of tracks from the collision: %d", tracks.size());
    LOGF(info, "Number of muons from the collision: %d", muons.size());*/

    TrigAssocRatio = muons.size()/tracks.size();

    countFilled++;
    
    EventStats.get<TH1>(HIST("VtxZ"))->Fill(event.posZ());
    EventStats.get<TH1>(HIST("TrigAssocRatio"))->Fill(TrigAssocRatio);

    

    for (auto& track : tracks)
    {
      UncorrelatedTracks.get<TH1>(HIST("BarrelPt"))->Fill(track.pt());
      UncorrelatedTracks.get<TH1>(HIST("BarrelEta"))->Fill(track.eta());
    }
    for (auto& muon : muons)
    {
      UncorrelatedTracks.get<TH1>(HIST("MuonPt"))->Fill(muon.pt());
      UncorrelatedTracks.get<TH1>(HIST("MuonEta"))->Fill(muon.eta());
    }

    constexpr static int pairTypeEMu = VarManager::kElectronMuon;
    
    for (auto& [track, muon] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks, muons))) 
    {
      //LOGF(info,"Track: %f", track.phi());
      //LOGF(info,"Muon: %f", muon.phi());
      float dPhi = ComputeDeltaPhi(track.phi(),muon.phi());
      //float dPhi = track.phi() - muon.phi();

      //LOGF(info,"dPhi: %f",dPhi);
      //LOGF(info,"dPhi Alt: %f",dPhiAlt);

      if(dPhi<0) {dPhi = dPhi + 2.*M_PI;}
      else if(dPhi>2*M_PI) {dPhi = dPhi - 2.*M_PI;}
      

      VarManager::FillPair<pairTypeEMu, gkTrackFillMap>(track, muon, fValues);

      PhiCorrelations.get<TH1>(HIST("PhiCorrelation_All"))->Fill(dPhi);
      PhiCorrelations.get<TH1>(HIST("PhiCorrelation_Sub"))->Fill(dPhi);

      //LOGF(info, "------------");
      //LOGF(info, "Track Sign: %d", track.sign());
      //LOGF(info, "Muon Sign: %d", muon.sign());


      if (track.sign() * muon.sign() < 0) {
        //LOGF(info, "Pair Type: Opp");
        PhiCorrelations.get<TH1>(HIST("PhiCorrelation_OS"))->Fill(dPhi);
        if(track.sign()<0){ PhiCorrelations.get<TH1>(HIST("PhiCorrelation_MP"))->Fill(dPhi);}
        if(muon.sign()<0){ PhiCorrelations.get<TH1>(HIST("PhiCorrelation_PM"))->Fill(dPhi);}
      }

      else{
        //LOGF(info, "Pair Type: Like");
        PhiCorrelations.get<TH1>(HIST("PhiCorrelation_LS"))->Fill(dPhi);
        if(track.sign()<0){ PhiCorrelations.get<TH1>(HIST("PhiCorrelation_PP"))->Fill(dPhi);}
        if(track.sign()>0){ PhiCorrelations.get<TH1>(HIST("PhiCorrelation_MM"))->Fill(dPhi);}
      }
      
    }

  }       
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    //adaptAnalysisTask<EventSelection>(cfgc),
    adaptAnalysisTask<BarrelSelection>(cfgc),
    adaptAnalysisTask<MuonSelection>(cfgc),
    //adaptAnalysisTask<EventMixing>(cfgc),
    //adaptAnalysisTask<AnalysisEMuHists>(cfgc),
    adaptAnalysisTask<AnalysisLight>(cfgc)};
  
}

void DefineHistograms(HistogramManager* histMan, TString histClasses, Configurable<std::string> configVar)
{
  std::unique_ptr<TObjArray> objArray(histClasses.Tokenize(";"));
  for (Int_t iclass = 0; iclass < objArray->GetEntries(); ++iclass) 
  {
    TString classStr = objArray->At(iclass)->GetName();
    histMan->AddHistClass(classStr.Data());

    TString histName = configVar.value;

    if (classStr.Contains("Event")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "event", histName);
    }

    if (classStr.Contains("Track") && !classStr.Contains("Pairs")) 
    {
      if (classStr.Contains("Barrel")) 
      {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", histName);
        if (classStr.Contains("PIDCalibElectron")) {
          dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "postcalib_electron");
        }
        if (classStr.Contains("PIDCalibPion")) {
          dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "postcalib_pion");
        }
        if (classStr.Contains("PIDCalibProton")) {
          dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", "postcalib_proton");
        }
      }

      if (classStr.Contains("Muon")) {
        dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "track", histName);
      }
    }

    if (classStr.Contains("Pairs")) {
      dqhistograms::DefineHistograms(histMan, objArray->At(iclass)->GetName(), "pair", histName);
    }

  }
}

/*
struct AnalysisLight
{
  //Filter filterEventSelected = aod::extrainfo::isEventSelected == 1;

  HistogramManager* fHistMan;
  std::vector<std::vector<TString>> fTrackMuonHistNames;
  int count = 0;
  int countFilled = 0;

  float* fValues = new float[VarManager::kNVars];
  int nBins = 100;
  HistogramRegistry UncorrelatedTracks{
    "UncorrelatedTracks",
    {

      {"BarrelPt", "Barrel Tracks", {HistType::kTH1F, {{nBins, 0, 50,"p_{T}^{e} (GeV)"}}}},
      {"MuonPt", "Muon Tracks", {HistType::kTH1F, {{nBins, 0, 30,"p_{T}^{#mu} (GeV)"}}}},
    }
    };

  HistogramRegistry BarrelMuon_All{
    "BarrelMuon_All",
    {
      {"BarrelPt_All", "p_{T}^{e} (All Pairs)", {HistType::kTH1F, {{nBins, 0, 100,"p_{T}^{e} (GeV)"}}}},
      {"MuonPt_All", "p_{T}^{#mu} (All Pairs)", {HistType::kTH1F, {{nBins, 0, 100,"p_{T}^{#mu} (GeV)"}}}},
      {"PhiCorrelation_All", "#Delta #phi (All Pairs)", {HistType::kTH1F, {{nBins, 0, 2*M_PI,"#Delta #phi"}}}},
    }
    };

  HistogramRegistry BarrelMuon_LS{
    "BarrelMuon_LS",
    {
      {"BarrelPt_LS", "p_{T}^{e} (Like-Sign Pairs)", {HistType::kTH1F, {{nBins, 0, 100,"p_{T}^{e} (GeV)"}}}},
      {"MuonPt_LS", "p_{T}^{#mu} (Like-Sign Pairs)", {HistType::kTH1F, {{nBins, 0, 100,"p_{T}^{#mu} (GeV)"}}}},
      {"PhiCorrelation_LS", "#Delta #phi (Like-Sign Pairs)", {HistType::kTH1F, {{nBins, 0, 2*M_PI,"#Delta #phi"}}}},
    }
    };

  HistogramRegistry BarrelMuon_OS{
    "BarrelMuon_OS",
    {
      {"BarrelPt_OS", "p_{T}^{e} (Opposite-Sign Pairs)", {HistType::kTH1F, {{nBins, 0, 100,"p_{T}^{e} (GeV)"}}}},
      {"MuonPt_OS", "p_{T}^{#mu} (Opposite-Sign Pairs)", {HistType::kTH1F, {{nBins, 0, 100,"p_{T}^{#mu} (GeV)"}}}},
      {"PhiCorrelation_OS", "#Delta #phi (Opposite-Sign Pairs)", {HistType::kTH1F, {{nBins, 0, 2*M_PI,"#Delta #phi"}}}},
    }
    };


  
  void process(MyEvents::iterator const& event, MyBarrelTracks const& tracks, MyMuonTracks const& muons)
  {
    //LOGF(info, "Event Number: %d", count);
    //count++;

    if (tracks.size()==0 || muons.size()==0) {
      return;
    }
    
    LOGF(info, "Event Number With Contributions: %d", countFilled);
    LOGF(info, "Number of tracks from the collision: %d", tracks.size());
    LOGF(info, "Number of muons from the collision: %d", muons.size());
    countFilled++;
    
    

    for (auto& track : tracks){UncorrelatedTracks.get<TH1>(HIST("BarrelPt"))->Fill(track.pt());}
    for (auto& muon : muons){UncorrelatedTracks.get<TH1>(HIST("MuonPt"))->Fill(muon.pt());}


    
    
    

    constexpr static int pairTypeEMu = VarManager::kElectronMuon;
    
    for (auto& [track, muon] : combinations(o2::soa::CombinationsFullIndexPolicy(tracks, muons))) 
    {

      //float dPhi = ComputeDeltaPhi(track.phi(),muon.phi());
      float dPhi = track.phi() - muon.phi();
      if(dPhi<0) {dPhi = dPhi + 2.*M_PI;}
      else if(dPhi>2*M_PI) {dPhi = dPhi - 2.*M_PI;}

      VarManager::FillPair<pairTypeEMu, gkTrackFillMap>(track, muon, fValues);

      BarrelMuon_All.get<TH1>(HIST("BarrelPt_All"))->Fill(track.pt());
      BarrelMuon_All.get<TH1>(HIST("MuonPt_All"))->Fill(muon.pt());
      BarrelMuon_All.get<TH1>(HIST("PhiCorrelation_All"))->Fill(dPhi);

      if (track.sign() * muon.sign() < 0) {
        BarrelMuon_OS.get<TH1>(HIST("BarrelPt_OS"))->Fill(track.pt());
        BarrelMuon_OS.get<TH1>(HIST("MuonPt_OS"))->Fill(muon.pt());
        BarrelMuon_OS.get<TH1>(HIST("PhiCorrelation_OS"))->Fill(dPhi);
      }

      else{
        BarrelMuon_LS.get<TH1>(HIST("BarrelPt_LS"))->Fill(track.pt());
        BarrelMuon_LS.get<TH1>(HIST("MuonPt_LS"))->Fill(muon.pt());
        BarrelMuon_LS.get<TH1>(HIST("PhiCorrelation_LS"))->Fill(dPhi);
      }
      
    }

  }       
};
*/