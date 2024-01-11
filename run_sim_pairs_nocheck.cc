#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

using namespace Pythia8;
//script to save information about the candidate e-mu pairs in all hardQCD processes

int main() {
    // Turn SoftQCD on/off
    bool softQCD = false;

    // Pythia object
    Pythia pythia;

    // ROOT file for histograms
    TFile* outFile = new TFile("sim_tuples_pairs_nocheck_hardQCDall.root", "RECREATE");

    // pTHat bins
    int nBins;
    const double* binEdges;
    const int* binWidths;
    if (softQCD) {
        nBins = 8;
        static const double tempArray[9] = {0.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 40.0, 60.0};
        //static const int tempArrayW[6] = {16,10,10,14,20,30};

        binEdges = &tempArray[0];
        //binWidths = &tempArrayW[0];
    } else {
        nBins = 7;
        static const double tempArray[8] = {5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 40.0, 60.0};

        binEdges = &tempArray[0];
    }


    // Histograms
    // Total Cross Section
    TH1F *hardPt = new TH1F("HardQCD:All","Process Total Cross-Section;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 55, 5.0, 60.0);
    TH1F *hardPtPart = new TH1F("hardQCD_part","", 55, 5.0, 60.0);

    // HF Cross Sections

    vector<TNtuple*> emuTuples(nBins);
    

    vector<double> binLuminocity(nBins); // luminocity from generated process sigma to calculate cross-sections
    vector<double> binEvCount(nBins);
    vector<double> weightSums(nBins);
    vector<double> sigmaGens(nBins);

    //decayFlag: does lepton come from a heavy quark decay
    //semilepFlag: does lepton come from semileptonic decay
    //pairDecFlag: does lepton have pair partner that comes from heavy quark decay
    //pairSemilepFlag:  does lepton have pair partner that comes from semileptonic decay
    //pairEtaFlag: if lepton has pair partner, is it in the correct eta range

    for (int i = 0; i < nBins; ++i) {                        
        emuTuples[i] = new TNtuple("emu_pairs", "emu_pairs", "idQ:idE:ptE:etaE:phiE:thetaE:idMu:ptMu:etaMu:phiMu:thetaMu:phiDiff");        

    }

    // Number of events to generate per bin.
    int N_events = 5000000;

    int code = 0;


    int genEvents;

    // run events for each ptHat bin 
    for (int iBin = 0; iBin < nBins; ++iBin) {
        if (softQCD && iBin < 1) {
            pythia.readString("HardQCD:all = off");
            pythia.readString("HardQCD:hardccbar = off");
            pythia.readString("HardQCD:hardbbbar = off");
            pythia.readString("SoftQCD:nonDiffractive = on");
            //genEvents = N_events*10;
        } else {
            // set pythia initialization variables
            genEvents = N_events;
            pythia.readString("HardQCD:all = on");
            //pythia.readString("HardQCD:hardccbar = on");
            //pythia.readString("HardQCD:hardbbbar = on");
            pythia.readString("SoftQCD:nonDiffractive = off");
        }

        pythia.readString("Beams:eCM = 13600.");
        pythia.readString("Tune:pp = 14");
        // pythia.readString("411:onMode=off");
        // pythia.readString("411:onIfAny=13");
        pythia.settings.parm("PhaseSpace:pTHatMin", binEdges[iBin]);
        pythia.settings.parm("PhaseSpace:pTHatMax", binEdges[iBin + 1]);
        pythia.init();

        int eventCount = 0;

        hardPtPart->Reset();
        cout<<"--------------------New Bin--------------------"<<endl;

        for (int iEvent = 0; iEvent < genEvents; ++iEvent) {
            //cout<<"\nEvent No: "<<iEvent+1<<endl;

            if (!pythia.next()) {continue;}

            double pTHat  = pythia.info.pTHat();

            if(pythia.info.isNonDiffractive()){code = pythia.info.codeSub();}
            else{code = pythia.info.code();}

            if (softQCD && iBin < 1 && pythia.info.isNonDiffractive()
            && pTHat > binEdges[iBin+1]) continue;

            if (pTHat < binEdges[iBin]) continue;

            hardPtPart->Fill(pTHat);
            eventCount++;

            std::vector<vector<double>> elecs;
            std::vector<vector<double>> muons;

            bool foundE = false;
            bool foundMu = false;
            
            int quarkID = 0;
            //cout << "====START OF NEW EVENT====" << endl;
            
            //if(121<=code && code<=124) //only checks quarks if they come from a hardest process which produces a heavy flavor pair
            //{
                for (int i = 0; i < pythia.event.size(); ++i) {
                    //cout <<"---------------------------------------- New Event"<<endl;
                    
                    int particleStatus = pythia.event[i].status();
                    int particleID = pythia.event[i].id();
                    double particlePt = pythia.event[i].pT();
                    double particleEta = pythia.event[i].eta();
                    double particlePhi = pythia.event[i].phi();
                    double particleTheta = pythia.event[i].theta();

                    if(particleStatus==-23)
                    {
                        quarkID = particleID;
                    }

                    if(particleStatus>0 && std::abs(particleID)==11)
                    {             
                        foundE = true;
                        elecs.push_back({(double)particleID,particlePt,particleEta,particlePhi,particleTheta});
                    }

                     if(particleStatus>0 && std::abs(particleID)==13)
                    {             
                        foundMu = true;
                        muons.push_back({(double)particleID,particlePt,particleEta,particlePhi,particleTheta});
                    }
                }
                
                if(foundE==true && foundMu==true){
                    for(vector<double> elec: elecs){
                        for(vector<double> muon: muons){
                            if(elec[0]*muon[0] < 1){

                                double dPhiEmu = elec[3] - muon[3];
                                if(dPhiEmu<0) {dPhiEmu = dPhiEmu + 2.*M_PI;}
                                else if(dPhiEmu>2*M_PI) {dPhiEmu = dPhiEmu - 2.*M_PI;}

                                
                                emuTuples[iBin]->Fill(quarkID,elec[0],elec[1],elec[2],elec[3],elec[4],muon[0],muon[1],muon[2],muon[3],muon[4],dPhiEmu);
                            }
                        }
                    }
                }
            //}
 
        }

        // cross-section for the bin
        //double luminocity_hard = N_events/(pythia.info.sigmaGen()*pow(10,9));
        
        double luminocity_hard = (pythia.info.weightSum())/(pythia.info.sigmaGen()*pow(10,9));

        binLuminocity[iBin] = luminocity_hard;
        binEvCount[iBin] = eventCount;

        weightSums[iBin] = pythia.info.weightSum();
        sigmaGens[iBin] = pythia.info.sigmaGen();

        hardPtPart->Scale(1/luminocity_hard,"width");

        // add to final distribution
        hardPt->Add(hardPtPart);
    }

    for (int i = 0; i < nBins; ++i) 
    {   
        emuTuples[i]->Write(Form("emu%d", i));
    }

    outFile->WriteObject(&binLuminocity, "luminocities");
    outFile->WriteObject(&binEvCount, "eventCounts");
    outFile->WriteObject(&weightSums, "weightSums");
    outFile->WriteObject(&sigmaGens, "sigmaGens");
    
    // Total Cross Section
    TCanvas *canvasTotal = new TCanvas("total_sigma","total_sigma");
    gPad->SetLogy();

    hardPt->SetLineColor(1);
    hardPt->Draw();

    canvasTotal->Write();



    delete outFile;

    return 0;
}
