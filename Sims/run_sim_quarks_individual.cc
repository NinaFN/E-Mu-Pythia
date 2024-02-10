#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

using namespace Pythia8;

//script to save information about individual heavy quarks and their decay products

int main() {
    // Turn SoftQCD on/off
    bool softQCD = false;

    // Pythia object
    Pythia pythia;

    // ROOT file for histograms
    TFile* outFile = new TFile("sim_tuples_136_semilep_quarks_individual.root", "RECREATE");

    // pTHat bins
    int nBins;
    const double* binEdges;
    const int* binWidths;
    if (softQCD) {
        nBins = 8;
        static const double tempArray[9] = {0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0, 120.0};
        //static const int tempArrayW[6] = {16,10,10,14,20,30};

        binEdges = &tempArray[0];
        //binWidths = &tempArrayW[0];
    } else {
        nBins = 7;
        static const double tempArray[8] = {5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 40.0, 60.0};

        binEdges = &tempArray[0];
        //binWidths = &tempArrayW[0];
    }


    // Histograms
    // Total Cross Section
    TH1F *hardPt = new TH1F("HardQCD:All","Process Total Cross-Section;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 50, 5.0, 60.0);
    TH1F *hardPtPart = new TH1F("hardQCD_part","", 50, 5.0, 60.0);

    // HF Cross Sections
    vector<TNtuple*> qeTuples(nBins);
    vector<TNtuple*> qmTuples(nBins);
    

    vector<double> binLuminocity(nBins); // luminocity from generated process sigma to calculate cross-sections
    vector<double> binEvCount(nBins);
    vector<double> weightSums(nBins);
    vector<double> sigmaGens(nBins);

    for (int i = 0; i < nBins; ++i) {

        qmTuples[i] = new TNtuple("quark_muon", "quark_muon", "decayFlag:semilepFlag:ptHat:idQ:etaQ:phiQ:thetaQ:idL:ptL:etaL:phiL:thetaL");        
        qeTuples[i] = new TNtuple("quark_elec", "quark_elec", "decayFlag:semilepFlag:ptHat:idQ:etaQ:phiQ:thetaQ:idL:ptL:etaL:phiL:thetaL");        

    }

    // Number of events to generate per bin.
    int N_events = 2000000;
    

    int decayMap = 0;
    //int decayMap2 = 0;

    int code = 0;


    int genEvents;

    // run events for each ptHat bin 
    for (int iBin = 0; iBin < nBins; ++iBin) {
        if (softQCD && iBin < 1) {
            pythia.readString("HardQCD:all = off");
            pythia.readString("HardQCD:hardccbar = off");
            pythia.readString("HardQCD:hardbbbar = off");
            pythia.readString("SoftQCD:nonDiffractive = on");
            genEvents = N_events*30;
        } else {
            // set pythia initialization variables
            genEvents = N_events;
            pythia.readString("HardQCD:all = off");
            pythia.readString("HardQCD:hardccbar = on");
            pythia.readString("HardQCD:hardbbbar = on");
            pythia.readString("SoftQCD:nonDiffractive = off");

           /*  pythia.readString("411:onIfAny=11 13 ");
            pythia.readString("421:onIfAny=11 13 ");
            pythia.readString("413:onIfAny=11 13 ");
            pythia.readString("423:onIfAny=11 13 ");
            pythia.readString("415:onIfAny=11 13 ");
            pythia.readString("425:onIfAny=11 13 ");
            pythia.readString("431:onIfAny=11 13 ");
            pythia.readString("433:onIfAny=11 13 ");
            pythia.readString("435:onIfAny=11 13 ");

            pythia.readString("511:onIfAny=11 13 ");
            pythia.readString("521:onIfAny=11 13 ");
            pythia.readString("513:onIfAny=11 13 ");
            pythia.readString("523:onIfAny=11 13 ");
            pythia.readString("515:onIfAny=11 13 ");
            pythia.readString("525:onIfAny=11 13 ");
            pythia.readString("531:onIfAny=11 13 ");
            pythia.readString("533:onIfAny=11 13 ");
            pythia.readString("535:onIfAny=11 13 ");
            pythia.readString("541:onIfAny=11 13 ");
            pythia.readString("543:onIfAny=11 13 ");
            pythia.readString("545:onIfAny=11 13 "); */
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

            int decayElec = -1;
            int decayMuon = -1;
            int barrelHit = -1;
            int forwardHit = -1;
            
            int semilepE = -1;
            int semilepMu = -1;

            int elecID = 0;
            double elecPt = 0;
            double elecEta = 0;
            double elecPhi = 0;
            double elecTheta = 0;

            int muonID = 0;
            double muonPt = 0;
            double muonEta = 0;
            double muonPhi = 0;
            double muonTheta = 0;

            int quarkCount = 0;

            
            //cout << "====START OF NEW EVENT====" << endl;
            
            if(121<=code && code<=124) //only checks quarks if they come from a hardest process which produces a heavy flavor pair
            {
                for (int i = 0; i < pythia.event.size(); ++i) {
                    //cout <<"---------------------------------------- New Event"<<endl;
                    
                    int particleStatus = pythia.event[i].status();
                    
                    if (particleStatus==-23)
                    {
                        quarkCount++;

                        int particleID = pythia.event[i].id();
                        double particlePt = pythia.event[i].pT();
                        double particleEta = pythia.event[i].eta();
                        double particlePhi = pythia.event[i].phi();
                        double particleTheta = pythia.event[i].theta();

                        std::vector<int> childList = pythia.event[i].daughterListRecursive();
                        std::sort(childList.begin(),childList.end());  
                        childList.erase( unique( childList.begin(), childList.end() ), childList.end() );


                        for(int child: childList)
                        {
                            if(std::abs(pythia.event[child].id())==11 && pythia.event[child].status()>0)
                            {
                                decayElec = 1;
                            
                                elecID = pythia.event[child].id();
                                elecPt = pythia.event[child].pT();
                                elecEta = pythia.event[child].eta();
                                elecPhi = pythia.event[child].phi();
                                elecTheta = pythia.event[child].theta();

                                
                                for(int mother: pythia.event[child].motherList())
                                {
                                    if((std::abs(pythia.event[mother].id())>=411 && std::abs(pythia.event[mother].id()<=435))||(std::abs(pythia.event[mother].id())>=511 && std::abs(pythia.event[mother].id()<=545)))
                                    {
                                        semilepE=1;
                                        break;
                                    }
                                }
                                
                                //int parentMes = pythia.event[pythia.event[child].mother1()].id();
                                //if((parentMes>=411 && parentMes<=435)||(parentMes>=511 && parentMes<=545))
                                //{semilepE = 1;}

                                if(semilepE==1){break;}

                                /* if(elecEta<=0.9 && elecEta>=-0.9)
                                {
                                    barrelHit=1;
                                    break; //prioritize storing info from leptons that do decay to correct eta range
                                } */
                            }
                        }

                        for(int child: childList)
                        {
                            if(std::abs(pythia.event[child].id())==13 && pythia.event[child].status()>0)
                            {
                                decayMuon = 1;
                            
                                muonID = pythia.event[child].id();
                                muonPt = pythia.event[child].pT();
                                muonEta = pythia.event[child].eta();
                                muonPhi = pythia.event[child].phi();
                                muonTheta = pythia.event[child].theta();

                                for(int mother: pythia.event[child].motherList())
                                {
                                    if((std::abs(pythia.event[mother].id())>=411 && std::abs(pythia.event[mother].id()<=435))||(std::abs(pythia.event[mother].id())>=511 && std::abs(pythia.event[mother].id()<=545)))
                                    {
                                        semilepMu=1;
                                        break;
                                    }
                                }

                                if(semilepMu==1){break;}

                                /* if(muonEta>=-4.0 && muonEta<=-2.5)
                                {
                                    forwardHit=1;
                                    break; //prioritize storing info from leptons that do decay to correct eta range
                                } */
                            }
                            
                        }
                        
                        //pythia.event.list(true);
                        
                        qeTuples[iBin]->Fill(decayElec,semilepE,pTHat,particleID,particleEta,particlePhi,particleTheta,elecID,elecPt,elecEta,elecPhi,elecTheta);
                        qmTuples[iBin]->Fill(decayMuon,semilepMu,pTHat,particleID,particleEta,particlePhi,particleTheta,muonID,muonPt,muonEta,muonPhi,muonTheta);

                    }

                    if(quarkCount>2){break;}
                }
            }
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
        qeTuples[i]->Write(Form("qe%d", i));
        qmTuples[i]->Write(Form("qm%d", i));
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
