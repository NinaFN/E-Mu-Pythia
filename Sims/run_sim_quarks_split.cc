#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

//script to save information about each quark in the pair for later comparison of Q and Qbar kinematic props

using namespace Pythia8;


int main() {
    // Turn SoftQCD on/off
    bool softQCD = false;

    // Pythia object
    Pythia pythia;

    // ROOT file for histograms
    TFile* outFile = new TFile("sim_tuples_quarks_split.root", "RECREATE");

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
        
        //static const int tempArrayW[5] = {10,10,14,20,30};

        binEdges = &tempArray[0];
        //binWidths = &tempArrayW[0];
    }


    // Histograms
    // Total Cross Section
    TH1F *hardPt = new TH1F("HardQCD:All","Process Total Cross-Section;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 50, 5.0, 60.0);
    TH1F *hardPtPart = new TH1F("hardQCD_part","", 50, 5.0, 60.0);

    // HF Cross Sections
    vector<TNtuple*> qpTuples(nBins);
    

    vector<double> binLuminocity(nBins); // luminocity from generated process sigma to calculate cross-sections
    vector<double> binEvCount(nBins);
    vector<double> weightSums(nBins);
    vector<double> sigmaGens(nBins);

    for (int i = 0; i < nBins; ++i) {

        qpTuples[i] = new TNtuple("quark_pair", "quark_pair", "ptHat:ptQ1:ptQ2:idQ1:idQ2:yQ1:yQ2:etaQ1:etaQ2:phiQ1:phiQ2:thetaQ1:thetaQ2");        

    }

    // Number of events to generate per bin.
    int N_events = 100000;
    

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
        }

        pythia.readString("Beams:eCM = 13700.");
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

            
            //cout << "====START OF NEW EVENT====" << endl;
            
            if(121<=code && code<=124) //only checks quarks if they come from a hardest process which produces a heavy flavor pair
            {
                for (int i = 0; i < pythia.event.size(); ++i) {
                    int particleStatus = pythia.event[i].status();
                    
                    if (particleStatus==-23)
                    {
                        //cout <<"---------------------------------------- New -23 Particle"<<endl;

                        int q1ID = pythia.event[i].id();
                        double q1Pt = pythia.event[i].pT();
                        double q1Y = pythia.event[i].y();
                        double q1Eta = pythia.event[i].eta();
                        double q1Phi = pythia.event[i].phi();
                        double q1Theta = pythia.event[i].theta();

                        int q2ID = pythia.event[i+1].id();
                        double q2Pt = pythia.event[i+1].pT();
                        double q2Y = pythia.event[i+1].y();
                        double q2Eta = pythia.event[i+1].eta();
                        double q2Phi = pythia.event[i+1].phi();
                        double q2Theta = pythia.event[i+1].theta();

                        
                        
                        //ptHat:ptQ1:ptQ2:idQ1:idQ2:yQ1:yQ2:etaQ1:etaQ2:phiQ1:phiQ2:thetaQ1:thetaQ2
                        qpTuples[iBin]->Fill(pTHat,q1Pt,q2Pt,q1ID,q2ID,q1Y,q2Y,q1Eta,q2Eta,q1Phi,q2Phi,q1Theta,q2Theta);
                        break;
                    }

                    //if(quarkCount>2){break;}
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
        qpTuples[i]->Write(Form("qp%d", i));
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
