#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

using namespace Pythia8;

//script to save detailed information of the true e-mu pairs

int main() {
    // Turn SoftQCD on/off
    bool softQCD = false;

    // Pythia object
    Pythia pythia;

    // ROOT file for histograms
    TFile* outFile = new TFile("sim_tuples_136_semilep_paircheck.root", "RECREATE");

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
    }


    // Histograms
    // Total Cross Section
    TH1F *hardPt = new TH1F("HardQCD:All","Process Total Cross-Section;#hat{p}_{T} (GeV/c);#frac{d#sigma}{dp_{T}} (pb/GeV/c)", 55, 5.0, 60.0);
    TH1F *hardPtPart = new TH1F("hardQCD_part","", 55, 5.0, 60.0);

    // HF Cross Sections
    vector<TNtuple*> qeTuples(nBins);
    vector<TNtuple*> qmTuples(nBins);
    vector<TNtuple*> qmTuplesSemilep(nBins);
    

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
        qmTuples[i] = new TNtuple("quark_muon", "quark_muon", "decayFlag:pairDecFlag:pairEta:pairPt:pairPhi:ptHat:idQ:ptL:etaQ:etaL:phiQ:phiL:dPhi");        
        qmTuplesSemilep[i] = new TNtuple("quark_muon_semilep", "quark_muon_semilep", "semilepFlag:pairSemilepFlag:pairEta:pairPt:pairPhi:ptHat:idQ:ptL:etaQ:etaL:phiQ:phiL:dPhi");        

        qeTuples[i] = new TNtuple("quark_elec", "quark_elec", "decayFlag:semilepFlag:pairDecFlag:pairSemilepFlag:pairEta:pairPt:pairPhi:ptHat:idQ:ptL:etaQ:etaL:phiQ:phiL");        

    }

    // Number of events to generate per bin.
    int N_events = 2000000;

    int code = 0;


    int genEvents;

    // run events for each ptHat bin 
    for (int iBin = 0; iBin < nBins; ++iBin) {
        if (softQCD && iBin < 1) {
            pythia.readString("HardQCD:all = off");
            pythia.readString("HardQCD:hardccbar = off");
            pythia.readString("HardQCD:hardbbbar = off");
            pythia.readString("SoftQCD:nonDiffractive = on");
            genEvents = N_events*10;
        } else {
            // set pythia initialization variables
            genEvents = N_events;
            pythia.readString("HardQCD:all = off");
            pythia.readString("HardQCD:hardccbar = on");
            pythia.readString("HardQCD:hardbbbar = on");
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

            std::vector<vector<double>> quarkElecPos;
            std::vector<vector<double>> quarkMuonPos;

            int muonEtaFound = -1;
            int elecEtaFound = -1;
            
            //cout << "====START OF NEW EVENT====" << endl;
           
            
            
            if(121<=code && code<=124) //only checks quarks if they come from a hardest process which produces a heavy flavor pair
            {
                for (int i = 0; i < pythia.event.size(); ++i) {
                    //cout <<"---------------------------------------- New Event"<<endl;
                    
                    int particleStatus = pythia.event[i].status();
                    int particleID = pythia.event[i].id();
                    double particlePt = pythia.event[i].pT();
                    double particleEta = pythia.event[i].eta();
                    double particlePhi = pythia.event[i].phi();
                    double particleTheta = pythia.event[i].theta();
                    
                    if (particleStatus==-23) //quarks coming from hardest process, our relevant parent particles
                    {
                        std::vector<int> childList = pythia.event[i].daughterListRecursive();
                        std::sort(childList.begin(),childList.end());  
                        childList.erase( unique( childList.begin(), childList.end() ), childList.end() );

                        for(int child: childList)
                        {

                            if(std::abs(pythia.event[child].id())==11 && pythia.event[child].status()>0){
                                //save location of child in event record, along with particle information about the parent quark for storage
                                quarkElecPos.push_back({(double)child, (double)pythia.event[child].id(), (double)particleID, particlePt, particleEta, particlePhi, particleTheta, pythia.event[child].eta(),pythia.event[child].pT(),pythia.event[child].phi()});
                            } 

                            if(std::abs(pythia.event[child].id())==13 && pythia.event[child].status()>0){
                                quarkMuonPos.push_back({(double)child, (double)pythia.event[child].id(), (double)particleID, particlePt, particleEta, particlePhi, particleTheta, pythia.event[child].eta(),pythia.event[child].pT(),pythia.event[child].phi()});
                            }
                        }   
                    }

/*                     if(particleStatus>0 && std::abs(particleID)==11)
                    {                       
                        int decayFlag = -1;
                        int semilepFlag = -1;
                        int pairDecFlag = -1;
                        int pairSemilepFlag = -1;
                        double pairEtaFlag = -1;
                        double pairEta = 10; //outside of range of alice detector
                        double pairPt = 0;
                        double pairPhi = 0;

                        int parentID = 0;
                        double parentPt = 0;
                        double parentEta = 0;
                        double parentPhi = 0;
                        double parentTheta = 0;

                        

                        for(vector<double> child: quarkElecPos)
                        {
                            if((int)child[0]==i) //checks if found final state electron is in the list of heavy quark decay particles
                            {
                                decayFlag = 1;
                                
                                for(int mother: pythia.event[i].motherList())
                                {
                                    cout<<"mother: "<<pythia.event[mother].id()<<endl;
                                    if((std::abs(pythia.event[mother].id())>=411 && std::abs(pythia.event[mother].id()<=435))||(std::abs(pythia.event[mother].id())>=511 && std::abs(pythia.event[mother].id()<=545)))
                                    {
                                        semilepFlag=1;
                                        cout<<"!!! E"<<endl;
                                        break;
                                    }
                                }

                                //int parentMes = pythia.event[pythia.event[i].mother1()].id();
                                //if((parentMes>=411 && parentMes<=435)||(parentMes>=511 && parentMes<=545)){semilepFlag=1;}

                                //save parent information from vector
                                parentID = (int)child[2];
                                parentPt = (int)child[3];
                                parentEta = child[4];
                                parentPhi = child[5];
                                parentTheta = child[6];

                                if (!quarkMuonPos.empty()) //only checks for pair partner if there are decay muons present
                                {
                                    for(vector<double> muon: quarkMuonPos)
                                    {
                                        //check to make sure that the muon is oppositely signed to the electron
                                        //further check to make sure that the muon's parent quark is oppositely signed to the electron's parent quark
                                        //this ensures that the electron and muon did not come from the same parent (back to back decay)

                                        if((int)muon[1]*particleID<1 && (int)muon[2]*parentID<1)
                                        {
                                            pairDecFlag=1;
                                            cout<<"!!!"<<endl;
                                            cout<<"electron: "<<i<<endl;
                                            
                                            pythia.event.list();

                                            for(int mother: pythia.event[(int)muon[0]].motherList())
                                            {
                                                cout<<"muon mother: "<<pythia.event[mother].id()<<endl;
                                                if((std::abs(pythia.event[mother].id())>=411 && std::abs(pythia.event[mother].id()<=435))||(std::abs(pythia.event[mother].id())>=511 && std::abs(pythia.event[mother].id()<=545)))
                                                {
                                                    pairSemilepFlag=1;
                                                    cout<<"!!! M"<<endl;
                                                    break;
                                                }
                                            }



                                            //int muonParentMes = pythia.event[pythia.event[(int)muon[0]].mother1()].id();
                                            //if((muonParentMes>=411 && muonParentMes<=435)||(muonParentMes>=511 && muonParentMes<=545)){pairSemilepFlag=1;}
                                            
                                            pairEta = muon[7];
                                            pairPt = muon[8];
                                            pairPhi = muon[9];
                                            if(pairEta<=-2.5 && pairEta>= -4.0){pairEtaFlag = 1;}
                                        }
                                    }
                                }
                                break; //don't keep searching through decay child list after electron has been found
                            }
                        }
                        //decayFlag:semilepFlag:pairDecFlag:pairSemilepFlag:pairEta:pairPt:pairPhi:ptHat:idQ:ptL:etaQ:etaL:phiQ:phiL
                        qeTuples[iBin]->Fill(decayFlag,semilepFlag,pairDecFlag,pairSemilepFlag,pairEta,pairPt,pairPhi,pTHat,parentID,particlePt,parentEta,particleEta,parentPhi,particlePhi);
                    } */

 
                    
                    
                    if(particleStatus>0 && std::abs(particleID)==13)
                    {                       
                        int decayFlag = -1;
                        int semilepFlag = -1;
                        int pairDecFlag = -1;
                        int pairSemilepFlag = -1;
                        int pairEtaFlag = -1;
                        double pairEta = 10; //outside of range of alice detector
                        double pairPt = 0;
                        double pairPhi = 0;
                        double dPhiEmu = -10;

                        int parentID = 0;
                        double parentPt = 0;
                        double parentEta = 0;
                        double parentPhi = 0;
                        double parentTheta = 0;

                        

                        for(vector<double> child: quarkMuonPos)
                        {
                            if((int)child[0]==i) //checks if found final state muon is in the list of heavy quark decay particles
                            {
                                decayFlag = 1;

                                for(int mother: pythia.event[i].motherList())
                                {
                                    //cout<<"mother: "<<pythia.event[mother].id()<<endl;
                                    if((std::abs(pythia.event[mother].id())>=411 && std::abs(pythia.event[mother].id()<=435))||(std::abs(pythia.event[mother].id())>=511 && std::abs(pythia.event[mother].id()<=545)))
                                    {
                                        semilepFlag=1;
                                       // cout<<"!!! M"<<endl;
                                        break;
                                    }
                                }
                                //save parent information from vector
                                parentID = (int)child[2];
                                parentPt = (int)child[3];
                                parentEta = child[4];
                                parentPhi = child[5];
                                parentTheta = child[6];

                                if (!quarkElecPos.empty()) //only checks for pair partner if there are decay electrons present
                                {
                                    for(vector<double> elec: quarkElecPos)
                                    {
                                        if(elec[1]*particleID<1 && elec[2]*parentID<1)
                                        {
                                            pairDecFlag=1;

                                            for(int mother: pythia.event[(int)elec[0]].motherList())
                                            {
                                               // cout<<"elec mother: "<<pythia.event[mother].id()<<endl;
                                                if((std::abs(pythia.event[mother].id())>=411 && std::abs(pythia.event[mother].id()<=435))||(std::abs(pythia.event[mother].id())>=511 && std::abs(pythia.event[mother].id()<=545)))
                                                {
                                                    pairSemilepFlag=1;
                                                   // cout<<"!!! E"<<endl;
                                                    break;
                                                }
                                            }

                                            pairEta = elec[7];
                                            pairPt = elec[8];
                                            pairPhi = elec[9];
                                            if(pairEta<=0.9 && pairEta>=-0.9){pairEtaFlag = 1;}

                                            dPhiEmu = pairPhi - particlePhi;
                                            if(dPhiEmu<0) {dPhiEmu = dPhiEmu + 2.*M_PI;}
                                            else if(dPhiEmu>2*M_PI) {dPhiEmu = dPhiEmu - 2.*M_PI;}

                                            if(pairSemilepFlag==1){break;}
                                        }
                                    }
                                }
                                if(semilepFlag==1){break;} //don't keep searching through decay child list after muon has been found
                            }
                        }
                        //decayFlag:semilepFlag:pairDecFlag:pairSemilepFlag:pairEta:pairPt:pairPhi:ptHat:idQ:ptL:etaQ:etaL:phiQ:phiL
                        qmTuples[iBin]->Fill(decayFlag,pairDecFlag,pairEta,pairPt,pairPhi,pTHat,parentID,particlePt,parentEta,particleEta,parentPhi,particlePhi,dPhiEmu);                    
                        qmTuplesSemilep[iBin]->Fill(semilepFlag,pairSemilepFlag,pairEta,pairPt,pairPhi,pTHat,parentID,particlePt,parentEta,particleEta,parentPhi,particlePhi,dPhiEmu);                    
                        
                    }

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
        //qeTuples[i]->Write(Form("qe%d", i));
        qmTuples[i]->Write(Form("qm%d", i));
        qmTuplesSemilep[i]->Write(Form("qms%d", i));
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
