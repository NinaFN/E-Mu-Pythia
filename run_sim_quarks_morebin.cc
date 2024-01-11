#include "Pythia8/Pythia.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

using namespace Pythia8;

//script to save event categorization information based on the decays of the heavy quark pair

int main() {
    // Turn SoftQCD on/off
    bool softQCD = false;

    // Pythia object
    Pythia pythia;

    // ROOT file for histograms
    TFile* outFile = new TFile("sim_tuples_136_semilep_quarks_pairprod.root", "RECREATE");

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
    vector<TNtuple*> quarkTuples(nBins);
    vector<TNtuple*> pairTuples(nBins);
    

    vector<double> binLuminocity(nBins); // luminocity from generated process sigma to calculate cross-sections
    vector<double> binEvCount(nBins);
    vector<double> weightSums(nBins);
    vector<double> sigmaGens(nBins);

    for (int i = 0; i < nBins; ++i) {
        //Decay Map:
            // 0: quark pair does not decay to RL or SL
            // 1: a quark decays to mu
            // 2: a quark decays to e
            // 3: decays to both e and mu happen
            // 4: quark pair decays to both e and mu, each lepton comes from seperate quark
            // 5: quark pair decays to both e and mu, each lepton comes from seperate quark, opp sign (matches quark sign)
            // 6: spicy 5 where both decays are semileptonic ("real pairs")

        quarkTuples[i] = new TNtuple("quarks", "quarks", "event:id:pt1:ptHat:eta1:eta2:phi1:phi2:decayMap");

    }

    // Number of events to generate per bin.
    int N_events = 2000000;
    

    int decayMap = 0;
    int decayMapPair = 0;
    int decayMap1 = 0;
    int decayMap2 = 0;

    int code = 0;
    int genEvents;


    // run events for each ptHat bin 
    for (int iBin = 0; iBin < nBins; ++iBin) {
        if (softQCD && iBin < 1) {
            pythia.readString("HardQCD:all = off");
            pythia.readString("HardQCD:hardccbar = off");
            pythia.readString("HardQCD:hardbbbar = off");
            pythia.readString("SoftQCD:nonDiffractive = on");
            genEvents = N_events*5;
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

            //eta, charge, semilep>
            std::vector<std::tuple<double, double, bool, int>> elec_q1;
            std::vector<std::tuple<double, double, bool, int>> elec_q2;
            std::vector<std::tuple<double, double, bool, int>> muon_q1;
            std::vector<std::tuple<double, double, bool, int>> muon_q2;

            double pTHat  = pythia.info.pTHat();

            bool decayElec_q1 = false;
            bool decayElec_q2 = false;
            bool decayMuon_q1 = false;
            bool decayMuon_q2 = false;

            bool barrelHit = false;
            bool forwardHit = false;

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
                    //cout <<"---------------------------------------- New Event"<<endl;
                    decayMap = 0;
                    decayMapPair = 0;
                    decayMap1 = 0;
                    decayMap2 = 0;
                    //decayMap2 = 0;
                    
                    int particleID = pythia.event[i].id();
                    int particleStatus = pythia.event[i].status();
                    double particlePt = pythia.event[i].pT();
                    double particleEta = pythia.event[i].eta();
                    
                    

                    if (particleStatus==-23)
                    {
                        int chargeQ1 = pythia.event[i].charge(); 
                        int chargeQ2 = pythia.event[i+1].charge(); 


                        std::vector<int> childList1 = pythia.event[i].daughterListRecursive();
                        std::sort(childList1.begin(),childList1.end());  
                        childList1.erase( unique( childList1.begin(), childList1.end() ), childList1.end() );

                        for(int child: childList1)
                        {
                            if(std::abs(pythia.event[child].id())==11 && pythia.event[child].status()>0)
                            {   
                                decayElec_q1 = true;

                                double eta = pythia.event[child].eta();
                                double charge = pythia.event[child].charge();

                                bool semilep = false;
                                //int parentMes1 = pythia.event[pythia.event[child].mother1()].id();
                                //int parentMes2 = pythia.event[pythia.event[child].mother2()].id();

                                for(int mother: pythia.event[child].motherList())
                                {
                                    if((std::abs(pythia.event[mother].id())>=411 && std::abs(pythia.event[mother].id()<=435))||(std::abs(pythia.event[mother].id())>=511 && std::abs(pythia.event[mother].id()<=545)))
                                    {
                                        semilep=true;
                                        break;
                                    }
                                }
                                //if((parentMes1>=411 && parentMes1<=435)||(parentMes1>=511 && parentMes1<=545)){semilep=true;}
                                //else if((parentMes2>=411 && parentMes2<=435)||(parentMes2>=511 && parentMes2<=545)){semilep=true;}

                                elec_q1.push_back(std::tuple<double,double,bool,int>{eta,charge,semilep,chargeQ1});
        
                                //cout<<"child particle found: "<<pythia.event[child].id()<<endl;
                            }

                            if(std::abs(pythia.event[child].id())==13 && pythia.event[child].status()>0)
                            {   
                                decayMuon_q1 = true;

                                double eta = pythia.event[child].eta();
                                double charge = pythia.event[child].charge();

                                bool semilep = false;
                                for(int mother: pythia.event[child].motherList())
                                {
                                    if((std::abs(pythia.event[mother].id())>=411 && std::abs(pythia.event[mother].id()<=435))||(std::abs(pythia.event[mother].id())>=511 && std::abs(pythia.event[mother].id()<=545)))
                                    {
                                        semilep=true;
                                        break;
                                    }
                                }

                                muon_q1.push_back(std::tuple<double,double,bool,int>{eta,charge,semilep,chargeQ1});


                                //cout<<"child particle found: "<<pythia.event[child].id()<<endl;
                            }
                        }
                        
                        std::vector<int> childList2 = pythia.event[i+1].daughterListRecursive();
                        std::sort(childList2.begin(),childList2.end());  
                        childList2.erase( unique( childList2.begin(), childList2.end() ), childList2.end() );

                        for(int child: childList2)
                        {
                            if(std::abs(pythia.event[child].id())==11 && pythia.event[child].status()>0)
                            {
                                decayElec_q2 = true;

                                double eta = pythia.event[child].eta();
                                double charge = pythia.event[child].charge();


                                bool semilep = false;
                                for(int mother: pythia.event[child].motherList())
                                {
                                    if((std::abs(pythia.event[mother].id())>=411 && std::abs(pythia.event[mother].id()<=435))||(std::abs(pythia.event[mother].id())>=511 && std::abs(pythia.event[mother].id()<=545)))
                                    {
                                        semilep=true;
                                        break;
                                    }
                                }

                                elec_q2.push_back(std::tuple<double,double,bool,int>{eta,charge,semilep,chargeQ2});
                                
                                //cout<<"child particle found: "<<pythia.event[child].id()<<endl;
                            }

                            if(std::abs(pythia.event[child].id())==13 && pythia.event[child].status()>0)
                            {
                                decayMuon_q2 = true;

                                double eta = pythia.event[child].eta();
                                double charge = pythia.event[child].charge();


                                bool semilep = false;
                                for(int mother: pythia.event[child].motherList())
                                {
                                    if((std::abs(pythia.event[mother].id())>=411 && std::abs(pythia.event[mother].id()<=435))||(std::abs(pythia.event[mother].id())>=511 && std::abs(pythia.event[mother].id()<=545)))
                                    {
                                        semilep=true;
                                        break;
                                    }
                                }


                                muon_q2.push_back(std::tuple<double,double,bool,int>{eta,charge,semilep,chargeQ2});
                                
                                //cout<<"child particle found: "<<pythia.event[child].id()<<endl;
                            }
                        }

                        if(decayMuon_q1==true || decayMuon_q2==true){decayMap = 1;} //one quark decays to muon
                        if(decayElec_q1==true || decayElec_q2==true){decayMap = 2;} //one quark decays to electron
                        if((decayElec_q1==true || decayElec_q2==true) && (decayMuon_q1==true || decayMuon_q2==true)){decayMap = 3;} //both electron and muon come from hf decay
                        
                        
                        if((decayElec_q1==true && decayMuon_q2==true))
                        {
                            //cout<<"yes1"<<endl;
                            decayMap1=4; //back to back decay
                            //etaFlag = -1;
                            
                            for (std::tuple<double, double, bool, int> elec: elec_q1){
                                for (std::tuple<double, double, bool, int> muon: muon_q2){
                                    if(std::get<1>(elec)*std::get<1>(muon) == -1){

                                        decayMap1=5;

                                        //if((std::get<0>(elec)<=0.9 && std::get<0>(elec)>=-0.9)&&(std::get<0>(muon)<=-2.5 && std::get<0>(muon)>=-4))
                                        //{etaFlag = 1;}

                                        if(std::get<2>(elec)==true && std::get<2>(muon)==true){
                                        //if(std::get<1>(elec)*std::get<3>(elec)==1 && std::get<1>(muon)*std::get<3>(muon)==1){
                                            decayMap1=6;
                                            //cout<<"!!!"<<endl;
                                            break; //prioritize saving semilep pairs
                                        }
                                    }
                                }
                                //cout<<decayMap1<<endl;
                                if(decayMap1==6){break;}
                            }         
                        }
                        //cout<<endl;
                
                        if((decayElec_q2==true && decayMuon_q1==true))
                        {
                            //cout<<"yes2"<<endl;
                            decayMap2=4; //back to back decay
                            //etaFlag = -1;
                            
                            for (std::tuple<double, double, bool, int> elec: elec_q2){
                                for (std::tuple<double, double, bool, int> muon: muon_q1){
                                    if(std::get<1>(elec)*std::get<1>(muon) == -1){
                                        decayMap2=5;

                                        //if((std::get<0>(elec)<=0.9 && std::get<0>(elec)>=-0.9)&&(std::get<0>(muon)<=-2.5 && std::get<0>(muon)>=-4))
                                        //{etaFlag = 1;}

                                        //if(std::get<2>(elec)==true && std::get<2>(muon)==true){
                                        if(std::get<1>(elec)*std::get<3>(elec)==1 && std::get<1>(muon)*std::get<3>(muon)==1){
                                            decayMap2=6;
                                            ///cout<<"!!!"<<endl;

                                            break; //prioritize saving semilep pairs
                                        }
                                    }
                                }
                                //cout<<decayMap2<<endl;
                                if(decayMap2==6){break;}
                            }         
                        }
                        
                        
                        //if(decayMap1>0 && decayMap2>0){cout<<"both"<<endl;}

                        if(decayMap1>=decayMap2){decayMapPair=decayMap1;}
                        else{decayMapPair=decayMap2;}

                        if(decayMapPair>decayMap){decayMap=decayMapPair;}

                        //cout<<endl<<decayMap<<endl;
                            
                        //event:id:pt1:ptHat:eta1:eta2:phi1:phi2:decayMap        
                        quarkTuples[iBin]->Fill(iEvent,particleID,particlePt,pTHat,particleEta,pythia.event[i+1].eta(),pythia.event[i].phi(),pythia.event[i+1].phi(),decayMap);                    
                        break;
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
        quarkTuples[i]->Write(Form("quark%d", i));
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
