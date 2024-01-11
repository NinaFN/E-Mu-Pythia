#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

//script to calculate and plot the ratios used to choose the CB parent quark eta boundaries in the approximation method

void read_sim_eta_range_opt_e() {
    TFile *f = TFile::Open("sim_tuples_136_semilep_quarks_individual.root","READ");
    
    
    vector<double> *binLuminocity;
    vector<double> *binEvCounts;
    vector<double> *weightSums;
    vector<double> *sigmaGens;

    vector<TNtuple*> qeTuples(6); //should be count of bins but 7 is easier
    vector<TNtuple*> qmTuples(6);

    f->GetObject("luminocities",binLuminocity);
    f->GetObject("eventCounts",binEvCounts);
    f->GetObject("weightSums",weightSums);
    f->GetObject("sigmaGens",sigmaGens);

    int binCount = 0;

    for(std::vector<double>::iterator it = binLuminocity->begin(); it != binLuminocity->end(); ++it){
        qeTuples[binCount] = (TNtuple*)f->Get(Form("qe%d", binCount));
        qmTuples[binCount] = (TNtuple*)f->Get(Form("qm%d", binCount));
        binCount++;
    }

    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle("Fraction of e decays into central barrel;|#eta_{Q}| #leq x;Ratio");
    auto ratioGraph = new TGraph();
    ratioGraph->SetTitle("#frac{N_{Q decaying to e with |#eta_{e}| #leq 0.9}}{N_{all Q decaying to e}}");
    auto yieldGraph = new TGraph();
    yieldGraph->SetTitle("#frac{N_{Q decaying to e}}{N_{Q total}}");


    double sample = 12;
    double edge = 0.7;

    cout<<"test 1"<<endl;
    for (int i=0; i<=sample; i++){

        double qeDecCount=0;
        double qeEtaCount=0;
        double qeYieldCount=0;
        double qeTotalCount=0;
        double qeRatio;
        double qeYield;

        //cout<<edge<<": ";

        std::ostringstream streamDec;
        streamDec << "etaQ>=-" << edge<<" && etaQ<="<<edge<<" && semilepFlag==1";
        std::string cutDec = streamDec.str();

        std::ostringstream streamEta;
        streamEta << "etaQ>=-" << edge<<" && etaQ<="<<edge<<" && semilepFlag==1 && etaL>=-0.9 && etaL<=0.9";
        std::string cutEta = streamEta.str();

        std::ostringstream streamYield;
        streamYield << "etaQ>=-" << edge<<" && etaQ<="<<edge;
        std::string cutYield = streamYield.str();


        //cout<<cutDec<<endl<<cutEta<<endl;
        binCount = 0;
        cout<<"test 2"<<endl;

        for(std::vector<double>::iterator it = binLuminocity->begin(); it != binLuminocity->end(); ++it){
            double qeDecCountBin = qeTuples[binCount]->GetEntries(cutDec.c_str());
            double qeEtaCountBin = qeTuples[binCount]->GetEntries(cutEta.c_str());
            //double qeYieldCountBin = qeTuples[binCount]->GetEntries(cutYield.c_str());
            //double qeTotalCountBin = qeTuples[binCount]->GetEntries();

            qeDecCount+=qeDecCountBin;
            qeEtaCount+=qeEtaCountBin;
            //qeYieldCount+=qeYieldCountBin;
            //qeTotalCount+=qeTotalCountBin;
            //cout<<binCount<<"\t"<<qeDecCountBin<<"\t"<<qeEtaCountBin<<endl;
            binCount++;
            //cout<<"test 3"<<endl;
        }

        cout<<"test 4"<<endl;
        qeRatio = qeEtaCount/qeDecCount;
        //qeYield = qeYieldCount/qeTotalCount;

        ratioGraph->AddPoint(edge, qeRatio);
        //yieldGraph->AddPoint(edge, qeYield);

        edge+=0.1;
        cout<<qeRatio<<endl;
        cout<<"test 5"<<endl;
    }
    
    cout<<"test 6"<<endl;
    ////Plotting
    TFile *outf =  new TFile("sim_plot_eta_range_opt_e_136.root", "RECREATE");
    ratioGraph->Write();
    cout<<"test 6.5"<<endl;

    TCanvas *canvas = new TCanvas("canvas","canvas");
    cout<<"test 7"<<endl;

    ratioGraph->SetMarkerStyle(20);
        cout<<"test 8"<<endl;

    ratioGraph->SetMarkerColor(kRed);
    mg->Add(ratioGraph);

    cout<<"test 9"<<endl;

    /* yieldGraph->SetMarkerStyle(21);
    yieldGraph->SetMarkerColor(kBlue);
    mg->Add(yieldGraph); */

    mg->Draw("AP");
    //ratioGraph->Draw("AP");
    canvas->BuildLegend();

    canvas->Write();
    
    delete outf;

}