#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

//script to calculate and plot the ratios used to choose the FWD parent quark eta boundaries in the approximation method


void read_sim_eta_range_opt_mu() {
    TFile *f = TFile::Open("sim_tuples_136_semilep_quarks_individual.root","READ");
    
    
    vector<double> *binLuminocity;
    vector<double> *binEvCounts;
    vector<double> *weightSums;
    vector<double> *sigmaGens;

    vector<TNtuple*> qmTuples(5); //should be count of bins but 7 is easier

    f->GetObject("luminocities",binLuminocity);
    f->GetObject("eventCounts",binEvCounts);
    f->GetObject("weightSums",weightSums);
    f->GetObject("sigmaGens",sigmaGens);

    int binCount = 0;

    for(std::vector<double>::iterator it = binLuminocity->begin(); it != binLuminocity->end(); ++it){
        qmTuples[binCount] = (TNtuple*)f->Get(Form("qm%d", binCount));
        binCount++;
    }

    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle("Fraction of #mu decays into foward region;#eta_{Q} #leq x;Ratio");
    auto ratioGraph = new TGraph();
    ratioGraph->SetTitle("#frac{N_{Q decaying to #mu with -4.0 #leq #eta_{#mu} #leq -2.5}}{N_{all Q decaying to #mu}}");
    auto yieldGraph = new TGraph();
    yieldGraph->SetTitle("#frac{N_{Q decaying to #mu}}{N_{Q total}}");


    double sample = 12;
    double edge1 = -2.3;
    double edge2 = -3.8;


    for (int i=0; i<=sample; i++){

        double qmDecCount=0;
        double qmEtaCount=0;
        double qmYieldCount=0;
        double qmTotalCount=0;
        double qmRatio;
        double qmYield;

        cout<<edge1<<": ";

        std::ostringstream streamDec;
        streamDec << "etaQ>=" << edge2<<" && etaQ<="<<edge1<<" && semilepFlag==1";
        std::string cutDec = streamDec.str();

        std::ostringstream streamEta;
        streamEta << "etaQ>=" << edge2<<" && etaQ<="<<edge1<<" && semilepFlag==1 && etaL>=-4 && etaL<=-2.5";
        std::string cutEta = streamEta.str();

        std::ostringstream streamYield;
        streamYield << "etaQ>=" << edge2<<" && etaQ<="<<edge1;
        std::string cutYield = streamYield.str();


        cout<<cutDec<<endl<<cutEta<<endl;
        binCount = 0;

        for(std::vector<double>::iterator it = binLuminocity->begin(); it != binLuminocity->end(); ++it){
            double qmDecCountBin = qmTuples[binCount]->GetEntries(cutDec.c_str());
            double qmEtaCountBin = qmTuples[binCount]->GetEntries(cutEta.c_str());
            double qmYieldCountBin = qmTuples[binCount]->GetEntries(cutYield.c_str());
            double qmTotalCountBin = qmTuples[binCount]->GetEntries();

            qmDecCount+=qmDecCountBin;
            qmEtaCount+=qmEtaCountBin;
            qmYieldCount+=qmYieldCountBin;
            qmTotalCount+=qmTotalCountBin;
            cout<<binCount<<"\t"<<qmDecCountBin<<"\t"<<qmEtaCountBin<<endl;
            binCount++;
        }

        
        qmRatio = qmEtaCount/qmDecCount;
        qmYield = qmYieldCount/qmTotalCount;

        ratioGraph->AddPoint(edge1, qmRatio);
        yieldGraph->AddPoint(edge1, qmYield);

        edge1+=0.1;
        edge2 = edge2 - 0.1;
        //cout<<qmRatio<<endl;
    }
    

    ////Plotting
    TFile *outf =  new TFile("sim_plot_eta_range_opt_mu_136.root", "RECREATE");
    TCanvas *canvas = new TCanvas("canvas","canvas");

    ratioGraph->SetMarkerStyle(20);
    ratioGraph->SetMarkerColor(kRed);
    mg->Add(ratioGraph);

    /* yieldGraph->SetMarkerStyle(21);
    yieldGraph->SetMarkerColor(kBlue);
    mg->Add(yieldGraph); */

    mg->Draw("AP");
    canvas->BuildLegend();

    canvas->Write();
    
    delete outf;

}