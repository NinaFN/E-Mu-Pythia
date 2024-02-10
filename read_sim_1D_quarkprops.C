#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

//script to plot the kinematic distributions of the parent quarks in 1D


void read_sim_1D_quarkprops() {
    TFile *f = TFile::Open("Data/sim_tuples_136_semilep_quarks_individual.root","READ");
    
    
    vector<double> *binLuminocity;
    vector<double> *binEvCounts;
    vector<double> *weightSums;
    vector<double> *sigmaGens;

    vector<TNtuple*> qTuples(8); //should be count of bins but 7 is easier
 

    f->GetObject("luminocities",binLuminocity);
    f->GetObject("eventCounts",binEvCounts);
    f->GetObject("weightSums",weightSums);
    f->GetObject("sigmaGens",sigmaGens);

    int binCount = 0;

    //Create Hists
    TH1F *qPt = new TH1F("q_full","Quark Transverse Momentum;#hat{p}_{T} (GeV/c);#frac{d#sigma}{d#hat{p}_{T}} (pb/GeV/c)", 50, 14.0, 60);
    TH1F *qPtPart = new TH1F("q_pt_part","", 50, 14.0, 60);
   
    TH1F *qEta = new TH1F("q_eta","Quark Pseudorapidity;#eta_{Q};#frac{d#sigma}{d#eta_{Q}}", 50, -6.0, 6.0);
    TH1F *qEtaCut = new TH1F("q_eta_cut","Quark Pseudorapidity (#hat{p}_{T}>10 GeV);#eta_{Q};#frac{d#sigma}{d#eta_{Q}}", 50, -6.0, 6.0);
    TH1F *qEtaPart = new TH1F("q_eta_part","", 50, -6.0, 6.0);
    
    TH2F *EtaPtPart = new TH2F("eta_pt_part","", 100, -6.0, 6.0, 100, 5.0, 10);   
    TH2F *EtaPt = new TH2F("eta_pt","#eta_{Q} vs #hat{p}_{T};#eta_{Q};#hat{p}_{T}", 100, -6.0, 6.0, 100, 5.0, 10);


    TH1F *qPhi = new TH1F("q_phi","Quark Azimuthal Angle;#phi_{Q};#frac{d#sigma}{d#phi_{Q}}", 50, -M_PI, M_PI);
    TH1F *qPhiPart = new TH1F("q_phi_part","", 50, -M_PI, M_PI);

    TH1F *qTheta = new TH1F("q_theta","Quark Polar Angle;#theta_{Q};#frac{d#sigma}{d#theta_{Q}}", 50, 0.0, M_PI);
    TH1F *qThetaPart = new TH1F("q_theta_part","", 50, 0.0, M_PI);

    for(std::vector<double>::iterator it = binLuminocity->begin(); it != binLuminocity->end(); ++it){
        
        //Fill Tuples
        qTuples[binCount] = (TNtuple*)f->Get(Form("qe%d", binCount));

        ////Fill Histograms

        //Electron Decays
      
        qPtPart->Reset();
        qTuples[binCount]->Draw("ptHat>>q_pt_part","ptHat>=14");
        qPtPart->Scale(1/(*it),"width");
        qPt->Add(qPtPart);

        qEtaPart->Reset();
        qTuples[binCount]->Draw("etaQ>>q_eta_part","ptHat>=14");
        qEtaPart->Scale(1/(*it),"width");
        qEta->Add(qEtaPart);

        qEtaPart->Reset();
        qTuples[binCount]->Draw("etaQ>>q_eta_part","ptHat>10");
        qEtaPart->Scale(1/(*it),"width");
        qEtaCut->Add(qEtaPart);

        qPhiPart->Reset();
        qTuples[binCount]->Draw("phiQ>>q_phi_part","ptHat>=14");
        qPhiPart->Scale(1/(*it),"width");
        qPhi->Add(qPhiPart);

        qThetaPart->Reset();
        qTuples[binCount]->Draw("thetaQ>>q_theta_part","ptHat>=14");
        qThetaPart->Scale(1/(*it),"width");
        qTheta->Add(qThetaPart);

        EtaPtPart->Reset();
        qTuples[binCount]->Draw("ptHat:etaQ>>eta_pt_part","ptHat>=14");
        EtaPtPart->Scale(1/(*it),"width");
        EtaPt->Add(EtaPtPart);
        

        binCount++;
    }

    ////Plotting
    // qs
    TFile *outf =  new TFile("Hists/1D_quarkprops_xcut.root", "RECREATE");
    TCanvas *canvasq = new TCanvas("q_sigma","q_sigma");
    canvasq->Divide(1,3);
    
    canvasq->cd(1);
    qPt->SetLineColor(kBlack);
    qPt->SetStats(0);
    qPt->Draw();

    canvasq->cd(2);
    qPhi->SetLineColor(kBlack);
    qPhi->SetStats(0);
    qPhi->SetAxisRange(1e6,10e6,"Y");
    qPhi->Draw();
    
    canvasq->cd(3);
    qEta->SetLineColor(kBlack);
    qEta->SetStats(0);
    qEta->Draw();


    /* canvasq->cd(4);
    qTheta->SetLineColor(kBlack);
    qTheta->SetStats(0);
    qTheta->Draw(); */


    canvasq->Write();

    ////

    TCanvas *canvasqeta = new TCanvas("qeta_sigma","qeta_sigma");
    canvasqeta->Divide(2,1);
    
    
    canvasqeta->cd(1);
    //qEtaCut->Fit("gaus");
    qEtaCut->SetLineColor(kBlack);
    qEtaCut->SetStats(0);
    qEtaCut->Draw();

    canvasqeta->cd(2);
    //gPad->SetLogz();
    EtaPt->SetStats(0);
    //qq2DPhi->GetZaxis()->SetRangeUser(0, 1e7);
    EtaPt->Draw("colz");


    canvasqeta->Write();

    
    delete outf;

}