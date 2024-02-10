#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

//script to compare the kinematic variables of parent quarks and child leptons


void read_sim_1D_quarklep() {
    TFile *f = TFile::Open("Data/sim_tuples_136_semilep_quarks_individual.root","READ");
    
    
    vector<double> *binLuminocity;
    vector<double> *binEvCounts;
    vector<double> *weightSums;
    vector<double> *sigmaGens;

    vector<TNtuple*> qeTuples(8); //should be count of bins but 7 is easier
    vector<TNtuple*> qmTuples(8);

    f->GetObject("luminocities",binLuminocity);
    f->GetObject("eventCounts",binEvCounts);
    f->GetObject("weightSums",weightSums);
    f->GetObject("sigmaGens",sigmaGens);

    int binCount = 0;

    //Create Hists

    TH1F *qe1DPartEta = new TH1F("qe_part_eta","", 50, -8.0, 8.0);
    TH1F *qe1DEta = new TH1F("qe_eta","Difference between quark and electron #eta;#eta_{Q} - #eta_{e};", 50, -8.0, 8.0);

    TH1F *qe1DPartPt = new TH1F("qe_part_pt","", 50, 5.0, 60.0);
    TH1F *qe1DPt = new TH1F("qe_pt","p_{T}^{Q} vs p_{T}^{e};p_{T}^{Q};p_{T}^{e}", 50, 5.0, 60.0);
    
    TH1F *qe1DPartPhi = new TH1F("qe_part_phi","", 50, -2*M_PI, 2*M_PI);   
    TH1F *qe1DPhi = new TH1F("qe_phi","Difference between quark and electron #phi;#phi_{Q} - #phi_{e};", 50, -2*M_PI, 2*M_PI);


    TH1F *qe1DPartFF = new TH1F("qe_part_ff","", 50, 0, 1);   
    TH1F *qe1DFF = new TH1F("qe_ff","Fragmentation Function: Q #rightarrow e;p_{T}^{e}/p_{T}^{Q};", 50, 0, 1);
    

    ///

    TH1F *qm1DPartEta = new TH1F("qm_part_eta","", 50, -12.0, 12.0);
    TH1F *qm1DEta = new TH1F("qm_eta","Difference between quark and muon #eta;#eta_{Q} - #eta_{#mu};", 50, -12.0, 12.0);

    TH1F *qm1DPartPt = new TH1F("qm_part_pt","", 50, 5.0, 60.0);
    TH1F *qm1DPt = new TH1F("qm_pt","p_{T}^{Q} vs p_{T}^{#mu};p_{T}^{Q};p_{T}^{#mu}", 50, 5.0, 60.0);
    
    TH1F *qm1DPartPhi = new TH1F("qm_part_phi","", 50, -2*M_PI, 2*M_PI);
    TH1F *qm1DPhi = new TH1F("qm_phi","Difference between quark and muon #phi;#phi_{Q} - #phi_{#mu};", 50, -2*M_PI, 2*M_PI);

    TH1F *qm1DPartFF = new TH1F("qm_part_ff","", 50, 0, 1);   
    TH1F *qm1DFF = new TH1F("qm_ff","Fragmentation Function: Q #rightarrow #mu;p_{T}^{#mu}/p_{T}^{Q};", 50, 0, 1);

    for(std::vector<double>::iterator it = binLuminocity->begin(); it != binLuminocity->end(); ++it){
        
        //Fill Tuples
        qeTuples[binCount] = (TNtuple*)f->Get(Form("qe%d", binCount));
        qmTuples[binCount] = (TNtuple*)f->Get(Form("qm%d", binCount));

        ////Fill Histograms

        //Electron Decays

        qe1DPartEta->Reset();
        qeTuples[binCount]->Draw("etaL-etaQ>>qe_part_eta","semilepFlag==1 && ((idQ*idL)<0) && ptHat>=14");
        qe1DPartEta->Scale(1/(*it),"width");
        qe1DEta->Add(qe1DPartEta);

       
        qe1DPartPhi->Reset();
        qeTuples[binCount]->Draw("phiL-phiQ>>qe_part_phi","semilepFlag==1 && ((idQ*idL)<0) && ptHat>=14"); //(((phiL-phiQ)<-2.5 && (phiL-phiQ)>-3.7) || ((phiL-phiQ)>2.5 && (phiL-phiQ)<3.7))
        qe1DPartPhi->Scale(1/(*it),"width");
        qe1DPhi->Add(qe1DPartPhi);

    
        qe1DPartFF->Reset();
        qeTuples[binCount]->Draw("ptL/ptHat>>qe_part_ff","semilepFlag==1 && ((idQ*idL)<0) && ptHat>=14");
        qe1DPartFF->Scale(1/(*it),"width");
        qe1DFF->Add(qe1DPartFF);

        // qm

        qm1DPartEta->Reset();
        qmTuples[binCount]->Draw("etaL-etaQ>>qm_part_eta","semilepFlag==1 && ((idQ*idL)<0) && ptHat>=14");
        qm1DPartEta->Scale(1/(*it),"width");
        qm1DEta->Add(qm1DPartEta);

    
        qm1DPartPhi->Reset();
        qmTuples[binCount]->Draw("phiL-phiQ>>qm_part_phi","semilepFlag==1 && ((idQ*idL)<0) && ptHat>=14");
        qm1DPartPhi->Scale(1/(*it),"width");
        qm1DPhi->Add(qm1DPartPhi);

        qm1DPartFF->Reset();
        qmTuples[binCount]->Draw("ptL/ptHat>>qm_part_ff","semilepFlag==1 && ((idQ*idL)<0) && ptHat>=14");
        qm1DPartFF->Scale(1/(*it),"width");
        qm1DFF->Add(qm1DPartFF);


        binCount++;
    }


    ////Plotting
    // qes
    TFile *outf =  new TFile("Hists/1D_quarklep_xcut.root", "RECREATE");

    TCanvas *canvasQEAll = new TCanvas("qe_all","qe_all");
    canvasQEAll->Divide(3,1);

    canvasQEAll->cd(1);
    qe1DFF->SetStats(0);
    qe1DFF->Draw();

    canvasQEAll->cd(2);
    qe1DPhi->SetStats(0);
    qe1DPhi->Draw();

    canvasQEAll->cd(3);
    qe1DEta->SetStats(0);
    qe1DEta->Draw();
    


    canvasQEAll->Write();

    //

    // qm

    TCanvas *canvasQMAll = new TCanvas("qm_all","qm_all");
    canvasQMAll->Divide(3,1);



    canvasQMAll->cd(1);
    qm1DFF->SetStats(0);
    qm1DFF->Draw();

    canvasQMAll->cd(2);
    qm1DPhi->SetStats(0);
    qm1DPhi->Draw();

    canvasQMAll->cd(3);
    qm1DEta->SetStats(0);
    qm1DEta->Draw();

    canvasQMAll->Write();



    TCanvas *canvasff = new TCanvas("ff","ff");
    canvasff->Divide(2,1);

    canvasff->cd(1);
    qe1DFF->Draw();

    canvasff->cd(2);
    qm1DFF->Draw();
    
    canvasff->Write();

    delete outf;

}