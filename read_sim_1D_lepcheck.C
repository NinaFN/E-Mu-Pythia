#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

//script to check the amount of central barrel electrons/forward region muons coming from quarks within the chosen eta approximation


void read_sim_1D_lepcheck() {
    TFile *f = TFile::Open("sim_tuples_136_semilep_quarks_individual.root","READ");
    
    
    vector<double> *binLuminocity;
    vector<double> *binEvCounts;
    vector<double> *weightSums;
    vector<double> *sigmaGens;

    vector<TNtuple*> qeTuples(7); //should be count of bins but 7 is easier
    vector<TNtuple*> qmTuples(7);

    f->GetObject("luminocities",binLuminocity);
    f->GetObject("eventCounts",binEvCounts);
    f->GetObject("weightSums",weightSums);
    f->GetObject("sigmaGens",sigmaGens);

    int binCount = 0;

    //Create Hists

    //hists of electrons entirely based on kinematics cuts
    

    TH1F *qeAll = new TH1F("qe_all","Semileptonic Q #rightarrow e decays in the central barrel;p_{T}^{e};#frac{d#sigma}{dp_{T}^{#mu}} (pb/GeV/c)", 20, 0.0, 40.0);
    TH1F *qePart = new TH1F("qe_part","", 20, 0.0, 40.0);
    TH1F *qeParentRange = new TH1F("qe_range",";p_{T}^{e};Ratio", 20, 0.0, 40.0);
    TH1F *etaParte = new TH1F("eta_parte","",  50, -6, 6);
    TH1F *qeParentEta = new TH1F("qe_eta","Electron Parent Quark Pseudorapidity;#eta_{Q};#frac{d#sigma}{d#eta_{Q}}", 50, -6, 6);

    TH1F *qmAll = new TH1F("qm_all","Semileptonic Q #rightarrow #mu decays in the forward region;p_{T}^{#mu};#frac{d#sigma}{dp_{T}^{#mu}} (pb/GeV/c)", 20, 0.0, 40.0);
    TH1F *qmPart = new TH1F("qm_part","", 20, 0.0, 40.0);
    TH1F *qmParentRange = new TH1F("qm_range",";p_{T}^{#mu};Ratio", 20, 0.0, 40.0);
    TH1F *etaPartm = new TH1F("eta_partm","",  50, -6, 6);
    TH1F *qmParentEta = new TH1F("qm_eta","Muon Parent Quark Pseudorapidity;#eta_{Q};#frac{d#sigma}{d#eta_{Q}}", 50, -6, 6);

    

    for(std::vector<double>::iterator it = binLuminocity->begin(); it != binLuminocity->end(); ++it){
        
        //Fill Tuples
        qeTuples[binCount] = (TNtuple*)f->Get(Form("qe%d", binCount));
        qmTuples[binCount] = (TNtuple*)f->Get(Form("qm%d", binCount));

        ////Fill Histograms

        //Electron Decays
      
        qePart->Reset();
        qeTuples[binCount]->Draw("ptL>>qe_part","semilepFlag==1 && etaL>=-0.9 && etaL<=0.9");// && etaQ<=1.6 && etaQ>=-1.6
        qePart->Scale(1/(*it),"width");
        qeAll->Add(qePart); 

        qePart->Reset();
        qeTuples[binCount]->Draw("ptL>>qe_part","semilepFlag==1 && etaQ<=1.6 && etaQ>=-1.6 && etaL>=-0.9 && etaL<=0.9");
        qePart->Scale(1/(*it),"width");
        qeParentRange->Add(qePart);

        etaParte->Reset();
        qeTuples[binCount]->Draw("etaQ>>eta_parte","etaL>=-0.9 && etaL<=0.9 && semilepFlag==1");
        etaParte->Scale(1/(*it),"width");
        qeParentEta->Add(etaParte);


        // qm

        qmPart->Reset();
        qmTuples[binCount]->Draw("ptL>>qm_part","semilepFlag==1 && etaL>=-4 && etaL<=-2.5");// && etaQ<=1.6 && etaQ>=-1.6
        qmPart->Scale(1/(*it),"width");
        qmAll->Add(qmPart); 

        qmPart->Reset();
        qmTuples[binCount]->Draw("ptL>>qm_part","semilepFlag==1 && etaQ<=-2.1 && etaQ>=-4.4 && etaL>=-4 && etaL<=-2.5");
        qmPart->Scale(1/(*it),"width");
        qmParentRange->Add(qmPart);

        etaPartm->Reset();
        qmTuples[binCount]->Draw("etaQ>>eta_partm","etaL>=-4 && etaL<=-2.5 && semilepFlag==1");
        etaPartm->Scale(1/(*it),"width");
        qmParentEta->Add(etaPartm);

        binCount++;
    }



    ////Plotting
    // qes
    TFile *outf =  new TFile("A_lepcheck.root", "RECREATE");

    TCanvas *canvasQEAll = new TCanvas("qe_all","qe_all");

    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetLogy();
    pad1->SetBottomMargin(0);
    qeAll->GetYaxis()->SetTitleOffset(1);
    pad1->Draw();
    pad1->cd(); 


    qeAll->SetLineColor(kBlack);
    qeAll->SetStats(0);
    qeAll->Draw("SAME");

    qeParentRange->SetLineColor(kRed);
    qeParentRange->SetStats(0);
    qeParentRange->DrawCopy("SAME");


    auto labelqeEta = new TLatex();
    labelqeEta->DrawLatex(0.0, 0.0, "-0.9 < #eta_{e} < 0.9");
    labelqeEta->Draw("SAME"); 

    auto legendqe = new TLegend();
    legendqe->AddEntry(qeAll,"All Q #rightarrow e decays","l");
    legendqe->AddEntry(qeParentRange,"Q #rightarrow e decays with |#eta_{Q}| #leq 1.6","l");
    legendqe->Draw("SAME");

    canvasQEAll->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);

     qeParentRange->GetXaxis()->SetTitleSize(13);
     qeParentRange->GetXaxis()->SetTitleOffset(1);
     qeParentRange->GetXaxis()->SetLabelSize(10);

     qeParentRange->GetYaxis()->SetTitleSize(13);
     qeParentRange->GetYaxis()->SetTitleOffset(1);
     qeParentRange->GetYaxis()->SetLabelSize(10);

    pad2->Draw("SAME");
    pad2->cd(); 

     qeParentRange->Divide(qeAll);
     qeParentRange->SetStats(0);
     qeParentRange->SetAxisRange(0.75,1.1,"Y");
     qeParentRange->DrawCopy();

   

    canvasQEAll->Write();

    // qm

   TCanvas *canvasQMAll = new TCanvas("qm_all","qm_all");

    TPad *pad1m = new TPad("pad1m", "pad1m", 0, 0.3, 1, 1.0);
    pad1m->SetLogy();
    pad1m->SetBottomMargin(0);
    qmAll->GetYaxis()->SetTitleOffset(1);
    pad1m->Draw();
    pad1m->cd(); 


    qmAll->SetLineColor(kBlack);
    qmAll->SetStats(0);
    qmAll->Draw("SAME");

    qmParentRange->SetLineColor(kRed);
    qmParentRange->SetStats(0);
    qmParentRange->DrawCopy("SAME"); 


    auto labelqmEta = new TLatex();
    labelqmEta->DrawLatex(0.0, 0.0, "-4 < #eta_{#mu} < -2.5");
    labelqmEta->Draw("SAME"); 

    auto legendqm = new TLegend();
    legendqm->AddEntry(qmAll,"All Q #rightarrow #mu decays","l");
    legendqm->AddEntry(qmParentRange,"Q #rightarrow #mu decays with -4.4 #leq #eta_{Q} #leq -2.1","l");
    legendqm->Draw("SAME");

    canvasQMAll->cd();
    TPad *pad2m = new TPad("pad2m", "pad2m", 0, 0.05, 1, 0.3);
    pad2m->SetTopMargin(0);
    pad2m->SetBottomMargin(0.2);

     qmParentRange->GetXaxis()->SetTitleSize(13);
     qmParentRange->GetXaxis()->SetTitleOffset(1);
     qmParentRange->GetXaxis()->SetLabelSize(10);

     qmParentRange->GetYaxis()->SetTitleSize(13);
     qmParentRange->GetYaxis()->SetTitleOffset(1);
     qmParentRange->GetYaxis()->SetLabelSize(10);

    pad2m->Draw("SAME");
    pad2m->cd(); 

     qmParentRange->Divide(qmAll);
     qmParentRange->SetStats(0);
     qmParentRange->SetAxisRange(0.55,1.1,"Y");
     qmParentRange->DrawCopy();

   

    canvasQMAll->Write();



    TCanvas *canvasE = new TCanvas("e_kine","e_kine");
    canvasE->Divide(1,2);
    
    
    canvasE->cd(1);
    qeParentEta->SetLineColor(kBlack);
    qeParentEta->SetStats(0);
    qeParentEta->Draw();
    auto labeleEta = new TLatex();
    labeleEta->DrawLatex(0.0, 0.0, "-0.9 < #eta_{e} < 0.9");
    labeleEta->Draw("SAME");

    canvasE->cd(2);
    qmParentEta->SetLineColor(kBlack);
    qmParentEta->SetStats(0);
    qmParentEta->Draw();
    auto labelmEta = new TLatex();
    labelmEta->DrawLatex(0.0, 0.0, "-4 < #eta_{#mu} < -2.5");
    labelmEta->Draw("SAME"); 


    canvasE->Write();
    
    delete outf;

}