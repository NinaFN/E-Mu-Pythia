#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

//script to plot the 2D distributions of e vs mu for the candidate pair samples


void read_sim_2D_emu_nocheck() {
    TFile *f = TFile::Open("Data/sim_tuples_pairs_nocheck_hardQCDall.root","READ");
    
    
    vector<double> *binLuminocity;
    vector<double> *binEvCounts;
    vector<double> *weightSums;
    vector<double> *sigmaGens;

    vector<TNtuple*> emuTuples(7); //should be count of bins but 7 is easier

    f->GetObject("luminocities",binLuminocity);
    f->GetObject("eventCounts",binEvCounts);
    f->GetObject("weightSums",weightSums);
    f->GetObject("sigmaGens",sigmaGens);

    int binCount = 0;

    //Create Hists

    TH2F *emu2DPartEta = new TH2F("emu_part_eta","", 50, -8.0, 8.0, 50, -8.0, 8.0);
    TH2F *emu2DPartEtaSmall = new TH2F("emu_part_eta_small","", 20, -0.9, 0.9, 20, -4, -2.5);
    TH2F *emu2DEta = new TH2F("emu_eta","#eta_{e} vs #eta_{#mu};#eta_{e};#eta_{#mu}", 50, -8.0, 8.0, 50, -8.0, 8.0);
    TH2F *emu2DEtaCut = new TH2F("emu_eta_cut","#eta_{e} vs #eta_{#mu};#eta_{e};#eta_{#mu}", 50, -8.0, 8.0, 50, -8.0, 8.0);
    TH2F *emu2DEtaCut1 = new TH2F("emu_eta_cut1","#eta_{e} vs #eta_{#mu};#eta_{e};#eta_{#mu}", 20, -0.9, 0.9, 20, -4, -2.5);
    TH2F *emu2DEtaCut2 = new TH2F("emu_eta_cut2","#eta_{e} vs #eta_{#mu};#eta_{e};#eta_{#mu}", 20, -0.9, 0.9, 20, -4, -2.5);


    /* TH2F *emu2DPartEta = new TH2F("emu_part_eta","", 50, -8.0, 8.0, 50, -8.0, 8.0);
    TH2F *emu2DEta = new TH2F("emu_eta","#eta_{e} vs #eta_{#mu};#eta_{e};#eta_{#mu}", 50, -8.0, 8.0, 50, -8.0, 8.0);
    TH2F *emu2DEtaCut = new TH2F("emu_eta_cut","#eta_{e} vs #eta_{#mu};#eta_{e};#eta_{#mu}", 50, -8.0, 8.0, 50, -8.0, 8.0);
    TH2F *emu2DEtaCut1 = new TH2F("emu_eta_cut1","#eta_{e} vs #eta_{#mu};#eta_{e};#eta_{#mu}", 50, -8.0, 8.0, 50, -8.0, 8.0);
    TH2F *emu2DEtaCut2 = new TH2F("emu_eta_cut2","#eta_{e} vs #eta_{#mu};#eta_{e};#eta_{#mu}", 50, -8.0, 8.0, 50, -8.0, 8.0);
 */
    TH2F *emu2DPartPt = new TH2F("emu_part_pt","", 40, 0, 40, 40, 0, 40.0);
    TH2F *emu2DPt = new TH2F("emu_pt","p_{T}^{e} vs p_{T}^{#mu};p_{T}^{e};p_{T}^{#mu}", 40, 0, 40, 40, 0, 40.0);
    TH2F *emu2DPtCut = new TH2F("emu_pt_cut","p_{T}^{e} vs p_{T}^{#mu};p_{T}^{e};p_{T}^{#mu}", 40, 0, 40, 40, 0, 40.0);
    TH2F *emu2DPtCut1 = new TH2F("emu_pt_cut1","p_{T}^{e} vs p_{T}^{#mu};p_{T}^{e};p_{T}^{#mu}", 40, 0, 40, 40, 0, 40.0);
    TH2F *emu2DPtCut2 = new TH2F("emu_pt_cut2","p_{T}^{e} vs p_{T}^{#mu};p_{T}^{e};p_{T}^{#mu}", 40, 0, 40, 40, 0, 40.0);


    TH2F *emu2DPartPhi = new TH2F("emu_part_phi","", 50, -M_PI, M_PI, 50, -M_PI, M_PI);   
    TH2F *emu2DPhi = new TH2F("emu_phi","#phi_{e} vs #phi_{#mu};#phi_{e};#phi_{#mu}", 50, -M_PI, M_PI, 50, -M_PI, M_PI);
    TH2F *emu2DPhiCut = new TH2F("emu_phi_cut","#phi_{e} vs #phi_{#mu};#phi_{e};#phi_{#mu}", 50, -M_PI, M_PI, 50, -M_PI, M_PI);
    TH2F *emu2DPhiCut1 = new TH2F("emu_phi_cut1","#phi_{e} vs #phi_{#mu};#phi_{e};#phi_{#mu}", 50, -M_PI, M_PI, 50, -M_PI, M_PI);
    TH2F *emu2DPhiCut2 = new TH2F("emu_phi_cut2","#phi_{e} vs #phi_{#mu};#phi_{e};#phi_{#mu}", 50, -M_PI, M_PI, 50, -M_PI, M_PI);

    TH2F *emu2DPartTheta = new TH2F("emu_part_theta","", 50, -0.5, M_PI, 50, -0.5, M_PI);   
    TH2F *emu2DTheta = new TH2F("emu_theta","#theta_{e} vs #theta_{#mu};#theta_{e};#theta_{#mu}", 50, -0.5, M_PI, 50, -0.5, M_PI);
    TH2F *emu2DThetaCut = new TH2F("emu_theta_cut","#theta_{e} vs #theta_{#mu};#theta_{e};#theta_{#mu}", 50, -0.5, M_PI, 50, -0.5, M_PI);

    double pairFoundTot = 0;
    double pairPtTot = 0;
    double pairTot = 0;


    for(std::vector<double>::iterator it = binLuminocity->begin(); it != binLuminocity->end(); ++it){
        
        //Fill Tuples
        emuTuples[binCount] = (TNtuple*)f->Get(Form("emu%d", binCount));

        ////Fill Histograms

        //Electron Decays

        emu2DPartEta->Reset();
        emuTuples[binCount]->Draw("etaMu:etaE>>emu_part_eta");
        emu2DPartEta->Scale(1/(*it),"width");
        emu2DEta->Add(emu2DPartEta);

        emu2DPartPt->Reset();
        emuTuples[binCount]->Draw("ptMu:ptE>>emu_part_pt");
        emu2DPartPt->Scale(1/(*it),"width");
        emu2DPt->Add(emu2DPartPt);

        emu2DPartPhi->Reset();
        emuTuples[binCount]->Draw("phiMu:phiE>>emu_part_phi");
        emu2DPartPhi->Scale(1/(*it),"width");
        emu2DPhi->Add(emu2DPartPhi);

/////

        emu2DPartEta->Reset();
        emuTuples[binCount]->Draw("etaMu:etaE>>emu_part_eta","ptE>=3 && ptMu>=5");
        emu2DPartEta->Scale(1/(*it),"width");
        emu2DEtaCut->Add(emu2DPartEta);

        emu2DPartPt->Reset();
        emuTuples[binCount]->Draw("ptMu:ptE>>emu_part_pt","ptE>=3 && ptMu>=5");
        emu2DPartPt->Scale(1/(*it),"width");
        emu2DPtCut->Add(emu2DPartPt);

        emu2DPartPhi->Reset();
        emuTuples[binCount]->Draw("phiMu:phiE>>emu_part_phi","ptE>=3 && ptMu>=5");
        emu2DPartPhi->Scale(1/(*it),"width");
        emu2DPhiCut->Add(emu2DPartPhi);

/////

        emu2DPartEtaSmall->Reset();
        emuTuples[binCount]->Draw("etaMu:etaE>>emu_part_eta_small","etaE>=-0.9 && etaE<=0.9 && etaMu>=-4 && etaMu<=-2.5");
        emu2DPartEta->Scale(1/(*it),"width");
        emu2DEtaCut1->Add(emu2DPartEtaSmall);

        emu2DPartPt->Reset();
        emuTuples[binCount]->Draw("ptMu:ptE>>emu_part_pt","etaE>=-0.9 && etaE<=0.9 && etaMu>=-4 && etaMu<=-2.5");
        emu2DPartPt->Scale(1/(*it),"width");
        emu2DPtCut1->Add(emu2DPartPt);

        emu2DPartPhi->Reset();
        emuTuples[binCount]->Draw("phiMu:phiE>>emu_part_phi","etaE>=-0.9 && etaE<=0.9 && etaMu>=-4 && etaMu<=-2.5");
        emu2DPartPhi->Scale(1/(*it),"width");
        emu2DPhiCut1->Add(emu2DPartPhi);

/////

        emu2DPartEtaSmall->Reset();
        emuTuples[binCount]->Draw("etaMu:etaE>>emu_part_eta_small","ptE>=3 && ptMu>=5 && etaE>=-0.9 && etaE<=0.9 && etaMu>=-4 && etaMu<=-2.5");
        emu2DPartEta->Scale(1/(*it),"width");
        emu2DEtaCut2->Add(emu2DPartEtaSmall);

        emu2DPartPt->Reset();
        emuTuples[binCount]->Draw("ptMu:ptE>>emu_part_pt","ptE>=3 && ptMu>=5 && etaE>=-0.9 && etaE<=0.9 && etaMu>=-4 && etaMu<=-2.5");
        emu2DPartPt->Scale(1/(*it),"width");
        emu2DPtCut2->Add(emu2DPartPt);

        emu2DPartPhi->Reset();
        emuTuples[binCount]->Draw("phiMu:phiE>>emu_part_phi","ptE>=3 && ptMu>=5 && etaE>=-0.9 && etaE<=0.9 && etaMu>=-4 && etaMu<=-2.5");
        emu2DPartPhi->Scale(1/(*it),"width");
        emu2DPhiCut2->Add(emu2DPartPhi);






        binCount++;
    }


    ////Plotting
    // emus
    TFile *outf =  new TFile("Hists/2D_emu_nocheck_hardQCDall_e3m5.root", "RECREATE");

    TCanvas *canvasEMUAll = new TCanvas("emu_all","emu_all");
    canvasEMUAll->Divide(3,1);


    canvasEMUAll->cd(1);
    //gPad->SetLogz();
    emu2DEta->SetStats(0);
    //emu2DEta->GetZaxis()->SetRangeUser(0, 1e7);
    emu2DEta->Draw("colz");

/*     TBox *Corr1 = new TBox(-0.9,-4.5,0.9,-2.5);
    Corr1->SetFillColor(kBlack);
    Corr1->SetFillStyle(3444);
    Corr1->Draw(); */

    canvasEMUAll->cd(2);
    gPad->SetLogz();
    emu2DPt->SetStats(0);
    //emu2DPt->GetZaxis()->SetRangeUser(0, 1e7);
    emu2DPt->Draw("colz");

    canvasEMUAll->cd(3);
    //gPad->SetLogz();
    emu2DPhi->SetStats(0);
   // emu2DPhi->GetZaxis()->SetRangeUser(0, 1e7);
    emu2DPhi->Draw("colz");

/*     canvasEMUAll->cd(4);
    gPad->SetLogz();
    emu2DTheta->SetStats(0);
    emu2DTheta->GetZaxis()->SetRangeUser(0, 1e7);
    emu2DTheta->Draw("colz"); */

    canvasEMUAll->Write();


    ////


     TCanvas *canvasEMUCut = new TCanvas("emu_cut_pt","emu_cut_pt");
    canvasEMUCut->Divide(3,1);


    canvasEMUCut->cd(1);
    //gPad->SetLogz();
    emu2DEtaCut->SetStats(0);
    //emu2DEta->GetZaxis()->SetRangeUser(0, 1e7);
    emu2DEtaCut->Draw("colz");

    canvasEMUCut->cd(2);
    gPad->SetLogz();
    emu2DPtCut->SetStats(0);
    //emu2DPt->GetZaxis()->SetRangeUser(0, 1e7);
    emu2DPtCut->Draw("colz");

    canvasEMUCut->cd(3);
    //gPad->SetLogz();
    emu2DPhiCut->SetStats(0);
   // emu2DPhi->GetZaxis()->SetRangeUser(0, 1e7);
    emu2DPhiCut->Draw("colz");

    canvasEMUCut->Write();


        ////


     ////


    TCanvas *canvasEMUCut1 = new TCanvas("emu_cut_eta","emu_cut_eta");
    canvasEMUCut1->Divide(3,1);


    canvasEMUCut1->cd(1);
    //gPad->SetLogz();
    emu2DEtaCut1->SetStats(0);
    //emu2DEta->GetZaxis()->SetRangeUser(0, 1e7);
    emu2DEtaCut1->Draw("colz");

    canvasEMUCut1->cd(2);
    gPad->SetLogz();
    emu2DPtCut1->SetStats(0);
    //emu2DPt->GetZaxis()->SetRangeUser(0, 1e7);
    emu2DPtCut1->Draw("colz");

    canvasEMUCut1->cd(3);
    //gPad->SetLogz();
    emu2DPhiCut1->SetStats(0);
   // emu2DPhi->GetZaxis()->SetRangeUser(0, 1e7);
    emu2DPhiCut1->Draw("colz");

    canvasEMUCut1->Write();


        ////



    TCanvas *canvasEMUCut2 = new TCanvas("emu_cut_all","emu_cut_all");
    canvasEMUCut2->Divide(3,1);


    canvasEMUCut2->cd(1);
    //gPad->SetLogz();
    emu2DEtaCut2->SetStats(0);
    //emu2DEta->GetZaxis()->SetRangeUser(0, 1e7);
    emu2DEtaCut2->Draw("colz");

    canvasEMUCut2->cd(2);
    gPad->SetLogz();
    emu2DPtCut2->SetStats(0);
    //emu2DPt->GetZaxis()->SetRangeUser(0, 1e7);
    emu2DPtCut2->Draw("colz");

    canvasEMUCut2->cd(3);
    //gPad->SetLogz();
    emu2DPhiCut2->SetStats(0);
   // emu2DPhi->GetZaxis()->SetRangeUser(0, 1e7);
    emu2DPhiCut2->Draw("colz");

    canvasEMUCut2->Write();


//


    delete outf;

}