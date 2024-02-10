#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

//script to compare the kinematic variables of parent quarks and child leptons


void read_sim_2D_quarklep() {
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

    TH2F *qe2DPartEta = new TH2F("qe_part_eta","", 50, -6.0, 6.0, 50, -6.0, 6.0);
    TH2F *qe2DEta = new TH2F("qe_eta","#eta_{Q} vs #eta_{e};#eta_{Q};#eta_{e}", 50, -6.0, 6.0, 50, -6.0, 6.0);
    
    TH2F *qe2DPartEtaCut = new TH2F("qe_part_etacut","", 50, -0.9, 0.9, 50, -1.6, 1.6);
    TH2F *qe2DEtaCut = new TH2F("qe_etacut","#eta_{Q} vs #eta_{e};#eta_{e};#eta_{Q}", 50, -0.9, 0.9, 50, -1.6, 1.6);
    

    TH2F *qe2DPartPt = new TH2F("qe_part_pt","", 50, 14.0, 60.0, 50, 0.0, 40.0);
    TH2F *qe2DPt = new TH2F("qe_pt","p_{T}^{Q} vs p_{T}^{e};p_{T}^{Q};p_{T}^{e}", 50, 14.0, 60.0, 50, 0.0, 40.0);
    
    TH2F *qe2DPartPhi = new TH2F("qe_part_phi","", 50, -M_PI, M_PI, 50, -M_PI, M_PI);   
    TH2F *qe2DPhi = new TH2F("qe_phi","#phi_{Q} vs #phi_{e};#phi_{Q};#phi_{e}", 50, -M_PI, M_PI, 50, -M_PI, M_PI);

    TH2F *qe2DPartTheta = new TH2F("qe_part_theta","", 50, 0, M_PI, 50, 0, M_PI);   
    TH2F *qe2DTheta = new TH2F("qe_theta","#theta_{Q} vs #theta_{e};#theta_{Q};#theta_{e}", 50, 0, M_PI, 50, 0, M_PI);

    TH2F *qe2DPartThetaCalc = new TH2F("qe_part_thetacalc","", 50, 0, M_PI, 50, 0, M_PI);   
    TH2F *qe2DThetaCalc = new TH2F("qe_thetacalc","#theta_{Q} vs #theta_{e};#theta_{Q};#theta_{e}", 50, 0, M_PI, 50, 0, M_PI);

    ///

    TH2F *qm2DPartEta = new TH2F("qm_part_eta","", 50, -6.0, 6.0, 50, -6.0, 6.0);
    TH2F *qm2DEta = new TH2F("qm_eta","#eta_{Q} vs #eta_{#mu};#eta_{Q};#eta_{#mu}", 50, -6.0, 6.0, 50, -6.0, 6.0);
    
    TH2F *qm2DPartEtaCut = new TH2F("qm_part_etacut","", 50, -4, -2.5, 50, -4.4, -2.1);
    TH2F *qm2DEtaCut = new TH2F("qm_etacut","#eta_{Q} vs #eta_{#mu};#eta_{#mu};#eta_{Q}", 50, -4, -2.5, 50, -4.4, -2.1);
    

    TH2F *qm2DPartPt = new TH2F("qm_part_pt","", 50, 14.0, 60.0, 50, 0.0, 40.0);
    TH2F *qm2DPt = new TH2F("qm_pt","p_{T}^{Q} vs p_{T}^{#mu};p_{T}^{Q};p_{T}^{#mu}", 50, 14.0, 60.0, 50, 0.0, 40.0);
    
    TH2F *qm2DPartPhi = new TH2F("qm_part_phi","", 50, -M_PI, M_PI, 50, -M_PI, M_PI);
    TH2F *qm2DPhi = new TH2F("qm_phi","#phi_{Q} vs #phi_{#mu};#phi_{Q};#phi_{#mu}", 50, -M_PI, M_PI, 50, -M_PI, M_PI);

    TH2F *qm2DPartTheta = new TH2F("qm_part_theta","", 50, 0, M_PI, 50, 0, M_PI);
    TH2F *qm2DTheta = new TH2F("qm_theta","#theta_{Q} vs #theta_{#mu};#theta_{Q};#theta_{#mu}", 50, 0, M_PI, 50, 0, M_PI);


    TH2F *EtaPhiDiffPart = new TH2F("eta_phi_diff_part","", 50, -2*M_PI, 2*M_PI, 50, -6.0, 6.0);
    TH2F *qeEtaPhiDiff = new TH2F("qe_eta_phi_diff","#Delta #eta vs #Delta #phi;#Delta #phi;#Delta #eta",  50, -2*M_PI, 2*M_PI, 50, -6.0, 6.0);
    TH2F *qmEtaPhiDiff = new TH2F("qm_eta_phi_diff","#Delta #eta vs #Delta #phi;#Delta #phi;#Delta #eta",  50, -2*M_PI, 2*M_PI, 50, -6.0, 6.0);



    for(std::vector<double>::iterator it = binLuminocity->begin(); it != binLuminocity->end(); ++it){
        
        //Fill Tuples
        qeTuples[binCount] = (TNtuple*)f->Get(Form("qe%d", binCount));
        qmTuples[binCount] = (TNtuple*)f->Get(Form("qm%d", binCount));

        ////Fill Histograms

        //Electron Decays

        qe2DPartEta->Reset();
        qeTuples[binCount]->Draw("etaL:etaQ>>qe_part_eta","semilepFlag==1 && ((idQ*idL)<0) && ptHat>=14");
        qe2DPartEta->Scale(1/(*it),"width");
        qe2DEta->Add(qe2DPartEta);

        qe2DPartEtaCut->Reset();
        qeTuples[binCount]->Draw("etaQ:etaL>>qe_part_etacut","semilepFlag==1 && ((idQ*idL)<0) && ptHat>=14");
        qe2DPartEtaCut->Scale(1/(*it),"width");
        qe2DEtaCut->Add(qe2DPartEtaCut);

        qe2DPartPt->Reset();
        qeTuples[binCount]->Draw("ptL:ptHat>>qe_part_pt","semilepFlag==1 && ((idQ*idL)<0) && ptHat>=14");
        qe2DPartPt->Scale(1/(*it),"width");
        qe2DPt->Add(qe2DPartPt);

        qe2DPartPhi->Reset();
        qeTuples[binCount]->Draw("phiL:phiQ>>qe_part_phi","semilepFlag==1 && ((idQ*idL)<0) && ptHat>=14");
        qe2DPartPhi->Scale(1/(*it),"width");
        qe2DPhi->Add(qe2DPartPhi);

        qe2DPartTheta->Reset();
        qeTuples[binCount]->Draw("thetaL:thetaQ>>qe_part_theta","semilepFlag==1 && ((idQ*idL)<0) && ptHat>=14");
        qe2DPartTheta->Scale(1/(*it),"width");
        qe2DTheta->Add(qe2DPartTheta);

        qe2DPartThetaCalc->Reset();
        qeTuples[binCount]->Draw("thetaL+thetaQ:thetaQ>>qe_part_thetacalc","semilepFlag==1 && ((idQ*idL)<0) && ptHat>=14");
        qe2DPartTheta->Scale(1/(*it),"width");
        qe2DThetaCalc->Add(qe2DPartThetaCalc);

        // qm

        qm2DPartEta->Reset();
        qmTuples[binCount]->Draw("etaL:etaQ>>qm_part_eta","semilepFlag==1 && ((idQ*idL)<0) && ptHat>=14");
        qm2DPartEta->Scale(1/(*it),"width");
        qm2DEta->Add(qm2DPartEta);

        qm2DPartEtaCut->Reset();
        qmTuples[binCount]->Draw("etaQ:etaL>>qm_part_etacut","semilepFlag==1 && ((idQ*idL)<0) && ptHat>=14");
        qm2DPartEtaCut->Scale(1/(*it),"width");
        qm2DEtaCut->Add(qm2DPartEtaCut);

        qm2DPartPt->Reset();
        qmTuples[binCount]->Draw("ptL:ptHat>>qm_part_pt","semilepFlag==1 && ((idQ*idL)<0) && ptHat>=14");
        qm2DPartPt->Scale(1/(*it),"width");
        qm2DPt->Add(qm2DPartPt);

        qm2DPartPhi->Reset();
        qmTuples[binCount]->Draw("phiL:phiQ>>qm_part_phi","semilepFlag==1 && ((idQ*idL)<0) && ptHat>=14");
        qm2DPartPhi->Scale(1/(*it),"width");
        qm2DPhi->Add(qm2DPartPhi);

        qm2DPartTheta->Reset();
        qmTuples[binCount]->Draw("thetaL:thetaQ>>qm_part_theta","semilepFlag==1 && ((idQ*idL)<0) && ptHat>=14");
        qm2DPartTheta->Scale(1/(*it),"width");
        qm2DTheta->Add(qm2DPartTheta);


        EtaPhiDiffPart->Reset();
        qmTuples[binCount]->Draw("etaQ-etaL:phiQ-phiL>>eta_phi_diff_part","semilepFlag==1 && ((idQ*idL)<0) && ptHat>=14");
        EtaPhiDiffPart->Scale(1/(*it),"width");
        qmEtaPhiDiff->Add(EtaPhiDiffPart);

        EtaPhiDiffPart->Reset();
        qeTuples[binCount]->Draw("etaQ-etaL:phiQ-phiL>>eta_phi_diff_part","semilepFlag==1 && ((idQ*idL)<0) && ptHat>=14");
        EtaPhiDiffPart->Scale(1/(*it),"width");
        qeEtaPhiDiff->Add(EtaPhiDiffPart);


        binCount++;
    }


    ////Plotting
    // qes
    TFile *outf =  new TFile("Hists/2D_quarklep_xcut.root", "RECREATE");

    TCanvas *canvasQEAll = new TCanvas("qe_all","qe_all");
    canvasQEAll->Divide(3,1);


    canvasQEAll->cd(1);
    //gPad->SetLogz();
    qe2DEta->SetStats(0);
    //qe2DEta->GetZaxis()->SetRangeUser(0, 4e5);
    qe2DEta->Draw("colz");

    canvasQEAll->cd(2);
    gPad->SetLogz();
    qe2DPt->SetStats(0);
    //qe2DPt->GetZaxis()->SetRangeUser(0, 4e5);
    qe2DPt->Draw("colz");

    /* TLine *l = new TLine(14,60,14,60);
    l->Draw(); */

    canvasQEAll->cd(3);
    //gPad->SetLogz();
    qe2DPhi->SetStats(0);
    //qe2DPhi->GetZaxis()->SetRangeUser(0, 4e5);
    qe2DPhi->Draw("colz");

    /* canvasQEAll->cd(4);
    gPad->SetLogz();
    qe2DTheta->SetStats(0);
    //qe2DTheta->GetZaxis()->SetRangeUser(0, 4e5);
    qe2DTheta->Draw("colz"); */

    canvasQEAll->Write();

    //

    TCanvas *canvasTheta = new TCanvas("theta","theta");
    canvasTheta->Divide(1,1);

    canvasQEAll->cd(4);
    //gPad->SetLogz();
    qe2DTheta->SetStats(0);
    //qe2DTheta->GetZaxis()->SetRangeUser(0, 4e5);
    qe2DTheta->Draw("colz");



    // qm

    TCanvas *canvasQMAll = new TCanvas("qm_all","qm_all");
    canvasQMAll->Divide(3,1);

    canvasQMAll->cd(1);
    //gPad->SetLogz();
    qm2DEta->SetStats(0);
    //qm2DEta->GetZaxis()->SetRangeUser(0, 4e5);
    qm2DEta->Draw("colz");

    canvasQMAll->cd(2);
    gPad->SetLogz();
    qm2DPt->SetStats(0);
    //qm2DPt->GetZaxis()->SetRangeUser(0, 4e5);
    qm2DPt->Draw("colz");
    /* TLine *lm = new TLine(14,60,14,60);
    lm->Draw(); */

    canvasQMAll->cd(3);
    //gPad->SetLogz();
    qm2DPhi->SetStats(0);
    //qm2DPhi->GetZaxis()->SetRangeUser(0, 4e5);
    qm2DPhi->Draw("colz");

    /* canvasQMAll->cd(4);
    gPad->SetLogz();
    qm2DTheta->SetStats(0);
    //qm2DTheta->GetZaxis()->SetRangeUser(0, 4e5);
    qm2DTheta->Draw("colz"); */

    canvasQMAll->Write();
    
    

    TCanvas *canvasEtaTest = new TCanvas("eta","eta");
    canvasEtaTest->Divide(2,1);

    canvasEtaTest->cd(1);
    gPad->SetLogz();
    qe2DEtaCut->SetStats(0);
    //qe2DEtaCut->GetZaxis()->SetRangeUser(0, 13e6);
    qe2DEtaCut->Draw("colz");

    canvasEtaTest->cd(2);
    gPad->SetLogz();
    qm2DEtaCut->SetStats(0);
    //qm2DEtaCut->GetZaxis()->SetRangeUser(0, 13e6);
    qm2DEtaCut->Draw("colz");


    canvasEtaTest->Write();


    TCanvas *canv = new TCanvas("canv","canv");
    canv->Divide(2,1);

    canv->cd(1);
    qeEtaPhiDiff->SetStats(0);
    qeEtaPhiDiff->Draw("colz");

    canv->cd(2);
    qmEtaPhiDiff->SetStats(0);
    qmEtaPhiDiff->Draw("colz");

    canv->Write();


    delete outf;

}