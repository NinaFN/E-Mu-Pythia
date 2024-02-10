#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

//script to plot the 2D distributions of e vs mu for the true pair samples


void read_sim_2D_emu() {
    TFile *f = TFile::Open("Data/sim_tuples_136_semilep_paircheck_signs.root","READ");
    
    
    vector<double> *binLuminocity;
    vector<double> *binEvCounts;
    vector<double> *weightSums;
    vector<double> *sigmaGens;

    vector<TNtuple*> qeTuples(7); //should be count of bins but 7 is easier
    vector<TNtuple*> qmTuples(7);
    vector<TNtuple*> qmsTuples(7);

    f->GetObject("luminocities",binLuminocity);
    f->GetObject("eventCounts",binEvCounts);
    f->GetObject("weightSums",weightSums);
    f->GetObject("sigmaGens",sigmaGens);

    int binCount = 0;

    //Create Hists

    TH2F *qe2DPartEta = new TH2F("qe_part_eta","", 50, -6.0, 6.0, 50, -6.0, 6.0);
    TH2F *qe2DPartEtaSmall = new TH2F("qe_part_eta_small","", 20, -0.9, 0.9, 20, -4, -2.5);
    TH2F *qe2DEta = new TH2F("qe_eta","#eta_{e} vs #eta_{#mu};#eta_{e};#eta_{#mu}", 50, -6.0, 6.0, 50, -6.0, 6.0);
    TH2F *qe2DEtaCut = new TH2F("qe_eta_cut","#eta_{e} vs #eta_{#mu};#eta_{e};#eta_{#mu}", 50, -6.0, 6.0, 50, -6.0, 6.0);
    TH2F *qe2DEtaCut1 = new TH2F("qe_eta_cut1","#eta_{e} vs #eta_{#mu};#eta_{e};#eta_{#mu}", 20, -0.9, 0.9, 20, -4, -2.5);
    TH2F *qe2DEtaCut2 = new TH2F("qe_eta_cut2","#eta_{e} vs #eta_{#mu};#eta_{e};#eta_{#mu}", 20, -0.9, 0.9, 20, -4, -2.5);

    TH2F *qe2DPartPt = new TH2F("qe_part_pt","", 40, 0, 40, 40, 0, 40.0);
    TH2F *qe2DPt = new TH2F("qe_pt","p_{T}^{e} vs p_{T}^{#mu};p_{T}^{e};p_{T}^{#mu}", 40, 0, 40, 40, 0, 40.0);
    TH2F *qe2DPtCut = new TH2F("qe_pt_cut","p_{T}^{e} vs p_{T}^{#mu};p_{T}^{e};p_{T}^{#mu}", 40, 0, 40, 40, 0, 40.0);
    TH2F *qe2DPtCut1 = new TH2F("qe_pt_cut1","p_{T}^{e} vs p_{T}^{#mu};p_{T}^{e};p_{T}^{#mu}", 40, 0, 40, 40, 0, 40.0);
    TH2F *qe2DPtCut2 = new TH2F("qe_pt_cut2","p_{T}^{e} vs p_{T}^{#mu};p_{T}^{e};p_{T}^{#mu}", 40, 0, 40, 40, 0, 40.0);


    TH2F *qe2DPartPhi = new TH2F("qe_part_phi","", 50, -M_PI, M_PI, 50, -M_PI, M_PI);   
    TH2F *qe2DPhi = new TH2F("qe_phi","#phi_{e} vs #phi_{#mu};#phi_{e};#phi_{#mu}", 50, -M_PI, M_PI, 50, -M_PI, M_PI);
    TH2F *qe2DPhiCut = new TH2F("qe_phi_cut","#phi_{e} vs #phi_{#mu};#phi_{e};#phi_{#mu}", 50, -M_PI, M_PI, 50, -M_PI, M_PI);
    TH2F *qe2DPhiCut1 = new TH2F("qe_phi_cut1","#phi_{e} vs #phi_{#mu};#phi_{e};#phi_{#mu}", 50, -M_PI, M_PI, 50, -M_PI, M_PI);
    TH2F *qe2DPhiCut2 = new TH2F("qe_phi_cut2","#phi_{e} vs #phi_{#mu};#phi_{e};#phi_{#mu}", 50, -M_PI, M_PI, 50, -M_PI, M_PI);

    TH2F *qe2DPartTheta = new TH2F("qe_part_theta","", 50, -0.5, M_PI, 50, -0.5, M_PI);   
    TH2F *qe2DTheta = new TH2F("qe_theta","#theta_{e} vs #theta_{#mu};#theta_{e};#theta_{#mu}", 50, -0.5, M_PI, 50, -0.5, M_PI);
    TH2F *qe2DThetaCut = new TH2F("qe_theta_cut","#theta_{e} vs #theta_{#mu};#theta_{e};#theta_{#mu}", 50, -0.5, M_PI, 50, -0.5, M_PI);

    double pairFoundTot = 0;
    double pairPtTot = 0;
    double pairTot = 0;

    TH2F *EtaPhiPart = new TH2F("eta_phi_part","", 50, -M_PI, M_PI, 50, -6.0, 6.0);
    TH2F *qeEtaPhi = new TH2F("qe_eta_phi","#eta_{e} vs #phi_{e};#phi_{e};#eta_{e}",  50, -M_PI, M_PI, 50, -6.0, 6.0);
    TH2F *qmEtaPhi = new TH2F("qm_eta_phi","#eta_{#mu} vs #phi_{#mu};#phi_{#mu};#eta_{#mu}",  50, -M_PI, M_PI, 50, -6.0, 6.0);

    TH2F *EtaPtPart = new TH2F("eta_pt_part","", 50, -6.0, 6.0, 50, 0.0, 10);   
    TH2F *qeEtaPt = new TH2F("qe_eta_pt","#eta_{e} vs p_{T}^{e};#eta_{e};p_{T}^{e}", 50, -6.0, 6.0, 50, 0.0, 10);
    TH2F *qmEtaPt = new TH2F("qm_eta_pt","#eta_{#mu} vs  p_{T}^{#mu};#eta_{#mu};p_{T}^{#mu}", 50, -6.0, 6.0, 50, 0.0, 10);

    TH2F *EtaPhiDiffPart = new TH2F("eta_phi_diff_part","", 50, -2*M_PI, 2*M_PI, 50, -6.0, 6.0);
    TH2F *qeEtaPhiDiff = new TH2F("qe_eta_phi_diff","#Delta #eta vs #Delta #phi;#Delta #phi;#Delta #eta",  50, -2*M_PI, 2*M_PI, 50, -6.0, 6.0);
    TH2F *qmEtaPhiDiff = new TH2F("qm_eta_phi_diff","#Delta #eta vs #Delta #phi;#Delta #phi;#Delta #eta",  50, -2*M_PI, 2*M_PI, 50, -6.0, 6.0);


    for(std::vector<double>::iterator it = binLuminocity->begin(); it != binLuminocity->end(); ++it){
        
        //Fill Tuples
        qmTuples[binCount] = (TNtuple*)f->Get(Form("qm%d", binCount));
        qmsTuples[binCount] = (TNtuple*)f->Get(Form("qms%d", binCount));
 
    
        ////Fill Histograms

        //Electron Decays

        qe2DPartEta->Reset();
        qmsTuples[binCount]->Draw("etaL:pairEta>>qe_part_eta","semilepFlag==1 && pairSemilepFlag==1 && ptHat>=14");
        qe2DPartEta->Scale(1/(*it),"width");
        qe2DEta->Add(qe2DPartEta);

        qe2DPartPt->Reset();
        qmsTuples[binCount]->Draw("ptL:pairPt>>qe_part_pt","semilepFlag==1 && pairSemilepFlag==1 && ptHat>=14");
        qe2DPartPt->Scale(1/(*it),"width");
        qe2DPt->Add(qe2DPartPt);

        qe2DPartPhi->Reset();
        qmsTuples[binCount]->Draw("phiL:pairPhi>>qe_part_phi","semilepFlag==1 && pairSemilepFlag==1 && ptHat>=14");
        qe2DPartPhi->Scale(1/(*it),"width");
        qe2DPhi->Add(qe2DPartPhi);

/////

        qe2DPartEta->Reset();
        qmsTuples[binCount]->Draw("etaL:pairEta>>qe_part_eta","semilepFlag==1 && pairSemilepFlag==1 && ptL>=3 && pairPt>=1 && ptHat>=14");
        qe2DPartEta->Scale(1/(*it),"width");
        qe2DEtaCut->Add(qe2DPartEta);

        qe2DPartPt->Reset();
        qmsTuples[binCount]->Draw("ptL:pairPt>>qe_part_pt","semilepFlag==1 && pairSemilepFlag==1 && ptL>=3 && pairPt>=1 && ptHat>=14");
        qe2DPartPt->Scale(1/(*it),"width");
        qe2DPtCut->Add(qe2DPartPt);

        qe2DPartPhi->Reset();
        qmsTuples[binCount]->Draw("phiL:pairPhi>>qe_part_phi","semilepFlag==1 && pairSemilepFlag==1 && ptL>=3 && pairPt>=1 && ptHat>=14");
        qe2DPartPhi->Scale(1/(*it),"width");
        qe2DPhiCut->Add(qe2DPartPhi);

/////

        qe2DPartEtaSmall->Reset();
        qmsTuples[binCount]->Draw("etaL:pairEta>>qe_part_eta_small","semilepFlag==1 && pairSemilepFlag==1 && pairEta>=-0.9 && pairEta<=0.9 && etaL>=-4 && etaL<=-2.5 && ptHat>=14");
        qe2DPartEtaSmall->Scale(1/(*it),"width");
        qe2DEtaCut1->Add(qe2DPartEtaSmall);

        qe2DPartPt->Reset();
        qmsTuples[binCount]->Draw("ptL:pairPt>>qe_part_pt","semilepFlag==1 && pairSemilepFlag==1 && pairEta>=-0.9 && pairEta<=0.9 && etaL>=-4 && etaL<=-2.5 && ptHat>=14");
        qe2DPartPt->Scale(1/(*it),"width");
        qe2DPtCut1->Add(qe2DPartPt);

        qe2DPartPhi->Reset();
        qmsTuples[binCount]->Draw("phiL:pairPhi>>qe_part_phi","semilepFlag==1 && pairSemilepFlag==1 && pairEta>=-0.9 && pairEta<=0.9 && etaL>=-4 && etaL<=-2.5 && ptHat>=14");
        qe2DPartPhi->Scale(1/(*it),"width");
        qe2DPhiCut1->Add(qe2DPartPhi);

/////

        qe2DPartEtaSmall->Reset();
        qmsTuples[binCount]->Draw("etaL:pairEta>>qe_part_eta_small","semilepFlag==1 && pairSemilepFlag==1 && ptL>=3 && pairPt>=1 && pairEta>=-0.9 && pairEta<=0.9 && etaL>=-4 && etaL<=-2.5 && ptHat>=14");
        qe2DPartEtaSmall->Scale(1/(*it),"width");
        qe2DEtaCut2->Add(qe2DPartEtaSmall);

        qe2DPartPt->Reset();
        qmsTuples[binCount]->Draw("ptL:pairPt>>qe_part_pt","semilepFlag==1 && pairSemilepFlag==1 && ptL>=3 && pairPt>=1 && pairEta>=-0.9 && pairEta<=0.9 && etaL>=-4 && etaL<=-2.5 && ptHat>=14");
        qe2DPartPt->Scale(1/(*it),"width");
        qe2DPtCut2->Add(qe2DPartPt);

        qe2DPartPhi->Reset();
        qmsTuples[binCount]->Draw("phiL:pairPhi>>qe_part_phi","semilepFlag==1 && pairSemilepFlag==1 && ptL>=3 && pairPt>=1 && pairEta>=-0.9 && pairEta<=0.9 && etaL>=-4 && etaL<=-2.5 && ptHat>=14");
        qe2DPartPhi->Scale(1/(*it),"width");
        qe2DPhiCut2->Add(qe2DPartPhi);

///////

        EtaPhiPart->Reset();
        qmsTuples[binCount]->Draw("pairEta:pairPhi>>eta_phi_part","semilepFlag==1 && pairSemilepFlag==1 && ptHat>=14");
        EtaPhiPart->Scale(1/(*it),"width");
        qeEtaPhi->Add(EtaPhiPart);

        EtaPhiPart->Reset();
        qmsTuples[binCount]->Draw("etaL:phiL>>eta_phi_part","semilepFlag==1 && pairSemilepFlag==1 && ptHat>=14");
        EtaPhiPart->Scale(1/(*it),"width");
        qmEtaPhi->Add(EtaPhiPart);

        EtaPhiDiffPart->Reset();
        qmsTuples[binCount]->Draw("etaL-pairEta:phiL-pairPhi>>eta_phi_diff_part","semilepFlag==1 && pairSemilepFlag==1 && ptHat>=14");
        EtaPhiDiffPart->Scale(1/(*it),"width");
        qmEtaPhiDiff->Add(EtaPhiDiffPart);

        EtaPtPart->Reset();
        qmsTuples[binCount]->Draw("pairPt:pairEta>>eta_pt_part","semilepFlag==1 && pairSemilepFlag==1 && ptHat>=14");
        EtaPtPart->Scale(1/(*it),"width");
        qeEtaPt->Add(EtaPtPart);

        EtaPtPart->Reset();
        qmsTuples[binCount]->Draw("ptL:etaL>>eta_pt_part","semilepFlag==1 && pairSemilepFlag==1 && ptHat>=14");
        EtaPtPart->Scale(1/(*it),"width");
        qmEtaPt->Add(EtaPtPart);

   



        binCount++;
    }
    cout<<endl;
    cout<<"Total pairs: "<<pairTot<<endl;
    cout<<"Found pairs: "<<pairFoundTot<<endl;


    ////Plotting
    // qes
    TFile *outf =  new TFile("Hists/2D_emu_e1m3_xcut.root", "RECREATE");

    TCanvas *canvasQEAll = new TCanvas("qe_all","qe_all");
    canvasQEAll->Divide(3,1);


    canvasQEAll->cd(1);
    //gPad->SetLogz();
    qe2DEta->SetStats(0);
    //qe2DEta->GetZaxis()->SetRangeUser(0, 1e7);
    qe2DEta->Draw("colz");
    TBox *Corr1 = new TBox(-0.9,-4.5,0.9,-2.5);
    Corr1->SetFillColor(kBlack);
    Corr1->SetFillStyle(3444);
    //Corr1->Draw(); 
    

    canvasQEAll->cd(2);
    gPad->SetLogz();
    qe2DPt->SetStats(0);
    //qe2DPt->GetZaxis()->SetRangeUser(0, 1e7);
    qe2DPt->Draw("colz");

    canvasQEAll->cd(3);
    //gPad->SetLogz();
    qe2DPhi->SetStats(0);
   // qe2DPhi->GetZaxis()->SetRangeUser(0, 1e7);
    qe2DPhi->Draw("colz");

/*     canvasQEAll->cd(4);
    gPad->SetLogz();
    qe2DTheta->SetStats(0);
    qe2DTheta->GetZaxis()->SetRangeUser(0, 1e7);
    qe2DTheta->Draw("colz"); */

    canvasQEAll->Write();


    ////


     TCanvas *canvasQECut = new TCanvas("qe_cut_pt","qe_cut_pt");
    canvasQECut->Divide(3,1);


    canvasQECut->cd(1);
    //gPad->SetLogz();
    qe2DEtaCut->SetStats(0);
    //qe2DEta->GetZaxis()->SetRangeUser(0, 1e7);
    qe2DEtaCut->Draw("colz");

    canvasQECut->cd(2);
    gPad->SetLogz();
    qe2DPtCut->SetStats(0);
    //qe2DPt->GetZaxis()->SetRangeUser(0, 1e7);
    qe2DPtCut->Draw("colz");

    canvasQECut->cd(3);
    //gPad->SetLogz();
    qe2DPhiCut->SetStats(0);
   // qe2DPhi->GetZaxis()->SetRangeUser(0, 1e7);
    qe2DPhiCut->Draw("colz");

    canvasQECut->Write();


        ////


     ////


    TCanvas *canvasQECut1 = new TCanvas("qe_cut_eta","qe_cut_eta");
    canvasQECut1->Divide(3,1);


    canvasQECut1->cd(1);
    //gPad->SetLogz();
    qe2DEtaCut1->SetStats(0);
    //qe2DEta->GetZaxis()->SetRangeUser(0, 1e7);
    qe2DEtaCut1->Draw("colz");

    canvasQECut1->cd(2);
    gPad->SetLogz();
    qe2DPtCut1->SetStats(0);
    //qe2DPt->GetZaxis()->SetRangeUser(0, 1e7);
    qe2DPtCut1->Draw("colz");

    canvasQECut1->cd(3);
    //gPad->SetLogz();
    qe2DPhiCut1->SetStats(0);
   // qe2DPhi->GetZaxis()->SetRangeUser(0, 1e7);
    qe2DPhiCut1->Draw("colz");

    canvasQECut1->Write();


        ////



    TCanvas *canvasQECut2 = new TCanvas("qe_cut_all","qe_cut_all");
    canvasQECut2->Divide(3,1);


    canvasQECut2->cd(1);
    //gPad->SetLogz();
    qe2DEtaCut2->SetStats(0);
    //qe2DEta->GetZaxis()->SetRangeUser(0, 1e7);
    qe2DEtaCut2->Draw("colz");

    canvasQECut2->cd(2);
    gPad->SetLogz();
    qe2DPtCut2->SetStats(0);
    //qe2DPt->GetZaxis()->SetRangeUser(0, 1e7);
    qe2DPtCut2->Draw("colz");

    canvasQECut2->cd(3);
    //gPad->SetLogz();
    qe2DPhiCut2->SetStats(0);
   // qe2DPhi->GetZaxis()->SetRangeUser(0, 1e7);
    qe2DPhiCut2->Draw("colz");

    canvasQECut2->Write();


//

     TCanvas *canvasQEAngles = new TCanvas("qe_angles","qe_angles");
    canvasQEAngles->Divide(2,2);


    canvasQEAngles->cd(1);
    qeEtaPhi->SetStats(0);
    qeEtaPhi->Draw("colz");

    
    canvasQEAngles->cd(2);
    qmEtaPhi->SetStats(0);
    qmEtaPhi->Draw("colz");


    canvasQEAngles->cd(3);
    qeEtaPt->SetStats(0);
    qeEtaPt->Draw("colz");

    canvasQEAngles->cd(4);
    qmEtaPt->SetStats(0);
    qmEtaPt->Draw("colz");

    canvasQEAngles->Write();

    TCanvas *canv = new TCanvas("canv","canv");

    qmEtaPhiDiff->SetStats(0);
    qmEtaPhiDiff->Draw("colz");

    canv->Write();


    delete outf;

}