#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

//script to compare the azimuthal difference between e-mu pairs in real and candidate samples at various cuts 


void read_sim_pairs_compare_phi() {
    TFile *f = TFile::Open("sim_tuples_pairs_nocheck_hardQCD_calc.root","READ");
    TFile *f2 = TFile::Open("sim_tuples_136_semilep_paircheck.root","READ");

    
    vector<double> *binLuminocity;
    vector<double> *binEvCounts;
    vector<double> *weightSums;
    vector<double> *sigmaGens;

    vector<double> *binLuminocity2;
    vector<double> *binEvCounts2;
    vector<double> *weightSums2;
    vector<double> *sigmaGens2;

    vector<TNtuple*> emuTuples(7); //should be count of bins but 7 is easier
    vector<TNtuple*> qmsTuples(7);

    f->GetObject("luminocities",binLuminocity);
    f->GetObject("eventCounts",binEvCounts);
    f->GetObject("weightSums",weightSums);
    f->GetObject("sigmaGens",sigmaGens);

    f2->GetObject("luminocities",binLuminocity2);
    f2->GetObject("eventCounts",binEvCounts2);
    f2->GetObject("weightSums",weightSums2);
    f2->GetObject("sigmaGens",sigmaGens2);

    int binCount = 0;
    int binCount2 = 0;

    //Create Hists
    TH1F *emuPart = new TH1F("emu_part","", 19, 0, 2*M_PI);
    TH1F *emuPhiDiff = new TH1F("emu_phi_diff",";#Delta #phi;#frac{d#sigma}{d#Delta#phi}", 19, 0, 2*M_PI);
    TH1F *emuPhiDiffCutLowEx = new TH1F("emu_phi_diff_cut_lowex",";#Delta #phi;#frac{d#sigma}{d#Delta#phi}", 19, 0, 2*M_PI);
    TH1F *emuPhiDiffCutLow = new TH1F("emu_phi_diff_cut_low",";#Delta #phi;#frac{d#sigma}{d#Delta#phi}", 19, 0, 2*M_PI);
    TH1F *emuPhiDiffCutHigh = new TH1F("emu_phi_diff_cut_high",";#Delta #phi;#frac{d#sigma}{d#Delta#phi}", 19, 0, 2*M_PI);

    TH1F *emuPartSL = new TH1F("SLemu_part","", 19, 0, 2*M_PI);
    TH1F *emuPhiDiffSL = new TH1F("SLemu_phi_diff",";#Delta #phi;", 19, 0, 2*M_PI);
    TH1F *emuPhiDiffCutLowExSL = new TH1F("SLemu_phi_diff_cut_lowex",";#Delta #phi;", 19, 0, 2*M_PI);
    TH1F *emuPhiDiffCutLowSL = new TH1F("SLemu_phi_diff_cut_low",";#Delta #phi;", 19, 0, 2*M_PI);
    TH1F *emuPhiDiffCutHighSL = new TH1F("SLemu_phi_diff_cut_high",";#Delta #phi;", 19, 0, 2*M_PI);
   

    for(std::vector<double>::iterator it = binLuminocity->begin(); it != binLuminocity->end(); ++it){
        cout<<*it<<endl;
        //Fill Tuples
        emuTuples[binCount] = (TNtuple*)f->Get(Form("emu%d", binCount));
        ////Fill Histograms

        //Electron Decays
    

        emuPart->Reset();
        emuTuples[binCount]->Draw("phiDiff>>emu_part");
        emuPart->Scale(1/(*it), "width");
        emuPhiDiff->Add(emuPart); 

        emuPart->Reset();
        emuTuples[binCount]->Draw("phiDiff>>emu_part","ptE>0.5 && ptMu>1");// && etaMu<-2.5 && etaMu>-4
        emuPart->Scale(1/(*it), "width");
        emuPhiDiffCutLowEx->Add(emuPart);

        emuPart->Reset();
        emuTuples[binCount]->Draw("phiDiff>>emu_part","ptE>1 && ptMu>3");// && etaMu<-2.5 && etaMu>-4
        emuPart->Scale(1/(*it), "width");
        emuPhiDiffCutLow->Add(emuPart);

        emuPart->Reset();
        emuTuples[binCount]->Draw("phiDiff>>emu_part","ptE>3 && ptMu>5");// && etaMu<-2.5 && etaMu>-4
        emuPart->Scale(1/(*it), "width");
        emuPhiDiffCutHigh->Add(emuPart);


        binCount++;
    }

    for(std::vector<double>::iterator it = binLuminocity2->begin(); it != binLuminocity2->end(); ++it){
        
        //Fill Tuples
        qmsTuples[binCount2] = (TNtuple*)f2->Get(Form("qms%d", binCount2));

        ////Fill Histograms

        //Electron Decays
    

        emuPartSL->Reset();
        qmsTuples[binCount2]->Draw("dPhi>>SLemu_part","semilepFlag==1 && pairSemilepFlag==1");
        emuPartSL->Scale(1/(*it), "width");
        emuPhiDiffSL->Add(emuPartSL); 

        emuPartSL->Reset();
        qmsTuples[binCount2]->Draw("dPhi>>SLemu_part","pairPt>0.5 && ptL>1 && semilepFlag==1 && pairSemilepFlag==1");// && etaL<-2.5 && etaL>-4
        emuPartSL->Scale(1/(*it), "width");
        emuPhiDiffCutLowExSL->Add(emuPartSL);

        emuPartSL->Reset();
        qmsTuples[binCount2]->Draw("dPhi>>SLemu_part","pairPt>1 && ptL>3 && semilepFlag==1 && pairSemilepFlag==1");// && etaL<-2.5 && etaL>-4
        emuPartSL->Scale(1/(*it), "width");
        emuPhiDiffCutLowSL->Add(emuPartSL);

        emuPartSL->Reset();
        qmsTuples[binCount2]->Draw("dPhi>>SLemu_part","pairPt>3 && ptL>5 && semilepFlag==1 && pairSemilepFlag==1");// && etaL<-2.5 && etaL>-4
        emuPartSL->Scale(1/(*it), "width");
        emuPhiDiffCutHighSL->Add(emuPartSL);
    


        binCount2++;
    }

    


    ////Plotting
    // emus
    TFile *outf =  new TFile("A_phi_compare.root", "RECREATE");
    TCanvas *canvasCutsNC = new TCanvas("emu_cutsnc","emu_cutsnc");
    canvasCutsNC->Divide(1,3);

    canvasCutsNC->cd(1);
    emuPhiDiffSL->SetLineColor(kBlack);
    emuPhiDiffSL->SetStats(0);
    emuPhiDiffSL->Draw(); 

    canvasCutsNC->cd(2);
    emuPhiDiffCutLowSL->SetLineColor(kBlack);
    emuPhiDiffCutLowSL->SetStats(0);
    emuPhiDiffCutLowSL->Draw(); 

    auto label1nc = new TLatex();
    label1nc->DrawLatex(0.0, 0.0, "p_{T}^{e} > 1 GeV");
    label1nc->DrawLatex(0.0, 0.0, "p_{T}^{#mu} > 3 GeV");
    label1nc->Draw("SAME"); 

    canvasCutsNC->cd(3);
    emuPhiDiffCutHighSL->SetLineColor(kBlack);
    emuPhiDiffCutHighSL->SetStats(0);
    emuPhiDiffCutHighSL->Draw(); 

    auto label2nc = new TLatex();
    label2nc->DrawLatex(0.0, 0.0, "p_{T}^{e} > 3 GeV");
    label2nc->DrawLatex(0.0, 0.0, "p_{T}^{#mu} > 5 GeV");
    label2nc->Draw("SAME"); 

    canvasCutsNC->Write();
    ////

    TCanvas *canvasCuts = new TCanvas("emu_cuts","emu_cuts");
    canvasCuts->Divide(1,3);

    canvasCuts->cd(1);
    emuPhiDiff->SetLineColor(kBlue);
    emuPhiDiff->SetAxisRange(0,1700e6,"Y");
    emuPhiDiff->SetStats(0);
    emuPhiDiff->Draw();
 
    emuPhiDiffSL->SetLineColor(kBlack);
    emuPhiDiffSL->SetStats(0);
    emuPhiDiffSL->Draw("SAME"); 

    auto legend0 = new TLegend();
    legend0->AddEntry(emuPhiDiff,"All e-#mu pair candidates","l");
    legend0->AddEntry(emuPhiDiffSL,"True e-#mu pairs","l");
    legend0->Draw("SAME");

/*     canvasCuts->cd(2);
    emuPhiDiffCutLowEx->SetLineColor(kMagenta);
    emuPhiDiffCutLowEx->SetAxisRange(0,3300e3,"Y");
    emuPhiDiffCutLowEx->SetStats(0);
    emuPhiDiffCutLowEx->Draw();

    emuPhiDiffCutLowExSL->SetLineColor(kBlack);
    emuPhiDiffCutLowExSL->SetStats(0);
    emuPhiDiffCutLowExSL->Draw("SAME");

    auto label1ex = new TLatex();
    label1ex->DrawLatex(0.0, 0.0, "p_{T}^{e} > 0.5 GeV");
    label1ex->DrawLatex(0.0, 0.0, "p_{T}^{#mu} > 1 GeV");
    label1ex->Draw("SAME"); 

    auto legend1ex = new TLegend();
    legend1ex->AddEntry(emuPhiDiffCutLowEx,"All e-#mu pair candidates","l");
    legend1ex->AddEntry(emuPhiDiffCutLowExSL,"True e-#mu pairs","l");
    legend1ex->Draw("SAME"); */


    canvasCuts->cd(2);
    emuPhiDiffCutLow->SetLineColor(kGreen);
    emuPhiDiffCutLow->SetAxisRange(0,3000e3,"Y");
    emuPhiDiffCutLow->SetStats(0);
    emuPhiDiffCutLow->Draw();

     emuPhiDiffCutLowSL->SetLineColor(kBlack);
    emuPhiDiffCutLowSL->SetStats(0);
    emuPhiDiffCutLowSL->Draw("SAME"); 

    auto label1 = new TLatex();
    label1->DrawLatex(0.0, 0.0, "p_{T}^{e} > 1 GeV");
    label1->DrawLatex(0.0, 0.0, "p_{T}^{#mu} > 3 GeV");
    label1->Draw("SAME"); 

    auto legend1 = new TLegend();
    legend1->AddEntry(emuPhiDiffCutLow,"All e-#mu pair candidates","l");
    legend1->AddEntry(emuPhiDiffCutLowSL,"True e-#mu pairs","l");
    legend1->Draw("SAME");



    canvasCuts->cd(3);
    emuPhiDiffCutHigh->SetLineColor(kRed);
   emuPhiDiffCutHigh->SetAxisRange(0,250e3,"Y");
    emuPhiDiffCutHigh->SetStats(0);
    emuPhiDiffCutHigh->Draw();

    emuPhiDiffCutHighSL->SetLineColor(kBlack);
    emuPhiDiffCutHighSL->SetStats(0);
    emuPhiDiffCutHighSL->Draw("SAME"); 

    auto label2 = new TLatex();
    label2->DrawLatex(0.0, 0.0, "p_{T}^{e} > 3 GeV");
    label2->DrawLatex(0.0, 0.0, "p_{T}^{#mu} > 5 GeV");
    label2->Draw("SAME"); 

    auto legend2 = new TLegend();
    legend2->AddEntry(emuPhiDiffCutHigh,"All e-#mu pair candidates","l");
    legend2->AddEntry(emuPhiDiffCutHighSL,"True e-#mu pairs","l");
    legend2->Draw("SAME");

    

    canvasCuts->Write();


    delete outf;

}