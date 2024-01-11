#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

//script to plot the reduced cross section and subsequent luminosity 

void read_sim_quarkpairs_calc() {
    TFile *f = TFile::Open("sim_tuples_136_semilep_quarks_pairprod.root","READ");
    TFile *f_events = TFile::Open("A_quarkpairs_full_semilep.root","READ");
    TFile *f_quarks = TFile::Open("A_cuts_all_20_semilep_e1m3.root","READ");

    
    vector<double> *binLuminocity;
    vector<double> *binEvCounts;

    vector<TNtuple*> quarkTuples(7); //should be count of bins but 7 is easier

    f->GetObject("luminocities",binLuminocity);
    f->GetObject("eventCounts",binEvCounts);


    int binCount = 0;

    //Create Hists
    TH1F *quarkPtTotal2 = new TH1F("hf_full_2","Events with Q#bar{Q} pairs detectable as e-#mu pairs;#hat{p}_{T} (GeV/c);#frac{d#sigma}{d#hat{p}_{T}} (pb/GeV/c)", 20, 5.0, 60.0);
    TH1F *quarkPtPart = new TH1F("quark_pt_part","", 20, 5.0, 60.0);
    TH1F *quarkPtLum = new TH1F("hf_lum","Events with Q#bar{Q} pairs detectable as e-#mu pairs;#hat{p}_{T} (GeV/c);#frac{N_{ev}}{year}", 20, 5.0, 60.0);


    TH1F* fullEvents = (TH1F*) f_events->Get("hf_full");
    TH1F* ratioEvents = (TH1F*) f_events->Get("quark_pt_6");
    
    TH1F* ratioQE = (TH1F*) f_quarks->Get("qe_pt_high");
    TH1F* ratioQM = (TH1F*) f_quarks->Get("qm_pt_high");

    quarkPtTotal2->Add(fullEvents);    
    quarkPtTotal2->Multiply(ratioEvents); //scale for num of e-mu events
    quarkPtTotal2->Multiply(ratioQE); //scale for electrons in range
    quarkPtTotal2->Multiply(ratioQM); //scale for muons in range
    quarkPtTotal2->Scale(0.11, "width"); //scale for quark parents in eta range


    ////Plotting
    // Quarks
    TFile *outf =  new TFile("A_reduced_cs_e1m3.root", "RECREATE");

    TCanvas *canvasQuark2 = new TCanvas("Quark_sigma1","Quark_sigma1");
    //gPad->SetLogy();

    //auto* f1  = new TF1("f1","crystalball",-5,5);
    //f1->SetParameters(1, 0, 1, 2, 0.5);
    //quarkPtTotal2->Fit(f1);
    quarkPtTotal2->SetLineColor(kRed);
    quarkPtTotal2->SetStats(0);
    quarkPtTotal2->Draw();
    
    auto label0 = new TLatex();
    label0->DrawLatex(0.0, 0.0, "p_{T}^{#mu} #geq 3 GeV");
    label0->DrawLatex(0.0, 0.0, "p_{T}^{e} #geq 1 GeV");
    label0->Draw("SAME"); 

    canvasQuark2->Write();

    //


    TCanvas *canvasQuarkLum = new TCanvas("Quark_lum","Quark_lum");
    //gPad->SetLogy();

    quarkPtLum->Add(fullEvents);    
    quarkPtLum->Multiply(ratioEvents); //scale for num of e-mu events
    quarkPtLum->Multiply(ratioQE); //scale for electrons in range
    quarkPtLum->Multiply(ratioQM); //scale for muons in range
    quarkPtLum->Scale(0.11, "width"); //scale for quark parents in eta range
    quarkPtLum->Scale(0.032); //scale for num of e-mu events

    quarkPtLum->SetLineColor(kTeal-1);
    quarkPtLum->SetStats(0);
    quarkPtLum->Draw();

    Double_t error_lum;
    Double_t integral_lum = quarkPtLum->IntegralAndError(0,30,error_lum,"width");
    cout<<"Total Events (Integral): "<<integral_lum<<" +- "<<error_lum<<endl;

    auto label = new TLatex();
    label->DrawLatex(0.0, 0.0, "Integrated N_{ev}/year:");
    label->DrawLatex(0.0, 0.0, "2301 #pm 13");
    label->DrawLatex(0.0, 0.0, "p_{T}^{#mu} #geq 3 GeV");
    label->DrawLatex(0.0, 0.0, "p_{T}^{e} #geq 1 GeV");
    label->Draw("SAME"); 

    
    canvasQuarkLum->Write();

    delete outf; 

}