#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

//script to plot the reduced cross section and subsequent luminosity 

void read_sim_quarkpairs_calc_xcut() {
    TFile *f = TFile::Open("Data/sim_tuples_136_semilep_quarks_pairprod.root","READ");
    TFile *f_events = TFile::Open("Hists/quarkpairs_full_xcut.root","READ");
    TFile *f_quarks = TFile::Open("Hists/cuts_all_e3m5_xcut.root","READ");

    
    vector<double> *binLuminocity;
    vector<double> *binEvCounts;

    vector<TNtuple*> quarkTuples(7); //should be count of bins but 7 is easier

    f->GetObject("luminocities",binLuminocity);
    f->GetObject("eventCounts",binEvCounts);


    int binCount = 0;

    //Create Hists
    TH1F *quarkPtTotal2 = new TH1F("hf_full_2","Events with Q#bar{Q} pairs detectable as e-#mu pairs;#hat{p}_{T} (GeV/c);#frac{d#sigma}{d#hat{p}_{T}} (pb/GeV/c)", 20, 14.0, 60.0);
    TH1F *quarkPtPart = new TH1F("quark_pt_part","", 20, 14.0, 60.0);
    TH1F *quarkPtLum = new TH1F("hf_lum","Events with Q#bar{Q} pairs detectable as e-#mu pairs;#hat{p}_{T} (GeV/c);#frac{N_{ev}}{year}", 20, 14.0, 60.0);


    TH1F* fullEvents = (TH1F*) f_events->Get("hf_full");
    TH1F* ratioEvents = (TH1F*) f_events->Get("quark_pt_6");
    
    TH1F* ratioQE = (TH1F*) f_quarks->Get("qe_pt_high");
    TH1F* ratioQM = (TH1F*) f_quarks->Get("qm_pt_high");

    quarkPtTotal2->Add(fullEvents);    
    quarkPtTotal2->Multiply(ratioEvents); //scale for num of e-mu events
    quarkPtTotal2->Multiply(ratioQE); //scale for electrons in range
    quarkPtTotal2->Multiply(ratioQM); //scale for muons in range
    quarkPtTotal2->Scale(0.093); //scale for quark parents in eta range
    //quarkPtTotal2->Scale(0.11, "width");

  /*   cout<<"!!! "<<fullEvents->GetBinContent(0)*ratioEvents->GetBinContent(0)*ratioQE->GetBinContent(0)*ratioQM->GetBinContent(0)<<endl;
    double sumcs = 0;
    for(int i = 0; i <= 21; i++)
    {
        sumcs+=quarkPtTotal2->GetBinContent(i);
        cout<<quarkPtTotal2->GetBinContent(i)<<endl;
    }

    cout<<"Sum: "<<sumcs<<endl; */

    ////Plotting
    // Quarks
    TFile *outf =  new TFile("Hists/reduced_cs_e3m5_xcut.root", "RECREATE");

    TCanvas *canvasQuark2 = new TCanvas("Quark_sigma1","Quark_sigma1");
    quarkPtTotal2->SetLineColor(kRed);
    quarkPtTotal2->SetStats(0);
    quarkPtTotal2->Draw();
    
    auto label0 = new TLatex();
    label0->DrawLatex(0.5, 0.5, "p_{T}^{#mu} #geq 5 GeV");
    label0->DrawLatex(0.5, 0.5, "p_{T}^{e} #geq 3 GeV");
    label0->Draw("SAME"); 

    canvasQuark2->Write();

    //


    TCanvas *canvasQuarkLum = new TCanvas("Quark_lum","Quark_lum");
    //gPad->SetLogy();

    quarkPtLum->Add(fullEvents);    
    quarkPtLum->Multiply(ratioEvents); //scale for num of e-mu events
    quarkPtLum->Multiply(ratioQE); //scale for electrons in range
    quarkPtLum->Multiply(ratioQM); //scale for muons in range
    quarkPtLum->Scale(0.093); //scale for quark parents in eta range
    quarkPtLum->Scale(16.9); //scale for num of e-mu events

    quarkPtLum->SetLineColor(kTeal-1);
    quarkPtLum->SetStats(0);
    quarkPtLum->Draw();

    Double_t error_lum;
    Double_t integral_lum = quarkPtLum->IntegralAndError(0,20,error_lum,"width");
    cout<<"Total Events (Integral): "<<integral_lum<<" +- "<<error_lum<<endl;

    /* double sum = 0;
    for(int i = 0; i <= 21; i++)
    {
        sum+=quarkPtLum->GetBinContent(i);
        cout<<quarkPtLum->GetBinContent(i)<<endl;
    }

    cout<<"Sum: "<<sum<<endl; */

    auto label = new TLatex();
    label->DrawLatex(3.0, 5.0, "Integrated N_{ev}/year:");
    //label->DrawLatex(3.0, 5.0, "1557 #pm 13");
    label->DrawLatex(3.0, 5.0, "142.6 #pm 1.2");
    label->DrawLatex(3.0, 5.0, "p_{T}^{#mu} #geq 5 GeV");
    label->DrawLatex(3.0, 5.0, "p_{T}^{e} #geq 3 GeV");
    label->Draw("SAME"); 

    
    canvasQuarkLum->Write();

    delete outf; 

}