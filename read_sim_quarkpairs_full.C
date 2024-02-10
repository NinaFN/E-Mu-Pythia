#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

//script to show the event classifications which describe how the heavy quark pair in each simulation decay

void read_sim_quarkpairs_full() {
    TFile *f = TFile::Open("Data/sim_tuples_136_semilep_quarks_pairprod.root","READ");
    
    
    
    vector<double> *binLuminocity;
    vector<double> *sigmaGen;
    vector<double> *binEvCounts;
    vector<double> sigmaErr = {7.6e-5, 3.7e-5, 1.5e-5, 3.5e-6, 1.7e-6, 3.1e-7, 1.4e-7};
    vector<TNtuple*> quarkTuples(7); //should be count of bins but 7 is easier

    f->GetObject("luminocities",binLuminocity);
    f->GetObject("eventCounts",binEvCounts);
    f->GetObject("sigmaGens",sigmaGen);

    int binCount = 0;

    //Create Hists
    TH1F *quarkPtTotal = new TH1F("hf_full","Events that produce c#bar{c} or b#bar{b} pairs in Pythia HardQCD processes;#hat{p}_{T} (GeV/c);#frac{d#sigma}{d#hat{p}_{T}} (pb/GeV/c)", 20, 14.0, 60);
    TH1F *quarkPtPart = new TH1F("quark_pt_part","", 20, 14.0, 60);
    TH1F *quarkPt1 = new TH1F("quark_pt_1",";#hat{p}_{T} (GeV/c);Ratio", 20, 14.0, 60);
    TH1F *quarkPt2 = new TH1F("quark_pt_2",";#hat{p}_{T} (GeV/c);Ratio", 20, 14.0, 60);
    TH1F *quarkPt3 = new TH1F("quark_pt_3","", 20, 14.0, 60);
    TH1F *quarkPt4 = new TH1F("quark_pt_4","", 20, 14.0, 60);
    TH1F *quarkPt5 = new TH1F("quark_pt_5","", 20, 14.0, 60);
    TH1F *quarkPt6 = new TH1F("quark_pt_6","", 20, 14.0, 60);
    TH1F *quarkPt7 = new TH1F("quark_pt_7","", 20, 14.0, 60);
    
    double pairFoundTot = 0;
    double pairTot = 0;

    for(std::vector<double>::iterator it = binLuminocity->begin(); it != binLuminocity->end(); ++it){
        
        //Fill Tuples
        quarkTuples[binCount] = (TNtuple*)f->Get(Form("quark%d", binCount));



        ////Fill Histograms

      
        quarkPtPart->Reset();
        quarkTuples[binCount]->Draw("ptHat>>quark_pt_part");
        quarkPtPart->Scale(1/(*it), "width");
        quarkPtTotal->Add(quarkPtPart);


        quarkPtPart->Reset();
        quarkTuples[binCount]->Draw("ptHat>>quark_pt_part","decayMap>=1 && decayMap!=2");
        quarkPtPart->Scale(1/(*it), "width");
        quarkPt1->Add(quarkPtPart);

        quarkPtPart->Reset();
        quarkTuples[binCount]->Draw("ptHat>>quark_pt_part","decayMap>=2");
        quarkPtPart->Scale(1/(*it), "width");
        quarkPt2->Add(quarkPtPart);

        quarkPtPart->Reset();
        quarkTuples[binCount]->Draw("ptHat>>quark_pt_part","decayMap>=3");
        quarkPtPart->Scale(1/(*it), "width");
        quarkPt3->Add(quarkPtPart);

        quarkPtPart->Reset();
        quarkTuples[binCount]->Draw("ptHat>>quark_pt_part","decayMap>=4");
        quarkPtPart->Scale(1/(*it), "width");
        quarkPt4->Add(quarkPtPart);

        quarkPtPart->Reset();
        quarkTuples[binCount]->Draw("ptHat>>quark_pt_part","decayMap>=5");
        quarkPtPart->Scale(1/(*it), "width");
        quarkPt5->Add(quarkPtPart);

        quarkPtPart->Reset();
        quarkTuples[binCount]->Draw("ptHat>>quark_pt_part","decayMap>=6");
        quarkPtPart->Scale(1/(*it), "width");
        quarkPt6->Add(quarkPtPart);

        quarkPtPart->Reset();
        quarkTuples[binCount]->Draw("ptHat>>quark_pt_part","decayMap>=7");
        quarkPtPart->Scale(1/(*it), "width");
        quarkPt7->Add(quarkPtPart);

        binCount++;
    }



    
    ////Plotting
    // Quarks
    TFile *outf =  new TFile("Hists/quarkpairs_full_xcut.root", "RECREATE");
    TCanvas *canvasQuark = new TCanvas("Quark_sigma","Quark_sigma");
    
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetLogy();
    pad1->SetBottomMargin(0);
    quarkPtTotal->GetYaxis()->SetTitleOffset(1);
    pad1->Draw();
    pad1->cd(); 
    quarkPtTotal->GetYaxis()->SetTitleOffset(1.2);

    quarkPtTotal->SetLineColor(kBlack);
    //quarkPtTotal->SetAxisRange(0,10e8,"Y");
    quarkPtTotal->SetStats(0);
    quarkPtTotal->Draw();

    quarkPt2->SetLineColor(kRed);
    quarkPt2->SetStats(0);
    quarkPt2->DrawCopy("SAME");

    quarkPt1->SetLineColor(kOrange);
    quarkPt1->SetStats(0);
    quarkPt1->DrawCopy("SAME");

    quarkPt3->SetLineColor(kGreen);
    quarkPt3->SetStats(0);
    quarkPt3->DrawCopy("SAME");

    quarkPt4->SetLineColor(kAzure+8);
    quarkPt4->SetStats(0);
    quarkPt4->DrawCopy("SAME");

    quarkPt5->SetLineColor(kViolet+2);
    quarkPt5->SetStats(0);
    quarkPt5->DrawCopy("SAME");

    quarkPt6->SetLineColor(kViolet);
    quarkPt6->SetStats(0);
    quarkPt6->DrawCopy("SAME");
 

    auto legendQuark = new TLegend();
    legendQuark->AddEntry(quarkPtTotal,"All","l");
    legendQuark->AddEntry(quarkPt2,"Single quark decay to #mu","l");
    legendQuark->AddEntry(quarkPt1,"Single quark decay to e","l");
    legendQuark->AddEntry(quarkPt3,"Pair decay to e and #mu","l");
    legendQuark->AddEntry(quarkPt4,"Back to back decay to e and #mu","l");
    legendQuark->AddEntry(quarkPt5,"Back to back decay to oppositely signed e and #mu","l");
    legendQuark->AddEntry(quarkPt6,"Semileptonic decay to oppositely signed e and #mu","l");
    legendQuark->Draw("SAME");

    canvasQuark->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);

    quarkPt2->GetXaxis()->SetTitleSize(13);
    quarkPt2->GetXaxis()->SetTitleOffset(1);
    quarkPt2->GetXaxis()->SetLabelSize(10);
    quarkPt2->GetYaxis()->SetTitleSize(13);
    quarkPt2->GetYaxis()->SetTitleOffset(1.6);
    quarkPt2->GetYaxis()->SetLabelSize(10);

    pad2->Draw("SAME");
    pad2->cd(); 

    quarkPt2->Divide(quarkPtTotal);
    quarkPt2->SetStats(0);
    quarkPt2->SetAxisRange(0,0.65,"Y");
    quarkPt2->DrawCopy();

    quarkPt1->Divide(quarkPtTotal);
    quarkPt1->SetStats(0);
    quarkPt1->DrawCopy("SAME");

    quarkPt3->Divide(quarkPtTotal);
    quarkPt3->SetStats(0);
    quarkPt3->DrawCopy("SAME");


    quarkPt4->Divide(quarkPtTotal);
    quarkPt4->SetStats(0);
    quarkPt4->DrawCopy("SAME");

    quarkPt5->Divide(quarkPtTotal);
    quarkPt5->SetStats(0);
    quarkPt5->DrawCopy("SAME");

    quarkPt6->Divide(quarkPtTotal);
    quarkPt6->SetStats(0);
    quarkPt6->DrawCopy("SAME");


    vector<double> ratio_events;
    vector<double> ratio_events_err;

    double intCalcQ5 = 0;
    double errCalcQ5_sq = 0;
    for(int j = 0; j != quarkPt6->GetSize(); ++j){
        cout<<quarkPt6->GetAt(j)<<"+-"<<sqrt(quarkPt6->GetSumw2()->GetAt(j))<<endl;

        ratio_events.push_back(quarkPt6->GetAt(j));
        ratio_events_err.push_back(quarkPt6->GetSumw2()->GetAt(j));

        intCalcQ5 += quarkPt6->GetAt(j);
        errCalcQ5_sq += quarkPt6->GetSumw2()->GetAt(j);
    } 
    cout<<intCalcQ5<<endl;
    cout<<intCalcQ5/20<<"+-"<<sqrt(errCalcQ5_sq)/20;
    
    outf->WriteObject(&ratio_events, "ratio_events");
    outf->WriteObject(&ratio_events_err, "ratio_events_err");

    quarkPtTotal->Write();
    quarkPt6->Write();

    

    canvasQuark->Write();

    delete outf;

}