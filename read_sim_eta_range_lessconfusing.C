#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

//script to calculate and plot the ratios used to choose the CB parent quark eta boundaries in the approximation method

void read_sim_eta_range_lessconfusing() {
    TFile *f = TFile::Open("Data/sim_tuples_136_semilep_quarks_individual.root","READ");
    
    
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
    TH1F *qeAll = new TH1F("qe_dec","Fraction of Q #rightarrow e decays into CB;#eta_{Q};#frac{d#sigma}{d#eta_{Q}} (pb/GeV/c)", 40, -6, 6);
    TH1F *qePart = new TH1F("qe_part","", 40, -6, 6);
    TH1F *qeDec = new TH1F("qe_det",";#eta_{Q};Ratio", 40, -6, 6);

    TH1F *qmAll = new TH1F("qm_dec","Fraction of Q #rightarrow #mu decays into FWD;#eta_{Q};#frac{d#sigma}{d#eta_{Q}} (pb/GeV/c)", 40, -6, 6);
    TH1F *qmPart = new TH1F("qm_part","", 40, -6, 6);
    TH1F *qmDec = new TH1F("qm_det",";#eta_{Q};Ratio", 40, -6, 6);



    for(std::vector<double>::iterator it = binLuminocity->begin(); it != binLuminocity->end(); ++it){
        
        //Line Tuples
        qeTuples[binCount] = (TNtuple*)f->Get(Form("qe%d", binCount));
        qmTuples[binCount] = (TNtuple*)f->Get(Form("qm%d", binCount));

        ////Line Histograms

        //Electron Decays
      
         qePart->Reset();
        qeTuples[binCount]->Draw("etaQ>>qe_part","semilepFlag==1 && ((idQ*idL)<0)");// &&  
        qePart->Scale(1/(*it),"width");
        qeAll->Add(qePart); 


        qePart->Reset();
        qeTuples[binCount]->Draw("etaQ>>qe_part","semilepFlag==1 && ((idQ*idL)<0) && etaL<=0.9 && etaL>=-0.9");
        qePart->Scale(1/(*it),"width");
        qeDec->Add(qePart);

        // qm

         qmPart->Reset();
        qmTuples[binCount]->Draw("etaQ>>qm_part","semilepFlag==1 && ((idQ*idL)<0)"); // && etaQ<=-2.1 && etaQ>=-4.4
        qmPart->Scale(1/(*it),"width");
        qmAll->Add(qmPart); 

        qmPart->Reset();
        qmTuples[binCount]->Draw("etaQ>>qm_part","semilepFlag==1 && ((idQ*idL)<0) && etaL<=-2 && etaL>=-4.5");
        qmPart->Scale(1/(*it),"width");
        qmDec->Add(qmPart);



        binCount++;
    }

    ////Plotting
    // qes
    TFile *outf =  new TFile("Hists/eta_range_lessconfusing.root", "RECREATE");

    TCanvas *canvasQEAll = new TCanvas("qe_all","qe_all");

   /*  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0);
    qeAll->GetYaxis()->SetTitleOffset(1);
    pad1->Draw();
    pad1->cd();  */


    qeAll->SetLineColor(kBlack);
    qeAll->SetAxisRange(0,13e6,"Y");
    qeAll->SetStats(0);
    qeAll->Draw("SAME");

    qeDec->SetLineColor(kRed);
    qeDec->SetStats(0);
    qeDec->DrawCopy("SAME");

    TLine *qeCorr1 = new TLine(-1.6,0,-1.6,13e6);
    qeCorr1->SetLineColor(kBlack);
    qeCorr1->SetLineStyle(9);
    qeCorr1->Draw();

    TLine *qeCorr2 = new TLine(1.6,0,1.6,13e6);
    qeCorr2->SetLineColor(kBlack);
    qeCorr2->SetLineStyle(9);
    qeCorr2->Draw();

    TLine *qeCor1 = new TLine(-0.9,0,-0.9,13e6);
    qeCor1->SetLineColor(kBlack);
    qeCor1->SetLineStyle(2);
    qeCor1->Draw();

    TLine *qeCor2 = new TLine(0.9,0,0.9,13e6);
    qeCor2->SetLineColor(kBlack);
    qeCor2->SetLineStyle(2);
    qeCor2->Draw();

    Double_t error_qe;
    Double_t integral_qe = qeAll->IntegralAndError(9,14,error_qe);
    cout<<"Total (Integral): "<<integral_qe<<" +- "<<error_qe<<endl;
    
    Double_t error_qelim;
    Double_t integral_qelim = qeDec->IntegralAndError(9,14,error_qelim);
    cout<<"Limited (Integral): "<<integral_qelim<<" +- "<<error_qelim<<endl;



    auto legendqe = new TLegend();
    legendqe->AddEntry(qeAll,"All Q #rightarrow e decays","l");
    legendqe->AddEntry(qeDec,"All Q #rightarrow e decays with |#eta_{e}| #leq 0.9","l");
    legendqe->Draw("SAME");

    /* canvasQEAll->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);

    qeDec->GetXaxis()->SetTitleSize(13);
    qeDec->GetXaxis()->SetTitleOffset(1);
    qeDec->GetXaxis()->SetLabelSize(10);

    qeDec->GetYaxis()->SetTitleSize(13);
    qeDec->GetYaxis()->SetTitleOffset(1);
    qeDec->GetYaxis()->SetLabelSize(10);

    pad2->Draw("SAME");
    pad2->cd(); 

    qeDec->Divide(qeAll);
    qeDec->SetStats(0);
    qeDec->DrawCopy(); */

    


    canvasQEAll->Write();

    // qm

   TCanvas *canvasQMAll = new TCanvas("qm_all","qm_all");

    /* TPad *pad1qm = new TPad("pad1qm", "pad1qm", 0, 0.3, 1, 1.0);
    pad1qm->SetBottomMargin(0);
    qmAll->GetYaxis()->SetTitleOffset(1);
    pad1qm->Draw();
    pad1qm->cd(); 
 */
    qmAll->SetLineColor(kBlack);
    qmAll->SetAxisRange(0,11e6,"Y");
    qmAll->SetStats(0);
    qmAll->Draw("SAME");

    qmDec->SetLineColor(kRed);
    qmDec->SetStats(0);
    qmDec->DrawCopy("SAME");

    
    TLine *qmCorr1 = new TLine(-4.4,0,-4.4,11e6);
    qmCorr1->SetLineColor(kBlack);
    qmCorr1->SetLineStyle(9);
    qmCorr1->Draw();

    TLine *qmCorr2 = new TLine(-2.1,0,-2.1,11e6);
    qmCorr2->SetLineColor(kBlack);
    qmCorr2->SetLineStyle(9);
    qmCorr2->Draw();


    TLine *qmCor1 = new TLine(-4,0,-4,11e6);
    qmCor1->SetLineColor(kBlack);
    qmCor1->SetLineStyle(2);
    qmCor1->Draw();

    TLine *qmCor2 = new TLine(-2.5,0,-2.5,11e6);
    qmCor2->SetLineColor(kBlack);
    qmCor2->SetLineStyle(2);
    qmCor2->Draw();
    
    


    auto legendqm = new TLegend();
    legendqm->AddEntry(qmAll,"All Q #rightarrow #mu decays","l");
    legendqm->AddEntry(qmDec,"All Q #rightarrow #mu decays with -4.0 #leq #eta_{#mu} #leq -2.5","l");
    legendqm->Draw("SAME");

    /* canvasQMAll->cd();
    TPad *pad2qm = new TPad("pad2qm", "pad2qm", 0, 0.05, 1, 0.3);
    pad2qm->SetTopMargin(0);
    pad2qm->SetBottomMargin(0.2);

    qmDec->GetXaxis()->SetTitleSize(13);
    qmDec->GetXaxis()->SetTitleOffset(1);
    qmDec->GetXaxis()->SetLabelSize(10);

    qmDec->GetYaxis()->SetTitleSize(13);
    qmDec->GetYaxis()->SetTitleOffset(1);
    qmDec->GetYaxis()->SetLabelSize(10);

    pad2qm->Draw("SAME");
    pad2qm->cd(); 

    qmDec->Divide(qmAll);
    qmDec->SetStats(0);
    qmDec->DrawCopy(); */

    canvasQMAll->Write();

    delete outf;

}