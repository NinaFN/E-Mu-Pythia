#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

void read_sim_1D_cuts_all() {
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
    TH1F *qeAll = new TH1F("qe_dec","Individual Q #rightarrow e decay breakdown;#hat{p_{T}};#frac{d#sigma}{d#hat{p_{T}}} (pb/GeV/c)", 20, 5.0, 60.0);
    TH1F *qePart = new TH1F("qe_part","", 20, 5.0, 60.0);
    TH1F *qeDec = new TH1F("qe_det",";#hat{p_{T}};Ratio", 20, 5.0, 60.0);
    TH1F *qeEta = new TH1F("qe_pt_low",";#hat{p_{T}};Ratio", 20, 5.0, 60.0);
    TH1F *qePt = new TH1F("qe_pt_high",";#hat{p_{T}};Ratio", 20, 5.0, 60.0);

    TH1F *qmAll = new TH1F("qm_dec","Individual Q #rightarrow #mu decay breakdown;#hat{p_{T}};#frac{d#sigma}{d#hat{p_{T}}} (pb/GeV/c)", 20, 5.0, 60.0);
    TH1F *qmPart = new TH1F("qm_part","", 20, 5.0, 60.0);
    TH1F *qmDec = new TH1F("qm_det",";#hat{p_{T}};Ratio", 20, 5.0, 60.0);
    TH1F *qmEta = new TH1F("qm_pt_low",";#hat{p_{T}};Ratio", 20, 5.0, 60.0);
    TH1F *qmPt = new TH1F("qm_pt_high",";#hat{p_{T}};Ratio", 20, 5.0, 60.0);



    for(std::vector<double>::iterator it = binLuminocity->begin(); it != binLuminocity->end(); ++it){
        
        //Fill Tuples
        qeTuples[binCount] = (TNtuple*)f->Get(Form("qe%d", binCount));
        qmTuples[binCount] = (TNtuple*)f->Get(Form("qm%d", binCount));

        ////Fill Histograms

        //Electron Decays
      
         qePart->Reset();
        qeTuples[binCount]->Draw("ptHat>>qe_part","etaQ<=1.6 && etaQ>=-1.6");// && etaQ<=1.6 && etaQ>=-1.6
        qePart->Scale(1/(*it),"width");
        qeAll->Add(qePart); 

       /*  qePart->Reset();
        qeTuples[binCount]->Draw("ptHat>>qe_part");
        qePart->Scale(1/(*it),"width");
        qeAll->Add(qePart); */


        qePart->Reset();
        qeTuples[binCount]->Draw("ptHat>>qe_part","semilepFlag==1 && etaQ<=1.6 && etaQ>=-1.6");
        qePart->Scale(1/(*it),"width");
        qeDec->Add(qePart);

        qePart->Reset();
        qeTuples[binCount]->Draw("ptHat>>qe_part","etaL>=-0.9 && etaL<=0.9 && semilepFlag==1 && etaQ<=1.6 && etaQ>=-1.6");
        qePart->Scale(1/(*it),"width");
        qeEta->Add(qePart);

        qePart->Reset();
        qeTuples[binCount]->Draw("ptHat>>qe_part","ptL>=3 && ptL<=20 && semilepFlag==1 && etaL>=-0.9 && etaL<=0.9 && etaQ<=1.6 && etaQ>=-1.6");
        qePart->Scale(1/(*it),"width");
        qePt->Add(qePart);

        // qm

         qmPart->Reset();
        qmTuples[binCount]->Draw("ptHat>>qm_part","etaQ<=-2.1 && etaQ>=-4.4"); // && etaQ<=-2.1 && etaQ>=-4.4
        qmPart->Scale(1/(*it),"width");
        qmAll->Add(qmPart); 

        /* qmPart->Reset();
        qmTuples[binCount]->Draw("ptHat>>qm_part"); // && etaQ<=-2.1 && etaQ>=-4.4
        qmPart->Scale(1/(*it),"width");
        qmAll->Add(qmPart); */

        qmPart->Reset();
        qmTuples[binCount]->Draw("ptHat>>qm_part","semilepFlag==1 && etaQ<=-2.1 && etaQ>=-4.4");
        qmPart->Scale(1/(*it),"width");
        qmDec->Add(qmPart);

        qmPart->Reset();
        qmTuples[binCount]->Draw("ptHat>>qm_part","etaL>=-4 && etaL<=-2.5 && semilepFlag==1  && etaQ<=-2.1 && etaQ>=-4.4");
        qmPart->Scale(1/(*it),"width");
        qmEta->Add(qmPart);

        qmPart->Reset();
        qmTuples[binCount]->Draw("ptHat>>qm_part","ptL>=5 && semilepFlag==1 && etaL>=-4 && etaL<=-2.5 && etaQ<=-2.1 && etaQ>=-4.4");
        qmPart->Scale(1/(*it),"width");
        qmPt->Add(qmPart);

        binCount++;
    }



/*     cout<<"Electrons"<<endl;
    Double_t error_tot_qe;
    Double_t error_found_qe;

    Double_t integral_tot_qe = qeAll->IntegralAndError(0,30,error_tot_qe,"width");
    Double_t integral_found_qe = qePt->IntegralAndError(0,30,error_found_qe,"width");

    cout<<"Total pairs (Integral): "<<integral_tot_qe<<" +- "<<error_tot_qe<<endl;
    cout<<"Found pairs (Integral): "<<integral_found_qe<<" +- "<<error_found_qe<<endl;
    
    Double_t unc_qe =  (integral_found_qe/integral_tot_qe) * pow( (pow((error_found_qe/integral_found_qe),2.0) + pow((error_tot_qe/integral_tot_qe),2.0)) , (0.5));
    cout<<"Ratio (Integral): "<<integral_found_qe/integral_tot_qe<<" +- "<<unc_qe<<endl<<endl;

    cout<<"Muons"<<endl;
    Double_t error_tot_qm;
    Double_t error_found_qm;

    Double_t integral_tot_qm = qmAll->IntegralAndError(0,30,error_tot_qm,"width");
    Double_t integral_found_qm = qmPt->IntegralAndError(0,30,error_found_qm,"width");

    cout<<"Total pairs (Integral): "<<integral_tot_qm<<" +- "<<error_tot_qm<<endl;
    cout<<"Found pairs (Integral): "<<integral_found_qm<<" +- "<<error_found_qm<<endl;
    
    Double_t unc_qm =  (integral_found_qm/integral_tot_qm) * pow( (pow((error_found_qm/integral_found_qm),2.0) + pow((error_tot_qm/integral_tot_qm),2.0)) , (0.5));
    cout<<"Ratio (Integral): "<<integral_found_qm/integral_tot_qm<<" +- "<<unc_qm<<endl;
 */
    ////Plotting
    // qes
    TFile *outf =  new TFile("A_cuts_all_20_semilep_e3m5_test.root", "RECREATE");

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

    qeDec->SetLineColor(kRed);
    qeDec->SetStats(0);
    qeDec->DrawCopy("SAME");

    qeEta->SetLineColor(kViolet);
    qeEta->SetStats(0);
    qeEta->DrawCopy("SAME");

    qePt->SetLineColor(kBlue);
    qePt->SetStats(0);
    qePt->DrawCopy("SAME");

     auto labelqeEta = new TLatex();
    labelqeEta->DrawLatex(0.0, 0.0, "-1.6 #leq #eta_{Q} #leq 1.6");
    labelqeEta->Draw("SAME"); 

    auto legendqe = new TLegend();
    legendqe->AddEntry(qeAll,"All quarks","l");
    legendqe->AddEntry(qeDec,"Quark decays to e","l");
    legendqe->AddEntry(qeEta,"Quark decays to e with -0.9 #leq #eta_{e} #leq 0.9","l");
    legendqe->AddEntry(qePt,"Quark decays to e with -0.9 #leq #eta_{e} #leq 0.9 and p_{T}^{e} #geq 3 Gev","l");
    legendqe->Draw("SAME");

    canvasQEAll->cd();
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
    qeDec->SetAxisRange(0,0.5,"Y");
    qeDec->DrawCopy();

    qeEta->Divide(qeAll);
    qeEta->SetStats(0);
    qeEta->DrawCopy("SAME");

    qePt->Divide(qeAll);
    qePt->SetStats(0);
    qePt->DrawCopy("SAME");

    vector<double> ratio_qe_pt;
    vector<double> ratio_qe_pt_err;
    double intCalcQE = 0;
    for(int j = 0; j != qePt->GetSize(); ++j){
        //cout<<qePt->GetAt(j)<<"+-"<<sqrt(qePt->GetSumw2()->GetAt(j))<<endl;

        ratio_qe_pt.push_back(qePt->GetAt(j));
        ratio_qe_pt_err.push_back(qePt->GetSumw2()->GetAt(j));

        intCalcQE += qePt->GetAt(j);
    } 
    //cout<<intCalcQE<<endl;
    //cout<<intCalcQE/30<<endl;
    
    outf->WriteObject(&ratio_qe_pt, "ratio_qe_pt");
    outf->WriteObject(&ratio_qe_pt_err, "ratio_qe_pt_err");

    qePt->Write();

    canvasQEAll->Write();

    // qm

   TCanvas *canvasQMAll = new TCanvas("qm_all","qm_all");

    TPad *pad1qm = new TPad("pad1qm", "pad1qm", 0, 0.3, 1, 1.0);
    pad1qm->SetLogy();
    pad1qm->SetBottomMargin(0);
    qmAll->GetYaxis()->SetTitleOffset(1);
    pad1qm->Draw();
    pad1qm->cd(); 

    qmAll->SetLineColor(kBlack);
    qmAll->SetStats(0);
    qmAll->Draw("SAME");

    qmDec->SetLineColor(kRed);
    qmDec->SetStats(0);
    qmDec->DrawCopy("SAME");

    qmEta->SetLineColor(kViolet);
    qmEta->SetStats(0);
    qmEta->DrawCopy("SAME");

    qmPt->SetLineColor(kBlue);
    qmPt->SetStats(0);
    qmPt->DrawCopy("SAME");

    auto labelqmEta = new TLatex();
    labelqmEta->DrawLatex(0.0, 0.0, "-4.4 #leq #eta_{Q} #leq -2.1");
    labelqmEta->Draw("SAME"); 

    auto legendqm = new TLegend();
    legendqm->AddEntry(qmAll,"All quarks","l");
    legendqm->AddEntry(qmDec,"Quark decays to #mu","l");
    legendqm->AddEntry(qmEta,"Quark decays to #mu with -4.0 #leq #eta_{#mu} #leq -2.5","l");
    legendqm->AddEntry(qmPt,"Quark decays to #mu with -4.0 #leq #eta_{#mu} #leq -2.5 and p_{T}^{#mu} #geq 5 GeV","l");
    legendqm->Draw("SAME");

    canvasQMAll->cd();
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
    qmDec->SetAxisRange(0,0.5,"Y");
    qmDec->DrawCopy();

    qmEta->Divide(qmAll);
    qmEta->SetStats(0);
    qmEta->DrawCopy("SAME");

    qmPt->Divide(qmAll);
    qmPt->SetStats(0);
    qmPt->DrawCopy("SAME");

    vector<double> ratio_qm_pt;
    vector<double> ratio_qm_pt_err;
    double intCalcQM = 0;
    for(int j = 0; j != qmPt->GetSize(); ++j){
        //cout<<qmPt->GetAt(j)<<"+-"<<sqrt(qmPt->GetSumw2()->GetAt(j))<<endl;

        ratio_qm_pt.push_back(qmPt->GetAt(j));
        ratio_qm_pt_err.push_back(qmPt->GetSumw2()->GetAt(j));

        intCalcQM += qmPt->GetAt(j);
    } 
    //cout<<intCalcQM<<endl;
    //cout<<intCalcQM/30<<endl;
    
    outf->WriteObject(&ratio_qm_pt, "ratio_qm_pt");
    outf->WriteObject(&ratio_qm_pt_err, "ratio_qm_pt_err");

    qmPt->Write();

    canvasQMAll->Write();
    
    delete outf;

}