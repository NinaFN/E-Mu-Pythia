#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

//script to compare the kinematic variables of the two quarks in a pair


void read_sim_2D_qqbar() {
    TFile *f = TFile::Open("sim_tuples_quarks_split_old.root","READ");
    
    
    vector<double> *binLuminocity;
    vector<double> *binEvCounts;
    vector<double> *weightSums;
    vector<double> *sigmaGens;

    vector<TNtuple*> qqTuples(6); //should be count of bins but 7 is easier

    f->GetObject("luminocities",binLuminocity);
    f->GetObject("eventCounts",binEvCounts);
    f->GetObject("weightSums",weightSums);
    f->GetObject("sigmaGens",sigmaGens);

    int binCount = 0;

    //Create Hists

    TH2F *qq2DPartEta = new TH2F("qq_part_eta","", 50, -8.0, 8.0, 50, -8.0, 8.0);
    TH2F *qq2DEta = new TH2F("qq_eta","#eta_{Q1} vs #eta_{Q2};#eta_{Q1};#eta_{Q2}", 50, -8.0, 8.0, 50, -8.0, 8.0);
    
    TH2F *qq2DPartPt = new TH2F("qq_part_pt","", 50, 5, 60, 50, 5, 60);
    TH2F *qq2DPt = new TH2F("qq_pt","p_{T}^{Q1} vs p_{T}^{Q2};p_{T}^{Q1};p_{T}^{Q2}", 50, 5, 60, 50, 5, 60);
    //TH2F *qq2DPt = new TH2F("qq_pt","y_{Q1} vs y_{Q2};y_{Q1};y_{Q2}", 50, -8.0, 8.0, 50, -8.0, 8.0);

    TH1F *qq2DPartDiffY1 = new TH1F("qq_part_diffy1","", 50, -15.0, 15.0);
    TH1F *qq2DPartDiffY2 = new TH1F("qq_part_diffy2","", 50, -15.0, 15.0);
    TH1F *qq2DDiffY = new TH1F("qq_diffy","Difference in quark pair rapidity;#Delta y;N", 50, -15.0, 15.0);


    TH2F *qq2DPartPhi = new TH2F("qq_part_phi","", 50, -M_PI, M_PI, 50, -M_PI, M_PI);   
    TH2F *qq2DPhi = new TH2F("qq_phi","#phi_{Q1} vs #phi_{Q2};#phi_{Q1};#phi_{Q2}", 50, -M_PI, M_PI, 50, -M_PI, M_PI);

    TH2F *qq2DPartTheta = new TH2F("qq_part_theta","", 50, -0, M_PI, 50, -0, M_PI);   
    TH2F *qq2DTheta = new TH2F("qq_theta","#theta_{Q1} vs #theta_{Q2};#theta_{Q1};#theta_{Q2}", 50, -0, M_PI, 50, -0, M_PI);


    TH1F *qq1DPart = new TH1F("qq_part","",50,5,60);
    TH1F *qq1DFull = new TH1F("qq_full",";#hat{p_{T}};#frac{d#sigma}{dp_{T}} (pb/GeV/c)",50,5,60);
    TH1F *qq1DObs = new TH1F("qq_obs",";#hat{p_{T}};#frac{d#sigma}{dp_{T}} (pb/GeV/c)",50,5,60);


    double pairFoundTot = 0;
    double pairTot = 0;

    for(std::vector<double>::iterator it = binLuminocity->begin(); it != binLuminocity->end(); ++it){
        
        //Fill Tuples
        qqTuples[binCount] = (TNtuple*)f->Get(Form("qp%d", binCount));

        //calculations
        double pairFoundBin = qqTuples[binCount]->GetEntries("((etaQ1<=1.6 && etaQ1>=-1.6 && etaQ2<=-2.1 && etaQ2>=-4.4) || (etaQ2<=1.6 && etaQ2>=-1.6 && etaQ1<=-2.1 && etaQ1>=-4.4))");
        pairFoundTot+=pairFoundBin;

        double pairBin = qqTuples[binCount]->GetEntries();
        pairTot+=pairBin;

        cout<<"Bin "<<binCount<<endl;
        cout<<"Total pairs: "<<pairBin<<endl;
        cout<<"Found pairs: "<<pairFoundBin<<endl;

        ////Fill Histograms

        //Electron Decays

        qq2DPartEta->Reset();
        qqTuples[binCount]->Draw("etaQ1:etaQ2>>qq_part_eta");
        qq2DPartEta->Scale(1/(*it),"width");
        qq2DEta->Add(qq2DPartEta);

        /* qq2DPartPt->Reset();
        qqTuples[binCount]->Draw("yQ1:yQ2>>qq_part_pt");
        qq2DPartPt->Scale(1/(*it),"width");
        qq2DPt->Add(qq2DPartPt); */

        qq2DPartPt->Reset();
        qqTuples[binCount]->Draw("ptQ1:ptQ2>>qq_part_pt");
        qq2DPartPt->Scale(1/(*it),"width");
        qq2DPt->Add(qq2DPartPt);

        qq2DPartPhi->Reset();
        qqTuples[binCount]->Draw("phiQ1:phiQ2>>qq_part_phi");
        qq2DPartPhi->Scale(1/(*it),"width");
        qq2DPhi->Add(qq2DPartPhi);

        qq2DPartTheta->Reset();
        qqTuples[binCount]->Draw("thetaQ1:thetaQ2>>qq_part_theta");
        qq2DPartTheta->Scale(1/(*it),"width");
        qq2DTheta->Add(qq2DPartTheta);

        //

        qq1DPart->Reset();
        qqTuples[binCount]->Draw("ptHat>>qq_part");
        qq1DPart->Scale(1/(*it),"width");
        qq1DFull->Add(qq1DPart);

        qq1DPart->Reset();
        qqTuples[binCount]->Draw("ptHat>>qq_part","((etaQ1<=1.6 && etaQ1>=-1.6 && etaQ2<=-2.1 && etaQ2>=-4.4) || (etaQ2<=1.6 && etaQ2>=-1.6 && etaQ1<=-2.1 && etaQ1>=-4.4))");
        qq1DPart->Scale(1/(*it),"width");
        qq1DObs->Add(qq1DPart);

        //

        qq2DPartDiffY1->Reset();
        qq2DPartDiffY2->Reset();
        qqTuples[binCount]->Draw("yQ1>>qq_part_diffy1");
        qqTuples[binCount]->Draw("yQ2>>qq_part_diffy2");
        qq2DPartDiffY1->Add(qq2DPartDiffY2,-1);
        qq2DPartDiffY1->Scale(1/(*it),"width");
        qq2DDiffY->Add(qq2DPartDiffY1);

        binCount++;
    }
    cout<<endl;
    cout<<"Total pairs: "<<pairTot<<endl;
    cout<<"Found pairs: "<<pairFoundTot<<endl;

    ////Plotting
    // qqs
    TFile *outf =  new TFile("sim_hists_2D_qqbar_hardQCD.root", "RECREATE");

    TCanvas *canvasQQAll = new TCanvas("qq_all","qq_all");
    canvasQQAll->Divide(2,2);


    canvasQQAll->cd(1);
    //gPad->SetLogz();
    qq2DEta->SetStats(0);
    qq2DEta->GetZaxis()->SetRangeUser(0, 1e7);
    qq2DEta->Draw("colz");

    canvasQQAll->cd(2);
    //gPad->SetLogz();
    qq2DPt->SetStats(0);
    qq2DPt->GetZaxis()->SetRangeUser(0, 1e7);
    qq2DPt->Draw("colz");

    canvasQQAll->cd(3);
    //gPad->SetLogz();
    qq2DPhi->SetStats(0);
    qq2DPhi->GetZaxis()->SetRangeUser(0, 1e7);
    qq2DPhi->Draw("colz");

    canvasQQAll->cd(4);
    //gPad->SetLogz();
    qq2DTheta->SetStats(0);
    qq2DTheta->GetZaxis()->SetRangeUser(0, 1e7);
    qq2DTheta->Draw("colz");

    canvasQQAll->Write();


//

    TCanvas *canvasQQRegions = new TCanvas("qq_reg","qq_reg");
    
    ////gPad->SetLogz();
    qq2DEta->SetStats(0);
    qq2DEta->GetZaxis()->SetRangeUser(0, 1e7);
    qq2DEta->Draw("colz");

    TBox *Corr1 = new TBox(-1.6,-4.4,1.6,-2.1);
    TBox *Corr2 = new TBox(-4.4,-1.6,-2.1,1.6);

    Corr1->SetFillColor(kWhite);
    Corr2->SetFillColor(kWhite);

    Corr1->SetFillStyle(3444);
    Corr2->SetFillStyle(3444);

    Corr1->Draw();
    Corr2->Draw();

    canvasQQRegions->Write();

    //

    TCanvas *canvas1D = new TCanvas("qq","qq");
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetLogy();
    pad1->SetBottomMargin(0);
    qq1DFull->GetYaxis()->SetTitleOffset(1);
    pad1->Draw();
    pad1->cd(); 

    qq1DFull->SetLineColor(kBlack);
    qq1DFull->SetStats(0);
    qq1DFull->Draw("SAME");

    qq1DObs->SetLineColor(kGreen);
    qq1DObs->SetStats(0);
    qq1DObs->DrawCopy("SAME");

    canvas1D->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);

    qq1DObs->GetXaxis()->SetTitleSize(13);
    qq1DObs->GetXaxis()->SetTitleOffset(1);
    qq1DObs->GetXaxis()->SetLabelSize(10);

    qq1DObs->GetYaxis()->SetTitleSize(13);
    qq1DObs->GetYaxis()->SetTitleOffset(1);
    qq1DObs->GetYaxis()->SetLabelSize(10);

    pad2->Draw("SAME");
    pad2->cd(); 

    qq1DObs->Divide(qq1DFull);
    qq1DObs->SetStats(0);
    qq1DObs->SetAxisRange(0,0.06,"Y");
    qq1DObs->DrawCopy();

    canvas1D->Write();



//

    TCanvas *canvasQQy = new TCanvas("qq_y","qq_y");
    canvasQQy->Divide(2,2);

    canvasQQy->cd(1);
    //gPad->SetLogz();
    qq2DPt->SetStats(0);
    qq2DPt->GetZaxis()->SetRangeUser(0, 1e7);
    qq2DPt->Draw("colz");

    canvasQQy->cd(2);
    //gPad->SetLogz();
    qq2DDiffY->SetStats(0);
    qq2DDiffY->Draw("colz");

    canvasQQy->Write();


    delete outf;

}