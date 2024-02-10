#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

void read_Run3() {
    TFile *f = TFile::Open("AnalysisResults_e1m3.root","READ");
    TFile *outf =  new TFile("hists_e1m3.root", "RECREATE");

    TH1F *dPhiAll;
    TH1F *dPhiOS;
    TH1F *dPhiLS;

    TH1F *ptE;
    TH1F *ptMu;

    f->GetObject("analysis-light/PhiCorrelation_All",dPhiAll);
    f->GetObject("analysis-light/PhiCorrelation_OS",dPhiOS);
    f->GetObject("analysis-light/PhiCorrelation_LS",dPhiLS);
    f->GetObject("analysis-light/BarrelPt",ptE);
    f->GetObject("analysis-light/MuonPt",ptMu);

    TCanvas *all = new TCanvas("all","all");
    all->Divide(1,3);

    all->cd(1);
    dPhiAll->Draw();
    all->cd(2);
    dPhiOS->Draw();
    all->cd(3);
    dPhiLS->Draw();

    all->Write();

    TCanvas *pt = new TCanvas("pt","pt");
    pt->Divide(1,2);

    pt->cd(1);
    gPad->SetLogy();
    ptE->SetTitle("Barrel Track Transverse Momenta");
    ptE->Draw();

    pt->cd(2);
    gPad->SetLogy();
    ptMu->SetTitle("Forward Track Transverse Momenta");
    ptMu->Draw();


    pt->Write();
    

    delete outf;
}