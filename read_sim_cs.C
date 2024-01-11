#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <TNtuple.h>
#include <cmath>

//script to read the pythia event information (lum, cs, etc) saved into the root data files 

void read_sim_cs() {
    TFile *f = TFile::Open("sim_tuples_pairs_nocheck_mb.root","READ");
    
    vector<double> *binLuminocity;
    vector<double> *binEvCounts;
    vector<double> *weightSums;
    vector<double> *sigmaGens;
    //vector<double> *sigmaErr;

    vector<TNtuple*> qeTuples(8); //should be count of bins but 7 is easier
    vector<TNtuple*> qmTuples(8);

    f->GetObject("luminocities",binLuminocity);
    f->GetObject("eventCounts",binEvCounts);
    f->GetObject("weightSums",weightSums);
    f->GetObject("sigmaGens",sigmaGens);
   // f->GetObject("sigmaErr",sigmaErr);

    int binCount = 0;

    //Create Hists

    cout<<"Bin Luminocity:"<<endl;
    for(std::vector<double>::iterator it = binLuminocity->begin(); it != binLuminocity->end(); ++it){
        
        //Fill Tuples
        cout<<*it<<"\t";
        binCount++;
    }
    cout<<endl;

    cout<<"Bin Ev Counts:"<<endl;
    for(std::vector<double>::iterator it = binEvCounts->begin(); it != binEvCounts->end(); ++it){
        
        //Fill Tuples
        cout<<*it<<"\t";
        binCount++;
    }
    cout<<endl;

    cout<<"Weight Sums:"<<endl;
    for(std::vector<double>::iterator it = weightSums->begin(); it != weightSums->end(); ++it){
        
        //Fill Tuples
        cout<<*it<<"\t";
        binCount++;
    }
    cout<<endl;

    cout<<"Sigma Gens:"<<endl;
    for(std::vector<double>::iterator it = sigmaGens->begin(); it != sigmaGens->end(); ++it){
        
        //Fill Tuples
        cout<<*it<<"\t";
        binCount++;
    }
    cout<<endl;

    cout<<"Sigma Gens:"<<endl;
    for(std::vector<double>::iterator it = sigmaGens->begin(); it != sigmaGens->end(); ++it){
        
        //Fill Tuples
        cout<<*it<<"\t";
        binCount++;
    }
    cout<<endl;

/*     cout<<"Sigma Errs:"<<endl;
    for(std::vector<double>::iterator it = sigmaErr->begin(); it != sigmaErr->end(); ++it){
        
        //Fill Tuples
        cout<<*it<<"\t";
        binCount++;
    } */
    cout<<endl; 

 
}