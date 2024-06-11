#include "TH2F.h"
#include <iostream>
#include <cstdio>
#include <TApplication.h>
#include <Eigen/Dense>
#include "ransac.hpp"
using namespace std;
#include <fstream>
#include "TFile.h"
#include "TTree.h"





   UInt_t          framenr;
   UChar_t         chipnr;
   UShort_t        xpix;
   UShort_t        ypix;
   UShort_t        ToT;
   UShort_t        ToA;
   UChar_t         FToA;
   Long64_t        GToA;
   UShort_t        SpidrTime;

   TBranch        *b_framenr;   //!
   TBranch        *b_chipnr;   //!
   TBranch        *b_xpix;   //!
   TBranch        *b_ypix;   //!
   TBranch        *b_ToT;   //!
   TBranch        *b_ToA;   //!
   TBranch        *b_FtoA;   //!
   TBranch        *b_GToA;   //!
   TBranch        *b_SpidrTime;
int main(int argc, char** argv) {
	TApplication theApp("App",&argc, argv);
    TFile *file = TFile::Open("../../FRIB/CF4_OK/CF4_50torr_OK/50torr.root");

    if (!file || file->IsZombie()) {
        std::cout << "Error: cannot open the ROOT file" << std::endl;
        return 1;
    }

    TTree *tree = (TTree*)file->Get("t2");

        if (!tree) {
        std::cout << "Error: cannot find the TTree in the ROOT file" << std::endl;
        return 1;
    }
    
   tree->SetBranchAddress("framenr", &framenr, &b_framenr);
   tree->SetBranchAddress("chipnr", &chipnr, &b_chipnr);
   tree->SetBranchAddress("xpix", &xpix, &b_xpix);
   tree->SetBranchAddress("ypix", &ypix, &b_ypix);
   tree->SetBranchAddress("ToT", &ToT, &b_ToT);
   tree->SetBranchAddress("ToA", &ToA, &b_ToA);
   tree->SetBranchAddress("FToA", &FToA, &b_FtoA);
   tree->SetBranchAddress("GToA", &GToA, &b_GToA);
   tree->SetBranchAddress("SpidrTime", &SpidrTime, &b_SpidrTime);

    Long64_t nentries = tree->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    nb = tree->GetEntry(0);   nbytes += nb;
    UInt_t framenrBuff = framenr;
  
    vector<vector<double>> dataset;

for (Long64_t jentry=1; jentry<nentries;jentry++) {
    Long64_t ientry = tree->LoadTree(jentry);
      
    if (ientry < 0) break;
    
    nb = tree->GetEntry(jentry);   nbytes += nb;
      
    UShort_t xpixc = xpix - 258;
    if(framenrBuff == framenr){      
        vector<double> row = {xpixc, ypix, ToA, ToT, framenr};
        dataset.push_back(row);
        
    }
    else{   
     
       
        Ransac ransac(dataset);
        ransac.execute();
        ransac.displayResult();

        // Clear the dataset
       dataset.clear();
        dataset.shrink_to_fit(); // Deallocate memory
        framenrBuff = framenr;
    }
   
}  


if (!dataset.empty()) {
    Ransac ransac(dataset);
    ransac.execute();
    ransac.displayResult();
}

    theApp.Run();
	return 0;
}



