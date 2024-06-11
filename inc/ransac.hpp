#pragma once
#include<string>
#include<vector>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TColor.h>
#include <TROOT.h>
#include "TChain.h"

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "Eigen/Dense"
#include <algorithm>
#include <utility>



using namespace std;
using namespace Eigen;

struct line{
    Vector3d pos_vec;
    Vector3d Dir_vec;
    double points_avg;
    double inlier_ratio;
};



struct Pixel
{
	Double_t x;
	Double_t y;
	Int_t cl;
	Double_t ToA;
	Double_t ToT;
};

struct Cluster
{
std::vector<Pixel> pixels;
Int_t ID;
};


class Ransac

{

private:

	string fileName;

	//ransacResult result;

public:

	vector < vector <double > > dataset;

    std::vector<Pixel> Pixels;

	std::vector<Pixel> AV;

	std::vector<Cluster> Clusters;

	//std::vector<int> labels_;

	//std::vector<int> normalizedLabels_;

	//std::vector<outlierScore>outlierScores_;

	//std::vector <double> membershipProbabilities_;

	uint32_t noisyPoints_;

	uint32_t numClusters_;
    
	std::map<Int_t,Pixel> pixelMap;

    std::vector<double> distances;
    std::vector<double> angles;
	
	//vector<vector<double>> dataset;

	double dist;
	double angle;

	Ransac(vector<vector<double>> dataset) {

		this->dataset = dataset;

	}

	string getFileName();
			   
	//int loadCsv(int numberOfValues, bool skipHeader=false);

	void execute();

	void displayResult();


};

