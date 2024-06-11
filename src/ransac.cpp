#include "ransac.hpp"
#include <iostream>
#include <filesystem>
#include <chrono>
#include <thread>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <set>
#include <map>
#include <cstdint>
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TPolyLine3D.h>
#include <TF1.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TStyle.h>
#include <TColor.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TAxis3D.h>
#include "TChain.h"
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "Eigen/Dense"
#include <utility>
#include <random>
#include <Fit/Fitter.h>
#include <Math/Vector3D.h>
#include <cassert>
#include <Math/Functor.h>
#include <TF2.h>
#include <TMath.h>

using namespace std;
using namespace ROOT::Math;

double calculate_dist(Vector3d P, Vector3d Q, Vector3d R)
{
    Vector3d num_vec,den_vec;
    double numerator, denominator, result;
   num_vec = (R - P).cross(R-Q);
   numerator = num_vec.norm();
   den_vec = Q -P;
   denominator = den_vec.norm();
   result = numerator/denominator;

   return result;
}

int num_of_tracks = 0;

vector<double> diff_ToA;
  vector<double>x_inliers;
  vector<double>y_inliers;
  vector<double>z_inliers;
  vector<double>x_outliers;
  vector<double>y_outliers;
  vector<double>z_outliers;
  vector<double>x_inliers_aux;  
  vector<double>y_inliers_aux;
  vector<double>z_inliers_aux;
  vector<double>x_outliers_aux;
  vector<double>y_outliers_aux;
  vector<double>z_outliers_aux;


  vector<double>x_line;
  vector<double>y_line;
  vector<double>z_line;

vector<double> x_centroid;
vector<double> y_centroid;
vector<double> ToT_values;
TH2F* track = new TH2F("track","Track; xpix ; ypix",256,0.,255.,256,0.,255.);

line line_fitting(double,double, vector<Vector3d>& points);

void linef(double t, const double *p, double &x, double &y, double &z) {
   x = p[0] + p[1]*t;
   y = p[2] + p[3]*t;
   z = t;
}

bool first = true;

struct SumDistance2 {
   TGraph2D *fGraph;
 
   SumDistance2(TGraph2D *g) : fGraph(g) {}
 
   double distance2(double x,double y,double z, const double *p) {
      XYZVector xp(x,y,z);
      XYZVector x0(p[0], p[2], 0. );
      XYZVector x1(p[0] + p[1], p[2] + p[3], 1. );
      XYZVector u = (x1-x0).Unit();
      double d2 = ((xp-x0).Cross(u)).Mag2();
      return d2;
   }
double operator() (const double *par) {
      assert(fGraph != nullptr);
      double * x = fGraph->GetX();
      double * y = fGraph->GetY();
      double * z = fGraph->GetZ();
      int npoints = fGraph->GetN();
      double sum = 0;
      for (int i  = 0; i < npoints; ++i) {
         double d = distance2(x[i],y[i],z[i],par);
         sum += d;
      }
      if (first) {
         //std::cout << "Total Initial distance square = " << sum << std::endl;
      }
      first = false;
      return sum;
   }
 
};


 vector<double>num_values;
 vector<vector<double>> dataset_aux;
 
bool th = true;
float ransac_ratio = 0.1;
double ransac_thresh;
vector<Vector3d>points;
vector<Vector3d>Direc_vec;
vector<Vector3d>Posit_vec;
vector<Vector3d>inliers;
Vector3d Pfinal,Qfinal;



TPolyLine3D *l1 = new TPolyLine3D();
void Ransac::execute() {
    
 if(th){
cout << "Threshold: ";
cin >> ransac_thresh;
th = false;
}

srand(time(NULL));
//set the sample data

//double m = 0.05;
float random = 0.5f;
float noise = 0.f;

const static int q = 15;
const static float c1 = (1 << q) - 1;
const static float c2 = ((int)(c1 / 3)) + 1;
const static float c3 = 1.f / c1;

random = ((float)rand()/(float)RAND_MAX + 1);
noise = (2.f * ((random * c2) + (random * c2) + (random * c2)) - 3.f * (c2 - 1.f)) * c3;
float noise_scale = 2.0;





for (uint32_t i = 0; i < dataset.size(); i++){
Vector3d Pt;
Pt[0] = dataset[i][0];
Pt[1] = dataset[i][1];
Pt[2] = dataset[i][2];
points.push_back(Pt);
}




line new_line = line_fitting(ransac_ratio,ransac_thresh,points);
Vector3d Dir = new_line.Dir_vec;
Vector3d Pos = new_line.pos_vec;

}
bool niter = true;

 double ransac_iter;
int num = 0;
line line_fitting(double ransac_ratio, double ransac_thresh,vector<Vector3d>& points)
{
    
   if(niter){
    std::cout << "Number of iterations: ";
    std::cin >> ransac_iter;
    niter = false;
     }
 line l;
    int n_samples = points.size();

std::random_device rd;
std::mt19937 g(rd());
double ratio = 0.0;

  vector<int>indices;



for(int i = 0;i<n_samples;i++)
    indices.push_back(i);


double pt1,pt2;
Vector3d P,Q, dir_vec;
for (int it = 0;it<ransac_iter;it++)
{
   
    int p=0;
        
    vector<vector<double>>maybe_pts;
    vector<vector<double>>test_pts;
    shuffle(indices.begin(),indices.end(),g);
    for(int j=0;j<indices.size();j++)
    {
        vector<double>maybe;
        vector<double>test;
        if(j<2)
        {
            int idx1 = indices.at(j);
            while(p<3)
            {
                pt1 = (points.at(idx1))(p);
               maybe.push_back(pt1);
               p++;


            }
            maybe_pts.push_back(maybe);

        }
        p = 0;
        if(j>= 2)
        {
            int idx2 = indices.at(j);
            while(p<3)
            {
                pt2 = (points.at(idx2))(p);
                test.push_back(pt2);
                p++;
            }
            test_pts.push_back(test);
        }
    }

    for(int i=0;i<2;i++)
    {
        if(i==0)
        {
            P(0) = maybe_pts.at(i).at(0);
            P(1) = maybe_pts.at(i).at(1);
            P(2) = maybe_pts.at(i).at(2);
        }
        if(i==1)
        {
            Q(0) = maybe_pts.at(i).at(0);
            Q(1) = maybe_pts.at(i).at(1);
            Q(2) = maybe_pts.at(i).at(2);
        }
    }
   dir_vec = Q - P;
  vector<double>out_dist;
  vector<double>points_dist;
  Vector3d test_points;
num = 0;
double dist;
double Sum = 0;
for(int id = 0;id<test_pts.size();id++)
{
   double x0 = test_pts.at(id).at(0);
   double y0 = test_pts.at(id).at(1);
   double z0 = test_pts.at(id).at(2);

   test_points(0) = x0;
   test_points(1) = y0;
   test_points(2) = z0;
   
  dist = calculate_dist(P,Q,test_points);
  points_dist.push_back(dist);

   if(dist <= ransac_thresh)
   {
       x_inliers_aux.push_back(x0);
       y_inliers_aux.push_back(y0);
       z_inliers_aux.push_back(z0);
       num++;


   }
   if(dist > ransac_thresh)
   {
       out_dist.push_back(dist);
       x_outliers_aux.push_back(x0);
       y_outliers_aux.push_back(y0);
       z_outliers_aux.push_back(z0);

   }
}
for(int i=0;i<points_dist.size();i++)
{
    Sum += points_dist.at(i);
}
double avg = Sum/n_samples;



num = num+2;

   if(num/(float)n_samples > ratio)
   {
       ratio = (num)/(float)n_samples;
       l.Dir_vec = dir_vec;
       l.pos_vec = P;
       l.inlier_ratio = ratio;
       l.points_avg = avg;
    
    x_inliers = x_inliers_aux;
    y_inliers = y_inliers_aux;
    z_inliers = z_inliers_aux;
    x_inliers_aux.clear();
    y_inliers_aux.clear();
    z_inliers_aux.clear();

    x_outliers = x_outliers_aux;
    y_outliers = y_outliers_aux;
    z_outliers = z_outliers_aux;
    x_outliers_aux.clear();
    y_outliers_aux.clear();
    z_outliers_aux.clear();

    Pfinal = P;
    Qfinal = Q;

   }

   else{
    x_inliers_aux.clear();
    y_inliers_aux.clear();
    z_inliers_aux.clear();
    x_outliers_aux.clear();
    y_outliers_aux.clear();
    z_outliers_aux.clear();
   }
}

if(x_inliers.size() > 90){
    Direc_vec.push_back(l.Dir_vec);
    Posit_vec.push_back(l.pos_vec);
    num_values.push_back(x_inliers.size() + 2);
    for (int i = 0; i < x_inliers.size(); i++){
     Vector3d Pt;
     Pt[0] = x_inliers[i];
     Pt[1] = y_inliers[i];
     Pt[2] = z_inliers[i];
     inliers.push_back(Pt);
    }

    inliers.push_back(Pfinal);
    inliers.push_back(Qfinal);

    for (int i = 0; i < x_inliers.size(); i++) {
     Vector3d point(x_inliers[i], y_inliers[i], z_inliers[i]);
     points.erase(std::remove(points.begin(), points.end(), point), points.end());
    }

     points.erase(std::remove(points.begin(), points.end(), Pfinal), points.end());
     points.erase(std::remove(points.begin(), points.end(), Qfinal), points.end());

    x_inliers.clear();
    y_inliers.clear();
    z_inliers.clear();


    line new_line = line_fitting(ransac_ratio, ransac_thresh, points);
}

return l;

}



void Ransac::displayResult() {
	cout << endl << endl;


    diff_ToA = z_inliers;


 std::sort(diff_ToA.begin(), diff_ToA.end());


  auto last2 = std::unique(diff_ToA.begin(), diff_ToA.end());
   
    diff_ToA.erase(last2, diff_ToA.end());


vector<Pixel> pixel_in;
for(int i = 0; i<diff_ToA.size(); i++){
   
    for (int j = 0; j < x_inliers.size(); j++)
    {
        for (int k = 0; k < dataset.size(); k++)
        {

            if(diff_ToA[i] == dataset[k][2] && x_inliers[j] == dataset[k][0] && y_inliers[j] == dataset[k][1])
            {
        
                Pixel p;
                p.x = x_inliers[j];
                p.y = y_inliers[j];
                p.cl = dataset[k][2];
                p.ToA = dataset[k][2];
                p.ToT = dataset[k][3];
                pixel_in.push_back(p);
            }

            }
        }
        
    }
    
vector<Pixel> A[diff_ToA.size()];

for(uint32_t i = 0; i < diff_ToA.size(); i++){
uint32_t index = diff_ToA[i];
	 auto px = std::copy_if(pixel_in.begin(), pixel_in.end(), std::back_inserter(A[i]), [index](Pixel p){
             return p.ToA == index;
			 });
}


if(num_values.size() == 0){
    cout<<"There are no tracks!!"<<endl;
   }
if(num_values.size() == 1){
    cout<<"There is 1 track!!"<<endl;
}
if(num_values.size() >= 2){

       cout<<"There are "<<num_values.size()<<" tracks!!"<<endl;
       cout<<endl;
       cout<<"Number of event: "<<dataset[0][4]<<endl;
       cout<<"Checking if it is a scattering event..."<<endl;
       cout<<endl;

  for (int i = 0; i < Posit_vec.size(); i++) {
    for (int j = i + 1; j < Posit_vec.size(); j++) {
        // Calculate the shortest distance between the two lines
        double P21[3] = {Posit_vec[j][0] - Posit_vec[i][0], Posit_vec[j][1] - Posit_vec[i][1], Posit_vec[j][2] - Posit_vec[i][2]};
        double cross[3] = {Direc_vec[i][1] * Direc_vec[j][2] - Direc_vec[i][2] * Direc_vec[j][1], Direc_vec[i][2] * Direc_vec[j][0] - Direc_vec[i][0] * Direc_vec[j][2], Direc_vec[i][0] * Direc_vec[j][1] - Direc_vec[i][1] * Direc_vec[j][0]};
        double dot = P21[0] * cross[0] + P21[1] * cross[1] + P21[2] * cross[2];
        double crossMagnitude = sqrt(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]);
        double dist = std::abs(dot) / crossMagnitude;
        distances.push_back(dist);
          //  if (dist <= 3){
        // Calculate the angle between the two direction vectors
        double dotProduct = Direc_vec[i][0] * Direc_vec[j][0] + Direc_vec[i][1] * Direc_vec[j][1] + Direc_vec[i][2] * Direc_vec[j][2];
        double magnitudeI = sqrt(Direc_vec[i][0] * Direc_vec[i][0] + Direc_vec[i][1] * Direc_vec[i][1] + Direc_vec[i][2] * Direc_vec[i][2]);
        double magnitudeJ = sqrt(Direc_vec[j][0] * Direc_vec[j][0] + Direc_vec[j][1] * Direc_vec[j][1] + Direc_vec[j][2] * Direc_vec[j][2]);
        double angle = acos(dotProduct / (magnitudeI * magnitudeJ));
        angle = angle * 180.0 / M_PI;
        angles.push_back(angle);


        std::cout << "Scattering event detected with an angle of "<<angle<<" degrees" << std::endl;
       //}
    }
}

}

num_values.clear();
Posit_vec.clear();
Direc_vec.clear();
inliers.clear();

}
