#ifndef constants_h
#define constants_h

const int Range_Bin = 3; 
const int Pt1_Bin = 3 ; 
const int Pt2_Bin = 5;
//cut on the track 
const float Track_Pt_LowCut = 0.15;
const float Track_Eta_Cut = 1.0;


//cut on the lambda 
const float Lambda_Rapidity_Cut = 1.0;
const float Lambda_Pt_lowCut    = 0.5; //GeV/c
const float Lambda_Pt_highCut   = 5.0; //GeV/c 

//lambda pt bin 
const float Lambda_Pt_LowerBin[Pt1_Bin] = {0.5,1.5,0.5};
const float Lambda_Pt_HigherBin[Pt1_Bin]= {1.5,5.0,5.0};

#endif 