#ifndef constants_h
#define constants_h

const int Range_Bin = 2; 
//cut on the track 
const float Track_Pt_LowCut = 0.15;
const float Track_Eta_Cut   = 1.0;


//cut on the lambda 
const float Lambda_Rapidity_Cut = 1.0;
const float Lambda_Pt_lowCut    = 0.5; //GeV/c
const float Lambda_Pt_highCut   = 5.0; //GeV/c
const float Lambda_mass_lowCut  = 1.112250;//GeV/c^2
const float Lambda_mass_highCut = 1.119833;//GeV/c^2


//MixEvent Limitation 
const float ptDiffLim = 0.1;
const float phiDiffLim = 0.1;
const float rapidityDiffLim = 0.1;

//physics constants 
double MASS_PROTON = 0.93827208816;//GeV/c^2
double MASS_PION   = 0.13957039;   //GeV/c^2
#endif
