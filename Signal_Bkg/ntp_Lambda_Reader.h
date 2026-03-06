#ifndef ntp_Lambda_Reader_h
#define ntp_Lambda_Reader_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <vector>
#include <string>

// Header file for the classes stored in the TTree if any.

class ntp_Lambda_Reader {
public :
   std::vector<std::string> InputFiles;

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           eventId;
   Float_t         Vz;
   Int_t           mNTrigs;
   Int_t           mTrigId[3];   //[mNTrigs]
   Int_t           mNTrks;
   Float_t         high_pt[6];   //[mNTrks]
   Float_t         high_phi[6];   //[mNTrks]
   Float_t         high_eta[6];   //[mNTrks]
   Int_t           NLambda;
   Int_t           p1_InEventID[5];   //[NLambda]
   Float_t         p1_pt[5];   //[NLambda]
   Float_t         p1_phi[5];   //[NLambda]
   Float_t         p1_eta[5];   //[NLambda]
   Float_t         p1_dca[5];   //[NLambda]
   Int_t           p1_ch[5];   //[NLambda]
   Int_t           p1_hasTOFinfo[5];   //[NLambda]
   Float_t         p1_dedx[5];   //[NLambda]
   Float_t         p1_beta[5];   //[NLambda]
   Int_t           p2_InEventID[5];   //[NLambda]
   Float_t         p2_pt[5];   //[NLambda]
   Float_t         p2_phi[5];   //[NLambda]
   Float_t         p2_eta[5];   //[NLambda]
   Float_t         p2_dca[5];   //[NLambda]
   Float_t         p2_ch[5];   //[NLambda]
   Int_t           p2_hasTOFinfo[5];   //[NLambda]
   Float_t         p2_dedx[5];   //[NLambda]
   Float_t         p2_beta[5];   //[NLambda]
   Int_t           pair_charge[5];   //[NLambda]
   Float_t         pair_DCAdaughters[5];   //[NLambda]
   Float_t         pair_theta[5];   //[NLambda]
   Float_t         pair_decayL[5];   //[NLambda]
   Float_t         pair_phi[5];   //[NLambda]
   Float_t         pair_eta[5];   //[NLambda]
   Float_t         pair_pt[5];   //[NLambda]
   Float_t         pair_mass[5];   //[NLambda]

   // List of branches
   TBranch        *b_eventId;   //!
   TBranch        *b_Vz;   //!
   TBranch        *b_mNTrigs;   //!
   TBranch        *b_mTrigId;   //!
   TBranch        *b_mNTrks;   //!
   TBranch        *b_high_pt;   //!
   TBranch        *b_high_phi;   //!
   TBranch        *b_high_eta;   //!
   TBranch        *b_NLambda;   //!
   TBranch        *b_p1_InEventID;   //!
   TBranch        *b_p1_pt;   //!
   TBranch        *b_p1_phi;   //!
   TBranch        *b_p1_eta;   //!
   TBranch        *b_p1_dca;   //!
   TBranch        *b_p1_ch;   //!
   TBranch        *b_p1_hasTOFinfo;   //!
   TBranch        *b_p1_dedx;   //!
   TBranch        *b_p1_beta;   //!
   TBranch        *b_p2_InEventID;   //!
   TBranch        *b_p2_pt;   //!
   TBranch        *b_p2_phi;   //!
   TBranch        *b_p2_eta;   //!
   TBranch        *b_p2_dca;   //!
   TBranch        *b_p2_ch;   //!
   TBranch        *b_p2_hasTOFinfo;   //!
   TBranch        *b_p2_dedx;   //!
   TBranch        *b_p2_beta;   //!
   TBranch        *b_pair_charge;   //!
   TBranch        *b_pair_DCAdaughters;   //!
   TBranch        *b_pair_theta;   //!
   TBranch        *b_pair_decayL;   //!
   TBranch        *b_pair_phi;   //!
   TBranch        *b_pair_eta;   //!
   TBranch        *b_pair_pt;   //!
   TBranch        *b_pair_mass;   //!

   ntp_Lambda_Reader(std::vector<std::string> Infiles);
   virtual ~ntp_Lambda_Reader();
   virtual void     Init(TTree *tree);
};




ntp_Lambda_Reader::ntp_Lambda_Reader(std::vector<std::string> Infiles) 
   : fChain(0), InputFiles(Infiles)
{

}

ntp_Lambda_Reader::~ntp_Lambda_Reader()
{
   
}



void ntp_Lambda_Reader::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventId", &eventId, &b_eventId);
   fChain->SetBranchAddress("Vz", &Vz, &b_Vz);
   fChain->SetBranchAddress("mNTrigs", &mNTrigs, &b_mNTrigs);
   fChain->SetBranchAddress("mTrigId", mTrigId, &b_mTrigId);
   fChain->SetBranchAddress("mNTrks", &mNTrks, &b_mNTrks);
   fChain->SetBranchAddress("high_pt", high_pt, &b_high_pt);
   fChain->SetBranchAddress("high_phi", high_phi, &b_high_phi);
   fChain->SetBranchAddress("high_eta", high_eta, &b_high_eta);
   fChain->SetBranchAddress("NLambda", &NLambda, &b_NLambda);
   fChain->SetBranchAddress("p1_InEventID", p1_InEventID, &b_p1_InEventID);
   fChain->SetBranchAddress("p1_pt", p1_pt, &b_p1_pt);
   fChain->SetBranchAddress("p1_phi", p1_phi, &b_p1_phi);
   fChain->SetBranchAddress("p1_eta", p1_eta, &b_p1_eta);
   fChain->SetBranchAddress("p1_dca", p1_dca, &b_p1_dca);
   fChain->SetBranchAddress("p1_ch", p1_ch, &b_p1_ch);
   fChain->SetBranchAddress("p1_hasTOFinfo", p1_hasTOFinfo, &b_p1_hasTOFinfo);
   fChain->SetBranchAddress("p1_dedx", p1_dedx, &b_p1_dedx);
   fChain->SetBranchAddress("p1_beta", p1_beta, &b_p1_beta);
   fChain->SetBranchAddress("p2_InEventID", p2_InEventID, &b_p2_InEventID);
   fChain->SetBranchAddress("p2_pt", p2_pt, &b_p2_pt);
   fChain->SetBranchAddress("p2_phi", p2_phi, &b_p2_phi);
   fChain->SetBranchAddress("p2_eta", p2_eta, &b_p2_eta);
   fChain->SetBranchAddress("p2_dca", p2_dca, &b_p2_dca);
   fChain->SetBranchAddress("p2_ch", p2_ch, &b_p2_ch);
   fChain->SetBranchAddress("p2_hasTOFinfo", p2_hasTOFinfo, &b_p2_hasTOFinfo);
   fChain->SetBranchAddress("p2_dedx", p2_dedx, &b_p2_dedx);
   fChain->SetBranchAddress("p2_beta", p2_beta, &b_p2_beta);
   fChain->SetBranchAddress("pair_charge", pair_charge, &b_pair_charge);
   fChain->SetBranchAddress("pair_DCAdaughters", pair_DCAdaughters, &b_pair_DCAdaughters);
   fChain->SetBranchAddress("pair_theta", pair_theta, &b_pair_theta);
   fChain->SetBranchAddress("pair_decayL", pair_decayL, &b_pair_decayL);
   fChain->SetBranchAddress("pair_phi", pair_phi, &b_pair_phi);
   fChain->SetBranchAddress("pair_eta", pair_eta, &b_pair_eta);
   fChain->SetBranchAddress("pair_pt", pair_pt, &b_pair_pt);
   fChain->SetBranchAddress("pair_mass", pair_mass, &b_pair_mass);
   
}

#endif



