//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Mar  6 21:58:59 2026 by ROOT version 6.24/06
// from TTree ntp_Lambda/Lambda TTree
// found on file: output_1.root
//////////////////////////////////////////////////////////

#ifndef ntp_Lambda_h
#define ntp_Lambda_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class ntp_Lambda {
public :
   std::vector<std::string> InPutFileList;
   std::string     OutPutFile;

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           eventId;
   Float_t         Vz;
   Float_t         lead_pt;
   Float_t         lead_phi;
   Float_t         lead_eta;
   Float_t         sublead_pt;
   Float_t         sublead_phi;
   Float_t         sublead_eta;
   Int_t           p1_InEventID;
   Float_t         p1_pt;
   Float_t         p1_phi;
   Float_t         p1_eta;
   Float_t         p1_dca;
   Int_t           p1_ch;
   Int_t           p1_hasTOFinfo;
   Int_t           p2_InEventID;
   Float_t         p2_pt;
   Float_t         p2_phi;
   Float_t         p2_eta;
   Float_t         p2_dca;
   Int_t           p2_hasTOFinfo;
   Int_t           pair_charge;
   Float_t         pair_DCAdaughters;
   Float_t         pair_theta;
   Float_t         pair_decayL;
   Float_t         pair_phi;
   Float_t         pair_eta;
   Float_t         pair_pt;
   Float_t         pair_mass;

   // List of branches
   TBranch        *b_eventId;   //!
   TBranch        *b_Vz;   //!
   TBranch        *b_lead_pt;   //!
   TBranch        *b_lead_phi;   //!
   TBranch        *b_lead_eta;   //!
   TBranch        *b_sublead_pt;   //!
   TBranch        *b_sublead_phi;   //!
   TBranch        *b_sublead_eta;   //!
   TBranch        *b_p1_InEventID;   //!
   TBranch        *b_p1_pt;   //!
   TBranch        *b_p1_phi;   //!
   TBranch        *b_p1_eta;   //!
   TBranch        *b_p1_dca;   //!
   TBranch        *b_p1_ch;   //!
   TBranch        *b_p1_hasTOFinfo;   //!
   TBranch        *b_p2_InEventID;   //!
   TBranch        *b_p2_pt;   //!
   TBranch        *b_p2_phi;   //!
   TBranch        *b_p2_eta;   //!
   TBranch        *b_p2_dca;   //!
   TBranch        *b_p2_hasTOFinfo;   //!
   TBranch        *b_pair_charge;   //!
   TBranch        *b_pair_DCAdaughters;   //!
   TBranch        *b_pair_theta;   //!
   TBranch        *b_pair_decayL;   //!
   TBranch        *b_pair_phi;   //!
   TBranch        *b_pair_eta;   //!
   TBranch        *b_pair_pt;   //!
   TBranch        *b_pair_mass;   //!

   ntp_Lambda(TTree *tree=0);
   virtual ~ntp_Lambda();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ntp_Lambda_cxx
ntp_Lambda::ntp_Lambda() : fChain(0) 
{

}

ntp_Lambda::~ntp_Lambda()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ntp_Lambda::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ntp_Lambda::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ntp_Lambda::Init(TTree *tree)
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
   fChain->SetBranchAddress("lead_pt", &lead_pt, &b_lead_pt);
   fChain->SetBranchAddress("lead_phi", &lead_phi, &b_lead_phi);
   fChain->SetBranchAddress("lead_eta", &lead_eta, &b_lead_eta);
   fChain->SetBranchAddress("sublead_pt", &sublead_pt, &b_sublead_pt);
   fChain->SetBranchAddress("sublead_phi", &sublead_phi, &b_sublead_phi);
   fChain->SetBranchAddress("sublead_eta", &sublead_eta, &b_sublead_eta);
   fChain->SetBranchAddress("p1_InEventID", &p1_InEventID, &b_p1_InEventID);
   fChain->SetBranchAddress("p1_pt", &p1_pt, &b_p1_pt);
   fChain->SetBranchAddress("p1_phi", &p1_phi, &b_p1_phi);
   fChain->SetBranchAddress("p1_eta", &p1_eta, &b_p1_eta);
   fChain->SetBranchAddress("p1_dca", &p1_dca, &b_p1_dca);
   fChain->SetBranchAddress("p1_ch", &p1_ch, &b_p1_ch);
   fChain->SetBranchAddress("p1_hasTOFinfo", &p1_hasTOFinfo, &b_p1_hasTOFinfo);
   fChain->SetBranchAddress("p2_InEventID", &p2_InEventID, &b_p2_InEventID);
   fChain->SetBranchAddress("p2_pt", &p2_pt, &b_p2_pt);
   fChain->SetBranchAddress("p2_phi", &p2_phi, &b_p2_phi);
   fChain->SetBranchAddress("p2_eta", &p2_eta, &b_p2_eta);
   fChain->SetBranchAddress("p2_dca", &p2_dca, &b_p2_dca);
   fChain->SetBranchAddress("p2_hasTOFinfo", &p2_hasTOFinfo, &b_p2_hasTOFinfo);
   fChain->SetBranchAddress("pair_charge", &pair_charge, &b_pair_charge);
   fChain->SetBranchAddress("pair_DCAdaughters", &pair_DCAdaughters, &b_pair_DCAdaughters);
   fChain->SetBranchAddress("pair_theta", &pair_theta, &b_pair_theta);
   fChain->SetBranchAddress("pair_decayL", &pair_decayL, &b_pair_decayL);
   fChain->SetBranchAddress("pair_phi", &pair_phi, &b_pair_phi);
   fChain->SetBranchAddress("pair_eta", &pair_eta, &b_pair_eta);
   fChain->SetBranchAddress("pair_pt", &pair_pt, &b_pair_pt);
   fChain->SetBranchAddress("pair_mass", &pair_mass, &b_pair_mass);
   Notify();
}

Bool_t ntp_Lambda::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ntp_Lambda::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ntp_Lambda::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ntp_Lambda_cxx
