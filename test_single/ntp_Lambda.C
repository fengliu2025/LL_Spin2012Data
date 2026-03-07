#define ntp_Lambda_cxx
#include "ntp_Lambda.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void ntp_Lambda::Loop()
{

   TH1D *h1D_Vz             = new TH1D("h1D_Vz","h1D_Vz",600,-20,100);
   TH1D *h1D_p1_pt          = new TH1D("h1D_p1_pt","h1D_p1_pt",100,-0.2,3);
   TH1D *h1D_p1_phi         = new TH1D("h1D_p1_phi","h1D_p1_phi",100,-2 *TMath::Pi(),2*TMath::Pi());
   TH1D *h1D_p1_eta         = new TH1D("h1D_p1_eta","h1D_p1_eta",100,-2,2);
   TH1D *h1D_p1_dca         = new TH1D("h1D_p1_dca","h1D_p1_dca",100,-1,10);
   TH1D *h1D_p1_ch          = new TH1D("h1D_p1_ch","h1D_p1_ch",3,-1.5,1.5);
   TH1D *h1D_p1_hasTOFinfo  = new TH1D("h1D_p1_hasTOFinfo","h1D_p1_hasTOFinfo",3,-1.5,1.5);
   TH1D *h1D_p2_pt          = new TH1D("h1D_p2_pt","h1D_p2_pt",100,-0.2,3);
   TH1D *h1D_p2_phi         = new TH1D("h1D_p2_phi","h1D_p2_phi",100,-2 *TMath::Pi(),2*TMath::Pi());
   TH1D *h1D_p2_eta         = new TH1D("h1D_p2_eta","h1D_p2_eta",100,-2,2);
   TH1D *h1D_p2_dca         = new TH1D("h1D_p2_dca","h1D_p2_dca",100,-1,10);
   TH1D *h1D_p2_hasTOFinfo  = new TH1D("h1D_p2_hasTOFinfo","h1D_p2_hasTOFinfo",3,-1.5,1.5);
   TH1D *h1D_pair_charge    = new TH1D("h1D_pair_charge","h1D_pair_charge",3,-1.5,1.5);
   TH1D *h1D_pair_DCAdaughters  = new TH1D("h1D_pair_DCAdaughters","h1D_pair_DCAdaughters",100,-1,10);
   TH1D *h1D_pair_theta     = new TH1D("h1D_pair_theta","h1D_pair_theta",100,-TMath::Pi(),TMath::Pi());
   TH1D *h1D_pair_decayL    = new TH1D("h1D_pair_decayL","h1D_pair_decayL",100,-1,30);
   TH1D *h1D_pair_phi       = new TH1D("h1D_pair_phi","h1D_pair_phi",100,-2 *TMath::Pi(),2*TMath::Pi());
   TH1D *h1D_pair_eta       = new TH1D("h1D_pair_eta","h1D_pair_eta",100,-2,2);
   TH1D *h1D_pair_pt        = new TH1D("h1D_pair_pt","h1D_pair_pt",100,-1,10);
   TH1D *h1D_pair_mass      = new TH1D("h1D_pair_mass","h1D_pair_mass",100,0.5,1.5);

   //------------------------------------iFile---------------------------
   for(int iFile = 0 ; iFile < InPutFileList.size(); iFile++){
      std::cout<<"iFile"<<iFile<<std::endl;
      // open the file 
      TFile *fin= TFile::Open( InPutFileList[iFile].c_str() );
      if(!fin){
         std::cout<<"Can not Open the File"<<std::endl;
         continue;
      }
      TTree *tree = (TTree *)fin->Get("ntp_Lambda");
      if(!tree){
         std::cout<<"Can not get the tree"<<std::endl;
         continue
      }
      Init(tree);

      Long64_t Nentries = fChain->GetEntriesFast();
      //------------------------------------iEntry---------------------------
      for(int iEntry = 0 ; iEntry < Nentries ; iEntry++){
         fChain->GetEntry(iEntry);

         h1D_Vz                   ->Fill(Vz);  
         h1D_p1_pt                ->Fill(p1_pt);  
         h1D_p1_phi               ->Fill(p1_phi );  
         h1D_p1_eta               ->Fill(p1_eta );  
         h1D_p1_dca               ->Fill(p1_dca );  
         h1D_p1_ch                ->Fill(p1_ch );  
         h1D_p1_hasTOFinfo        ->Fill(p1_hasTOFinfo );  
         h1D_p2_pt                ->Fill(p2_pt);  
         h1D_p2_phi               ->Fill(p2_phi  );  
         h1D_p2_eta               ->Fill(p2_eta );  
         h1D_p2_dca               ->Fill(p2_dca );  
         h1D_p2_hasTOFinfo        ->Fill(p2_hasTOFinfo);  
         h1D_pair_charge          ->Fill(pair_charge );  
         h1D_pair_DCAdaughters    ->Fill(pair_DCAdaughters);  
         h1D_pair_theta           ->Fill(pair_theta );  
         h1D_pair_decayL          ->Fill(pair_decayL);  
         h1D_pair_phi             ->Fill(pair_phi);  
         h1D_pair_eta             ->Fill(pair_eta);  
         h1D_pair_pt              ->Fill(pair_pt);  
         h1D_pair_mass            ->Fill(pair_mass );  

      }
      //------------------------------------iEntry---------------------------
   
      fin->Close();

   }
   //------------------------------------iFile---------------------------

   TFile *fout = TFile::Open(OutPutFile.c_str(),"RECREATE");
   h1D_Vz                   ->Write();
   h1D_p1_pt                ->Write();
   h1D_p1_phi               ->Write(); 
   h1D_p1_eta               ->Write();
   h1D_p1_dca               ->Write(); 
   h1D_p1_ch                ->Write();
   h1D_p1_hasTOFinfo        ->Write();
   h1D_p2_pt                ->Write();
   h1D_p2_phi               ->Write(); 
   h1D_p2_eta               ->Write(); 
   h1D_p2_dca               ->Write();
   h1D_p2_hasTOFinfo        ->Write();
   h1D_pair_charge          ->Write();
   h1D_pair_DCAdaughters    ->Write();
   h1D_pair_theta           ->Write();
   h1D_pair_decayL          ->Write();  
   h1D_pair_phi             ->Write();
   h1D_pair_eta             ->Write();  
   h1D_pair_pt              ->Write(); 
   h1D_pair_mass            ->Write();

   fout->Close();

   delete fout;
 

   
}
