#define ntp_Lambda_cxx
#include "ntp_Lambda.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void ntp_Lambda::Loop()
{
   //-----------------------------------------------------Define Histograms----------------------------------------------------
   TH1D * h1D_NLambda = new TH1D("h1D_NLambda","h1D_NLambda" ,5,0.5,5.5);
   //daughter particle 1
   TH1D * h1D_p1_pt   = new TH1D("h1D_p1_pt"  ,"h1D_p1_pt"   ,200,0,20 );
   TH1D * h1D_p1_phi  = new TH1D("h1D_p1_phi" ,"h1D_p1_phi"  ,200,-2*TMath::Pi(),2*TMath::Pi() );
   TH1D * h1D_p1_eta  = new TH1D("h1D_p1_eta" ,"h1D_p1_eta"  ,100,-5,5);
   //daughter particle 2 
   TH1D * h1D_p2_pt   = new TH1D("h1D_p2_pt"  ,"h1D_p2_pt"   ,200,0,20 );
   TH1D * h1D_p2_phi  = new TH1D("h1D_p2_phi" ,"h1D_p2_phi"  ,200,-2*TMath::Pi(),2*TMath::Pi());
   TH1D * h1D_p2_eta  = new TH1D("h1D_p2_eta" ,"h1D_p2_eta"  ,100,-5,5);

   //pair 
   TH1D *h1D_pair_phi  = new TH1D("h1D_pair_phi","h1D_pair_phi",100 , -2*TMath::Pi() , 2*TMath::Pi() );
   TH1D *h1D_pair_eta  = new TH1D("h1D_pair_eta","h1D_pair_eta",100 , -5, 5  );
   TH1D *h1D_pair_pt   = new TH1D("h1D_pair_pt" ,"h1D_pair_pt" ,100 ,  0, 10 );
   TH1D *h1D_pair_mass = new TH1D("h1D_pair_mass","h1D_pair_mass",100,1.0,1.3); 

   //-----------------------------------------------------Define Histograms----------------------------------------------------


   //*********************************************ENTER i_file Loop*********************************
   for(long int i_file=0 ; i_file < InputFiles.size() ; i_file++){
      if(i_file%100==0){
         std::cout << "i_file" << i_file << std::endl;
      }
      //Initialzie the tree
      //Open the file 
      TFile *fin =TFile::Open(InputFiles[i_file].c_str());
      if (fin==0){
         std::cout<<"Can't open file:"<<InputFiles[i_file]<<". Skip this file."<<std::endl;
         continue;
      }
      TTree *tree = (TTree *)fin->Get("ntp_Lambda");
      if(tree ==0){
         std::cout<<"Can't get the tree. Skip this file" <<std::endl;
         continue;
      }
      Init(tree);
      //Number of entries in this tree 
      long int N_Entry= fChain->GetEntries();

      //******************************************ENTER i_Event Loop*******************************
      for(long int i_Event = 0; i_Event < N_Entry ; i_Event++){
         fChain->GetEntry(i_Event);
         h1D_NLambda->Fill(NLambda);
         for(int i_lambda = 0; i_lambda < NLambda ; i_lambda++ ){
            h1D_p1_pt   ->Fill(p1_pt[i_lambda]) ;
            h1D_p1_eta  ->Fill(p1_eta[i_lambda]);
            h1D_p1_phi  ->Fill(p1_phi[i_lambda]); 

            h1D_p2_pt   ->Fill(p2_pt[i_lambda]);
            h1D_p2_eta  ->Fill(p2_eta[i_lambda]);
            h1D_p2_phi  ->Fill(p2_phi[i_lambda]);

            h1D_pair_phi ->Fill(pair_phi[i_lambda] );
            h1D_pair_pt  ->Fill(pair_pt[i_lambda]  );
            h1D_pair_eta ->Fill(pair_eta[i_lambda] );
            h1D_pair_mass->Fill(pair_mass[i_lambda]);
         }

      }
      //******************************************END i_Event Loop*******************************

      fin->Close();
   }
   //*********************************************END i_file Loop*********************************


   TFile *fout = TFile::Open(OutputFile.c_str(),"RECREATE");
   h1D_NLambda ->Write();
   h1D_p1_pt   ->Write();
   h1D_p1_phi  ->Write();
   h1D_p1_eta  ->Write();
  
   h1D_p2_pt   ->Write();
   h1D_p2_phi  ->Write();
   h1D_p2_eta  ->Write();

   //pair 
   h1D_pair_phi  ->Write();
   h1D_pair_eta  ->Write();
   h1D_pair_pt   ->Write();
   h1D_pair_mass ->Write();

   fout->Close();





}
