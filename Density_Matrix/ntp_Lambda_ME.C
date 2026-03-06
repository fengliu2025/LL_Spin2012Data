#define ntp_Lambda_ME_cxx
#include "ntp_Lambda_ME.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>



void ntp_Lambda_ME::Rotation(TVector3 a_hat,TVector3 b_hat, TVector c_hat, TLorentzVector *v){
   double x1 = v->Vect().Dot(a_hat); // x1 coordiante in the helicity frame a-b-c
   double x2 = v->Vect().Dot(b_hat);
   double x3 = v->Vect().Dot(c_hat);

   v->SetX(x1);
   v->SetY(x2);
   v->SetZ(x3);

}



void ntp_Lambda_ME::Loop()
{
   //-----------------------------------------------------Define Histograms----------------------------------------------------
    TH1D *h1d_B_a_1 = new TH1D("h1d_B_a_1","h1d_B_a_1");

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
         //Only use two Lambda events
         if(NLambda!=2) continue;
         if( pair_ch[0] == 1 || pair_ch[1] == 1 ) continue; //if any one of the two Lambda candidates is combina background, skip this event
         if( p1_ch[0] * p1_ch[1] != -1 ) continue; 
         if (!( (p1_ch[0] ==-1 && p1_ch[1] == +1) || (p1_ch[0] == +1 && p1_ch[1]==-1 ) ) ){
            std::cout<<"error with p1_ch[0] * p1_ch[1] !=-1" << std::endl;
            continue; 
         }
         //make some cut on the pion , proton and pair

         //construct the Lambda 
         TLorentzVector *Lambda     = new TLorentzVector ;
         TLorentzVector *proton     = new TLorentzVector ;
         TLorentzVector *pion_minus = new TLorentzVector ; 
         //construct the Lambda bar 
         TLorentzVector *Lambda_bar = new TLorentzVector ; 
         TLorentzVector *proton_bar = new TLorentzVector ;
         TLorentzVector *pion_plus  = new TLorentzVector ;

         int id_Lambda = 0 ; 
         int id_Lambda_bar = 1;

         if(p1_ch[0] == -1 && p1_ch[1] == +1  ){
            //The pair[1] is Lambda and pair[0] is Lambda_bar
            id_Lambda = 1; 
            id_Lambda_bar = 0;
         }

         Lambda->SetPtEtaPhiM(     pair_pt[id_Lambda], pair_eta[id_Lambda], pair_phi[id_Lambda], pair_mass[id_Lambda] ); 
         proton->SetPtEtaPhiM(     p1_pt[id_Lambda]  , p1_eta[id_Lambda]  , p1_phi[id_Lambda]  , PROTON_MASS          );
         pion_minus->SetPtEtaPhiM( p2_pt[id_Lambda]  , p2_eta[id_Lambda]  , p2_phi[id_Lambda]  , PION_MASS            );

         Lambda_bar->SetPtEtaPhiM( pair_pt[id_Lambda_bar], pair_eta[id_Lambda_bar], pair_phi[id_Lambda_bar], pair_mass[id_Lambda_bar]  );
         proton_bar->SetPtEtaPhiM( p1_pt[id_Lambda_bar]  , p1_eta[id_Lambda_bar]  ,  p1_phi[id_Lambda_bar] ,  PROTON_MASS              );
         pion_plus->SetPtEtaPhiM(  p2_pt[id_Lambda_bar]  , p2_eta[id_Lambda_bar]  ,  p2_phi[id_Lambda_bar] ,  PION_MASS                );

         TLorentzVector LL_bar = (*Lambda) + (*Lambda_bar);
         // construc the beta vector of LL_bar
         TVector3 beta_LL_bar = LL_bar.BoostVector();
         // Boost everything into the rest frame of LL_bar 
         LL_bar.Boost(  -beta_LL_bar      );

         Lambda->Boost( -beta_LL_bar      );
         proton->Boost( -beta_LL_bar      );
         pion_minus->Boost( -beta_LL_bar );

         Lambda_bar->Boost( -beta_LL_bar  );
         proton_bar->Boost( -beta_LL_bar  );
         pion_plus->Boost(  -beta_LL_bar  );
         // Rotate into the helicity frame

         //construct the helicity coordinate axis 
         TVector3 p_beam = (0,0,1);

         TVector3 c_hat = Lambda->Vect()/(Lambda->Vect().Mag());
         TVector3 b_hat = p_beam->Cross(c_hat).Unit();
         TVector3 a_hat = b_hat.Cross(c_hat);
         //Rotate Everything 
         Rotation(a_hat,b_hat,c_hat,Lambda);
         Rotation(a_hat,b_hat,c_hat,proton);
         Rotation(a_hat,b_hat,c_hat,pion_minus);

         Rotation(a_hat,b_hat,c_hat,Lambda_bar);
         Rotation(a_hat,b_hat,c_hat,proton_bar);
         Rotation(a_hat,b_hat,c_hat,pion_plus );

         // consruct the beta vector of Lambda and Lambda_bar 
         TVector3 beta_Lambda = Lambda.BoostVector();
         TVector3 beta_Lambda_bar = Lambda_bar.BoostVector();

         // Boost the daughter particles into the mother's rest frame, 
         proton->Boost(-beta_Lambda);
         pion_minus->Boost(-beta_Lambda);

         proton_bar->Boost(-beta_Lambda);
         pion_plus->Boost(-beta_Lambda);



         // calculate the theta_1(2), phi_1(2) 
         double theta_1 = proton.Theta();
         double phi_1   = proton.Phi();

         double theta_2 = proton_bar.Theta();
         double phi_2   = proton_bar.Phi();

         double B_a_1 = TMath::Sin(theta_1) * TMath::Cos(phi_1);
         double B_a_2 = TMath::Sin(theta_1) * TMath::Sin(phi_1);
         double B_a_3 = TMath::Cos(theta_1); 
         double B_b_1 = TMath::Sin(theta_2) * TMath::Cos(phi_2);
         double B_b_2 = TMath::Sin(theta_2) * TMath::Sin(phi_2);
         double B_b_3 = TMath::Cos(theta_2);

         double C_11 = TMath::Sin(theta_1) * TMath::Cos(phi_1) * TMath::Sin(theta_2) * TMath::Cos(phi_2);
         double C_12 = TMath::Sin(theta_1) * TMath::Cos(phi_1) * TMath::Sin(theta_2) * TMath::Sin(phi_2);
         double C_13 = TMath::Sin(theta_1) * TMath::Cos(phi_1) * TMath::Cos(theta_2);
         
         double C_21 = TMath::Sin(theta_1) * TMath::Sin(phi_1) * TMath::Sin(theta_2) * TMath::Cos(phi_2);
         double C_22 = TMath::Sin(theta_1) * TMath::Sin(phi_1) * TMath::Sin(theta_2) * TMath::Sin(phi_2);
         double C_23 = TMath::Sin(theta_1) * TMath::Sin(phi_1) * TMath::Cos(theta_2);
         
         double C_31 = TMath::Cos(theta_1) * TMath::Sin(theta_2) * TMath::Cos(phi_2);
         double C_32 = TMath::Cos(theta_1) * TMath::Sin(theta_2) * TMath::Sin(phi_2);
         double C_33 = TMath::Cos(theta_1) * TMath::Cos(theta_2);

         double cos_theta_star = TMath::Sin(theta_1) * TMath::Cos(phi_1) * TMath::Sin(theta_2) * TMath::Cos(phi_2) + TMath::Sin(theta_1) * TMath::Sin(phi_1) * TMath::Sin(theta_2) * TMath::Sin(phi_2) + TMath::Cos(theta_1) * TMath::Cos(theta_2);    

         // fill the histgram  
         h1D_B_a_1 -> Fill(B_a_1);
         h1D_B_a_2 -> Fill(B_a_2);
         h1D_B_a_3 -> Fill(B_a_3);

         h1D_B_b_1 -> Fill(B_b_1);
         h1D_B_b_2 -> Fill(B_b_2);
         h1D_B_b_3 -> Fill(B_b_3);

         h1D_C_11 -> Fill(C_11);
         h1D_C_12 -> Fill(C_12);
         h1D_C_13 -> Fill(C_13);
         h1D_C_21 -> Fill(C_21);
         h1D_C_22 -> Fill(C_22);
         h1D_C_23 -> Fill(C_23);
         h1D_C_31 -> Fill(C_31);
         h1D_C_32 -> Fill(C_32);
         h1D_C_33 -> Fill(C_33);

         h1D_cos_theta_star -> Fill(cos_theta_star);
         








         

      }
      //******************************************END i_Event Loop*******************************

      fin->Close();
   }
   //*********************************************END i_file Loop*********************************


   TFile *fout = TFile::Open(OutputFile.c_str(),"RECREATE");
  

   fout->Close();





}
