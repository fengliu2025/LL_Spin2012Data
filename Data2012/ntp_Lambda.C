#define ntp_Lambda_cxx
#include "ntp_Lambda.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "new_Tree.h"
void ntp_Lambda::Loop()
{
   //------------------------------------------ENTER FILE LOOP-----------------------------------
   for(int ifile = 0; ifile < InPutFileList.size(); ifile ++){
         //Open the .root file 
        TFile *fin = TFile::Open(InPutFileList[ifile].c_str(),"READ"); 
	     if(!fin){
            std::cout<<"Can not opne the file:"<<InPutFileList[ifile]<<std::endl;
            std::cout<<"Skip this root file"<<std::endl;
            continue;
         }

         TTree *Tree = (TTree *)fin->Get("ntp_Lambda"); 
         if(!Tree){
            std::cout<<"Can not open the TTree,skip this file"<<std::endl;
            continue;
         }
         //open the output
	      TFile *fout = new TFile(OutPutFileList[ifile].c_str(),"RECREATE");
         TTree *outputTree = new TTree("ntp_Lambda","ntp_Lambda");

         //Initialize the TTree 
         Init(Tree);
         int current_eventId = -1;
         int current_Nlambda = 0;
	     

         new_Tree *newTree = new new_Tree();
        
	
	      newTree->MakeNewTree(outputTree);
         
	      
	      newTree->ResetTree();
         Long64_t NEntries =  fChain->GetEntries();
         fChain->GetEntry(0);
         //-------------------------------------ENTER ENTRY LOOP--------------------------------
         for(int ientry=0; ientry < NEntries; ientry++ ){
          if(ientry%100000==0) cout<<"ientry="<<ientry<<std::endl;   
	
            newTree->p1_InEventID[current_Nlambda] = p1_InEventID;
            newTree->p1_pt[current_Nlambda]        = p1_pt; 
            newTree->p1_phi[current_Nlambda]       = p1_phi;
            newTree->p1_eta[current_Nlambda]       = p1_eta; 
            newTree->p1_dca[current_Nlambda]       = p1_dca;
            newTree->p1_ch[current_Nlambda]        = p1_ch;
            newTree->p1_hasTOFinfo[current_Nlambda]= p1_hasTOFinfo;

            newTree->p2_InEventID[current_Nlambda] = p2_InEventID;
            newTree->p2_pt[current_Nlambda]        = p2_pt;
            newTree->p2_phi[current_Nlambda]       = p2_phi;
            newTree->p2_eta[current_Nlambda]       = p2_eta;
            newTree->p2_dca[current_Nlambda]       = p2_dca;
            newTree->p2_hasTOFinfo[current_Nlambda]= p2_hasTOFinfo;

            newTree->pair_charge[current_Nlambda]  = pair_charge;
            newTree->pair_DCAdaughters[current_Nlambda]=pair_DCAdaughters;
            newTree->pair_theta[current_Nlambda]   = pair_theta;
            newTree->pair_decayL[current_Nlambda]  = pair_decayL;
            newTree->pair_phi[current_Nlambda]     = pair_phi;
            newTree->pair_eta[current_Nlambda]     = pair_eta;
            newTree->pair_pt[current_Nlambda]      = pair_pt;
            newTree->pair_mass[current_Nlambda]    = pair_mass;
       
  
            current_Nlambda++;
            if(ientry==NEntries-1){
               newTree->NLambda = current_Nlambda;
               newTree->Fill();
               break;
            }
            fChain->GetEntry(ientry+1);
            if(eventId!=current_eventId){
               newTree->NLambda = current_Nlambda;
               newTree->Fill();

               //prepare for next events 
               current_eventId = eventId;
               current_Nlambda =0;

               newTree->eventId = eventId;
               newTree->Vz      = Vz;
            }


         }
         //-------------------------------------END ENTRY LOOP--------------------------------

         //Write to OutPutFile 
         outputTree->Write();
         fout->Close(); 
	      fin->Close();
         delete fin;

   }
   //------------------------------------------END FILE LOOP  -----------------------------------
}
