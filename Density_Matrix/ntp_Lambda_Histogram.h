#ifndef ntp_Lambda_Histogram_h
#define ntp_Lambda_Histogram_h

//include the header file here 
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TMath.h>
#include <TLorentzVector.h>


#include <string>
#include <iostream>


#include "ntp_Lambda_Reader.h"
#include "ntp_Lambda_Calculator.h"

#include "constants.h"
class ntp_Lambda_Histogram {
public : 
	ntp_Lambda_Reader *Reader;
	ntp_Lambda_Calculator *Calculator;

	std::string OutPutFile;
	//Histograms for QA plot
	//QA plots. Basics Distribution of daughter particles.(After Selection).
	TH1D *h1D_p1_Pt;
	TH1D *h1D_p1_Eta;
	TH1D *h1D_p1_Phi;

	TH1D *h1D_p2_Pt;
	TH1D *h1D_p2_Eta;
	TH1D *h1D_p2_Phi;
	//QA plots. Basics Distribuion of Lambda Candidates.(After Selection).
	TH1D *h1D_pair_Pt;
	TH1D *h1D_pair_Eta;
	TH1D *h1D_pair_Phi;
	TH1D *h1D_pair_Mass;

	//Histograms for Pair Plots
	TH2D *h2D_L1Mass_L2Mass[Range_Bin][3];
	TH2D *h2D_L1Pt_L2Pt[Range_Bin][3];
	TH2D *h2D_L1Phi_L2Phi[Range_Bin][3];
	TH2D *h2D_L1Eta_L2Eta[Range_Bin][3];

	TH2D *h2D_PairMass_deltaR[3];
	TH1D *h1D_L1L2_PtDifference[3];

	//Histograms for density matrix 
	TH1D *h1D_B_a_1[Range_Bin][3];
	TH1D *h1D_B_a_2[Range_Bin][3];
	TH1D *h1D_B_a_3[Range_Bin][3];

	TH1D *h1D_B_b_1[Range_Bin][3];
	TH1D *h1D_B_b_2[Range_Bin][3];
	TH1D *h1D_B_b_3[Range_Bin][3];

	TH1D *h1D_C_11[Range_Bin][3];
	TH1D *h1D_C_12[Range_Bin][3];
	TH1D *h1D_C_13[Range_Bin][3];

	TH1D *h1D_C_21[Range_Bin][3];
	TH1D *h1D_C_22[Range_Bin][3];
	TH1D *h1D_C_23[Range_Bin][3];


	TH1D *h1D_C_31[Range_Bin][3];
	TH1D *h1D_C_32[Range_Bin][3];
	TH1D *h1D_C_33[Range_Bin][3];

	TH1D *h1D_cos_theta_star_Lab[Range_Bin][3];
	TH1D *h1D_cos_theta_star_PairRest[Range_Bin][3];

	//Trigger 
	TH1D *h1D_Trigger;
	TH1D *h1D_MB_Trigger;
	TH1D *h1D_HM_Trigger;

	TH1D *h1D_MB_HM_Trigger;

	ntp_Lambda_Histogram();
	ntp_Lambda_Histogram(ntp_Lambda_Reader *reader,ntp_Lambda_Calculator *calculator,std::string outPutFile);
	void InitHitogram();
	void Fill_QAplots(std::vector<int> GoodLambdaFlag);
	void Fill_PairPlots(TLorentzVector *v1, TLorentzVector*v2,int type_bin , int range_bin);
	void Fill_PairPlots(TLorentzVector *v1, TLorentzVector*v2,int type_bin , int range_bin,double weight);
	void Fill_DensityMatrix(int type_bin , int range_bin);
	void Fill_DensityMatrix(int type_bin , int range_bin, double weight) ;
	void Reset();

	void WriteAll();

};

ntp_Lambda_Histogram::ntp_Lambda_Histogram()
{}

ntp_Lambda_Histogram::ntp_Lambda_Histogram(ntp_Lambda_Reader *reader, ntp_Lambda_Calculator *calculator,std::string outPutFile)
:Reader(reader),Calculator(calculator),OutPutFile(outPutFile)
{

}

void ntp_Lambda_Histogram::InitHitogram(){
	h1D_p1_Pt      = new TH1D("h1D_p1_Pt","h1D_p1_Pt",100,0,10);
	h1D_p1_Eta     = new TH1D("h1D_p1_Eta","h1D_p1_Eta",100,-3,3); 
	h1D_p1_Phi     = new TH1D("h1D_p1_Phi","h1D_p1_Phi",100,-2*TMath::Pi(),2*TMath::Pi() ); 

	h1D_p2_Pt      = new TH1D("h1D_p2_Pt","h1D_p2_Pt",100,0,10);
	h1D_p2_Eta     = new TH1D("h1D_p2_Eta","h1D_p2_Eta",100,-3,3);
	h1D_p2_Phi     = new TH1D("h1D_p2_Phi","h1D_p2_Phi",100,-2*TMath::Pi(),2*TMath::Pi() ); 
	
	h1D_pair_Pt    = new TH1D("h1D_pair_Pt","h1D_pair_Pt",100,0,10);
	h1D_pair_Eta   = new TH1D("h1D_pair_Eta","h1D_pair_Eta",100,-3,3);
	h1D_pair_Phi   = new TH1D("h1D_pair_Phi","h1D_pair_Phi",100,-2*TMath::Pi(),2*TMath::Pi());
	h1D_pair_Mass  = new TH1D("h1D_pair_Mass","h1D_pair_Mass",120,1.07,1.20);


	for(int i_type = 0 ; i_type < 3 ;i_type++){
		h2D_PairMass_deltaR[i_type]   = new TH2D(Form("h2D_PairMass_deltaR_%d",i_type),Form("h2D_PairMass_deltaR_%d",i_type),400,2,6,100,0,4);
		h1D_L1L2_PtDifference[i_type] = new TH1D(Form("h1D_L1L2_PtDifference_%d",i_type),Form("h1D_L1L2_PtDifference_%d",i_type),100,-5,5   );

		for(int i_ran = 0 ; i_ran < Range_Bin ; i_ran++){
			h2D_L1Mass_L2Mass[i_ran][i_type] = new TH2D(Form("h2D_L1Mass_L2Mass_%d_%d",i_ran,i_type),Form("h2D_L1Mass_L2Mass_%d_%d",i_ran,i_type),200,1.07,1.20,200,1.07,1.20);
			h2D_L1Pt_L2Pt[i_ran][i_type]     = new TH2D(Form("h2D_L1Pt_L2Pt_%d_%d",i_ran,i_type),Form("h2D_L1Pt_L2Pt_%d_%d",i_ran,i_type),100,0,10,100,0,10);
			h2D_L1Phi_L2Phi[i_ran][i_type]   = new TH2D(Form("h2D_L1Phi_L2Phi_%d_%d",i_ran,i_type),Form("h2D_L1Phi_L2Phi_%d_%d",i_ran,i_type),100,-2*TMath::Pi(),2*TMath::Pi(),100,-2*TMath::Pi(),2*TMath::Pi());
			h2D_L1Eta_L2Eta[i_ran][i_type]   = new TH2D(Form("h2D_L1Eta_L2Eta_%d_%d",i_ran,i_type),Form("h2D_L1Eta_L2Eta_%d_%d",i_ran,i_type),100,-3,3,100,-3,3);

	
			h1D_B_a_1[i_ran][i_type] = new TH1D(Form("h1D_B_a_1_%d_%d",i_ran,i_type),Form("h1D_B_a_1_%d_%d",i_ran,i_type),100,-1,1);
			h1D_B_a_2[i_ran][i_type] = new TH1D(Form("h1D_B_a_2_%d_%d",i_ran,i_type),Form("h1D_B_a_2_%d_%d",i_ran,i_type),100,-1,1);
			h1D_B_a_3[i_ran][i_type] = new TH1D(Form("h1D_B_a_3_%d_%d",i_ran,i_type),Form("h1D_B_a_3_%d_%d",i_ran,i_type),100,-1,1);

			h1D_B_b_1[i_ran][i_type] = new TH1D(Form("h1D_B_b_1_%d_%d",i_ran,i_type),Form("h1D_B_b_1_%d_%d",i_ran,i_type),100,-1,1);
			h1D_B_b_2[i_ran][i_type] = new TH1D(Form("h1D_B_b_2_%d_%d",i_ran,i_type),Form("h1D_B_b_2_%d_%d",i_ran,i_type),100,-1,1);
			h1D_B_b_3[i_ran][i_type] = new TH1D(Form("h1D_B_b_3_%d_%d",i_ran,i_type),Form("h1D_B_b_3_%d_%d",i_ran,i_type),100,-1,1);

			h1D_C_11[i_ran][i_type]  = new TH1D(Form("h1D_C_11_%d_%d",i_ran,i_type),Form("h1D_C_11_%d_%d",i_ran,i_type),100,-1,1);
			h1D_C_12[i_ran][i_type]  = new TH1D(Form("h1D_C_12_%d_%d",i_ran,i_type),Form("h1D_C_12_%d_%d",i_ran,i_type),100,-1,1);
			h1D_C_13[i_ran][i_type]  = new TH1D(Form("h1D_C_13_%d_%d",i_ran,i_type),Form("h1D_C_13_%d_%d",i_ran,i_type),100,-1,1);

			h1D_C_21[i_ran][i_type]  = new TH1D(Form("h1D_C_21_%d_%d",i_ran,i_type),Form("h1D_C_21_%d_%d",i_ran,i_type),100,-1,1);
			h1D_C_22[i_ran][i_type]  = new TH1D(Form("h1D_C_22_%d_%d",i_ran,i_type),Form("h1D_C_22_%d_%d",i_ran,i_type),100,-1,1);
			h1D_C_23[i_ran][i_type]  = new TH1D(Form("h1D_C_23_%d_%d",i_ran,i_type),Form("h1D_C_23_%d_%d",i_ran,i_type),100,-1,1);


			h1D_C_31[i_ran][i_type]  = new TH1D(Form("h1D_C_31_%d_%d",i_ran,i_type),Form("h1D_C_31_%d_%d",i_ran,i_type),100,-1,1);
			h1D_C_32[i_ran][i_type]  = new TH1D(Form("h1D_C_32_%d_%d",i_ran,i_type),Form("h1D_C_32_%d_%d",i_ran,i_type),100,-1,1);
			h1D_C_33[i_ran][i_type]  = new TH1D(Form("h1D_C_33_%d_%d",i_ran,i_type),Form("h1D_C_33_%d_%d",i_ran,i_type),100,-1,1);

			h1D_cos_theta_star_Lab[i_ran][i_type] = new TH1D(Form("h1D_cos_theta_star_Lab_%d_%d",i_ran,i_type),Form("h1D_cos_theta_star_Lab_%d_%d",i_ran,i_type),100,-1,1);
			h1D_cos_theta_star_PairRest[i_ran][i_type] =new TH1D(Form("h1D_cos_theta_star_PairRest_%d_%d",i_ran,i_type),Form("h1D_cos_theta_star_PairRest_%d_%d",i_ran,i_type),100,-1,1);
		}
	}

	//
	h1D_Trigger       = new TH1D("h1D_Trigger","h1D_Trigger",6,0.5,6.5);
	h1D_MB_Trigger    = new TH1D("h1D_MB_Trigger","h1D_MB_Trigger",3,0.5,3.5);
	h1D_HM_Trigger    = new TH1D("h1D_HM_Trigger","h1D_HM_Trigger",2,0.5,2.5);

	h1D_MB_HM_Trigger = new TH1D("h1D_MB_HM_Trigger","h1D_MB_HM_Trigger",3,0.5,3.5);


}



void ntp_Lambda_Histogram::Fill_QAplots(std::vector<int> GoodLambdaFlag){
	for (unsigned int i=0; i < Reader->NLambda; i ++){
		if(GoodLambdaFlag[i]==0) continue;
		h1D_p1_Pt   ->Fill(Reader->p1_pt[i]);
		h1D_p1_Eta  ->Fill(Reader->p1_eta[i]);
		h1D_p1_Phi  ->Fill(Reader->p1_phi[i]);

		h1D_p2_Pt   ->Fill(Reader->p2_pt[i]);
		h1D_p2_Eta  ->Fill(Reader->p2_eta[i]);
		h1D_p2_Phi  ->Fill(Reader->p2_phi[i]);
	
		h1D_pair_Pt ->Fill(Reader->pair_pt[i]);
		h1D_pair_Eta->Fill(Reader->pair_eta[i]);
		h1D_pair_Phi->Fill(Reader->pair_phi[i]);
		h1D_pair_Mass->Fill(Reader->pair_mass[i]);
	}
	/*
	int MB_TriggerFlag = 0;
	int HM_TriggerFlag = 0; 


	for(unsigned int i =0 ; i < Reader->mNTrigs;i++){
		if( Reader->mTrigId[i] == 		910001 ){h1D_Trigger->Fill(1);h1D_MB_Trigger->Fill(1);MB_TriggerFlag=1;   }
		else if( Reader->mTrigId[i] ==  910003 ){h1D_Trigger->Fill(2);h1D_MB_Trigger->Fill(2);MB_TriggerFlag=1;   }
		else if( Reader->mTrigId[i] ==  910013 ){h1D_Trigger->Fill(3);h1D_MB_Trigger->Fill(3);MB_TriggerFlag=1;   }
		else if( Reader->mTrigId[i] ==  910802 ){h1D_Trigger->Fill(4);h1D_HM_Trigger->Fill(1);HM_TriggerFlag=1;   }
		else if( Reader->mTrigId[i] ==  910804 ){h1D_Trigger->Fill(5);h1D_HM_Trigger->Fill(2);HM_TriggerFlag=1;   }
		else {h1D_Trigger->Fill(6);}
	}

	if(MB_TriggerFlag == 1 && HM_TriggerFlag == 1 ) h1D_MB_HM_Trigger->Fill(3);
	if(MB_TriggerFlag == 1 && HM_TriggerFlag == 0 ) h1D_MB_HM_Trigger->Fill(1);
	if(MB_TriggerFlag == 0 && HM_TriggerFlag == 1 ) h1D_MB_HM_Trigger->Fill(2);
	*/
}



void ntp_Lambda_Histogram::Fill_PairPlots(TLorentzVector *v1, TLorentzVector *v2,int type_bin,int range_bin){
	h2D_L1Mass_L2Mass[range_bin][type_bin]->Fill(v1->M(),v2->M());          //h2D_L1Mass_L2Mass[Range_Bin-1][type_bin]->Fill(v1->M(),v2->M());
	h2D_L1Pt_L2Pt[range_bin][type_bin]->Fill(v1->Pt(),v2->Pt());			//h2D_L1Pt_L2Pt[Range_Bin-1][type_bin]->Fill(v1->Pt(),v2->Pt());					
	h2D_L1Phi_L2Phi[range_bin][type_bin]->Fill(v1->Phi(),v2->Phi());		//h2D_L1Phi_L2Phi[Range_Bin-1][type_bin]->Fill(v1->Phi(),v2->Phi());	
	h2D_L1Eta_L2Eta[range_bin][type_bin]->Fill(v1->Eta(),v2->Eta());		//h2D_L1Eta_L2Eta[Range_Bin-1][type_bin]->Fill(v1->Eta(),v2->Eta());

	TLorentzVector pair = (*v1)+(*v2);
	double delta_y = v1->Rapidity() - v2->Rapidity();
	double delta_phi = TMath::ACos( TMath::Cos( v1->Phi() - v2->Phi() ) );
	//double phi_1   = std::fmod(v1->Phi() + TMath::Pi()*2, TMath::Pi()*2 );
	//double phi_2   = std::fmod(v2->Phi() + TMath::Pi()*2, TMath::Pi()*2 ); 
	//double delta_phi = std::fmod(phi_1 - phi_2  + TMath::Pi()*2, TMath::Pi()*2  );
	//if( delta_phi>TMath::Pi() ) delta_phi= 2*TMath::Pi() - delta_phi;
	//std::cout<<delta_phi <<std::endl;
	//std::cout<<TMath::ACos( TMath::Cos( v1->Phi() - v2->Phi() ) ) <<std::endl;
	double deltaR = TMath::Sqrt( delta_phi * delta_phi + delta_y * delta_y  );



	h2D_PairMass_deltaR[type_bin]->Fill(pair.M(),deltaR);
	h1D_L1L2_PtDifference[type_bin]->Fill( v1->Pt() - v2->Pt() );
}

void ntp_Lambda_Histogram::Fill_PairPlots(TLorentzVector *v1, TLorentzVector *v2,int type_bin,int range_bin,double weight){
	h2D_L1Mass_L2Mass[range_bin][type_bin]->Fill(v1->M(),v2->M(),weight);          //h2D_L1Mass_L2Mass[Range_Bin-1][type_bin]->Fill(v1->M(),v2->M());
	h2D_L1Pt_L2Pt[range_bin][type_bin]->Fill(v1->Pt(),v2->Pt(),weight);			//h2D_L1Pt_L2Pt[Range_Bin-1][type_bin]->Fill(v1->Pt(),v2->Pt());					
	h2D_L1Phi_L2Phi[range_bin][type_bin]->Fill(v1->Phi(),v2->Phi(),weight);		//h2D_L1Phi_L2Phi[Range_Bin-1][type_bin]->Fill(v1->Phi(),v2->Phi());	
	h2D_L1Eta_L2Eta[range_bin][type_bin]->Fill(v1->Eta(),v2->Eta(),weight);		//h2D_L1Eta_L2Eta[Range_Bin-1][type_bin]->Fill(v1->Eta(),v2->Eta());

	TLorentzVector pair = (*v1)+(*v2);
	double delta_y = v1->Rapidity() - v2->Rapidity();
	double delta_phi = TMath::ACos( TMath::Cos( v1->Phi() - v2->Phi() ) );
	//double phi_1   = std::fmod(v1->Phi() + TMath::Pi()*2, TMath::Pi()*2 );
	//double phi_2   = std::fmod(v2->Phi() + TMath::Pi()*2, TMath::Pi()*2 ); 
	//double delta_phi = std::fmod(phi_1 - phi_2  + TMath::Pi()*2, TMath::Pi()*2  );
	//if( delta_phi>TMath::Pi() ) delta_phi= 2*TMath::Pi() - delta_phi;
	//std::cout<<delta_phi <<std::endl;
	//std::cout<<TMath::ACos( TMath::Cos( v1->Phi() - v2->Phi() ) ) <<std::endl;
	double deltaR = TMath::Sqrt( delta_phi * delta_phi + delta_y * delta_y  );



	h2D_PairMass_deltaR[type_bin]->Fill(pair.M(),deltaR,weight);
	h1D_L1L2_PtDifference[type_bin]->Fill( v1->Pt() - v2->Pt(),weight );
}




void ntp_Lambda_Histogram::Fill_DensityMatrix(int type_bin , int range_bin){
	h1D_B_a_1[range_bin][type_bin]->Fill( Calculator->B_a_1 );     // h1D_B_a_1[Range_Bin-1][type_bin]->Fill( Calculator->B_a_1 );
	h1D_B_a_2[range_bin][type_bin]->Fill( Calculator->B_a_2 );	   //	h1D_B_a_2[Range_Bin-1][type_bin]->Fill( Calculator->B_a_2 );
	h1D_B_a_3[range_bin][type_bin]->Fill( Calculator->B_a_3 );	   //	h1D_B_a_3[Range_Bin-1][type_bin]->Fill( Calculator->B_a_3 );

	h1D_B_b_1[range_bin][type_bin]->Fill( Calculator->B_b_1 );	   //	h1D_B_b_1[Range_Bin-1][type_bin]->Fill( Calculator->B_b_1 );
	h1D_B_b_2[range_bin][type_bin]->Fill( Calculator->B_b_2 );	   //	h1D_B_b_2[Range_Bin-1][type_bin]->Fill( Calculator->B_b_2 );
	h1D_B_b_3[range_bin][type_bin]->Fill( Calculator->B_b_3 );	   //	h1D_B_b_3[Range_Bin-1][type_bin]->Fill( Calculator->B_b_3 );

	h1D_C_11[range_bin][type_bin]->Fill( Calculator->C_11 );	   //	h1D_C_11[Range_Bin-1][type_bin]->Fill( Calculator->C_11 );
	h1D_C_12[range_bin][type_bin]->Fill( Calculator->C_12 );	   //	h1D_C_12[Range_Bin-1][type_bin]->Fill( Calculator->C_12 );
	h1D_C_13[range_bin][type_bin]->Fill( Calculator->C_13 );	   //	h1D_C_13[Range_Bin-1][type_bin]->Fill( Calculator->C_13 );

	h1D_C_21[range_bin][type_bin]->Fill( Calculator->C_21 );	   //	h1D_C_21[Range_Bin-1][type_bin]->Fill( Calculator->C_21 );
	h1D_C_22[range_bin][type_bin]->Fill( Calculator->C_22 );	   //	h1D_C_22[Range_Bin-1][type_bin]->Fill( Calculator->C_22 );
	h1D_C_23[range_bin][type_bin]->Fill( Calculator->C_23 );	   //	h1D_C_23[Range_Bin-1][type_bin]->Fill( Calculator->C_23 );


	h1D_C_31[range_bin][type_bin]->Fill( Calculator->C_31 );	   //	h1D_C_31[Range_Bin-1][type_bin]->Fill( Calculator->C_31 );
	h1D_C_32[range_bin][type_bin]->Fill( Calculator->C_32 );	   //	h1D_C_32[Range_Bin-1][type_bin]->Fill( Calculator->C_32 );
	h1D_C_33[range_bin][type_bin]->Fill( Calculator->C_33 );	   //	h1D_C_33[Range_Bin-1][type_bin]->Fill( Calculator->C_33 );

	h1D_cos_theta_star_Lab[range_bin][type_bin]->Fill(Calculator->cos_theta_star_Lab);				
	h1D_cos_theta_star_PairRest[range_bin][type_bin]->Fill(Calculator->cos_theta_star_PairRest);

	//h1D_cos_theta_star_Lab[Range_Bin-1][type_bin]->Fill(Calculator->cos_theta_star_Lab);
	//h1D_cos_theta_star_PairRest[Range_Bin-1][type_bin]->Fill(Calculator->cos_theta_star_PairRest);

	
}

void ntp_Lambda_Histogram::Fill_DensityMatrix(int type_bin , int range_bin,double weight){
	h1D_B_a_1[range_bin][type_bin]->Fill( Calculator->B_a_1 ,weight);      //h1D_B_a_1[Range_Bin-1][type_bin]->Fill( Calculator->B_a_1 ,weight);
	h1D_B_a_2[range_bin][type_bin]->Fill( Calculator->B_a_2 ,weight);		//h1D_B_a_2[Range_Bin-1][type_bin]->Fill( Calculator->B_a_2,weight );
	h1D_B_a_3[range_bin][type_bin]->Fill( Calculator->B_a_3 ,weight);		//h1D_B_a_3[Range_Bin-1][type_bin]->Fill( Calculator->B_a_3,weight );

	h1D_B_b_1[range_bin][type_bin]->Fill( Calculator->B_b_1 ,weight);		//h1D_B_b_1[Range_Bin-1][type_bin]->Fill( Calculator->B_b_1,weight );
	h1D_B_b_2[range_bin][type_bin]->Fill( Calculator->B_b_2 ,weight);		//h1D_B_b_2[Range_Bin-1][type_bin]->Fill( Calculator->B_b_2,weight );
	h1D_B_b_3[range_bin][type_bin]->Fill( Calculator->B_b_3 ,weight);		//h1D_B_b_3[Range_Bin-1][type_bin]->Fill( Calculator->B_b_3,weight );

	h1D_C_11[range_bin][type_bin]->Fill( Calculator->C_11 ,weight);		//h1D_C_11[Range_Bin-1][type_bin]->Fill( Calculator->C_11 ,weight);
	h1D_C_12[range_bin][type_bin]->Fill( Calculator->C_12 ,weight);		//h1D_C_12[Range_Bin-1][type_bin]->Fill( Calculator->C_12 ,weight);
	h1D_C_13[range_bin][type_bin]->Fill( Calculator->C_13 ,weight);		//h1D_C_13[Range_Bin-1][type_bin]->Fill( Calculator->C_13 ,weight);

	h1D_C_21[range_bin][type_bin]->Fill( Calculator->C_21 ,weight);		//h1D_C_21[Range_Bin-1][type_bin]->Fill( Calculator->C_21 ,weight);
	h1D_C_22[range_bin][type_bin]->Fill( Calculator->C_22 ,weight);		//h1D_C_22[Range_Bin-1][type_bin]->Fill( Calculator->C_22 ,weight);
	h1D_C_23[range_bin][type_bin]->Fill( Calculator->C_23 ,weight);		//h1D_C_23[Range_Bin-1][type_bin]->Fill( Calculator->C_23 ,weight);


	h1D_C_31[range_bin][type_bin]->Fill( Calculator->C_31 ,weight);		//h1D_C_31[Range_Bin-1][type_bin]->Fill( Calculator->C_31 ,weight);
	h1D_C_32[range_bin][type_bin]->Fill( Calculator->C_32 ,weight);		//h1D_C_32[Range_Bin-1][type_bin]->Fill( Calculator->C_32 ,weight);
	h1D_C_33[range_bin][type_bin]->Fill( Calculator->C_33 ,weight);		//h1D_C_33[Range_Bin-1][type_bin]->Fill( Calculator->C_33 ,weight);

	h1D_cos_theta_star_Lab[range_bin][type_bin]->Fill(Calculator->cos_theta_star_Lab,weight);				
	h1D_cos_theta_star_PairRest[range_bin][type_bin]->Fill(Calculator->cos_theta_star_PairRest,weight);

	//h1D_cos_theta_star_Lab[Range_Bin-1][type_bin]->Fill(Calculator->cos_theta_star_Lab,weight);
	//h1D_cos_theta_star_PairRest[Range_Bin-1][type_bin]->Fill(Calculator->cos_theta_star_PairRest,weight);

	
}


void ntp_Lambda_Histogram::Reset(){
	h1D_p1_Pt->Reset("ICES");
	h1D_p1_Eta->Reset("ICES");
	h1D_p1_Phi->Reset("ICES");

	h1D_p2_Pt->Reset("ICES");
	h1D_p2_Eta->Reset("ICES");
	h1D_p2_Phi->Reset("ICES");
	
	h1D_pair_Pt->Reset("ICES");
	h1D_pair_Eta->Reset("ICES");
	h1D_pair_Phi->Reset("ICES");
	h1D_pair_Mass->Reset("ICES");

	

	for(int i_r=0;i_r < Range_Bin ; i_r ++){
		for(int i_t=0;i_t < 3; i_t++){
			h2D_L1Mass_L2Mass[i_r][i_t]->Reset("ICES");
			h2D_L1Pt_L2Pt[i_r][i_t]->Reset("ICES");
			h2D_L1Phi_L2Phi[i_r][i_t]->Reset("ICES");
			h2D_L1Eta_L2Eta[i_r][i_t]->Reset("ICES");

			h1D_B_a_1[i_r][i_t]->Reset("ICES");
			h1D_B_a_2[i_r][i_t]->Reset("ICES");
			h1D_B_a_3[i_r][i_t]->Reset("ICES");

			h1D_B_b_1[i_r][i_t]->Reset("ICES");
			h1D_B_b_2[i_r][i_t]->Reset("ICES");
			h1D_B_b_3[i_r][i_t]->Reset("ICES");

			h1D_C_11[i_r][i_t]->Reset("ICES");
			h1D_C_12[i_r][i_t]->Reset("ICES");
			h1D_C_13[i_r][i_t]->Reset("ICES");

			h1D_C_21[i_r][i_t]->Reset("ICES");
			h1D_C_22[i_r][i_t]->Reset("ICES");
			h1D_C_23[i_r][i_t]->Reset("ICES");


			h1D_C_31[i_r][i_t]->Reset("ICES");
			h1D_C_32[i_r][i_t]->Reset("ICES");
			h1D_C_33[i_r][i_t]->Reset("ICES");

			h1D_cos_theta_star_Lab[i_r][i_t]->Reset("ICES");
			h1D_cos_theta_star_PairRest[i_r][i_t]->Reset("ICES");
		}
	}

	
	for(int i_t=0;i_t<3;i_t++){
		h2D_PairMass_deltaR[i_t]->Reset("ICES");
		h1D_L1L2_PtDifference[i_t]->Reset("ICES");
	}
}

void ntp_Lambda_Histogram::WriteAll(){
	TFile *fout= TFile::Open(OutPutFile.c_str(),"RECREATE");

	h1D_p1_Pt->Write();
	h1D_p1_Eta->Write();
	h1D_p1_Phi->Write();

	h1D_p2_Pt->Write();
	h1D_p2_Eta->Write();
	h1D_p2_Phi->Write();
	
	h1D_pair_Pt->Write();
	h1D_pair_Eta->Write();
	h1D_pair_Phi->Write();
	h1D_pair_Mass->Write();

	

	for(int i_r=0;i_r < Range_Bin ; i_r ++){
		for(int i_t=0;i_t < 3; i_t++){
			h2D_L1Mass_L2Mass[i_r][i_t]->Write();
			h2D_L1Pt_L2Pt[i_r][i_t]->Write();
			h2D_L1Phi_L2Phi[i_r][i_t]->Write();
			h2D_L1Eta_L2Eta[i_r][i_t]->Write();

			h1D_B_a_1[i_r][i_t]->Write();
			h1D_B_a_2[i_r][i_t]->Write();
			h1D_B_a_3[i_r][i_t]->Write();

			h1D_B_b_1[i_r][i_t]->Write();
			h1D_B_b_2[i_r][i_t]->Write();
			h1D_B_b_3[i_r][i_t]->Write();

			h1D_C_11[i_r][i_t]->Write();
			h1D_C_12[i_r][i_t]->Write();
			h1D_C_13[i_r][i_t]->Write();

			h1D_C_21[i_r][i_t]->Write();
			h1D_C_22[i_r][i_t]->Write();
			h1D_C_23[i_r][i_t]->Write();


			h1D_C_31[i_r][i_t]->Write();
			h1D_C_32[i_r][i_t]->Write();
			h1D_C_33[i_r][i_t]->Write();

			h1D_cos_theta_star_Lab[i_r][i_t]->Write();
			h1D_cos_theta_star_PairRest[i_r][i_t]->Write();
		}
	}

	
	for(int i_t=0;i_t<3;i_t++){
		h2D_PairMass_deltaR[i_t]->Write();
		h1D_L1L2_PtDifference[i_t]->Write();
	}
	h1D_Trigger      ->Write();
	h1D_MB_Trigger   ->Write();
	h1D_HM_Trigger   ->Write();

	h1D_MB_HM_Trigger ->Write();


	fout->Close();
	delete fout;
	

	
	

}






#endif 
