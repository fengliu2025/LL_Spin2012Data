#ifndef ntp_Lambda_Histogram_h
#define ntp_Lambda_Histogram_h

//include the header file here 
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TMath.h>

#include <vector>
#include <string>
#include <iostream>


#include "ntp_Lambda_Reader.h"
#include "constants.h"

class ntp_Lambda_Histogram {
public : 
	ntp_Lambda_Reader *Reader;

	std::string OutPutFile;
	//Histograms for QA plot
	//QA plots. Basics Distribution of daughter particles.(After Selection).
	TH1D *h1D_TrigId;
	TH1D *h1D_NTrig;

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

	TH1D *h1D_pair_DCAdaughters;
	TH1D *h1D_pair_theta;
	TH1D *h1D_pair_decayL;


	//Histograms for Pair Mass
	TH2D *h2D_USL_USLbar_MassDist[Pt2_Bin][Range_Bin];
	TH2D *h2D_USL_USL_MassDist[Pt2_Bin][Range_Bin];
	TH2D *h2D_USLbar_USLbar_MassDist[Pt2_Bin][Range_Bin];

	TH2D *h2D_USL_LS_MassDist[Pt2_Bin][Range_Bin];
	TH2D *h2D_USLbar_LS_MassDist[Pt2_Bin][Range_Bin];
	TH2D *h2D_LS_LS_MassDist[Pt2_Bin][Range_Bin];
	//Histograms for 1D-mass
	TH1D *h1D_USL_MassDist[Pt1_Bin];
	TH1D *h1D_USLbar_MassDist[Pt1_Bin];
	TH1D *h1D_LS_Plus_MassDist[Pt1_Bin];
	TH1D *h1D_LS_Minus_MassDist[Pt1_Bin];
	TH1D *h1D_LS_MassDist[Pt1_Bin];


	

	ntp_Lambda_Histogram();
	ntp_Lambda_Histogram(ntp_Lambda_Reader *reader,std::string OutPutFile);
	void InitHistogram();
	void Fill_QAplots();
	void Fill_Mass();
	void Reset();

	void WriteAll();

};

ntp_Lambda_Histogram::ntp_Lambda_Histogram()
{}

ntp_Lambda_Histogram::ntp_Lambda_Histogram(ntp_Lambda_Reader *reader,std::string outputFile)
:Reader(reader),OutPutFile(outputFile)
{

}

void ntp_Lambda_Histogram::InitHistogram(){
	h1D_TrigId	   = new TH1D("h1D_TrigId","h1D_TrigId",6,0.5,6.5);
	h1D_NTrig      = new TH1D("h1D_NTrig" ,"h1D_NTrig" ,6,0.5,6.5);

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

	h1D_pair_DCAdaughters = new TH1D("h1D_pair_DCAdaughters","h1D_pair_DCAdaughters", 100 , 0,5); 
	h1D_pair_theta        = new TH1D("h1D_pair_theta","h1D_pair_theta",1000,-TMath::Pi()/10,-TMath::Pi()/10);
	h1D_pair_decayL 	  = new TH1D("h1D_pair_decayL","h1D_pair_decayL",300,0,30); 
	for(int i_pt=0; i_pt < Pt2_Bin ; i_pt ++){
		for(int i_ran = 0; i_ran < Range_Bin ; i_ran ++){
			h2D_USL_USLbar_MassDist[i_pt][i_ran]     = new TH2D( Form("h2D_USL_USLbar_MassDist_%d_%d",i_pt,i_ran), Form("h2D_USL_USLbar_MassDist_%d_%d",i_pt,i_ran),120,1.07,1.20,120,1.07,1.20 );
			h2D_USL_USL_MassDist[i_pt][i_ran]        = new TH2D( Form("h2D_USL_USL_MassDist_%d_%d",i_pt,i_ran)	 , Form("h2D_USL_USL_MassDist_%d_%d",i_pt,i_ran)   ,120,1.07,1.20,120,1.07,1.20 );
			h2D_USLbar_USLbar_MassDist[i_pt][i_ran]  = new TH2D( Form("h2D_USLbar_USLbar_MassDist_%d_%d",i_pt,i_ran), Form("h2D_USLbar_USLbar_MassDist_%d_%d",i_pt,i_ran),120,1.07,1.20,120,1.07,1.20);

			h2D_USL_LS_MassDist[i_pt][i_ran]         = new TH2D( Form("h2D_USL_LS_MassDist_%d_%d",i_pt,i_ran) , Form("h2D_USL_LS_MassDist_%d_%d",i_pt,i_ran), 120, 1.07, 1.20, 120, 1.07 , 1.20 );
			h2D_USLbar_LS_MassDist[i_pt][i_ran]      = new TH2D( Form("h2D_USLbar_LS_MassDist_%d_%d",i_pt,i_ran),Form("h2D_USLbar_LS_MassDist_%d_%d",i_pt,i_ran),120, 1.07 ,1.20, 120, 1.07 ,1.20);
			h2D_LS_LS_MassDist[i_pt][i_ran]          = new TH2D( Form("h2D_LS_LS_MassDist_%d_%d",i_pt,i_ran), Form("h2D_LS_LS_MassDist_%d_%d",i_pt,i_ran), 120, 1.07 , 1.20, 120 , 1.07 ,1.20);
		}
	}
	
	for(int i_pt=0; i_pt < Pt1_Bin ; i_pt++){
		h1D_USL_MassDist[i_pt]      = new TH1D( Form("h1D_USL_MassDist_%d",i_pt) , Form("h1D_USL_MassDist_%d",i_pt), 120, 1.07 ,1.20);
		h1D_USLbar_MassDist[i_pt]   = new TH1D( Form("h1D_USLbar_MassDist_%d",i_pt), Form("h1D_USLbar_MassDist_%d",i_pt), 120, 1.07 ,1.20);
		h1D_LS_Plus_MassDist[i_pt]  = new TH1D( Form("h1D_LS_Plus_MassDist_%d",i_pt),Form("h1D_LS_Plus_MassDist_%d",i_pt), 120, 1.07, 1.20);
		h1D_LS_Minus_MassDist[i_pt] = new TH1D( Form("h1D_LS_Minus_MassDist_%d",i_pt),Form("h1D_LS_Minus_MassDist_%d",i_pt),120,1.07,1.20);
		h1D_LS_MassDist[i_pt]       = new TH1D( Form("h1D_LS_MassDist_%d",i_pt),Form("h1D_LS_MassDist_%d",i_pt),120,1.07,1.20);
	}
	
	


}



void ntp_Lambda_Histogram::Fill_QAplots(){
	h1D_NTrig  ->Fill(Reader->mNTrigs);
	for(int i =0 ;i < Reader->mNTrigs;i++){
		if(Reader->mTrigId[i] == 910001 ) 		h1D_TrigId->Fill(1);
		else if(Reader->mTrigId[i] == 910003) h1D_TrigId->Fill(2);
		else if(Reader->mTrigId[i] == 910013) h1D_TrigId->Fill(3);
		else if(Reader->mTrigId[i] == 910802) h1D_TrigId->Fill(4);
		else if(Reader->mTrigId[i] == 910804) h1D_TrigId->Fill(5);
		else h1D_TrigId->Fill(6);
	}

	
	for (unsigned int i=0; i < Reader->NLambda; i ++){
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

		h1D_pair_DCAdaughters ->Fill(Reader->pair_DCAdaughters[i]); 
		h1D_pair_theta        ->Fill(Reader->pair_theta[i]);
		h1D_pair_decayL 	  ->Fill(Reader->pair_decayL[i]); 
	}
}


void ntp_Lambda_Histogram::Fill_Mass(){

	//judge the type of lambdas(pair)
	int lambda_type[2] = {-1,-1};
	//lambda_type = 1: lambda 
	//lambda_type = 2: lambda_bar 
	//lambda_type = 3: LS + sign 
	//lambda_type = 4: LS - sign 

	for(int i=0;i <2; i++){
		if(Reader->pair_charge[i] == 0 && Reader->p1_ch[i] == +1 ) lambda_type[i] = 1; 
		if(Reader->pair_charge[i] == 0 && Reader->p1_ch[i] == -1 ) lambda_type[i] = 2;
		if(Reader->pair_charge[i] == 1 && Reader->p1_ch[i] == +1 ) lambda_type[i] = 3;
		if(Reader->pair_charge[i] == 1 && Reader->p1_ch[i] == -1 ) lambda_type[i] = 4;

	}
	//judge the bin of pt 
	int lambda_pt_idbin[2]={-1,-1};

	for(int i = 0 ; i < 2; i ++ ){
		for(int i_pt1=0; i_pt1 < (Pt1_Bin-1); i_pt1 ++ ){

			if( Reader->pair_pt[i] > Lambda_Pt_LowerBin[i_pt1] &&  Reader->pair_pt[i] < Lambda_Pt_HigherBin[i_pt1] ) lambda_pt_idbin[i] = i_pt1;
		
		}
	}

	//Fill the single pair histograms 
	//--------------------------------------------------------------------------------------------------------
	for(int i = 0 ; i <2 ;i++){
		if( lambda_type[i] == -1 ) std::cout<<"Something went wrong with the pair classification" <<std::endl;
		if( lambda_type[i] == 1  ) {
			h1D_USL_MassDist[lambda_pt_idbin[i]] ->Fill( Reader->pair_mass[i] );
			h1D_USL_MassDist[Pt1_Bin-1]->Fill( Reader->pair_mass[i] );
		}

		if( lambda_type[i] == 2  ) {
			h1D_USLbar_MassDist[lambda_pt_idbin[i]] ->Fill( Reader->pair_mass[i] );
			h1D_USLbar_MassDist[Pt1_Bin-1] ->Fill( Reader->pair_mass[i] );
		}

		if( lambda_type[i] == 3  ) {
			h1D_LS_Plus_MassDist[lambda_pt_idbin[i]] ->Fill( Reader->pair_mass[i] );
			h1D_LS_Plus_MassDist[Pt1_Bin-1] ->Fill( Reader->pair_mass[i] );
			h1D_LS_MassDist[lambda_pt_idbin[i]]->Fill( Reader->pair_mass[i] );
			h1D_LS_MassDist[Pt1_Bin-1]->Fill( Reader->pair_mass[i] );
		}

		if( lambda_type[i] == 4  ) {
			h1D_LS_Minus_MassDist[lambda_pt_idbin[i]]->Fill( Reader->pair_mass[i] );
			h1D_LS_Minus_MassDist[Pt1_Bin-1]->Fill( Reader->pair_mass[i] );
			h1D_LS_MassDist[lambda_pt_idbin[i]]->Fill( Reader->pair_mass[i] );
			h1D_LS_MassDist[Pt1_Bin-1]->Fill( Reader->pair_mass[i] );
		}

	}
	//-----------------------------------------------------------------------------------------------------------

	//judge the bin of range 
	//range_bin = 0 : short range 
	//range_bin = 1 : long range 
	//range_bin = 2 : all range  
	int range_bin = Range_Bin-1;
	TLorentzVector v1; 
	v1.SetPtEtaPhiM( Reader->pair_pt[0],Reader->pair_eta[0],Reader->pair_phi[0],Reader->pair_mass[0] );
	TLorentzVector v2; 
	v2.SetPtEtaPhiM( Reader->pair_pt[1],Reader->pair_eta[1],Reader->pair_phi[1],Reader->pair_mass[1] );
	double delta_y = TMath::Abs( v1.Rapidity() - v2.Rapidity() );

	double phi_1   = std::fmod(v1.Phi() + TMath::Pi()*2, TMath::Pi()*2 );
	double phi_2   = std::fmod(v2.Phi() + TMath::Pi()*2, TMath::Pi()*2 ); 
	double delta_phi = std::fmod( phi_1 - phi_2  + TMath::Pi()*2, TMath::Pi()*2  );
	if( delta_phi>TMath::Pi() ) delta_phi= 2*TMath::Pi() - delta_phi;
	//std::cout<<delta_phi <<"="<<TMath::ACos( TMath::Cos( v1.Phi() - v2.Phi() ) )<<std::endl;
	

	if( delta_y < 0.5 && delta_phi < TMath::Pi()/3. ) range_bin = 0;
	else range_bin = 1;


	
	int pair_type = 0;
	//pair_type = 1 : (US)L 	(US)L_bar
	//pair_type = 2 : (US)L 	(US)L
	//pair_type = 3 : (US)L_bar (US)L_bar
	//pair_type = 4 : (US)L 	(LS)
	//pair_tyoe = 5 : (US)L_bar (LS)
	//pair_tyoe = 6 : (LS)		(LS)

	int idx1 = 0;
	int idx2 = 1; 

	if(lambda_type[0] == 1 && lambda_type[1] == 1 ){idx1 = 0; idx2 = 1;pair_type = 2;}
	if(lambda_type[0] == 1 && lambda_type[1] == 2 ){idx1 = 0; idx2 = 1;pair_type = 1;}
	if(lambda_type[0] == 1 && lambda_type[1] > 2  ){idx1 = 0; idx2 = 1;pair_type = 4;}

	if(lambda_type[0] == 2 && lambda_type[1] == 1 ){idx1 = 1; idx2 = 0;pair_type = 1;}
	if(lambda_type[0] == 2 && lambda_type[1] == 2 ){idx1 = 0; idx2 = 1;pair_type = 3;}
	if(lambda_type[0] == 2 && lambda_type[1] > 2  ){idx1 = 0; idx2 = 1;pair_type = 5;}

	if(lambda_type[0] > 2  && lambda_type[1] == 1 ){idx1 = 1; idx2 = 0;pair_type = 4;}
	if(lambda_type[0] > 2  && lambda_type[1] == 2 ){idx1 = 1; idx2 = 0;pair_type = 5;}
	if(lambda_type[0] > 2  && lambda_type[1] >2   ){idx1 = 0; idx2 = 1;pair_type = 6;}

	//judege the bin of pt2
	// pt2_bin = 0 : 0 0 
	// pt2_bin = 1 : 0 1
	// pt2_bin = 2 : 1 0 
	// pt2_bin = 3 : 1 1
	// pt2_bin = 4 : all 
	int pt2_bin = Pt2_Bin - 1;
	if( lambda_pt_idbin[idx1] == 0 && lambda_pt_idbin[idx2] == 0  ) pt2_bin = 0;
	if( lambda_pt_idbin[idx1] == 0 && lambda_pt_idbin[idx2] == 1  ) pt2_bin = 1;
	if( lambda_pt_idbin[idx1] == 1 && lambda_pt_idbin[idx2] == 0  ) pt2_bin = 2;
	if( lambda_pt_idbin[idx1] == 1 && lambda_pt_idbin[idx2] == 1  ) pt2_bin = 3;

	//Fill the 2D-mass distribution 
	if(pair_type == 1 ){
		h2D_USL_USLbar_MassDist[pt2_bin][range_bin]->Fill(Reader->pair_mass[idx1],Reader->pair_mass[idx2]);
		h2D_USL_USLbar_MassDist[Pt2_Bin-1][range_bin]->Fill(Reader->pair_mass[idx1],Reader->pair_mass[idx2]);
		h2D_USL_USLbar_MassDist[pt2_bin][Range_Bin-1]->Fill(Reader->pair_mass[idx1],Reader->pair_mass[idx2]);
		h2D_USL_USLbar_MassDist[Pt2_Bin-1][Range_Bin-1]->Fill(Reader->pair_mass[idx1],Reader->pair_mass[idx2]);
	}
	if(pair_type == 2 ){
		h2D_USL_USL_MassDist[pt2_bin][range_bin]->Fill(Reader->pair_mass[idx1],Reader->pair_mass[idx2]);
		h2D_USL_USL_MassDist[Pt2_Bin-1][range_bin]->Fill(Reader->pair_mass[idx1],Reader->pair_mass[idx2]);
		h2D_USL_USL_MassDist[pt2_bin][Range_Bin-1]->Fill(Reader->pair_mass[idx1],Reader->pair_mass[idx2]);
		h2D_USL_USL_MassDist[Pt2_Bin-1][Range_Bin-1]->Fill(Reader->pair_mass[idx1],Reader->pair_mass[idx2]);
	}
	if(pair_type == 3 ){
		h2D_USLbar_USLbar_MassDist[pt2_bin][range_bin]->Fill(Reader->pair_mass[idx1],Reader->pair_mass[idx2]);
		h2D_USLbar_USLbar_MassDist[Pt2_Bin-1][range_bin]->Fill(Reader->pair_mass[idx1],Reader->pair_mass[idx2]);
		h2D_USLbar_USLbar_MassDist[pt2_bin][Range_Bin-1]->Fill(Reader->pair_mass[idx1],Reader->pair_mass[idx2]);
		h2D_USLbar_USLbar_MassDist[Pt2_Bin-1][Range_Bin-1]->Fill(Reader->pair_mass[idx1],Reader->pair_mass[idx2]);
	}
	if(pair_type == 4 ){
		h2D_USL_LS_MassDist[pt2_bin][range_bin]->Fill(Reader->pair_mass[idx1],Reader->pair_mass[idx2]);
		h2D_USL_LS_MassDist[Pt2_Bin-1][range_bin]->Fill(Reader->pair_mass[idx1],Reader->pair_mass[idx2]);
		h2D_USL_LS_MassDist[pt2_bin][Range_Bin-1]->Fill(Reader->pair_mass[idx1],Reader->pair_mass[idx2]);
		h2D_USL_LS_MassDist[Pt2_Bin-1][Range_Bin-1]->Fill(Reader->pair_mass[idx1],Reader->pair_mass[idx2]);
	}
	if(pair_type == 5){
		h2D_USLbar_LS_MassDist[pt2_bin][range_bin]->Fill(Reader->pair_mass[idx1],Reader->pair_mass[idx2]);
		h2D_USLbar_LS_MassDist[Pt2_Bin-1][range_bin]->Fill(Reader->pair_mass[idx1],Reader->pair_mass[idx2]);
		h2D_USLbar_LS_MassDist[pt2_bin][Range_Bin-1]->Fill(Reader->pair_mass[idx1],Reader->pair_mass[idx2]);
		h2D_USLbar_LS_MassDist[Pt2_Bin-1][Range_Bin-1]->Fill(Reader->pair_mass[idx1],Reader->pair_mass[idx2]);
	}
	if(pair_type == 6){
		h2D_LS_LS_MassDist[pt2_bin][range_bin]->Fill(Reader->pair_mass[idx1],Reader->pair_mass[idx2]);
		h2D_LS_LS_MassDist[Pt2_Bin-1][range_bin]->Fill(Reader->pair_mass[idx1],Reader->pair_mass[idx2]);
		h2D_LS_LS_MassDist[pt2_bin][Range_Bin-1]->Fill(Reader->pair_mass[idx1],Reader->pair_mass[idx2]);
		h2D_LS_LS_MassDist[Pt2_Bin-1][Range_Bin-1]->Fill(Reader->pair_mass[idx1],Reader->pair_mass[idx2]);

	}

}











void ntp_Lambda_Histogram::Reset(){
	h1D_TrigId->Reset("ICES");
	h1D_NTrig ->Reset("ICES");    

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


	h1D_pair_DCAdaughters ->Reset("ICES");
	h1D_pair_theta        ->Reset("ICES");
	h1D_pair_decayL 	  ->Reset("ICES");

	for(int i_pt=0; i_pt < Pt2_Bin ; i_pt ++){
		for(int i_ran = 0; i_ran < Range_Bin ; i_ran ++){
			h2D_USL_USLbar_MassDist[i_pt][i_ran]   ->Reset("ICES") ;
			h2D_USL_USL_MassDist[i_pt][i_ran]      ->Reset("ICES") ;
			h2D_USLbar_USLbar_MassDist[i_pt][i_ran]->Reset("ICES") ;

			h2D_USL_LS_MassDist[i_pt][i_ran]       ->Reset("ICES")  ;
			h2D_USLbar_LS_MassDist[i_pt][i_ran]    ->Reset("ICES")  ;
			h2D_LS_LS_MassDist[i_pt][i_ran]        ->Reset("ICES")  ;
		}
	}
	
	for(int i_pt=0; i_pt < Pt1_Bin ; i_pt++){
		h1D_USL_MassDist[i_pt]      ->Reset("ICES");
		h1D_USLbar_MassDist[i_pt]   ->Reset("ICES");
		h1D_LS_Plus_MassDist[i_pt]  ->Reset("ICES");
		h1D_LS_Minus_MassDist[i_pt] ->Reset("ICES");
		h1D_LS_MassDist[i_pt]       ->Reset("ICES");
	}

	
	
}

void ntp_Lambda_Histogram::WriteAll(){
	std::cout<<OutPutFile.c_str()<<std::endl;
	TFile *fout= TFile::Open(OutPutFile.c_str(),"RECREATE");

	h1D_TrigId->Write();
	h1D_NTrig ->Write();

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

	h1D_pair_DCAdaughters ->Write();
	h1D_pair_theta        ->Write();
	h1D_pair_decayL 	  ->Write();

	for(int i_pt=0; i_pt < Pt2_Bin ; i_pt ++){
		for(int i_ran = 0; i_ran < Range_Bin ; i_ran ++){
			h2D_USL_USLbar_MassDist[i_pt][i_ran]   ->Write() ;
			h2D_USL_USL_MassDist[i_pt][i_ran]      ->Write() ;
			h2D_USLbar_USLbar_MassDist[i_pt][i_ran]->Write() ;

			h2D_USL_LS_MassDist[i_pt][i_ran]       ->Write() ;
			h2D_USLbar_LS_MassDist[i_pt][i_ran]    ->Write() ;
			h2D_LS_LS_MassDist[i_pt][i_ran]        ->Write() ;
		}
	}
	
	for(int i_pt=0; i_pt < Pt1_Bin ; i_pt++){
		h1D_USL_MassDist[i_pt]      ->Write();
		h1D_USLbar_MassDist[i_pt]   ->Write();
		h1D_LS_Plus_MassDist[i_pt]  ->Write();
		h1D_LS_Minus_MassDist[i_pt] ->Write();
		h1D_LS_MassDist[i_pt]       ->Write();
	}

	

	fout->Close();
	delete fout;
	

	
	

}






#endif 