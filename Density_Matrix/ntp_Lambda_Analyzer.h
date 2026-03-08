#ifndef ntp_Lambda_Analyzer_h 
#define ntp_Lambda_Analyzer_h

//include the header files here 

#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TLorentzVector.h>

#include <vector>
#include <string>
#include <iostream>
#include <numeric>

#include "constants.h"
#include "ntp_Lambda_Calculator.h"
#include "ntp_Lambda_Histogram.h"
#include "ntp_Lambda_Reader.h"

class ntp_Lambda_Analyzer{
public:
	ntp_Lambda_Reader *SameEvent_Reader;
	ntp_Lambda_Reader *MixEvent_Reader;
	ntp_Lambda_Calculator *Calculator;
	ntp_Lambda_Histogram  *Histogramer; 


	ntp_Lambda_Analyzer();
	ntp_Lambda_Analyzer(ntp_Lambda_Reader *reader1, ntp_Lambda_Reader *reader2, ntp_Lambda_Calculator *calculator, ntp_Lambda_Histogram *hisogram );
	ntp_Lambda_Analyzer(ntp_Lambda_Reader *reader1, ntp_Lambda_Calculator *calculator, ntp_Lambda_Histogram *hisogram);

	~ntp_Lambda_Analyzer();
	

	int Pair_Type_Classifier(int idx1, int idx2);
	int Range_Type_Classifier(TLorentzVector *v1, TLorentzVector *v2);
	int isGoodTrigger();
	int Analyze_SEPair(int i_lambda,int j_lambda);
	void Analysis_SameEvent();
	void FindCounterparts(std::vector<TLorentzVector> *Lambda_counterpart,std::vector<TLorentzVector> *proton_counterpart,std::vector<TLorentzVector> *pion_counterpart,double pt, double rapidity, double phi, int p1Charge,int I_FILE);
	int Analyze_MEPair(int i_lambda,int j_lambda,int i_file);
	void Analysis_MixEvent();


};

ntp_Lambda_Analyzer::ntp_Lambda_Analyzer(){

}

ntp_Lambda_Analyzer::ntp_Lambda_Analyzer(ntp_Lambda_Reader *reader1, ntp_Lambda_Reader *reader2, ntp_Lambda_Calculator *calculator, ntp_Lambda_Histogram *hisogram )
: SameEvent_Reader(reader1),MixEvent_Reader(reader2),Calculator(calculator),Histogramer(hisogram)
{

}


ntp_Lambda_Analyzer::ntp_Lambda_Analyzer(ntp_Lambda_Reader *reader1, ntp_Lambda_Calculator *calculator, ntp_Lambda_Histogram *hisogram )
: SameEvent_Reader(reader1),MixEvent_Reader(0),Calculator(calculator),Histogramer(hisogram)
{

}



int ntp_Lambda_Analyzer::Pair_Type_Classifier(int idx1, int idx2){
	//pair type 
	//1: Lambda Lambda_bar           -1:Lambda_bar 
	//2: Lambda Lambda 			
	//3: Lambda_bar Lambda_bar 		 
	//4: Lambda LS 					 -4:LS Lambda
	//5: Lambda_bar LS 				 -5:LS Lambda_bar
	//6: LS LS 


	//1: Lambda Lambda_bar 
	if(SameEvent_Reader->pair_charge[idx1] == 0 && SameEvent_Reader->pair_charge[idx2] == 0 && SameEvent_Reader->p1_ch[idx1] == 1 && SameEvent_Reader->p1_ch[idx2] == -1 ) return  1 ;
	//-1:Lambda_bar 
	if(SameEvent_Reader->pair_charge[idx1] == 0 && SameEvent_Reader->pair_charge[idx2] == 0 && SameEvent_Reader->p1_ch[idx1] ==-1 && SameEvent_Reader->p1_ch[idx2] ==  1 ) return -1 ;
	//2: Lambda Lambda 	
	if(SameEvent_Reader->pair_charge[idx1] == 0 && SameEvent_Reader->pair_charge[idx2] == 0 && SameEvent_Reader->p1_ch[idx1] == 1 && SameEvent_Reader->p1_ch[idx2] ==  1 ) return  2 ;
	//3: Lambda_bar Lambda_bar 		
	if(SameEvent_Reader->pair_charge[idx1] == 0 && SameEvent_Reader->pair_charge[idx2] == 0 && SameEvent_Reader->p1_ch[idx1] ==-1 && SameEvent_Reader->p1_ch[idx2] == -1 ) return  3 ;
	//4: Lambda LS 	
	if(SameEvent_Reader->pair_charge[idx1] == 0 && SameEvent_Reader->pair_charge[idx2] == 1 && SameEvent_Reader->p1_ch[idx1] == 1  ) return  4 ;
	//-4:LS Lambda
	if(SameEvent_Reader->pair_charge[idx1] == 1 && SameEvent_Reader->pair_charge[idx2] == 0 && SameEvent_Reader->p1_ch[idx2] == 1  ) return -4 ;
	//5: Lambda_bar LS
	if(SameEvent_Reader->pair_charge[idx1] == 0 && SameEvent_Reader->pair_charge[idx2] == 1 && SameEvent_Reader->p1_ch[idx1] == -1 ) return  5 ;
	//-5:LS Lambda_bar
	if(SameEvent_Reader->pair_charge[idx1] == 1 && SameEvent_Reader->pair_charge[idx2] == 0 && SameEvent_Reader->p1_ch[idx2] == -1 ) return -5 ;
	//6:LS Lambda_bar
	if(SameEvent_Reader->pair_charge[idx1] == 1 && SameEvent_Reader->pair_charge[idx2] == 1) return 6 ;
	return 0;
	

}

int ntp_Lambda_Analyzer::Range_Type_Classifier(TLorentzVector *v1, TLorentzVector *v2){
	double delta_phi 		= TMath::ACos(TMath::Cos( v1->Phi()-v2->Phi() ) );
	double delta_rapidity 	= TMath::Abs( v1->Rapidity() - v2->Rapidity()   );

	//range_Type 
	//0 : short range 
	//1 : long  range 
	if(delta_phi < TMath::Pi()/3. && delta_rapidity < 0.5 ) return 0;
	if(delta_phi > TMath::Pi()/3. || delta_rapidity > 0.5 ) return 1; 
	return -1;  

}


int ntp_Lambda_Analyzer::Analyze_SEPair(int i_lambda,int j_lambda){
		if( SameEvent_Reader->p1_InEventID[i_lambda] == SameEvent_Reader->p1_InEventID[j_lambda] || 
	   		SameEvent_Reader->p1_InEventID[i_lambda] == SameEvent_Reader->p2_InEventID[j_lambda] || 
	   		SameEvent_Reader->p1_InEventID[i_lambda] == SameEvent_Reader->p2_InEventID[j_lambda] || 
	   		SameEvent_Reader->p1_InEventID[j_lambda] == SameEvent_Reader->p2_InEventID[i_lambda] || 
	   		SameEvent_Reader->p1_InEventID[j_lambda] == SameEvent_Reader->p2_InEventID[j_lambda] ||
	   		SameEvent_Reader->p2_InEventID[i_lambda] == SameEvent_Reader->p2_InEventID[j_lambda] ) return 0;
				
			

		//---------------------------Pair Type Classification-------------------------------
		int Pair_Type = 0 ; 
		int id_Lambda1= i_lambda ; 
		int id_Lambda2= j_lambda ;


		Pair_Type = Pair_Type_Classifier(i_lambda,j_lambda);

		if(Pair_Type == 1){
			id_Lambda1 = i_lambda;id_Lambda2 = j_lambda;
		}

		else if(Pair_Type == -1){
			id_Lambda1 = j_lambda;id_Lambda2 = i_lambda;
		}

		else if(Pair_Type == 2){
			id_Lambda1 = i_lambda;id_Lambda2 = j_lambda;
		}
		else if(Pair_Type == 3){
			id_Lambda1 = i_lambda;id_Lambda2 = j_lambda;
		}
		else{
			return 0;
			std::cout<<"Something went wrong with the pair type classifier"<<std::endl;
		}
			
		//---------------------------Create Same-Event Pairs-------------------------------
		TLorentzVector Lambda1; Lambda1.SetPtEtaPhiM( SameEvent_Reader->pair_pt[id_Lambda1], SameEvent_Reader->pair_eta[id_Lambda1], SameEvent_Reader->pair_phi[id_Lambda1], SameEvent_Reader->pair_mass[id_Lambda1]  );
		TLorentzVector proton1; proton1.SetPtEtaPhiM( SameEvent_Reader->p1_pt[id_Lambda1]  , SameEvent_Reader->p1_eta[id_Lambda1]  , SameEvent_Reader->p1_phi[id_Lambda1]  , MASS_PROTON                              );
		TLorentzVector pion1  ;   pion1.SetPtEtaPhiM( SameEvent_Reader->p2_pt[id_Lambda1]  , SameEvent_Reader->p2_eta[id_Lambda1]  , SameEvent_Reader->p2_phi[id_Lambda1]  , MASS_PION                                );

		TLorentzVector Lambda2; Lambda2.SetPtEtaPhiM( SameEvent_Reader->pair_pt[id_Lambda2], SameEvent_Reader->pair_eta[id_Lambda2], SameEvent_Reader->pair_phi[id_Lambda2], SameEvent_Reader->pair_mass[id_Lambda2]  );
		TLorentzVector proton2; proton2.SetPtEtaPhiM( SameEvent_Reader->p1_pt[id_Lambda2]  , SameEvent_Reader->p1_eta[id_Lambda2]  , SameEvent_Reader->p1_phi[id_Lambda2]  , MASS_PROTON                              );
		TLorentzVector pion2;     pion2.SetPtEtaPhiM( SameEvent_Reader->p2_pt[id_Lambda2]  , SameEvent_Reader->p2_eta[id_Lambda2]  , SameEvent_Reader->p2_phi[id_Lambda2]  , MASS_PION                                );


		//---------------------------Range Type Classification-------------------------------
		int Range_Type  = Range_Type_Classifier(&Lambda1,&Lambda2);
			
		//Fill the pair plots 
		Histogramer->Fill_PairPlots(&Lambda1,&Lambda2,TMath::Abs(Pair_Type)-1,Range_Type);

		//Calculate the Density Matrix
		Calculator->Reset(&Lambda1,&proton1,&pion1,&Lambda2,&proton2,&pion2);
		Calculator->Calculation();

		//Fill Density Matrix Histgrams 
		Histogramer->Fill_DensityMatrix(TMath::Abs(Pair_Type)-1,Range_Type);
		return 1;
}

int ntp_Lambda_Analyzer::isGoodTrigger(){
	int GoodTriggerFlag = 0;
	//for(int i=0;i < SameEvent_Reader->mNTrigs;i++){
		//if(SameEvent_Reader->mTrigId[i] == 910802 || SameEvent_Reader->mTrigId[i] == 910804 ) GoodTriggerFlag = 1; 
		//if(SameEvent_Reader->mTrigId[i] == 910001 || SameEvent_Reader->mTrigId[i] == 910003 || SameEvent_Reader->mTrigId[i] == 910013) GoodTriggerFlag = 1; 

	//}
	return GoodTriggerFlag;
}




void ntp_Lambda_Analyzer::Analysis_SameEvent(){

	//Start looping over all inputfiles inthe SameEvent_Reader
	unsigned long N_Inputfiles_SE = SameEvent_Reader->InputFiles.size();
	//-------------------------------Enter i_file loop---------------------------------
	for(unsigned long i_file = 0 ; i_file < N_Inputfiles_SE ; i_file ++){
		//Open the file 
		if(i_file%1==0) std::cout<<"i_file="<<i_file<<std::endl;
		TFile *fin=TFile::Open( SameEvent_Reader->InputFiles[i_file].c_str() ,"READ");
		if(fin==0){
			std::cout<<"Can not open the file:"<< SameEvent_Reader->InputFiles[i_file] <<", Skip this" << std::endl;
			continue;
		}
		TTree *tmp_tree = (TTree *)fin->Get("ntp_Lambda");
		if(tmp_tree==0){
			std::cout<<"Can not get the tree, skip this file."<<std::endl;
			continue;
		}

		SameEvent_Reader->Init(tmp_tree);
		//Get the number of entries in current tree
		Long64_t N_Events=SameEvent_Reader->fChain->GetEntries();

		//---------------------------Enter i_event loop----------------------------
		for(Long64_t i_event=0; i_event < N_Events ; i_event++ ){
			if(i_event % 10000==0) std::cout<<"i_event="<<i_event<<std::endl;
			SameEvent_Reader->fChain->GetEntry(i_event);
			//std::cout<<"NL"<<SameEvent_Reader->NLambda<<std::endl;
			
			
			//------------------------Make some selections on the events-----------------------------
			//if(isGoodTrigger() == 1 ) continue; //select on the triggers ,only events that has MB trigger
			//if(SameEvent_Reader->NLambda<2) continue; // current we only select on multi-Lambdas Events 
			std::vector<int> GoodLambdaFlag;
			//if(SameEvent_Reader->pair_charge[0] ==1 || SameEvent_Reader->pair_charge[1] ==1 ) continue;
			//std::cout<<"NL"<<SameEvent_Reader->NLambda<<std::endl;
			for(int i_lambda = 0; i_lambda<SameEvent_Reader->NLambda;i_lambda++){
				int isGoodLambda = 1; //1 is good lambda, 0 is not good lambad 
				//To judge if this lambda is US.
				if( SameEvent_Reader->pair_charge[i_lambda]==1 ){isGoodLambda = 0;}//
				//cut on the track
				if( SameEvent_Reader->p1_pt[i_lambda] < Track_Pt_LowCut || SameEvent_Reader->p2_pt[i_lambda] < Track_Pt_LowCut ) {isGoodLambda=0;}
				if( TMath::Abs(SameEvent_Reader->p1_eta[i_lambda]) > Track_Eta_Cut || TMath::Abs(SameEvent_Reader->p2_eta[i_lambda]) > Track_Eta_Cut  ) {isGoodLambda=0;}
				if( SameEvent_Reader->p1_dca[i_lambda] < 0.1 || SameEvent_Reader->p2_dca[i_lambda] < 0.3 ) {isGoodLambda=0;}

				//cut on the pair_DCA 
				if( SameEvent_Reader->pair_DCAdaughters[i_lambda] > 1.) {isGoodLambda=0;}
				//cut on the decay length
				if( SameEvent_Reader->pair_decayL[i_lambda] < 2. || SameEvent_Reader->pair_decayL[i_lambda] > 25.  ){isGoodLambda=0;}
				//cut on the cos_theta
				if( TMath::Cos(SameEvent_Reader->pair_theta[i_lambda]) < 0.996 ){isGoodLambda=0;}
				//cut on DCA of L to PV
				if ( ( TMath::Sin(SameEvent_Reader->pair_theta[i_lambda]) * SameEvent_Reader->pair_decayL[i_lambda] ) > 1.0  ){isGoodLambda=0;}
				TLorentzVector v;
				v.SetPtEtaPhiM(SameEvent_Reader->pair_pt[i_lambda],SameEvent_Reader->pair_eta[i_lambda],SameEvent_Reader->pair_phi[i_lambda],SameEvent_Reader->pair_mass[i_lambda]);
				//cut on the lambda_Candidates 
				if( v.Pt() < Lambda_Pt_lowCut || v.Pt() > Lambda_Pt_highCut || TMath::Abs( v.Rapidity() )  > Lambda_Rapidity_Cut ){isGoodLambda=0;}
				if( v.M()  < Lambda_mass_lowCut || v.M() > Lambda_mass_highCut ) {isGoodLambda=0;}

				GoodLambdaFlag.push_back(isGoodLambda);
			}

			
			int NGoodLambda = std::accumulate(GoodLambdaFlag.begin(), GoodLambdaFlag.end(), 0);
			//std::cout<<NGoodLambda<<std::endl;
			if(NGoodLambda < 2) continue;
			//std::cout<<NGoodLambda<<std::endl;

			//Fill Histograms of QA plot;
			Histogramer->Fill_QAplots(GoodLambdaFlag);//Notice that here, if two Lambda Candiates share some tracks, both Lambdas will be filled inthe QAplots.

			for(int i_lambda=0;i_lambda < SameEvent_Reader->NLambda;i_lambda++){
				if( GoodLambdaFlag[i_lambda] == 0 ) continue;
				for(int j_lambda=i_lambda+1;j_lambda < SameEvent_Reader->NLambda;j_lambda++){
					if( GoodLambdaFlag[j_lambda] == 0 ) continue;
					Analyze_SEPair(i_lambda,j_lambda);
				}
			}



		}
		//---------------------------End i_event loop----------------------------
		fin->Close();
		delete fin;
	}

	//-------------------------------Enter i_file loop---------------------------------

}




void ntp_Lambda_Analyzer::FindCounterparts(std::vector<TLorentzVector> *Lambda_counterpart,std::vector<TLorentzVector> *proton_counterpart,std::vector<TLorentzVector> *pion_counterpart,double pt, double rapidity, double phi, int p1Charge,int I_FILE){
	//Start looping over all inputfiles of Mix_Event_Reader
	unsigned long N_Inputfiles_ME = MixEvent_Reader->InputFiles.size();
	//--------------------------------Enter i_file loop------------------------------
	for(unsigned i_file =0;i_file < N_Inputfiles_ME ; i_file++){
		if(i_file!=I_FILE) continue;
		TFile *fin=TFile::Open(MixEvent_Reader->InputFiles[i_file].c_str(),"READ");
		if(fin==0){
			std::cout<<"Can not open the file:"<< SameEvent_Reader->InputFiles[i_file] <<", Skip this" << std::endl;
			continue;
		}
		TTree *tmp_tree = (TTree *)fin->Get("ntp_Lambda");
		if(tmp_tree==0){
			std::cout<<"Can not get the tree, skip this file."<<std::endl;
			continue;
		}
		MixEvent_Reader->Init(tmp_tree);
		//Get the number of entries in current trees 
		Long64_t N_Events=MixEvent_Reader->fChain->GetEntries();

		//-------------------------------Enter i_event loop----------------------------
		for(Long64_t i_event=0;i_event<N_Events;i_event++){
			//if(i_event%10000==0)std::cout<<"i_event"<<i_event<<std::endl;
			MixEvent_Reader->fChain->GetEntry(i_event);
			//------------------------Make some selections on the events-----------------------------
			if(MixEvent_Reader->NLambda!=1) continue;
			if(MixEvent_Reader->pair_charge[0]==1) continue; //must be the lambda/lambda_bar candidates
			if(MixEvent_Reader->p1_ch[0]!=p1Charge) continue;

			if( TMath::Abs(MixEvent_Reader->pair_pt[0] - pt  ) > ptDiffLim ) continue;
			if( TMath::ACos(TMath::Cos(MixEvent_Reader->pair_phi[0]-phi) ) > phiDiffLim ) continue;

			//cut on the tracks 
			if( TMath::Abs(MixEvent_Reader->p1_eta[0]) > Track_Eta_Cut || TMath::Abs(MixEvent_Reader->p2_eta[0]) > Track_Eta_Cut ) continue;
			if( MixEvent_Reader->p1_pt[0] < Track_Pt_LowCut || MixEvent_Reader->p2_pt[0] < Track_Pt_LowCut ) continue;


			//cut on the tracks dca
			if( MixEvent_Reader->p1_dca[0] < 0.1 || MixEvent_Reader->p2_dca[0] < 0.3 ) continue;

			//cut on the pair_DCA 
			if( MixEvent_Reader->pair_DCAdaughters[0] > 1.) continue;
			//cut on the decay length
			if( MixEvent_Reader->pair_decayL[0] < 2. || MixEvent_Reader->pair_decayL[0] > 25.  ) continue;
			//cut on the cos_theta
			if( TMath::Cos(MixEvent_Reader->pair_theta[0]) < 0.996 ) continue;
			//cut on DCA of L to PV
			if ( ( TMath::Sin(MixEvent_Reader->pair_theta[0]) * MixEvent_Reader->pair_decayL[0] ) > 1.0  )continue;



			//cut on the lambda candadates 
			TLorentzVector v;
			v.SetPtEtaPhiM(MixEvent_Reader->pair_pt[0],MixEvent_Reader->pair_eta[0],MixEvent_Reader->pair_phi[0],MixEvent_Reader->pair_mass[0]);
			if( v.Pt() < Lambda_Pt_lowCut || v.Pt() > Lambda_Pt_highCut ) continue;
			if( TMath::Abs( v.Rapidity() ) > Lambda_Rapidity_Cut ) continue;

			if( TMath::Abs(v.Rapidity() - rapidity ) > rapidityDiffLim   ) continue;
			//
			TLorentzVector lambda; lambda.SetPtEtaPhiM( MixEvent_Reader->pair_pt[0] , MixEvent_Reader->pair_eta[0] , MixEvent_Reader->pair_phi[0] , MixEvent_Reader->pair_mass[0] );
			TLorentzVector proton; proton.SetPtEtaPhiM( MixEvent_Reader->p1_pt[0]   , MixEvent_Reader->p1_eta[0]   , MixEvent_Reader->p1_phi[0]   , MASS_PROTON                   );
			TLorentzVector pion;   pion.SetPtEtaPhiM(   MixEvent_Reader->p2_pt[0]   , MixEvent_Reader->p2_eta[0]   , MixEvent_Reader->p2_phi[0]   , MASS_PION                     );

			Lambda_counterpart->push_back(lambda);
			proton_counterpart->push_back(proton);
			pion_counterpart->push_back(pion);


		}

		//-------------------------------End i_event loop----------------------------
		fin->Close();
		delete fin;
	}
	//--------------------------------End i_file loop------------------------------


}

int ntp_Lambda_Analyzer::Analyze_MEPair(int i_lambda,int j_lambda,int i_file){
		if(SameEvent_Reader->p1_InEventID[i_lambda] == SameEvent_Reader->p1_InEventID[j_lambda] || 
		   SameEvent_Reader->p1_InEventID[i_lambda] == SameEvent_Reader->p2_InEventID[i_lambda] || 
		   SameEvent_Reader->p1_InEventID[i_lambda] == SameEvent_Reader->p2_InEventID[j_lambda] || 
		   SameEvent_Reader->p1_InEventID[j_lambda] == SameEvent_Reader->p2_InEventID[i_lambda] || 
		   SameEvent_Reader->p1_InEventID[j_lambda] == SameEvent_Reader->p2_InEventID[j_lambda] ||
		   SameEvent_Reader->p2_InEventID[i_lambda] == SameEvent_Reader->p2_InEventID[j_lambda] ) return 0;

			

		//---------------------------Event Type Classification-------------------------------
		int Pair_Type = 0 ; 
		int id_Lambda1= i_lambda ; 
		int id_Lambda2= j_lambda ;

		

		Pair_Type = Pair_Type_Classifier(i_lambda,j_lambda);

		if(Pair_Type == 1){
			id_Lambda1 = i_lambda;id_Lambda2 = j_lambda;
		}

		else if(Pair_Type == -1){
			id_Lambda1 = j_lambda;id_Lambda2 = i_lambda;
		}

		else if(Pair_Type == 2){
			id_Lambda1 = i_lambda;id_Lambda2 = j_lambda;
		}
		else if(Pair_Type == 3){
			id_Lambda1 = i_lambda;id_Lambda2 = j_lambda;
		}
		else{
			return 0;
			std::cout<<"Something went wrong with the pair type classifier"<<std::endl;
		}
			
		//---------------------------Create Same-Event Pairs-------------------------------
		TLorentzVector Lambda1; Lambda1.SetPtEtaPhiM( SameEvent_Reader->pair_pt[id_Lambda1], SameEvent_Reader->pair_eta[id_Lambda1], SameEvent_Reader->pair_phi[id_Lambda1], SameEvent_Reader->pair_mass[id_Lambda1]  );
		TLorentzVector proton1; proton1.SetPtEtaPhiM( SameEvent_Reader->p1_pt[id_Lambda1]  , SameEvent_Reader->p1_eta[id_Lambda1]  , SameEvent_Reader->p1_phi[id_Lambda1]  , MASS_PROTON                              );
		TLorentzVector pion1  ;   pion1.SetPtEtaPhiM( SameEvent_Reader->p2_pt[id_Lambda1]  , SameEvent_Reader->p2_eta[id_Lambda1]  , SameEvent_Reader->p2_phi[id_Lambda1]  , MASS_PION                                );

		TLorentzVector Lambda2; Lambda2.SetPtEtaPhiM( SameEvent_Reader->pair_pt[id_Lambda2], SameEvent_Reader->pair_eta[id_Lambda2], SameEvent_Reader->pair_phi[id_Lambda2], SameEvent_Reader->pair_mass[id_Lambda2]  );
		TLorentzVector proton2; proton2.SetPtEtaPhiM( SameEvent_Reader->p1_pt[id_Lambda2]  , SameEvent_Reader->p1_eta[id_Lambda2]  , SameEvent_Reader->p1_phi[id_Lambda2]  , MASS_PROTON                              );
		TLorentzVector pion2;     pion2.SetPtEtaPhiM( SameEvent_Reader->p2_pt[id_Lambda2]  , SameEvent_Reader->p2_eta[id_Lambda2]  , SameEvent_Reader->p2_phi[id_Lambda2]  , MASS_PION                                );
		//---------------------------Range Type Classification-------------------------------
		int Range_Type  = Range_Type_Classifier(&Lambda1,&Lambda2);

		//----------------------------Create Mixed-Event Pairs-------------------------------
		std::vector<TLorentzVector> Lambda1_counterpart;
		std::vector<TLorentzVector> proton1_counterpart;
		std::vector<TLorentzVector> pion1_counterpart;
		std::vector<TLorentzVector> Lambda2_counterpart;
		std::vector<TLorentzVector> proton2_counterpart;
		std::vector<TLorentzVector> pion2_counterpart;
		FindCounterparts(&Lambda2_counterpart,&proton2_counterpart,&pion2_counterpart,Lambda2.Pt(),Lambda2.Rapidity(),Lambda2.Phi(),SameEvent_Reader->p1_ch[id_Lambda2],i_file );
		FindCounterparts(&Lambda1_counterpart,&proton1_counterpart,&pion1_counterpart,Lambda1.Pt(),Lambda1.Rapidity(),Lambda1.Phi(),SameEvent_Reader->p1_ch[id_Lambda1],i_file );

		for(int k_lambda = 0; k_lambda < Lambda2_counterpart.size();k_lambda++){
			//Fill the pair plots 
			Histogramer->Fill_PairPlots(&Lambda1,&Lambda2_counterpart[k_lambda],TMath::Abs(Pair_Type)-1,Range_Type, 0.5/Lambda2_counterpart.size() );
			//Calculate the Density Matrix
			Calculator->Reset(&Lambda1,&proton1,&pion1,&Lambda2_counterpart[k_lambda],&proton2_counterpart[k_lambda],&pion2_counterpart[k_lambda]);
			Calculator->Calculation();
			Histogramer->Fill_DensityMatrix(TMath::Abs(Pair_Type)-1,Range_Type,0.5/Lambda2_counterpart.size() );
		}

		for(int k_lambda = 0; k_lambda < Lambda1_counterpart.size();k_lambda++){
			//Fill the pair plots 
			Histogramer->Fill_PairPlots(&Lambda1_counterpart[k_lambda],&Lambda2,TMath::Abs(Pair_Type)-1,Range_Type, 0.5/Lambda1_counterpart.size() );
			//Calculate the Density Matrix
			Calculator->Reset(&Lambda1_counterpart[k_lambda],&proton1_counterpart[k_lambda],&pion1_counterpart[k_lambda],&Lambda2,&proton2,&pion2);
			Calculator->Calculation();
			Histogramer->Fill_DensityMatrix(TMath::Abs(Pair_Type)-1,Range_Type,0.5/Lambda1_counterpart.size() );
		}
		return 1;

}



void ntp_Lambda_Analyzer::Analysis_MixEvent(){
	//Start looping over all inputfiles of Same_Event_Reader 
	unsigned long N_Inputfiles_SE = SameEvent_Reader->InputFiles.size();
	//-------------------------------Enter i_file loop---------------------------------
	for(unsigned i_file = 0 ;i_file < N_Inputfiles_SE ;i_file++){
		//if(i_file%1==0) std::cout<<"i_file="<<i_file<<std::endl;
		std::cout<<"i_file="<<i_file<<std::endl;
		TFile *fin= TFile::Open(SameEvent_Reader->InputFiles[i_file].c_str(),"READ");
		if(fin==0){
			std::cout<<"Can not open the file:"<< SameEvent_Reader->InputFiles[i_file] <<", Skip this" << std::endl;
			continue;
		}
		TTree *tmp_tree = (TTree *)fin->Get("ntp_Lambda");
		if(tmp_tree==0){
			std::cout<<"Can not get the tree, skip this file."<<std::endl;
			continue;
		} 
		SameEvent_Reader->Init(tmp_tree);
		//Get the number of entries in current tree
		Long64_t N_Events=SameEvent_Reader->fChain->GetEntries();
		//---------------------------Ener i_event loop----------------------------
		for(Long64_t i_event=0;i_event< N_Events;i_event++){
			if(i_event%10000==0)std::cout<<"i_event"<<i_event<<std::endl;
			SameEvent_Reader->fChain->GetEntry(i_event);



			//------------------------Make some selections on the events-----------------------------
			//if(isGoodTrigger() == 1 ) continue; //select on the triggers 
			//if(SameEvent_Reader->NLambda < 3) continue;// current we only select on multi-Lambdas Events
			
			std::vector<int> GoodLambdaFlag;

			for(int i_lambda = 0; i_lambda<SameEvent_Reader->NLambda;i_lambda++){
				int isGoodLambda = 1; //1 is good lambda, 0 is not good lambad 
				if( SameEvent_Reader->pair_charge[i_lambda] == 1 ){ isGoodLambda =0; }
				if( SameEvent_Reader->p1_pt[i_lambda] < Track_Pt_LowCut || SameEvent_Reader->p2_pt[i_lambda] < Track_Pt_LowCut ) {isGoodLambda=0;}
				if( TMath::Abs(SameEvent_Reader->p1_eta[i_lambda]) > Track_Eta_Cut || TMath::Abs(SameEvent_Reader->p2_eta[i_lambda]) > Track_Eta_Cut  ) {isGoodLambda=0;}
				if( SameEvent_Reader->p1_dca[i_lambda] < 0.1 || SameEvent_Reader->p2_dca[i_lambda] < 0.3 ) {isGoodLambda=0;}

				//cut on the pair_DCA 
				if( SameEvent_Reader->pair_DCAdaughters[i_lambda] > 1.) {isGoodLambda=0;}
				//cut on the decay length
				if( SameEvent_Reader->pair_decayL[i_lambda] < 2. || SameEvent_Reader->pair_decayL[i_lambda] > 25.  ){isGoodLambda=0;}
				//cut on the cos_theta
				if( TMath::Cos(SameEvent_Reader->pair_theta[i_lambda]) < 0.996 ){isGoodLambda=0;}
				//cut on DCA of L to PV
				if ( ( TMath::Sin(SameEvent_Reader->pair_theta[i_lambda]) * SameEvent_Reader->pair_decayL[i_lambda] ) > 1.0  ){isGoodLambda=0;}


				TLorentzVector v;
				v.SetPtEtaPhiM(SameEvent_Reader->pair_pt[i_lambda],SameEvent_Reader->pair_eta[i_lambda],SameEvent_Reader->pair_phi[i_lambda],SameEvent_Reader->pair_mass[i_lambda]);
				//cut on the lambda_Candidates 
				if( v.Pt() < Lambda_Pt_lowCut || v.Pt() > Lambda_Pt_highCut || TMath::Abs( v.Rapidity() )  > Lambda_Rapidity_Cut ){isGoodLambda=0;}
				if( v.M()  < Lambda_mass_lowCut || v.M() > Lambda_mass_highCut ) {isGoodLambda=0;}
				GoodLambdaFlag.push_back(isGoodLambda);
			}
			
			int NGoodLambda = std::accumulate(GoodLambdaFlag.begin(), GoodLambdaFlag.end(), 0);
			if(NGoodLambda < 2) continue;

			//Fill Histograms of QA plot;
			Histogramer->Fill_QAplots(GoodLambdaFlag);


			for(int i_lambda=0;i_lambda < SameEvent_Reader->NLambda;i_lambda++){
				if( GoodLambdaFlag[i_lambda] == 0) continue;
				for(int j_lambda=i_lambda+1;j_lambda < SameEvent_Reader->NLambda;j_lambda++){
					if( GoodLambdaFlag[j_lambda] == 0 ) continue;
					Analyze_MEPair(i_lambda,j_lambda,i_file);
				}
			}

		}
		//---------------------------End i_event loop----------------------------
		fin->Close();
		delete fin;

	}
	//-------------------------------Enter i_file loop---------------------------------

}





#endif
