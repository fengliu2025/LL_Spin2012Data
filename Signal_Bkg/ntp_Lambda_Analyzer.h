#ifndef ntp_Lambda_Analyzer_h 
#define ntp_Lambda_Analyzer_h

//include the header files here 
#include <TFile.h>
#include <TLorentzVector.h>
#include <TMath.h>

#include <iostream>

#include "ntp_Lambda_Reader.h"
#include "ntp_Lambda_Histogram.h"

class ntp_Lambda_Analyzer{
public:
	ntp_Lambda_Reader *Reader;
	ntp_Lambda_Histogram  *Histogramer; 


	ntp_Lambda_Analyzer();
	ntp_Lambda_Analyzer(ntp_Lambda_Reader *reader1, ntp_Lambda_Histogram *hisogram);

	virtual ~ntp_Lambda_Analyzer();
	virtual void Analysis_MassSpectrum();


};

ntp_Lambda_Analyzer::~ntp_Lambda_Analyzer(){

}
ntp_Lambda_Analyzer::ntp_Lambda_Analyzer(){

}

ntp_Lambda_Analyzer::ntp_Lambda_Analyzer(ntp_Lambda_Reader *reader1, ntp_Lambda_Histogram *hisogram)
:Reader(reader1),Histogramer(hisogram){

}


void ntp_Lambda_Analyzer::Analysis_MassSpectrum(){


	//Start looping over all inputfiles inthe SameEvent_Reader
	unsigned long N_Inputfiles_SE = Reader->InputFiles.size();

	//-------------------------------Enter i_file loop---------------------------------
	for(unsigned long i_file = 0 ; i_file < N_Inputfiles_SE ; i_file ++){
		if(i_file%100==0)std::cout<<"i_file="<<i_file<<std::endl;
		//Open the file 
		TFile *fin=TFile::Open( (Reader->InputFiles[i_file]).c_str() ) ;
		
		if(fin==0){
			std::cout<<"Can not open the file:"<< Reader->InputFiles[i_file] <<", Skip this" << std::endl;
			continue;
		}
		TTree *tmp_tree = (TTree *)fin->Get("ntp_Lambda");
		if(tmp_tree==0){
			std::cout<<"Can not get the tree, skip this file."<<std::endl;
			continue;
		}

		Reader->Init(tmp_tree);
		//Get the number of entries in current tree
		Long64_t N_Events=Reader->fChain->GetEntries();

		//---------------------------Ener i_event loop----------------------------
		for(Long64_t i_event=0; i_event < N_Events ; i_event++ ){
			Reader->fChain->GetEntry(i_event);

			//Make some selections on the events
			if(Reader->NLambda!=2) continue; // currently we only select on two-Lambdas Events 
			int isGoodEvent = 1;// 1 is good event, 0 is not good event. 
			//Cuts on track and Lambda candidates
			for(int i_lambda=0 ; i_lambda < Reader->NLambda; i_lambda++){
				if(Reader->p1_pt[i_lambda] < Track_Pt_LowCut || Reader->p2_pt[i_lambda] < Track_Pt_LowCut ) {isGoodEvent=0;break;}
				if( TMath::Abs( Reader->p1_eta[i_lambda] ) > Track_Eta_Cut ||  TMath::Abs( Reader->p2_eta[i_lambda] ) > Track_Eta_Cut ){isGoodEvent=0;break;}

				TLorentzVector v;
				v.SetPtEtaPhiM(Reader->pair_pt[i_lambda],Reader->pair_eta[i_lambda],Reader->pair_phi[i_lambda],Reader->pair_mass[i_lambda]);
				if( v.Pt() < Lambda_Pt_lowCut  || v.Pt() > Lambda_Pt_highCut || TMath::Abs( v.Rapidity() ) > Lambda_Rapidity_Cut ) {isGoodEvent=0;break;} 

			}
			if(Reader->p1_InEventID[0] == Reader->p1_InEventID[1] || 
			   Reader->p1_InEventID[0] == Reader->p2_InEventID[0] || 
			   Reader->p1_InEventID[0] == Reader->p2_InEventID[1] || 
			   Reader->p1_InEventID[1] == Reader->p2_InEventID[0] || 
			   Reader->p1_InEventID[1] == Reader->p2_InEventID[1] ||
			   Reader->p2_InEventID[0] == Reader->p2_InEventID[1] ) continue;
			if( isGoodEvent ==0 ) continue;

			//Fill Histograms of QA plot;
			Histogramer->Fill_QAplots();

			//
			Histogramer->Fill_Mass();

		}
		//---------------------------End i_event loop----------------------------
		fin->Close();
		delete fin;

		

	}

	//-------------------------------Enter i_file loop---------------------------------



}

#endif
