//include the header files here 


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







void ntp_Lambda_Analyzer::Analysis_SameEvent(){

	//Start looping over all inputfiles inthe SameEvent_Reader
	unsigned long N_Inputfiles_SE = SameEvent_Reader->InputFiles.size();
	//-------------------------------Enter i_file loop---------------------------------
	for(unsigned long i_file = 0 ; i_file < N_Inputfiles_SE ; i_file ++){
		//Open the file 
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
		Long64_t N_Events SameEvent_Reader->fChain->GetEntries();

		//---------------------------Enter i_event loop----------------------------
		for(Long64_t i_event=0; i_event < N_Events ; i_event++ ){
			SameEvent_Reader->fChain->GetEntry(i_event);

			
			
			//------------------------Make some selections on the events-----------------------------
			if(SameEvent_Reader->NLambda!=2) continue; // current we only select on two-Lambdas Events 
			if(SameEvent_Reader->pair_ch[0] ==1 || SameEvent_Reader->pair_ch[1] ==1 ) continue; // if any one of the pair is background, pass the event. 
			int isGoodEvent = 1; //1 is good event, 0 is not good event 
			for(int i_lambda = 0; i_lambda<SameEvent_Reader->NLambda;i_lambda++){
				//cut on the track
				if( SameEvent_Reader->p1_pt[i_lambda] < Track_Pt_LowCut || SameEvent_Reader->p2_pt[i_lambda] < Track_Pt_LowCut ) {isGoodEvent=0;break;}
				if( TMath::Abs(SameEvent_Reader->p1_eta[i_lambda]) > Track_Eta_Cut || TMath::Abs(SameEvent_Reader->p2_eta[i_lambda]) > Track_Eta_Cut  ) {isGoodEvent=0;break;}
				
				TLorentzVector v;
				v.SetPtEtaPhiM(SameEvent_Reader->pair_pt[i_lambda],SameEvent_Reader->pair_eta[i_lambda],SameEvent_Reader->pair_phi[i_lambda],SameEvent_Reader->pair_mass[i_lambda]);
				//cut on the lambda_Candidates 
				if( v.Pt() < Lambda_Pt_lowCut || v.Pt() > Lambda_Pt_highCut || TMath::Abs( v.Rapidity() )  > Lambda_Rapidity_Cut ){isGoodEvent=0;break;}
				if( v.M()  < Lambda_mass_lowCut || v.M() > Lambda_mass_highCut ) {isGoodEvent=0;break;}
			}

			if(isGoodEvent==0) continue;
			if(SameEvent_Reader->p1_InEventID[0] == SameEvent_Reader->p1_InEventID[1] || 
			   SameEvent_Reader->p1_InEventID[0] == SameEvent_Reader->p2_InEventID[0] || 
			   SameEvent_Reader->p1_InEventID[0] == SameEvent_Reader->p2_InEventID[1] || 
			   SameEvent_Reader->p1_InEventID[1] == SameEvent_Reader->p2_InEventID[0] || 
			   SameEvent_Reader->p1_InEventID[1] == SameEvent_Reader->p2_InEventID[1] ||
			   SameEvent_Reader->p2_InEventID[0] == SameEvent_Reader->p2_InEventID[1] ) continue;
				
			//Fill Histograms of QA plot;
			Histogramer->FillQA();

			//---------------------------Event Type Classification-------------------------------
			int Pair_Type = 0 ; 
			int id_Lambda1= 0 ; 
			int id_Lambda2= 1 ;

			//Event_flag = 0: Lambda + Lambda_bar 
			//Event_flag = 1: Lambda + Lambda
			//Event_flag = 2: Lambda_bar + Lambda_bar 

			Pair_Type = Pair_Type_Classifier(0,1);

			if(Pair_Type == 1){
				id_Lambda1 = 0;id_Lambda2 = 1;
			}

			else if(Pair_Type == -1){
				id_Lambda1 = 1;id_Lambda2 = 0;
			}

			else if(Pair_Type == 2){
				id_Lambda1 = 0;id_Lambda2 = 1;
			}
			else if(Pair_Type == 3){
				id_Lambda1 = 0;id_Lambda2 = 1;
			}
			else{
				continue;
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
			Fill_PairPlots(&Lambda1,&Lambda2,TMath::Abs(Pair_Type)-1,Range_Type);

			//Calculate the Density Matrix
			Calculator->Reset(&Lambda1,&proton1,&pion1,&Lambda2,&proton2,&pion2);
			Calculator->Calculation();

			//Fill Density Matrix Histgrams 
			Histogramer->FillDensityMatrix(TMath::Abs(Pair_Type)-1,Range_Type);



		}
		//---------------------------End i_event loop----------------------------
		fin->Close();
		delete fin;
	}

	//-------------------------------Enter i_file loop---------------------------------

}




void ntp_Lambda_Analyzer::FindCounterparts(std::vector<TLorentzVector> *Lambda_counterpart,std::vector<TLorentzVector> *proton_counterpart,std::vector<TLorentzVector> *pion_counterpart,double pt, double rapidity, double phi, int p1Charge){
	//Start looping over all inputfiles of Mix_Event_Reader
	unsigned long N_Inputfiles_ME = MixEvent_Reader->InputFiles.size();
	//--------------------------------Enter i_file loop------------------------------
	for(unsigned i_file =0;i_file < N_Inputfiles_ME ; i_file++){
		TFile *fin=TFile::Open(MixEvent_Reader->InputFiles[i_file],c_str(),"READ");
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
		Long64_t N_Events MixEvent_Reader->fChain->GetEntries();

		//-------------------------------Enter i_event loop----------------------------
		for(Long64_t i_event=0;i_event<N_Events;i_event++){
			MixEvent_Reader->fChain->GetEntry(i_event);
			//------------------------Make some selections on the events-----------------------------
			if(MixEvent_Reader->NLambda!=1) continue;
			if(MixEvent_Reader->pair_charge[0]==1) continue; //must be the lambda/lambda_bar candidates
			if(MixEvent_Reader->p1_ch[0]!=p1Charge) continue;

			if( TMath::Abs(MixEvent_Reader->pair_pt[0] - pt  ) > ptDiffLim ) continue;
			if( TMath::ACos(MixEvent_Reader->pair_phi[0]-phi ) > phiDiffLim ) continue;

			//cut on the tracks 
			if( TMath::Abs(MixEvent_Reader->p1_eta[0]) > Track_Eta_Cut || TMath::Abs(MixEvent_Reader->p2_eta) > Track_Eta_Cut ) continue;
			if( MixEvent_Reader->p1_pt[0] < Track_Pt_LowCut || MixEvent_Reader->p2_pt[0] < Track_Pt_LowCut ) continue;
			//cut on the lambda candadates 
			TLorentzVector v;
			v.SetPtEtaPhiM(MixEvent_Reader->pair_pt[0],MixEvent_Reader->pair_eta[0],MixEvent_Reader->pair_phi[0],MixEvent_Reader->pair_mass[0]);
			if( v.Pt() < Lambda_Pt_lowCut || v.Pt() > Lambda_Pt_highCut ) continue;
			if( TMath::Abs( v.Rapidity() ) > Lambda_Rapidity_Cut ) continue;

			if( TMath::Abs(v.Rapidity - rapidity ) > rapidityDiffLim   ) continue;
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

void ntp_Lambda_Analyzer::Analysis_MixEvent(){
	//Start looping over all inputfiles of Same_Event_Reader 
	unsigned long N_Inputfiles_SE = SameEvent_Reader->InputFiles.size();
	//-------------------------------Enter i_file loop---------------------------------
	for(unsigned i_file = 0 ;i_file < N_Inputfiles_SE ;i_file++){
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
		Long64_t N_Events SameEvent_Reader->fChain->GetEntries();
		//---------------------------Ener i_event loop----------------------------
		for(Long64_t i_event=0;i_event< N_Events;i_event++){
			SameEvent_Reader->fChain->GetEntry(i_event);
			//------------------------Make some selections on the events-----------------------------
			if(SameEvent_Reader->NLambda!=2) continue;// current we only select on two-Lambdas Events, will extend to multi-lambda events 
			if(SameEvent_Reader->pair_charge[0] == 1 || SameEvent_Reader->pair_charge[1] == 1 ) continue;
			int isGoodEvent = 1;
			for(int i_lambda = 0; i_lambda<SameEvent_Reader->NLambda;i_lambda++){
				if( SameEvent_Reader->p1_pt[i_lambda] < Track_Pt_LowCut || SameEvent_Reader->p2_pt[i_lambda] < Track_Pt_LowCut ) {isGoodEvent=0;break;}
				if( TMath::Abs(SameEvent_Reader->p1_eta[i_lambda]) > Track_Eta_Cut || TMath::Abs(SameEvent_Reader->p2_eta[i_lambda]) > Track_Eta_Cut  ) {isGoodEvent=0;break;}
				
				TLorentzVector v;
				v.SetPtEtaPhiM(SameEvent_Reader->pair_pt[i_lambda],SameEvent_Reader->pair_eta[i_lambda],SameEvent_Reader->pair_phi[i_lambda],SameEvent_Reader->pair_mass[i_lambda]);
				//cut on the lambda_Candidates 
				if( v.Pt() < Lambda_Pt_lowCut || v.Pt() > Lambda_Pt_highCut || TMath::Abs( v.Rapidity() )  > Lambda_Rapidity_Cut ){isGoodEvent=0;break;}
				if( v.M()  < Lambda_mass_lowCut || v.M() > Lambda_mass_highCut ) {isGoodEvent=0;break;}
			}

			if(isGoodEvent ==0 ) continue;

			if(SameEvent_Reader->p1_InEventID[0] == SameEvent_Reader->p1_InEventID[1] || 
			   SameEvent_Reader->p1_InEventID[0] == SameEvent_Reader->p2_InEventID[0] || 
			   SameEvent_Reader->p1_InEventID[0] == SameEvent_Reader->p2_InEventID[1] || 
			   SameEvent_Reader->p1_InEventID[1] == SameEvent_Reader->p2_InEventID[0] || 
			   SameEvent_Reader->p1_InEventID[1] == SameEvent_Reader->p2_InEventID[1] ||
			   SameEvent_Reader->p2_InEventID[0] == SameEvent_Reader->p2_InEventID[1] ) continue;

			//Fill Histograms of QA plot;
			Histogramer->FillQA();

			//---------------------------Event Type Classification-------------------------------
			int Pair_Type = 0 ; 
			int id_Lambda1= 0 ; 
			int id_Lambda2= 1 ;

			//Event_flag = 0: Lambda + Lambda_bar 
			//Event_flag = 1: Lambda + Lambda
			//Event_flag = 2: Lambda_bar + Lambda_bar 

			Pair_Type = Pair_Type_Classifier(0,1);

			if(Pair_Type == 1){
				id_Lambda1 = 0;id_Lambda2 = 1;
			}

			else if(Pair_Type == -1){
				id_Lambda1 = 1;id_Lambda2 = 0;
			}

			else if(Pair_Type == 2){
				id_Lambda1 = 0;id_Lambda2 = 1;
			}
			else if(Pair_Type == 3){
				id_Lambda1 = 0;id_Lambda2 = 1;
			}
			else{
				continue;
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
			std::vector<TLorentzVector> pion1_conunterpart;
			std::vector<TLorentzVector> Lambda2_counterpart;
			std::vector<TLorentzVector> proton2_counterpart;
			std::vector<TLorentzVector> pion2_counterpart;
			FindCounterparts(&Lambda2_counterpart,&proton2_counterpart,pion2_counterpart,Lambda2.Pt(),Lambda2.Rapidity(),Lambda2.Phi(),SameEvent_Reader->p1_ch[id_Lambda1] );
			FindCounterparts(&Lambda1_counterpart,&proton1_counterpart,pion1_counterpart,Lambda1.Pt(),Lambda1.Rapidity(),Lambda1.Phi(),SameEvent_Reader->p1_ch[id_Lambda2] );

			for(int i_lambda = 0; i_lambda < Lambda2_counterpart.size();i_lambda++){
				//Fill the pair plots 
				Fill_PairPlots(&Lambda1,&Lambda2_counterpart[i_lambda],TMath::Abs(Pair_Type)-1,Range_Type, 0.5/Lambda2_counterpart.size() );
				//Calculate the Density Matrix
				Calculator->Reset(&Lambda1,&proton1,&pion1,&Lambda2_counterpart[i_lambda],&proton2_counterpart[i_lambda],&pion2_counterpart[i_lambda]);
				Calculator->Calculation();
				Histogramer->FillDensityMatrix(TMath::Abs(Pair_Type)-1,Range_Type,0.5/Lambda2_counterpart.size() );
			}

			for(int i_lambda = 0; i_lambda < Lambda1_counterpart.size();i_lambda++){
				//Fill the pair plots 
				Fill_PairPlots(&Lambda1_counterpart,&Lambda2[i_lambda],TMath::Abs(Pair_Type)-1,Range_Type, 0.5/Lambda1_counterpart.size() );
				//Calculate the Density Matrix
				Calculator->Reset(&Lambda1_counterpart[i_lambda],&proton1_counterpart[i_lambda],&pion1_counterpart[i_lambda],&Lambda2,&proton2,&pion2);
				Calculator->Calculation();
				Histogramer->FillDensityMatrix(TMath::Abs(Pair_Type)-1,Range_Type,0.5/Lambda1_counterpart.size() );
			}

		}
		//---------------------------End i_event loop----------------------------
		fin->Close();
		delete fin;

	}
	//-------------------------------Enter i_file loop---------------------------------

}







