

void SplitTree(){
	std::vector<std::string> InPutFileList;
	std::string OutPutFileAddress;

	//number of events per little tree 
	Long64_t NEntries_LT = 1000000;


	for(int i =0;i< 6;i++){
		InPutFileList.push_back( Form("/gpfs01/star/pwg/fliu/LL_Spin_Correlation/2012data/BigTree/Output_%d.root",i+1 ) );
	}
	OutPutFileAddress="/gpfs01/star/pwg/fliu/LL_Spin_Correlation/2012data/LittleTree/";




	//-----------------------------------Enter input file loop-----------------------------------
	for(int ifile=0;ifile<InPutFileList.size();ifile++){
		std::cout<<"ifile"<<ifile<<std::endl;
		TFile *fin=TFile::Open( InPutFileList[ifile].c_str() );
		if(!fin){
			std::cout<<"Can't Open the file:"<<InPutFileList[ifile]<<std::endl;
			continue;
		}

		//get the tree 
		TTree *InTree = (TTree*)fin->Get("ntp_Lambda");
		if(!InTree){
			std::cout<<"Can't get the tree:"<<std::endl;
			continue;
		}

		//Get the total entries of this tree
		Long64_t NEntries =  InTree->GetEntries();
		//calculate the number of little trees to be produced 
		int N_LT = NEntries/NEntries_LT + 1;

		for(int i_ltree = 0 ; i_ltree < N_LT;i_ltree++){
			//The output file for the little tree 
			std::string OutPutFile = OutPutFileAddress + Form("OutPutLittleTree_2012_OutPut_%d_%d.root",ifile+1,i_ltree);
			TFile *fout = TFile::Open(OutPutFile.c_str(),"RECREATE");
			TTree *LittleTree = InTree->CloneTree(0);
			Long64_t i_statr = i_ltree * NEntries_LT;
			Long64_t i_end   = (i_ltree+1) * NEntries_LT ; 

			if(i_ltree == N_LT-1){ i_end = NEntries; }

			for(Long64_t i_ELT = i_statr; i_ELT < i_end ; i_ELT++ ){
				InTree->GetEntry(i_ELT);
				LittleTree->Fill();

			}

			LittleTree->Write();
			fout->Close();

		}

	}









}
