

class new_Tree{
public:

   TTree          *tree;	

   Int_t           eventId;
   Float_t         Vz;
   //Int_t           mNTrigs;
   //Int_t           mTrigId[3];   //[mNTrigs]
   //Int_t           mNTrks;
   //Float_t         high_pt[6];   //[mNTrks]
   //Float_t         high_phi[6];   //[mNTrks]
   //Float_t         high_eta[6];   //[mNTrks]
   Int_t           NLambda;
   Int_t           p1_InEventID[1000];   //[NLambda]
   Float_t         p1_pt[1000];   //[NLambda]
   Float_t         p1_phi[1000];   //[NLambda]
   Float_t         p1_eta[1000];   //[NLambda]
   Float_t         p1_dca[1000];   //[NLambda]
   Int_t           p1_ch[1000];   //[NLambda]
   Int_t           p1_hasTOFinfo[1000];   //[NLambda]
  
   Int_t           p2_InEventID[1000];   //[NLambda]
   Float_t         p2_pt[1000];   //[NLambda]
   Float_t         p2_phi[1000];   //[NLambda]
   Float_t         p2_eta[1000];   //[NLambda]
   Float_t         p2_dca[1000];   //[NLambda]
   Int_t           p2_hasTOFinfo[1000];   //[NLambda]
  
   Int_t           pair_charge[1000];   //[NLambda]
   Float_t         pair_DCAdaughters[1000];   //[NLambda]
   Float_t         pair_theta[1000];   //[NLambda]
   Float_t         pair_decayL[1000];   //[NLambda]
   Float_t         pair_phi[1000];   //[NLambda]
   Float_t         pair_eta[1000];   //[NLambda]
   Float_t         pair_pt[1000];   //[NLambda]
   Float_t         pair_mass[1000];   //[NLambda]

   new_Tree();
   void MakeNewTree(TTree *OutPutTree);
   void ResetTree();
   void Fill();
   void WriteTree(std::string OutPutFile);
   
};

new_Tree::new_Tree(){
tree =  0;	
}

void new_Tree::MakeNewTree(TTree *OutPutTree){
	tree = OutPutTree;
	if(!tree){
		std::cout<<"the outpull tree doesn't exist"<<std::endl;
		return ;
	}
	
	tree->Branch("eventId", &eventId, "eventId/I");       //Int_t eventId
   tree->Branch("Vz", &Vz, "Vz/F"); //VzVzVPDmax
    
    /*
    tree->Branch("high_pt", high_pt, "high_pt[mNTrks]/F");  
    tree->Branch("high_phi", high_phi, "high_phi[mNTrks]/F");
    tree->Branch("high_eta", high_eta, "high_eta[mNTrks]/F");
	*/
    /*
    tree->Branch("lead_pt", &lead_pt, "lead_pt/F");               //Float_t p1_pt
    tree->Branch("lead_phi", &lead_phi, "lead_phi/F");             //Float_t p1_phi
    tree->Branch("lead_eta", &lead_eta, "lead_eta/F");             //Float_t p1_eta
    
    tree->Branch("sublead_pt", &sublead_pt, "sublead_pt/F");               //Float_t p1_pt
    tree->Branch("sublead_phi", &sublead_phi, "sublead_phi/F");             //Float_t p1_phi
    tree->Branch("sublead_eta", &sublead_eta, "sublead_eta/F");             //Float_t p1_eta
    */
    tree->Branch("NLambda", &NLambda, "NLambda/I");
    
    //proton
    tree->Branch("p1_InEventID", p1_InEventID, "p1_InEventID[NLambda]/I");               //Float_t p1_InEventID
    tree->Branch("p1_pt", p1_pt, "p1_pt[NLambda]/F");               //Float_t p1_pt
    tree->Branch("p1_phi", p1_phi, "p1_phi[NLambda]/F");             //Float_t p1_phi
    tree->Branch("p1_eta", p1_eta, "p1_eta[NLambda]/F");             //Float_t p1_eta
    tree->Branch("p1_dca", p1_dca, "p1_dca[NLambda]/F");               //Float_t p1_dca
    tree->Branch("p1_ch", p1_ch, "p1_ch[NLambda]/I");               //Float_t p1_ch
    tree->Branch("p1_hasTOFinfo", p1_hasTOFinfo, "p1_hasTOFinfo[NLambda]/I");   //Float_t p1_hasTOFinfo
    //tree->Branch("p1_dedx", p1_dedx, "p1_dedx[NLambda]/F");
    //tree->Branch("p1_beta", p1_beta, "p1_beta[NLambda]/F");

    //pion
    tree->Branch("p2_InEventID", p2_InEventID, "p2_InEventID[NLambda]/I");               //Float_t p2_InEventID
    tree->Branch("p2_pt", p2_pt, "p2_pt[NLambda]/F");               //Float_t p2_pt
    tree->Branch("p2_phi", p2_phi, "p2_phi[NLambda]/F");             //Float_t p2_phi
    tree->Branch("p2_eta", p2_eta, "p2_eta[NLambda]/F");             //Float_t p2_eta
    tree->Branch("p2_dca", p2_dca, "p2_dca[NLambda]/F");               //Float_t p2_dca
    //tree->Branch("p2_ch", p2_ch, "p2_ch[NLambda]/I");               //Float_t p2_ch
    tree->Branch("p2_hasTOFinfo", p2_hasTOFinfo, "p2_hasTOFinfo[NLambda]/I");   //Float_t p2_hasTOFinfo
    //tree->Branch("p2_dedx", p2_dedx, "p2_dedx[NLambda]/F");
    //tree->Branch("p2_beta", p2_beta, "p2_beta[NLambda]/F");

    //pair
    tree->Branch("pair_charge", pair_charge, "pair_charge[NLambda]/I");      //Int_t charge
    tree->Branch("pair_DCAdaughters", pair_DCAdaughters, "pair_DCAdaughters[NLambda]/F");      //Int_t pair_DCAdaughters
    tree->Branch("pair_theta", pair_theta, "pair_theta[NLambda]/F");   //Float_t pair_theta
    //tree->Branch("pair_cosThetaStar", pair_cosThetaStar, "pair_cosThetaStar[NLambda]/F");   //Float_t pair_cosThetaStar
    tree->Branch("pair_decayL", pair_decayL, "pair_decayL[NLambda]/F");  //Float_t pair_decayL
    tree->Branch("pair_phi", pair_phi, "pair_phi[NLambda]/F");      //Float_t pair_phi
    tree->Branch("pair_eta", pair_eta, "pair_eta[NLambda]/F");      //Float_t pair_eta
    tree->Branch("pair_pt", pair_pt, "pair_pt[NLambda]/F");        //Float_t pair_pt
    tree->Branch("pair_mass", pair_mass, "pair_mass[NLambda]/F");    //Float_t pair_mass
	//...
}

void new_Tree::ResetTree(){
	tree->Reset();
}

void new_Tree::Fill(){
	tree->Fill();
}	



void new_Tree::WriteTree(std::string OutPutFile){
	


}

