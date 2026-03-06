#ifndef ntp_Lambda_Calculator_h 
#define ntp_Lambda_Calculator_h
//include the header file hear 
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector3.h>


class ntp_Lambda_Calculator{
public : 
	TLorentzVector Lambda_1;
	TLorentzVector proton_1; 
	TLorentzVector pion_1; 

	TLorentzVector Lambda_2; 
	TLorentzVector proton_2;
	TLorentzVector pion_2;

	double theta_1;
	double phi_1;
	double theta_2; 
	double phi_2;

	double B_a_1; 
	double B_a_2;
	double B_a_3;

	double B_b_1;
	double B_b_2;
	double B_b_3;

	double C_11;
	double C_12;
	double C_13;

	double C_21;
	double C_22;
	double C_23;

	double C_31;
	double C_32;
	double C_33;
	double cos_theta_star_Lab;
	double cos_theta_star_PairRest;
	ntp_Lambda_Calculator();
	ntp_Lambda_Calculator(TLorentzVector *lambda1,TLorentzVector *proton1,TLorentzVector *pion1, TLorentzVector *lambda2,TLorentzVector *proton2, TLorentzVector *pion2);
	
	void Reset();
	void Reset(TLorentzVector *lambda1,TLorentzVector *proton1,TLorentzVector *pion1, TLorentzVector *lambda2,TLorentzVector *proton2, TLorentzVector *pion2);

	void Rotation(TVector3 a_hat,TVector3 b_hat, TVector3 c_hat, TLorentzVector *v);
	void Calculation();

};
ntp_Lambda_Calculator::ntp_Lambda_Calculator(){}
ntp_Lambda_Calculator::ntp_Lambda_Calculator(TLorentzVector *lambda1,TLorentzVector *proton1,TLorentzVector *pion1, TLorentzVector *lambda2,TLorentzVector *proton2, TLorentzVector *pion2)

{
	Lambda_1.SetPxPyPzE( lambda1->Px(),lambda1->Py(),lambda1->Pz(),lambda1->E()  );
	proton_1.SetPxPyPzE( proton1->Px(),proton1->Py(),proton1->Pz(),proton1->E()  );
	pion_1.SetPxPyPzE(     pion1->Px(),  pion1->Py(),  pion1->Pz(),  pion1->E()  );

	Lambda_2.SetPxPyPzE( lambda2->Px(),lambda2->Py(),lambda2->Pz(),lambda2->E()  );
	proton_2.SetPxPyPzE( proton2->Px(),proton2->Py(),proton2->Pz(),proton2->E()  );
	pion_2.SetPxPyPzE(     pion2->Px(),  pion2->Py(),  pion2->Pz(),  pion2->E()  );
}

void ntp_Lambda_Calculator::Reset(){
	theta_1   = 0.;
	phi_1     = 0.;
	theta_2   = 0.; 
	phi_2     = 0.;

	B_a_1     = 0.; 
	B_a_2     = 0.;
	B_a_3     = 0.;

	B_b_1     = 0.;
	B_b_2     = 0.;
	B_b_3     = 0.;

	C_11      = 0.;
	C_12      = 0.;
	C_13      = 0.;

	C_21      = 0.;
	C_22      = 0.;
	C_23      = 0.;

	C_31      = 0.;
	C_32      = 0.;
	C_33      = 0.;
	cos_theta_star_Lab= 0.;
	cos_theta_star_PairRest= 0.;
}	

void ntp_Lambda_Calculator::Reset(TLorentzVector *lambda1,TLorentzVector *proton1,TLorentzVector *pion1, TLorentzVector *lambda2,TLorentzVector *proton2, TLorentzVector *pion2)
{
	Lambda_1.SetPxPyPzE( lambda1->Px(),lambda1->Py(),lambda1->Pz(),lambda1->E()  );
	proton_1.SetPxPyPzE( proton1->Px(),proton1->Py(),proton1->Pz(),proton1->E()  );
	pion_1.SetPxPyPzE(     pion1->Px(),  pion1->Py(),  pion1->Pz(),  pion1->E()  );

	Lambda_2.SetPxPyPzE( lambda2->Px(),lambda2->Py(),lambda2->Pz(),lambda2->E()  );
	proton_2.SetPxPyPzE( proton2->Px(),proton2->Py(),proton2->Pz(),proton2->E()  );
	pion_2.SetPxPyPzE(     pion2->Px(),  pion2->Py(),  pion2->Pz(),  pion2->E()  );


	theta_1   = 0.;
	phi_1     = 0.;
	theta_2   = 0.; 
	phi_2     = 0.;

	B_a_1     = 0.; 
	B_a_2     = 0.;
	B_a_3     = 0.;

	B_b_1     = 0.;
	B_b_2     = 0.;
	B_b_3     = 0.;

	C_11      = 0.;
	C_12      = 0.;
	C_13      = 0.;

	C_21      = 0.;
	C_22      = 0.;
	C_23      = 0.;

	C_31      = 0.;
	C_32      = 0.;
	C_33      = 0.;
	cos_theta_star_Lab= 0.;
	cos_theta_star_PairRest= 0.;
}


void ntp_Lambda_Calculator::Rotation(TVector3 a_hat,TVector3 b_hat, TVector3 c_hat, TLorentzVector *v){
   double x1 = v->Vect().Dot(a_hat); // x1 coordiante in the helicity frame a-b-c
   double x2 = v->Vect().Dot(b_hat);
   double x3 = v->Vect().Dot(c_hat);

   v->SetX(x1);
   v->SetY(x2);
   v->SetZ(x3);

}



void ntp_Lambda_Calculator::Calculation(){
	 //calculate the cos_theta_star_lab first ()
	 TLorentzVector proton_1_copy = proton_1;
	 TLorentzVector proton_2_copy = proton_2;

	 proton_1_copy.Boost( -Lambda_1.BoostVector()  );
	 proton_2_copy.Boost( -Lambda_2.BoostVector()  );

	 cos_theta_star_PairRest = proton_1_copy.Vect().Dot( proton_2_copy.Vect() )/( proton_1_copy.Vect().Mag() * proton_2_copy.Vect().Mag()  );


	 TLorentzVector LL_bar = Lambda_1 + Lambda_2;
    // construc the beta vector of LL_bar
    TVector3 beta_LL_bar = LL_bar.BoostVector();
    // Boost everything into the rest frame of LL_bar 
    LL_bar.Boost(  -beta_LL_bar      );

    Lambda_1.Boost( -beta_LL_bar      );
    proton_1.Boost( -beta_LL_bar      );
    pion_1.Boost(   -beta_LL_bar      );

    Lambda_2.Boost( -beta_LL_bar  );
    proton_2.Boost( -beta_LL_bar  );
    pion_2. Boost(  -beta_LL_bar  );
    // Rotate into the helicity frame

    //construct the helicity coordinate axis 

    TVector3 p_beam = TVector3(0,0,1);
    TVector3 Lambda_1_Vect =  Lambda_1.Vect();
    TVector3 c_hat = Lambda_1_Vect.Unit();
    TVector3 b_hat = p_beam.Cross(c_hat).Unit();
    TVector3 a_hat = b_hat.Cross(c_hat);
    //Rotate Everything 
    Rotation(a_hat,b_hat,c_hat,&Lambda_1);
    Rotation(a_hat,b_hat,c_hat,&proton_1);
    Rotation(a_hat,b_hat,c_hat,&pion_1);

    Rotation(a_hat,b_hat,c_hat,&Lambda_2);
    Rotation(a_hat,b_hat,c_hat,&proton_2);
    Rotation(a_hat,b_hat,c_hat,&pion_2 );

    // consruct the beta vector of Lambda and Lambda_bar 
    TVector3 beta_Lambda_1 = Lambda_1.BoostVector();
    TVector3 beta_Lambda_2 = Lambda_2.BoostVector();

    // Boost the daughter particles into the mother's rest frame, 
    proton_1.Boost(-beta_Lambda_1);
    pion_1.Boost(-beta_Lambda_1);

    proton_2.Boost(-beta_Lambda_2);
    pion_2.Boost(-beta_Lambda_2);



    // calculate the theta_1(2), phi_1(2) 
    theta_1 = proton_1.Theta();
    phi_1   = proton_1.Phi();

    theta_2 = proton_2.Theta();
    phi_2   = proton_2.Phi();

    B_a_1 = TMath::Sin(theta_1) * TMath::Cos(phi_1);
    B_a_2 = TMath::Sin(theta_1) * TMath::Sin(phi_1);
    B_a_3 = TMath::Cos(theta_1); 
    B_b_1 = TMath::Sin(theta_2) * TMath::Cos(phi_2);
    B_b_2 = TMath::Sin(theta_2) * TMath::Sin(phi_2);
    B_b_3 = TMath::Cos(theta_2);

    C_11 = TMath::Sin(theta_1) * TMath::Cos(phi_1) * TMath::Sin(theta_2) * TMath::Cos(phi_2);
    C_12 = TMath::Sin(theta_1) * TMath::Cos(phi_1) * TMath::Sin(theta_2) * TMath::Sin(phi_2);
    C_13 = TMath::Sin(theta_1) * TMath::Cos(phi_1) * TMath::Cos(theta_2);
         
    C_21 = TMath::Sin(theta_1) * TMath::Sin(phi_1) * TMath::Sin(theta_2) * TMath::Cos(phi_2);
    C_22 = TMath::Sin(theta_1) * TMath::Sin(phi_1) * TMath::Sin(theta_2) * TMath::Sin(phi_2);
    C_23 = TMath::Sin(theta_1) * TMath::Sin(phi_1) * TMath::Cos(theta_2);
         
    C_31 = TMath::Cos(theta_1) * TMath::Sin(theta_2) * TMath::Cos(phi_2);
    C_32 = TMath::Cos(theta_1) * TMath::Sin(theta_2) * TMath::Sin(phi_2);
    C_33 = TMath::Cos(theta_1) * TMath::Cos(theta_2);

    cos_theta_star_Lab = TMath::Sin(theta_1) * TMath::Cos(phi_1) * TMath::Sin(theta_2) * TMath::Cos(phi_2) + TMath::Sin(theta_1) * TMath::Sin(phi_1) * TMath::Sin(theta_2) * TMath::Sin(phi_2) + TMath::Cos(theta_1) * TMath::Cos(theta_2);    


}


#endif