



#include "constants.h"




void plot(){
	TFile *fin = TFile::Open("QA.root");

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



    h1D_TrigId     = (TH1D*)fin->Get("h1D_TrigId");
    h1D_NTrig      = (TH1D*)fin->Get("h1D_NTrig");

    h1D_p1_Pt      = (TH1D*)fin->Get("h1D_p1_Pt");
    h1D_p1_Eta     = (TH1D*)fin->Get("h1D_p1_Eta"); 
    h1D_p1_Phi     = (TH1D*)fin->Get("h1D_p1_Phi"); 

    h1D_p2_Pt      = (TH1D*)fin->Get("h1D_p2_Pt");
    h1D_p2_Eta     = (TH1D*)fin->Get("h1D_p2_Eta");
    h1D_p2_Phi     = (TH1D*)fin->Get("h1D_p2_Phi"); 
    
    h1D_pair_Pt    = (TH1D*)fin->Get("h1D_pair_Pt");
    h1D_pair_Eta   = (TH1D*)fin->Get("h1D_pair_Eta");
    h1D_pair_Phi   = (TH1D*)fin->Get("h1D_pair_Phi");
    h1D_pair_Mass  = (TH1D*)fin->Get("h1D_pair_Mass");

    h1D_pair_DCAdaughters = (TH1D*)fin->Get("h1D_pair_DCAdaughters"); 
    h1D_pair_theta        = (TH1D*)fin->Get("h1D_pair_theta");
    h1D_pair_decayL       = (TH1D*)fin->Get("h1D_pair_decayL"); 
    for(int i_pt=0; i_pt < Pt2_Bin ; i_pt ++){
        for(int i_ran = 0; i_ran < Range_Bin ; i_ran ++){
            h2D_USL_USLbar_MassDist[i_pt][i_ran]     = (TH2D*)fin->Get( Form("h2D_USL_USLbar_MassDist_%d_%d",i_pt,i_ran));
            h2D_USL_USL_MassDist[i_pt][i_ran]        = (TH2D*)fin->Get( Form("h2D_USL_USL_MassDist_%d_%d",i_pt,i_ran)   );
            h2D_USLbar_USLbar_MassDist[i_pt][i_ran]  = (TH2D*)fin->Get( Form("h2D_USLbar_USLbar_MassDist_%d_%d",i_pt,i_ran));

            h2D_USL_LS_MassDist[i_pt][i_ran]         = (TH2D*)fin->Get( Form("h2D_USL_LS_MassDist_%d_%d",i_pt,i_ran)  );
            h2D_USLbar_LS_MassDist[i_pt][i_ran]      = (TH2D*)fin->Get( Form("h2D_USLbar_LS_MassDist_%d_%d",i_pt,i_ran));
            h2D_LS_LS_MassDist[i_pt][i_ran]          = (TH2D*)fin->Get( Form("h2D_LS_LS_MassDist_%d_%d",i_pt,i_ran)    );
        }
    }
    
    for(int i_pt=0; i_pt < Pt1_Bin ; i_pt++){
        h1D_USL_MassDist[i_pt]      = (TH1D*)fin->Get( Form("h1D_USL_MassDist_%d",i_pt) );
        h1D_USLbar_MassDist[i_pt]   = (TH1D*)fin->Get( Form("h1D_USLbar_MassDist_%d",i_pt) );
        h1D_LS_Plus_MassDist[i_pt]  = (TH1D*)fin->Get( Form("h1D_LS_Plus_MassDist_%d",i_pt) );
        h1D_LS_Minus_MassDist[i_pt] = (TH1D*)fin->Get( Form("h1D_LS_Minus_MassDist_%d",i_pt) );
        h1D_LS_MassDist[i_pt]       = (TH1D*)fin->Get( Form("h1D_LS_MassDist_%d",i_pt)  );
    }




  	//p1 
  	TCanvas *c1 = new TCanvas("c1","c1",800,600);
  	h1D_p1_Pt->GetXaxis()->SetRangeUser(0,6);
  	h1D_p1_Pt->GetXaxis()->SetTitle("p_{t}");
  	h1D_p1_Pt->Draw();
  	c1->SaveAs("plot/p1_pt.png");

  	TCanvas *c2 = new TCanvas("c2","c2",800,600);
  	h1D_p1_Phi->GetXaxis()->SetTitle("#phi");
  	h1D_p1_Phi->Draw();
  	c2->SaveAs("plot/p1_phi.png");

  	TCanvas *c3 = new TCanvas("c3","c3",800,600);
  	h1D_p1_Eta->GetXaxis()->SetTitle("#eta");
  	h1D_p1_Eta->Draw();
  	c3->SaveAs("plot/p1_eta.png");


  	//p2
  	TCanvas *c4 = new TCanvas("c4","c4",800,600);
  	h1D_p2_Pt->GetXaxis()->SetRangeUser(0,6);
  	h1D_p2_Pt->GetXaxis()->SetTitle("p_{t}");
  	h1D_p2_Pt->Draw();
  	c4->SaveAs("plot/p2_pt.png");

  	TCanvas *c5 = new TCanvas("c5","c5",800,600);
  	h1D_p2_Phi->GetXaxis()->SetTitle("#phi");
  	h1D_p2_Phi->Draw();
  	c5->SaveAs("plot/p2_phi.png");

  	TCanvas *c6 = new TCanvas("c6","c6",800,600);
  	h1D_p2_Eta->GetXaxis()->SetTitle("#eta");
  	h1D_p2_Eta->Draw();
  	c6->SaveAs("plot/p2_eta.png");

  	//pair 
  	TCanvas *c7 = new TCanvas("c7","c7",800,600);
  	h1D_pair_Pt->GetXaxis()->SetRangeUser(0,6);
  	h1D_pair_Pt->GetXaxis()->SetTitle("p_{t}");
  	h1D_pair_Pt->Draw();
  	c7->SaveAs("plot/pair_pt.png");

  	TCanvas *c8 = new TCanvas("c8","c8",800,600);
  	h1D_pair_Phi->GetXaxis()->SetTitle("#phi");
  	h1D_pair_Phi->Draw();
  	c8->SaveAs("plot/pair_phi.png");

  	TCanvas *c9 = new TCanvas("c9","c9",800,600);
  	h1D_pair_Eta->GetXaxis()->SetTitle("#eta");
  	h1D_pair_Eta->Draw();
  	c9->SaveAs("plot/pair_eta.png");

  	TCanvas *c10 = new TCanvas("c10","c10",800,600);
  	h1D_pair_Mass->GetXaxis()->SetTitle("mass");
  	h1D_pair_Mass->Draw();
  	c10->SaveAs("plot/pair_mass.png");
    /*
  	//number of lambda
  	TCanvas *c11 = new TCanvas("c11","c11",800,600);
  	h1D_NLambda->GetXaxis()->SetTitle("N_{#Lambda}");
  	gPad->SetLogy();
  	h1D_NLambda->Draw();
  	c11->SaveAs("plot/N_Lambda.png");
    */

    TCanvas *c12 = new TCanvas("c12","c12",800,600);
    h1D_TrigId->GetXaxis()->SetTitle("TrigId");
    gPad->SetLogy();
    h1D_TrigId->Draw();
    c12->SaveAs("plot/TrigId.png");



    TCanvas *c13 = new TCanvas("c13","c13",800,600);
    h1D_NTrig->GetXaxis()->SetTitle("NTrig");
    gPad->SetLogy();
    h1D_NTrig->Draw();
    c13->SaveAs("plot/NTrig.png");


    TCanvas *c14 = new TCanvas("c14","c14",800,600);
    h1D_pair_DCAdaughters->GetXaxis()->SetTitle("pair_DCAdaughters");
    gPad->SetLogy(0);
    h1D_pair_DCAdaughters->Draw();
    c14->SaveAs("plot/pair_DCAdaughters.png");


    TCanvas *c15 = new TCanvas("c15","c15",800,600);
    h1D_pair_theta->GetXaxis()->SetTitle("pair_theta");
    gPad->SetLogy(0);
    h1D_pair_theta->Draw();
    c15->SaveAs("plot/pair_theta.png");


    TCanvas *c16 = new TCanvas("c16","c16",800,600);
    h1D_pair_decayL->GetXaxis()->SetTitle("pair_decayL");
    gPad->SetLogy(0);
    h1D_pair_decayL->Draw();
    c16->SaveAs("plot/pair_decayL.png");





  	







}