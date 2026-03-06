








void plot(){
	TFile *fin = TFile::Open("QA.root");
	TH1D *h1D_NLambda = (TH1D*)fin->Get("h1D_NLambda");
   	//daughter particle 1
    TH1D *h1D_p1_pt   = (TH1D*)fin->Get("h1D_p1_pt" );
   	TH1D *h1D_p1_phi  = (TH1D*)fin->Get("h1D_p1_phi");
    TH1D *h1D_p1_eta  = (TH1D*)fin->Get("h1D_p1_eta");
   	//daughter particle 2 
   	TH1D *h1D_p2_pt   = (TH1D*)fin->Get("h1D_p2_pt" );
    TH1D *h1D_p2_phi  = (TH1D*)fin->Get("h1D_p2_phi");
  	TH1D *h1D_p2_eta  = (TH1D*)fin->Get("h1D_p2_eta");

   //pair 
   	TH1D *h1D_pair_phi  = (TH1D*)fin->Get("h1D_pair_phi"  );
   	TH1D *h1D_pair_eta  = (TH1D*)fin->Get("h1D_pair_eta"  );
   	TH1D *h1D_pair_pt   = (TH1D*)fin->Get("h1D_pair_pt"   );
  	TH1D *h1D_pair_mass = (TH1D*)fin->Get("h1D_pair_mass" ); 
  	//p1 
  	TCanvas *c1 = new TCanvas("c1","c1",800,600);
  	h1D_p1_pt->GetXaxis()->SetRangeUser(0,6);
  	h1D_p1_pt->GetXaxis()->SetTitle("p_{t}");
  	h1D_p1_pt->Draw();
  	c1->SaveAs("plot/p1_pt.png");

  	TCanvas *c2 = new TCanvas("c2","c2",800,600);
  	h1D_p1_phi->GetXaxis()->SetTitle("#phi");
  	h1D_p1_phi->Draw();
  	c2->SaveAs("plot/p1_phi.png");

  	TCanvas *c3 = new TCanvas("c3","c3",800,600);
  	h1D_p1_eta->GetXaxis()->SetTitle("#eta");
  	h1D_p1_eta->Draw();
  	c3->SaveAs("plot/p1_eta.png");


  	//p2
  	TCanvas *c4 = new TCanvas("c4","c4",800,600);
  	h1D_p2_pt->GetXaxis()->SetRangeUser(0,6);
  	h1D_p2_pt->GetXaxis()->SetTitle("p_{t}");
  	h1D_p2_pt->Draw();
  	c4->SaveAs("plot/p2_pt.png");

  	TCanvas *c5 = new TCanvas("c5","c5",800,600);
  	h1D_p2_phi->GetXaxis()->SetTitle("#phi");
  	h1D_p2_phi->Draw();
  	c5->SaveAs("plot/p2_phi.png");

  	TCanvas *c6 = new TCanvas("c6","c6",800,600);
  	h1D_p2_eta->GetXaxis()->SetTitle("#eta");
  	h1D_p2_eta->Draw();
  	c6->SaveAs("plot/p2_eta.png");

  	//pair 
  	TCanvas *c7 = new TCanvas("c7","c7",800,600);
  	h1D_pair_pt->GetXaxis()->SetRangeUser(0,6);
  	h1D_pair_pt->GetXaxis()->SetTitle("p_{t}");
  	h1D_pair_pt->Draw();
  	c7->SaveAs("plot/pair_pt.png");

  	TCanvas *c8 = new TCanvas("c8","c8",800,600);
  	h1D_pair_phi->GetXaxis()->SetTitle("#phi");
  	h1D_pair_phi->Draw();
  	c8->SaveAs("plot/pair_phi.png");

  	TCanvas *c9 = new TCanvas("c9","c9",800,600);
  	h1D_pair_eta->GetXaxis()->SetTitle("#eta");
  	h1D_pair_eta->Draw();
  	c9->SaveAs("plot/pair_eta.png");

  	TCanvas *c10 = new TCanvas("c10","c10",800,600);
  	h1D_pair_mass->GetXaxis()->SetTitle("mass");
  	h1D_pair_mass->Draw();
  	c10->SaveAs("plot/pair_mass.png");

  	//number of lambda
  	TCanvas *c11 = new TCanvas("c11","c11",800,600);
  	h1D_NLambda->GetXaxis()->SetTitle("N_{#Lambda}");
  	gPad->SetLogy();
  	h1D_NLambda->Draw();
  	c11->SaveAs("plot/N_Lambda.png");





  	







}