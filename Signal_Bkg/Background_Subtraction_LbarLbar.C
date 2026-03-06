#include "constants.h"

#include "Swap.C"

const double sig_low = 1.10;
const double sig_high = 1.13;

double Gauss(double x, double mean,double sigma){
    return TMath::Exp(-0.5*(x-mean)*(x-mean)/ (sigma*sigma) );
}

void plot(TH2D *h2D_SS_Subtracted,TH2D *h2D_SS,TH2D *h2D_SB_BS,std::vector<std::string> str_pair_type,std::vector<std::string> str_pt_bin, std::string str_range_bin,std::string output){
    h2D_SS_Subtracted->GetXaxis()->SetRangeUser(sig_low,sig_high);
    h2D_SS_Subtracted->GetYaxis()->SetRangeUser(sig_low,sig_high);

    h2D_SS->GetXaxis()->SetRangeUser(sig_low,sig_high);
    h2D_SS->GetYaxis()->SetRangeUser(sig_low,sig_high);

    h2D_SB_BS->GetXaxis()->SetRangeUser(sig_low,sig_high);
    h2D_SB_BS->GetYaxis()->SetRangeUser(sig_low,sig_high);
    TH1D *h1D_SS_X = h2D_SS->ProjectionX("h1D_SS_X");
    TH1D *h1D_SS_Y = h2D_SS->ProjectionY("h1D_SS_Y");

    TH1D *h1D_SB_BS_X = h2D_SB_BS->ProjectionX("h1D_SB_BS_X");
    TH1D *h1D_SB_BS_Y = h2D_SB_BS->ProjectionY("h1D_SB_BS_Y");

    TH1D *h1D_SS_Subtracted_X = h2D_SS_Subtracted->ProjectionX("h1D_SS_Subtracted_X");
    TH1D *h1D_SS_Subtracted_Y = h2D_SS_Subtracted->ProjectionY("h1D_SS_Subtracted_Y");


   
    /*
    h1D_SS_X->GetXaxis()->SetRangeUser(sig_low,sig_high);
    h1D_SB_BS_X->GetXaxis()->SetRangeUser(sig_low,sig_high);
    h1D_SS_Subtracted_X->GetXaxis()->SetRangeUser(sig_low,sig_high);

    h1D_SS_Y->GetXaxis()->SetRangeUser(sig_low,sig_high);
    h1D_SB_BS_Y->GetXaxis()->SetRangeUser(sig_low,sig_high);
    h1D_SS_Subtracted_Y->GetXaxis()->SetRangeUser(sig_low,sig_high);
    */

    TLegend* leg = new TLegend(0.15, 0.65, 0.3, 0.85);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);


    TLegend* leg2 = new TLegend(0.15, 0.65, 0.3, 0.85);
    leg2->SetBorderSize(0);
    leg2->SetFillColor(0);

    TLatex* tex = new TLatex();
    tex->SetTextSize(0.035);
    tex->SetNDC(true);

   


    //-------------------------------------------------------------------------------------------------------
    TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
    c1->Divide(2,2);
    c1->cd(1);
    
    h2D_SS->SetTitle("");
    h2D_SS->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");h2D_SS->GetXaxis()->CenterTitle();
    h2D_SS->GetYaxis()->SetTitle("Mass [GeV/c^{2}]");h2D_SS->GetYaxis()->CenterTitle();
    h2D_SS->GetZaxis()->SetTitle("Entries");         h2D_SS->GetZaxis()->CenterTitle();
    h2D_SS->Draw("surf");
    tex->DrawLatex(0.15,0.95,( str_pair_type[0] + "-" + str_pair_type[1] ).c_str() );
    tex->DrawLatex(0.15,0.89,"(US-US)");
    tex->DrawLatex(0.15,0.83,"Signal+Background"    );

    tex->DrawLatex(0.55,0.95, str_pt_bin[0].c_str() );
    tex->DrawLatex(0.55,0.89, str_pt_bin[1].c_str() );
    tex->DrawLatex(0.55,0.83, str_range_bin.c_str() );
    
    c1->cd(2);
    h2D_SB_BS->SetTitle("");
    h2D_SB_BS->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");h2D_SB_BS->GetXaxis()->CenterTitle();
    h2D_SB_BS->GetYaxis()->SetTitle("Mass [GeV/c^{2}]");h2D_SB_BS->GetYaxis()->CenterTitle();
    h2D_SB_BS->GetZaxis()->SetTitle("Entries");         h2D_SB_BS->GetZaxis()->CenterTitle(); 
    h2D_SB_BS->Draw("surf");
    tex->DrawLatex(0.15,0.95,( str_pair_type[0] + "-" + str_pair_type[1] ).c_str() );
    tex->DrawLatex(0.15,0.89,"(US-LS)");
    tex->DrawLatex(0.15,0.83,"Background"    );

    tex->DrawLatex(0.55,0.95, str_pt_bin[0].c_str() );
    tex->DrawLatex(0.55,0.89, str_pt_bin[1].c_str() );
    tex->DrawLatex(0.55,0.83, str_range_bin.c_str() );

    c1->cd(3);
    h1D_SS_X->SetMarkerStyle(20);
    h1D_SS_X->SetMarkerSize(1.0);
    h1D_SS_X->SetMarkerColor(kBlack);
    h1D_SS_X->SetLineColor(kBlack);
    h1D_SS_X->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");h1D_SS_X->GetXaxis()->CenterTitle();
    h1D_SS_X->GetYaxis()->SetTitle("Entries");         h1D_SS_X->GetYaxis()->CenterTitle();
    h1D_SS_X->SetTitle("");
    h1D_SS_X->Draw("E1");
    h1D_SB_BS_X->SetMarkerStyle(20);
    h1D_SB_BS_X->SetMarkerSize(1.0);
    h1D_SB_BS_X->SetMarkerColor(kRed);
    h1D_SB_BS_X->SetLineColor(kRed);
    h1D_SB_BS_X->SetTitle("");
    h1D_SB_BS_X->Draw("same E1");
    tex->DrawLatex(0.15,0.85, ("Projection to "+str_pair_type[0]  ).c_str() );

    tex->DrawLatex(0.70,0.75, str_pt_bin[0].c_str() );
    tex->DrawLatex(0.70,0.69, str_pt_bin[1].c_str() );
    tex->DrawLatex(0.70,0.63, str_range_bin.c_str() );

    leg->AddEntry(h1D_SS_X,"Signal + Background");
    leg->AddEntry(h1D_SB_BS_X,"Background");
    leg->Draw();


    c1->cd(4);
    h1D_SS_Y->SetMarkerStyle(20);
    h1D_SS_Y->SetMarkerSize(1.0);
    h1D_SS_Y->SetMarkerColor(kBlack);
    h1D_SS_Y->SetLineColor(kBlack);
    h1D_SS_Y->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");h1D_SS_Y->GetXaxis()->CenterTitle();
    h1D_SS_Y->GetYaxis()->SetTitle("Entries");         h1D_SS_Y->GetYaxis()->CenterTitle();
    h1D_SS_Y->SetTitle("");
    h1D_SS_Y->Draw("E1");
    h1D_SB_BS_Y->SetMarkerStyle(20);
    h1D_SB_BS_Y->SetMarkerSize(1.0);
    h1D_SB_BS_Y->SetMarkerColor(kRed);
    h1D_SB_BS_Y->SetLineColor(kRed);
    h1D_SB_BS_Y->SetTitle("");
    h1D_SB_BS_Y->Draw("same E1");
    tex->DrawLatex(0.15,0.85, ("Projection to "+str_pair_type[1]  ).c_str() );

    tex->DrawLatex(0.70,0.75, str_pt_bin[0].c_str() );
    tex->DrawLatex(0.70,0.69, str_pt_bin[1].c_str() );
    tex->DrawLatex(0.70,0.63, str_range_bin.c_str() );


    c1->SaveAs((output+"_c1.png").c_str());
    
    //------------------------------------------------------------------------------------------------------
    TCanvas *c2 = new TCanvas("c2","c2",1600,1200);
    c2->Divide(2,2);
    c2->cd(1);
    TF2 *f2 = new TF2("f2","[0]*exp(-0.5*((x-[1])/[2])**2 - 0.5 *((y-[3])/[4])**2)",1.1, 1.13, 1.1, 1.13);
    f2->SetParameters(h2D_SS_Subtracted->GetMaximum(),1.116,0.001,1.116,0.001);
    h2D_SS_Subtracted->Fit(f2,"N");
    double amplitude = f2->GetParameter(0);
    double mean_x    = f2->GetParameter(1);
    double sigma_x   = f2->GetParameter(2);
    double mean_y    = f2->GetParameter(3);
    double sigma_y   = f2->GetParameter(4);
    TF1* f1 = new TF1("f1","[0]*exp(-0.5*((x-[1])/[2])**2)",1.1, 1.13);
    double X_scale = 0 ; 
    double Y_scale = 0 ; 

    for(int i=1;i < h2D_SS->GetNbinsX();i++){
        X_scale = X_scale + Gauss(h2D_SS_Subtracted->GetYaxis()->GetBinCenter(i),mean_y,sigma_y  );
        Y_scale = Y_scale + Gauss(h2D_SS_Subtracted->GetXaxis()->GetBinCenter(i),mean_x,sigma_x  );
    }

    int fs_bin_x1 = h2D_SS_Subtracted->GetXaxis()->FindBin(mean_x-2*sigma_x);
    int fs_bin_x2 = h2D_SS_Subtracted->GetXaxis()->FindBin(mean_x+2*sigma_x);
    int fs_bin_y1 = h2D_SS_Subtracted->GetYaxis()->FindBin(mean_y-2*sigma_y);
    int fs_bin_y2 = h2D_SS_Subtracted->GetYaxis()->FindBin(mean_y+2*sigma_y);


    std::cout<<fs_bin_x1<<std::endl;
    std::cout<<fs_bin_x2<<std::endl;
    std::cout<<fs_bin_y1<<std::endl;
    std::cout<<fs_bin_y2<<std::endl;
    double fs_x1 = h2D_SS_Subtracted->GetXaxis()->GetBinLowEdge(fs_bin_x1);
    double fs_x2 = h2D_SS_Subtracted->GetXaxis()->GetBinUpEdge(fs_bin_x2);
    double fs_y1 = h2D_SS_Subtracted->GetYaxis()->GetBinLowEdge(fs_bin_y1);
    double fs_y2 = h2D_SS_Subtracted->GetYaxis()->GetBinUpEdge(fs_bin_y2);
    
     
     c2->cd(1);
     
     h2D_SS_Subtracted->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");h2D_SS_Subtracted->GetXaxis()->CenterTitle();
     h2D_SS_Subtracted->GetYaxis()->SetTitle("Mass [GeV/c^{2}]");h2D_SS_Subtracted->GetYaxis()->CenterTitle();
     h2D_SS_Subtracted->GetZaxis()->SetTitle("Entries");         h2D_SS_Subtracted->GetZaxis()->CenterTitle();
     h2D_SS_Subtracted->SetTitle("");
     h2D_SS_Subtracted->Draw("surf");
     tex->DrawLatex(0.15,0.95,( str_pair_type[0] + "-" + str_pair_type[1] ).c_str() );
     tex->DrawLatex(0.15,0.89,"Signal");
     tex->DrawLatex(0.70,0.95, str_pt_bin[0].c_str() );
     tex->DrawLatex(0.70,0.89, str_pt_bin[1].c_str() );
     tex->DrawLatex(0.70,0.83, str_range_bin.c_str() );


     
     c2->cd(2);
     f2->SetTitle("");
     f2->Draw("surf");
     tex->DrawLatex(0.15,0.95,( str_pair_type[0] + "-" + str_pair_type[1] ).c_str() );
     tex->DrawLatex(0.15,0.89,"2D Gaussian Fit");
     tex->DrawLatex(0.15,0.83,Form("Amplitude: %f",amplitude));
     tex->DrawLatex(0.15,0.77,Form("mean_x: %f",mean_x));
     tex->DrawLatex(0.15,0.71,Form("sigma_x: %f",sigma_x));
     tex->DrawLatex(0.15,0.65,Form("mean_y: %f",mean_y));
     tex->DrawLatex(0.15,0.59,Form("sigma_y: %f",sigma_y));
     tex->DrawLatex(0.70,0.95, str_pt_bin[0].c_str() );
     tex->DrawLatex(0.70,0.89, str_pt_bin[1].c_str() );
     tex->DrawLatex(0.70,0.83, str_range_bin.c_str() );
     tex->DrawLatex(0.70,0.77, Form("S/B = %.3f",h2D_SS_Subtracted->Integral(fs_bin_x1,fs_bin_x2,fs_bin_y1,fs_bin_y2)/h2D_SB_BS->Integral(fs_bin_x1,fs_bin_x2,fs_bin_y1,fs_bin_y2) ) );
     tex->DrawLatex(0.55,0.71, Form("[%f,%f]#times[%f,%f]",fs_x1,fs_x2,fs_y1,fs_y2));
     c2->cd(3);
     h1D_SS_Subtracted_X->GetXaxis()->SetTitle("Mass [GeV/c^{2}]"); h1D_SS_Subtracted_X->GetXaxis()->CenterTitle();
     h1D_SS_Subtracted_X->GetYaxis()->SetTitle("Entries");          h1D_SS_Subtracted_X->GetYaxis()->CenterTitle(); 
     h1D_SS_Subtracted_X->SetMarkerStyle(20);
     h1D_SS_Subtracted_X->SetMarkerSize(1.0);
     h1D_SS_Subtracted_X->SetMarkerColor(kBlack);
     h1D_SS_Subtracted_X->SetLineColor(kBlack);
     h1D_SS_Subtracted_X->SetTitle("");
     h1D_SS_Subtracted_X->Draw("E1");
     f1->SetParameters(amplitude * X_scale, mean_x,sigma_x);
     f1->Draw("same");
     leg2->AddEntry(h1D_SS_Subtracted_X,"Signal");
     leg2->AddEntry(f1,"Gaussian Fit");
     leg2->Draw();

     tex->DrawLatex(0.15,0.85,("Projection to "+str_pair_type[0]  ).c_str() );
     tex->DrawLatex(0.70,0.75, str_pt_bin[0].c_str() );
     tex->DrawLatex(0.70,0.69, str_pt_bin[1].c_str() );
     tex->DrawLatex(0.70,0.63, str_range_bin.c_str() );
     tex->DrawLatex(0.70,0.57, Form("S/B =%.3f",h1D_SS_Subtracted_X->Integral()/h1D_SB_BS_X->Integral() ) );




     c2->cd(4);
     h1D_SS_Subtracted_Y->GetXaxis()->SetTitle("Mass [GeV/c^{2}]"); h1D_SS_Subtracted_Y->GetXaxis()->CenterTitle();
     h1D_SS_Subtracted_Y->GetYaxis()->SetTitle("Entries");          h1D_SS_Subtracted_Y->GetYaxis()->CenterTitle(); 
     h1D_SS_Subtracted_Y->SetMarkerStyle(20);
     h1D_SS_Subtracted_Y->SetMarkerSize(1.0);
     h1D_SS_Subtracted_Y->SetMarkerColor(kBlack);
     h1D_SS_Subtracted_Y->SetLineColor(kBlack);
     h1D_SS_Subtracted_Y->SetTitle("");
     h1D_SS_Subtracted_Y->Draw("E1");
     f1->SetParameters(amplitude * Y_scale, mean_y,sigma_y);
     f1->Draw("same");

     tex->DrawLatex(0.15,0.85,("Projection to "+str_pair_type[1]  ).c_str() );
     tex->DrawLatex(0.70,0.75, str_pt_bin[0].c_str() );
     tex->DrawLatex(0.70,0.69, str_pt_bin[1].c_str() );
     tex->DrawLatex(0.70,0.63, str_range_bin.c_str() );
     tex->DrawLatex(0.70,0.57, Form("S/B =%.3f",h1D_SS_Subtracted_Y->Integral()/h1D_SB_BS_Y->Integral() ) );
    


     c2->SaveAs((output+"_c2.png").c_str());
    
}

void Background_Subtraction_2D(TH2D *h2D_SS,TH2D *h2D_SB,TH2D *h2D_BS, TH2D *h2D_BB,double p1,double p2,std::vector<std::string> str_pair_type,std::vector<std::string> str_pt_bin, std::string str_range_bin,std::string output){
    int s_bin1 = h2D_SS->GetXaxis()->FindBin(1.11);
    int s_bin2 = h2D_SS->GetXaxis()->FindBin(1.12);

    int bin1= h2D_SS->GetXaxis()->FindBin(1.08);
    int bin2= h2D_SS->GetXaxis()->FindBin(1.09);
    h2D_SB->Scale(1./h2D_SB->Integral()*p2);
    h2D_BS->Scale(1./h2D_BS->Integral()*p1);

    TH2D *h2D_SB_BS = (TH2D *)h2D_SB->Clone("h2D_SB_BS");
    h2D_SB_BS->Add(h2D_BS);

    double scale_factor_1 = h2D_SB_BS->Integral(bin1,bin2,bin1,bin2)/h2D_BB->Integral(bin1,bin2,bin1,bin2);
    h2D_BB->Scale(scale_factor_1);
    //h2D_SB_BS->Add(h2D_BB,-1);

    double scale_factor_2 = h2D_SS->Integral(s_bin1,s_bin2,bin1,bin2)/h2D_SB_BS->Integral(s_bin1,s_bin2,bin1,bin2);
    h2D_SB_BS->Scale(scale_factor_2);
    TH2D *h2D_SS_Subtracted = (TH2D *)h2D_SS->Clone("h2D_SS_Subtracted");
    h2D_SS_Subtracted->Add(h2D_SB_BS,-1);

    plot(h2D_SS_Subtracted, h2D_SS, h2D_SB_BS, str_pair_type, str_pt_bin, str_range_bin, output);
    
}


void Background_Subtraction_LbarLbar(){
    TFile *fin = TFile::Open("QA.root");
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

   double pl[3] = {1./(1.+4.639375),1./(1.+3.368993),1./(1.+4.099671)};
   double plbar[3] = {1./(1.+4.708089),1./(1.+3.220400),1./(1.+4.053884)};
   std::vector<std::string> str_pair_type;
   str_pair_type.push_back("#bar{#Lambda}");
   str_pair_type.push_back("#bar{#Lambda}");
   std::vector<std::string> str_pt_bin;
   str_pt_bin.push_back("0.5 < p_{1,t} < 5.0 GeV/c");
   str_pt_bin.push_back("0.5 < p_{2,t} < 5.0 GeV/c");
   std::string str_range_bin = "Full Range";
   std::string output="Background_Subtraction_LbarLbar/pt2_Bin_4_FullRange_LbarLbar";
   //-------------------------------------
   //lambda lambda_bar 0.5- 5.0 0.5-5.0 full range 
   Background_Subtraction_2D(h2D_USLbar_USLbar_MassDist[4][2],h2D_USLbar_LS_MassDist[4][2],SwapAxesManually(h2D_USLbar_LS_MassDist[4][2]),h2D_LS_LS_MassDist[4][2],plbar[2],plbar[2],str_pair_type,str_pt_bin, str_range_bin,output);
    //-------------------------------------

    //-------------------------------------
   str_range_bin = "Short Range";
   output="Background_Subtraction_LbarLbar/pt2_Bin_4_ShortRange_LbarLbar";
   //lambda lambda_bar 0.5- 5.0 0.5-5.0 short range 
   Background_Subtraction_2D(h2D_USLbar_USLbar_MassDist[4][0],h2D_USLbar_LS_MassDist[4][0],SwapAxesManually(h2D_USLbar_LS_MassDist[4][0]),h2D_LS_LS_MassDist[4][0],plbar[2],plbar[2],str_pair_type,str_pt_bin, str_range_bin,output);
    //-------------------------------------


    //-------------------------------------
   str_range_bin = "Long Range";
   output="Background_Subtraction_LbarLbar/pt2_Bin_4_LongRange_LbarLbar";
   //lambda lambda_bar 0.5- 5.0 0.5-5.0 long range 
   Background_Subtraction_2D(h2D_USLbar_USLbar_MassDist[4][1],h2D_USLbar_LS_MassDist[4][1],SwapAxesManually(h2D_USLbar_LS_MassDist[4][1]),h2D_LS_LS_MassDist[4][1],plbar[2],plbar[2],str_pair_type,str_pt_bin, str_range_bin,output);
    //-------------------------------------

    //-------------------------------------
   //lambda lambda_bar 0.5- 1.5 0.5-1.5 full range 
   str_range_bin = "Full Range";
   str_pt_bin[0] = "0.5 < p_{1,t} < 1.5 GeV/c";
   str_pt_bin[1] = "0.5 < p_{2,t} < 1.5 GeV/c";
   output="Background_Subtraction_LbarLbar/pt2_Bin_0_FullRange_LbarLbar";
   Background_Subtraction_2D(h2D_USLbar_USLbar_MassDist[0][2],h2D_USLbar_LS_MassDist[0][2],SwapAxesManually(h2D_USLbar_LS_MassDist[0][2]),h2D_LS_LS_MassDist[0][2],plbar[0],plbar[0],str_pair_type,str_pt_bin, str_range_bin,output);
   //-------------------------------------

   //-------------------------------------
   //lambda lambda_bar 0.5- 1.5 0.5-1.5 short range 
   str_range_bin = "Short Range";
   output="Background_Subtraction_LbarLbar/pt2_Bin_0_ShortRange_LbarLbar";
   Background_Subtraction_2D(h2D_USLbar_USLbar_MassDist[0][0],h2D_USLbar_LS_MassDist[0][0],SwapAxesManually(h2D_USLbar_LS_MassDist[0][0]),h2D_LS_LS_MassDist[0][0],plbar[0],plbar[0],str_pair_type,str_pt_bin, str_range_bin,output);
   //-------------------------------------

   //-------------------------------------
   //lambda lambda_bar 0.5- 1.5 0.5-1.5 long range 
   str_range_bin = "Long Range";
   output="Background_Subtraction_LbarLbar/pt2_Bin_0_LongRange_LbarLbar";
   Background_Subtraction_2D(h2D_USLbar_USLbar_MassDist[0][1],h2D_USLbar_LS_MassDist[0][1],SwapAxesManually(h2D_USLbar_LS_MassDist[0][1]),h2D_LS_LS_MassDist[0][1],plbar[0],plbar[0],str_pair_type,str_pt_bin, str_range_bin,output);
   //-------------------------------------


     //-------------------------------------
   //lambda lambda_bar 1.5- 5.0 1.5-5.0 full range 
   str_range_bin = "Full Range";
   str_pt_bin[0] = "1.5 < p_{1,t} < 5.0 GeV/c";
   str_pt_bin[1] = "1.5 < p_{2,t} < 5.0 GeV/c";
   output="Background_Subtraction_LbarLbar/pt2_Bin_3_FullRange_LbarLbar";
   Background_Subtraction_2D(h2D_USLbar_USLbar_MassDist[3][2],h2D_USLbar_LS_MassDist[3][2],SwapAxesManually(h2D_USLbar_LS_MassDist[3][2]),h2D_LS_LS_MassDist[3][2],plbar[1],plbar[1],str_pair_type,str_pt_bin, str_range_bin,output);
   //-------------------------------------

   //-------------------------------------
   //lambda lambda_bar 1.5- 5.0 1.5-5.0 short range 
   str_range_bin = "Short Range";
   output="Background_Subtraction_LbarLbar/pt2_Bin_3_ShortRange_LbarLbar";
   Background_Subtraction_2D(h2D_USLbar_USLbar_MassDist[3][0],h2D_USLbar_LS_MassDist[3][0],SwapAxesManually(h2D_USLbar_LS_MassDist[3][0]),h2D_LS_LS_MassDist[3][0],plbar[1],plbar[1],str_pair_type,str_pt_bin, str_range_bin,output);
   //-------------------------------------

   //-------------------------------------
   //lambda lambda_bar 1.5- 5.0 1.5-5.0 long range 
   str_range_bin = "Long Range";
   output="Background_Subtraction_LbarLbar/pt2_Bin_3_LongRange_LbarLbar";
   Background_Subtraction_2D(h2D_USLbar_USLbar_MassDist[3][1],h2D_USLbar_LS_MassDist[3][1],SwapAxesManually(h2D_USLbar_LS_MassDist[3][1]),h2D_LS_LS_MassDist[3][1],plbar[1],plbar[1],str_pair_type,str_pt_bin, str_range_bin,output);
   //-------------------------------------

     //-------------------------------------
   //lambda lambda_bar 1.5- 5.0 0.5-1.5 full range 
   str_range_bin = "Full Range";
   str_pt_bin[0] = "0.5 < p_{1,t} < 1.5 GeV/c";
   str_pt_bin[1] = "1.5 < p_{2,t} < 5.0 GeV/c";
   output="Background_Subtraction_LbarLbar/pt2_Bin_1_FullRange_LbarLbar";
   Background_Subtraction_2D(h2D_USLbar_USLbar_MassDist[1][2],h2D_USLbar_LS_MassDist[1][2],SwapAxesManually(h2D_USLbar_LS_MassDist[2][2]),h2D_LS_LS_MassDist[1][2],plbar[0],plbar[1],str_pair_type,str_pt_bin, str_range_bin,output);
   //-------------------------------------

   //-------------------------------------
   //lambda lambda_bar 1.5- 5.0 1.5-5.0 short range 
   str_range_bin = "Short Range";
   output="Background_Subtraction_LbarLbar/pt2_Bin_1_ShortRange_LbarLbar";
   Background_Subtraction_2D(h2D_USLbar_USLbar_MassDist[1][0],h2D_USLbar_LS_MassDist[1][0],SwapAxesManually(h2D_USLbar_LS_MassDist[2][0]),h2D_LS_LS_MassDist[1][0],plbar[0],plbar[1],str_pair_type,str_pt_bin, str_range_bin,output);
   //-------------------------------------

   //-------------------------------------
   //lambda lambda_bar 1.5- 5.0 1.5-5.0 long range 
   str_range_bin = "Long Range";
   output="Background_Subtraction_LbarLbar/pt2_Bin_1_LongRange_LbarLbar";
   Background_Subtraction_2D(h2D_USLbar_USLbar_MassDist[1][1],h2D_USLbar_LS_MassDist[1][1],SwapAxesManually(h2D_USLbar_LS_MassDist[2][1]),h2D_LS_LS_MassDist[1][1],plbar[0],plbar[1],str_pair_type,str_pt_bin, str_range_bin,output);
   //-------------------------------------


   //-------------------------------------
   //lambda lambda_bar 1.5- 5.0 0.5-1.5 full range 
   str_range_bin = "Full Range";
   str_pt_bin[0] = "1.5 < p_{1,t} < 5.0 GeV/c";
   str_pt_bin[1] = "0.5 < p_{2,t} < 1.5 GeV/c";
   output="Background_Subtraction_LbarLbar/pt2_Bin_2_FullRange_LbarLbar";
   Background_Subtraction_2D(h2D_USLbar_USLbar_MassDist[2][2],h2D_USLbar_LS_MassDist[2][2],SwapAxesManually(h2D_USLbar_LS_MassDist[1][2]),h2D_LS_LS_MassDist[2][2],plbar[1],plbar[0],str_pair_type,str_pt_bin, str_range_bin,output);
   //-------------------------------------

   //-------------------------------------
   //lambda lambda_bar 1.5- 5.0 1.5-5.0 short range 
   str_range_bin = "Short Range";
   output="Background_Subtraction_LbarLbar/pt2_Bin_2_ShortRange_LbarLbar";
   Background_Subtraction_2D(h2D_USLbar_USLbar_MassDist[2][0],h2D_USLbar_LS_MassDist[2][0],SwapAxesManually(h2D_USLbar_LS_MassDist[1][0]),h2D_LS_LS_MassDist[2][0],plbar[1],plbar[0],str_pair_type,str_pt_bin, str_range_bin,output);
   //-------------------------------------

   //-------------------------------------
   //lambda lambda_bar 1.5- 5.0 1.5-5.0 long range 
   str_range_bin = "Long Range";
   output="Background_Subtraction_LbarLbar/pt2_Bin_2_LongRange_LbarLbar";
   Background_Subtraction_2D(h2D_USLbar_USLbar_MassDist[2][1],h2D_USLbar_LS_MassDist[2][1],SwapAxesManually(h2D_USLbar_LS_MassDist[1][1]),h2D_LS_LS_MassDist[2][1],plbar[1],plbar[0],str_pair_type,str_pt_bin, str_range_bin,output);
   //-------------------------------------





}





