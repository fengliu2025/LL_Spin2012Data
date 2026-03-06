#include "constants.h"

#include "Swap.C"

void Background_Subtraction(TH1D *h_sigbkg,TH1D *h_background,std::string str_Pt_range,std::string str_lambda_flag,int pt_bin,int lambda_flag){
    //background range: 1.09 < M < 1.10;
    //signal range: 1.11 < M <1.12 
    int s_bin1 = h_sigbkg->FindBin(1.11);
    int s_bin2 = h_sigbkg->FindBin(1.12);

    int bin1= h_sigbkg->FindBin(1.09);
    int bin2= h_sigbkg->FindBin(1.10);


    double scale_factor = h_sigbkg->Integral(bin1,bin2)/ h_background->Integral(bin1,bin2);
    h_background->Scale(scale_factor);
    TH1D *h1D_signal = (TH1D*)h_sigbkg->Clone();
    h1D_signal->Add(h_background,-1.0);


    double S_over_SB_FullMass = h1D_signal->Integral()/h_sigbkg->Integral() ;
    double S_over_B_FullMass  = h1D_signal->Integral()/h_background->Integral();

    double S_over_SB_SigBand  = h1D_signal->Integral(s_bin1,s_bin2)/h_sigbkg->Integral(s_bin1,s_bin2);
    double S_over_B_SigBnad   = h1D_signal->Integral(s_bin1,s_bin2)/h_background->Integral(s_bin1,s_bin2);
    
    TCanvas *c1 = new TCanvas("c1","c1",1600,500);
    TLegend* leg = new TLegend(0.55, 0.8, 0.75, 0.9);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);

    TLatex* tex = new TLatex();
    tex->SetTextSize(0.035);
    tex->SetNDC(true);


    c1->Divide(2,1);
    c1->cd(1);
    h_sigbkg->SetLineColor(kRed);
    h_sigbkg->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
    h_sigbkg->GetXaxis()->CenterTitle();
    h_sigbkg->GetYaxis()->SetTitle("Entries");
    h_sigbkg->GetYaxis()->CenterTitle();
    h_sigbkg->Draw();


    h_background->SetLineColor(kBlue);
    h_background->Draw("same");

    h1D_signal->SetLineColor(kGreen);
    h1D_signal->Draw("same");

    leg->AddEntry(h_sigbkg,"Signal + Background");
    leg->AddEntry(h_background,"Background");
    leg->AddEntry(h1D_signal,"Signal");
    leg->Draw();

    tex->DrawLatex(0.15, 0.85, str_Pt_range.c_str());
    tex->DrawLatex(0.15, 0.75, str_lambda_flag.c_str());
    tex->DrawLatex(0.45, 0.55, Form("S/B=%f, full mass range.",S_over_B_FullMass));
    tex->DrawLatex(0.45, 0.45, Form("S/B=%f, signal band range.",S_over_B_SigBnad));

    c1->cd(2);

    double max_content = h_sigbkg->GetMaximum();
    TH1D *h_sigbkg_clone = (TH1D*)h_sigbkg->Clone();
    h_sigbkg_clone->GetYaxis()->SetRangeUser(0,max_content/20.);
    h_sigbkg_clone->Draw();
    h_background->Draw("same");
    h1D_signal->Draw("same");

    c1->SaveAs(Form("Mass_plot/pt_bin_%d_lambda_flag_%d.pdf",pt_bin,lambda_flag));
    

}

double Gauss(double x, double mean,double sigma){
    return TMath::Exp(-0.5*(x-mean)*(x-mean)/ (sigma*sigma) );
}

void Background_Subtraction_2D(TH2D *h2D_SS,TH2D *h2D_SB,TH2D *h2D_BS, TH2D *h2D_BB){
    int s_bin1 = h2D_SS->GetXaxis()->FindBin(1.11);
    int s_bin2 = h2D_SS->GetXaxis()->FindBin(1.12);

    int bin1= h2D_SS->GetXaxis()->FindBin(1.08);
    int bin2= h2D_SS->GetXaxis()->FindBin(1.09);
    h2D_SB->Scale(1./h2D_SB->Integral());
    h2D_BS->Scale(1./h2D_BS->Integral());

    TH2D *h2D_SB_BS = (TH2D *)h2D_SB->Clone();
    h2D_SB_BS->Add(h2D_BS);

    double scale_factor_1 = h2D_SB_BS->Integral(bin1,bin2,bin1,bin2)/h2D_BB->Integral(bin1,bin2,bin1,bin2);
    h2D_BB->Scale(scale_factor_1);
    h2D_SB_BS->Add(h2D_BB,-1);

    double scale_factor_2 = h2D_SS->Integral(s_bin1,s_bin2,bin1,bin2)/h2D_SB_BS->Integral(s_bin1,s_bin2,bin1,bin2);
    h2D_SB_BS->Scale(scale_factor_2);
    //h2D_SS->Add(h2D_SB_BS,-1);

    h2D_SS->GetXaxis()->SetRangeUser(1.1,1.13);
    h2D_SS->GetYaxis()->SetRangeUser(1.1,1.13);

    h2D_SB_BS->GetXaxis()->SetRangeUser(1.1,1.13);
    h2D_SB_BS->GetYaxis()->SetRangeUser(1.1,1.13);

    TH1D *h1D_SS_X = h2D_SS->ProjectionX();
    TH1D *h1D_SS_Y = h2D_SS->ProjectionY();
    TH1D *h1D_SB_BS_X = h2D_SB_BS->ProjectionX();
    TH1D *h1D_SB_BS_Y = h2D_SB_BS->ProjectionY();
    
    TF2* f2 = new TF2("f2", 
        "[0]*exp(-0.5*((x-[1])/[2])**2 - 0.5*((y-[3])/[4])**2)",
        1.1, 1.13, 1.1, 1.13);
    
    
    f2->SetParameters(
        h2D_SS->GetMaximum(),  
        1.116, 0.001,          
        1.116, 0.001           
    );


    
    //h2D_SS->Fit(f2);

    // get the result of fit 
    double amplitude = f2->GetParameter(0);
    double mean_x    = f2->GetParameter(1);
    double sigma_x   = f2->GetParameter(2);
    double mean_y    = f2->GetParameter(3);
    double sigma_y   = f2->GetParameter(4);

    std::cout << "=== 2D Gaussian Fit Results ===" << std::endl;
    std::cout << "Amplitude: " << amplitude << std::endl;
    std::cout << "Mean X: " << mean_x << " ± " << f2->GetParError(1) << std::endl;
    std::cout << "Sigma X: " << sigma_x << " ± " << f2->GetParError(2) << std::endl;
    std::cout << "Mean Y: " << mean_y << " ± " << f2->GetParError(3) << std::endl;
    std::cout << "Sigma Y: " << sigma_y << " ± " << f2->GetParError(4) << std::endl;
    std::cout << "Chi2/NDF: " << f2->GetChisquare() << "/" << f2->GetNDF() << " = " << f2->GetChisquare()/f2->GetNDF() << std::endl;

    double X_scale = 0 ; 
    double Y_scale = 0 ; 
    for(int i=1;i < h2D_SS->GetNbinsX();i++){
        X_scale = X_scale + Gauss(h2D_SS->GetYaxis()->GetBinCenter(i),mean_y,sigma_y  );
        Y_scale = Y_scale + Gauss(h2D_SS->GetXaxis()->GetBinCenter(i),mean_x,sigma_x  );
    }
    
    
    TF1* f1 = new TF1("f1","[0]*exp(-0.5*((x-[1])/[2])**2)",1.1, 1.13);
    
    //int fbin_x1 = h2D_SS->GetXaxis()->FindBin(mean_x-2*TMath::Abs(sigma_x));
    //int fbin_x2 = h2D_SS->GetXaxis()->FindBin(mean_x+2*TMath::Abs(sigma_x));
    //int fbin_y1 = h2D_SS->GetYaxis()->FindBin(mean_y-2*TMath::Abs(sigma_y));
    //int fbin_y2 = h2D_SS->GetYaxis()->FindBin(mean_y+2*TMath::Abs(sigma_y));
    
    //double NSignal = h2D_SS->Integral(fbin_x1,fbin_x2,fbin_y1,fbin_y2);
    //double NBackground = h2D_SB_BS->Integral(fbin_x1,fbin_x2,fbin_y1,fbin_y2);
    //std::cout<<"S/B="<<NSignal/NBackground<<std::endl;

    int fbin_x1 = h1D_SS_X->FindBin(mean_x-2*TMath::Abs(sigma_x));
    int fbin_x2 = h1D_SS_X->FindBin(mean_x+2*TMath::Abs(sigma_x));
    std::cout<<"fbin_x1="<<fbin_x1<<std::endl;
    std::cout<<"fbin_x12"<<fbin_x2<<std::endl;
    //double NSignal = h1D_SS_X->Integral(fbin_x1,fbin_x2);
    //double NBackground = h1D_SB_BS_X->Integral(fbin_x1,fbin_x2);
    double NSignal = h1D_SS_X->Integral();
    double NBackground = h1D_SB_BS_X->Integral();
    std::cout<<"S="<<NSignal<<std::endl;
    std::cout<<"B="<<NBackground<<std::endl;
    std::cout<<"S/B="<<NSignal/NBackground<<std::endl;
  

    

    TCanvas *c1 =new TCanvas("c2","c2",1200,600);
    c1->Divide(3,1);
    c1->cd(1);

    
    //f2->Draw("surf");
    h2D_SS->Draw("surf");

    c1->cd(2);
    //h1D_SS_X->Fit(f1);
    //f1->SetParameters(amplitude * X_scale,mean_x, sigma_x);
    //h1D_SS_X->GetYaxis()->SetRangeUser(-10000,120000);
    h1D_SS_X->Draw();
    f1->Draw("same");
    h1D_SB_BS_X->Draw("same");
    c1->cd(3);
    f1->SetParameters(amplitude * Y_scale,mean_y, sigma_y);
    //h1D_SS_Y->GetYaxis()->SetRangeUser(-10000,120000);
    h1D_SS_Y->Draw();
    f1->Draw("same");
    h1D_SB_BS_Y->Draw("same");

}




void Mass_plot(){
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

    int Lambda_Flag =0;
   /*
    for(int i = 0 ;i <Pt1_Bin ;i ++){
        std::string s1;
        s1 = Form("%.2f < Pt < %.2f GeV/c",Lambda_Pt_LowerBin[i],Lambda_Pt_HigherBin[i]);
        if(Lambda_Flag ==0 ){
            Background_Subtraction(h1D_USL_MassDist[i],h1D_LS_MassDist[i],s1,"#Lambda",i,Lambda_Flag);
        }
        if(Lambda_Flag == 1 ){
            Background_Subtraction(h1D_USLbar_MassDist[i],h1D_LS_MassDist[i],s1,"#bar{#Lambda}",i,Lambda_Flag);
        }


    }
    */
   
   


    

    
    //Background_Subtraction_2D(h2D_USL_USLbar_MassDist[4][1],h2D_USL_LS_MassDist[4][1],SwapAxesManually(h2D_USLbar_LS_MassDist[4][1]),h2D_LS_LS_MassDist[4][1]);
    Background_Subtraction_2D(h2D_USL_USL_MassDist[4][1],h2D_USL_LS_MassDist[4][1],SwapAxesManually(h2D_USL_LS_MassDist[4][1]),h2D_LS_LS_MassDist[4][1]);








}