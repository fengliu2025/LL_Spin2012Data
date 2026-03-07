//  Centrality binning:
//     Bin       Centrality (16)   Centrality (9)
//     0            75-80%            70-80%
//     1            70-75%            60-70%
//     2            65-70%            50-60%
//     3            60-65%            40-50%
//     4            55-60%            30-40%
//     5            50-55%            20-30%
//     6            45-50%            10-20%
//     7            40-45%             5-10%
//     8            35-40%             0- 5%
//     9            30-35%
//    10            25-30%
//    11            20-25%
//    12            15-20%
//    13            10-15%
//    14             5-10%
//    15             0- 5%
//centrality bins used (centrality 9):
// 0-10% (7+8), 10-40% (4+5+6), 40-80 (0+1+2+3), 0-80% (0-8)

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include"TH1.h"
#include"TH2.h"
#include"TH3.h"
#include"TF1.h"
#include"TF2.h"
#include"TF12.h"
#include"TMath.h"
#include"TCanvas.h"
#include"TFile.h"
#include"TLatex.h"
#include"TStyle.h"
#include"TPad.h"
#include"TLegend.h"
#include"TPaveText.h"
#include"TAxis.h"
#include"TTree.h"
#include"TFitResultPtr.h"
#include"TFitResult.h"
#include"TString.h"
#include"TLine.h"
#include"TChain.h"
#include"TLorentzVector.h"
#include"TGraphErrors.h"


using namespace std;

//const int nPtBins = 9;
//float const pT_bins[nPtBins+1] = { 0., 0.5, 1.,1.5, 2., 2.5, 3., 4., 5., 7.};

const int nPtBins = 8;
float const pT_bins[nPtBins+1] = { 0., 0.5, 1.,1.5, 2., 2.5, 3., 4., 5.};

const int nPtBins_corr = 2;
float const pT_bins_corr[nPtBins_corr+1] = { 0.5, 1.5, 5.};

const int nEtaBins = 3;
//float const eta_bins[nEtaBins+1] = { -1, -0.2, 0.2, 1 };
float const eta_bins[nEtaBins+1] = { -1, -0.4, 0.4, 1 };

const double pi_mass_PDG = 0.13957039; //p mass on GeV/c^2 from latest PDG
const double p_mass_PDG = 0.93827208816; //p mass on GeV/c^2 from latest PDG
const double L_mass_PDG = 1.115683; //mass in GeV/c^2 from latest PDG

const float L_y_cut = 1;
//0 - hybrid TOF for both daughters, 1 - strict TOF for pions, 2 - strict TOF for both pion and proton
//also check candidate tree - may have TOF requirement in production
//const int strictTOF_cut = 0;
const float L_cos_theta_cut = 0.996;
const float L_decayL_cut = 25;


bool Lambda_corr_2D_get_corr_Delta_R(const int cut_type = 0, const int energy = 510, const int year = 2017)
{
  //analyze stored Lambda pairs and save cos(theta*) histograms

  //2024 values
  const float L0_alpha = 0.747; //decay parameter of L0
  const float L0_alpha_relat_err = 0.009/L0_alpha; //relative error of decay parameter

  const float L0bar_alpha = -0.757; //decay paramteter of L0bar
  const float L0bar_alpha_relat_err = 0.004/fabs(L0bar_alpha); //relative error of decay paramteter


  //_______________________________________________________________________________________________________________________________________________

  //systematic error taken from short range in Delta y Delta phi
  TFile *inFile_sys_err_data = new TFile(Form("./output/Polarization/Polarization_Delta_y_Delta_phi.root", year), "read");

  TGraphErrors* pLLbar_short_work_sys = (TGraphErrors*) inFile_sys_err_data->Get("PolarizationGraph_delta_y_delta_phi_sys_err_average_0");


  float sys_err_tot_abs = pLLbar_short_work_sys->GetErrorX(1); //store systematic error from short-range


  //_______________________________________________________________________________________________________________________________________________


  //output file with polarization graphs
  TFile *out_file = new TFile(Form("./output/Polarization/%i/Polarization_Delta_R_bins.root", year), "recreate");


  TFile *LLbarOutFile; //input file with data histograms

  if(cut_type == 0)
  {
    //Run12

    LLbarOutFile = new TFile(Form("./output/ProdPlane/%i/ProdPlane_Lambda_ana_cuts.root", year), "read"); //analysis cuts    

    //--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  }
  else
  {
    cout<<"Wrong cut type"<<endl;

    return false;
  }


  //data histograms
  TH2F *L0_L0bar_cos_theta_star_vs_delta_R_scan_US_hist = (TH2F*)LLbarOutFile->Get("L0_L0bar_cos_theta_star_vs_delta_R_scan_US_hist");

  TH2F *L0_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_hist = (TH2F*)LLbarOutFile->Get("L0_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_hist");

  //--------------------------------------------------------

  TH2F *L0_L0_cos_theta_star_vs_delta_R_scan_US_hist = (TH2F*)LLbarOutFile->Get("L0_L0_cos_theta_star_vs_delta_R_scan_US_hist");

  TH2F *L0_L0_cos_theta_star_vs_delta_R_scan_US_LS_hist = (TH2F*)LLbarOutFile->Get("L0_L0_cos_theta_star_vs_delta_R_scan_US_LS_hist");

  //--------------------------------------------------------

  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_hist = (TH2F*)LLbarOutFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_hist");

  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_hist = (TH2F*)LLbarOutFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_hist");

  //--------------------------------------------------------

  //mixed event
  TH2F *L0_L0bar_cos_theta_star_vs_delta_R_scan_US_ME_hist = (TH2F*)LLbarOutFile->Get("L0_L0bar_cos_theta_star_vs_delta_R_scan_US_ME_hist");
  TH2F *L0_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_ME_hist = (TH2F*)LLbarOutFile->Get("L0_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_ME_hist");

  TH2F *L0_L0bar_cos_theta_star_vs_delta_R_scan_US_ME_alt_hist = (TH2F*)LLbarOutFile->Get("L0_L0bar_cos_theta_star_vs_delta_R_scan_US_ME_alt_hist");
  TH2F *L0_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_ME_alt_hist = (TH2F*)LLbarOutFile->Get("L0_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_ME_alt_hist");

  //--------------------------------------------------------

  TH2F *L0_L0_cos_theta_star_vs_delta_R_scan_US_ME_hist = (TH2F*)LLbarOutFile->Get("L0_L0_cos_theta_star_vs_delta_R_scan_US_ME_hist");
  TH2F *L0_L0_cos_theta_star_vs_delta_R_scan_US_LS_ME_hist = (TH2F*)LLbarOutFile->Get("L0_L0_cos_theta_star_vs_delta_R_scan_US_LS_ME_hist");

  TH2F *L0_L0_cos_theta_star_vs_delta_R_scan_US_ME_alt_hist = (TH2F*)LLbarOutFile->Get("L0_L0_cos_theta_star_vs_delta_R_scan_US_ME_alt_hist");
  TH2F *L0_L0_cos_theta_star_vs_delta_R_scan_US_LS_ME_alt_hist = (TH2F*)LLbarOutFile->Get("L0_L0_cos_theta_star_vs_delta_R_scan_US_LS_ME_alt_hist");

  //--------------------------------------------------------

  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_ME_hist = (TH2F*)LLbarOutFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_ME_hist");
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_ME_hist = (TH2F*)LLbarOutFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_ME_hist");

  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_ME_alt_hist = (TH2F*)LLbarOutFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_ME_alt_hist");
  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_ME_alt_hist = (TH2F*)LLbarOutFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_ME_alt_hist");

  //--------------------------------------------------------

  //TH1F *L0_L0bar_Delta_R_US_hist = (TH1F*)LLbarOutFile->Get("L0_L0bar_Delta_R_US_hist");
  //TH1F *L0_L0bar_Delta_R_LS_hist = (TH1F*)LLbarOutFile->Get("L0_L0bar_Delta_R_US_LS_hist");

  //________________________________________________________________________________________



  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLegendFont(42);

  //polarization graph for Delta y
  TGraphErrors *PolarizationGraph_Delta_R = new TGraphErrors(5);
  TGraphErrors *PolarizationGraph_Delta_R_sys_err = new TGraphErrors(5);
  TGraphErrors *PolarizationGraph_Delta_R_sys_err_average = new TGraphErrors(5);


  out_file->cd();

  for( unsigned int Delta_R_bin = 1; Delta_R_bin <= 5; Delta_R_bin++)
  {

    //int bin_min = 1+(Delta_R_bin-1)*10; //10 bin projections (0.5 wide bins)
    //int bin_max = 10+(Delta_R_bin-1)*10;

    int bin_min = 0; //10 bin projections (0.5 wide bins)
    int bin_max = 0;
/*
    if( Delta_R_bin == 1 )
    {
      bin_min = 1; //10 bin projections (0.5 wide bins)
      bin_max = 20;
    }
    else
    {
      bin_min = 21; //10 bin projections (0.5 wide bins)
      bin_max = 60;
    }
*/
    if( Delta_R_bin < 5 )
    {
      bin_min = 1+(Delta_R_bin-1)*10; //10 bin projections (0.5 wide bins)
      bin_max = 10+(Delta_R_bin-1)*10;
    }
    else
    {
      bin_min = 1+(Delta_R_bin-1)*10; //10 bin projections (0.5 wide bins)
      bin_max = 60;
    }


    float bin_edge_low = L0_L0bar_cos_theta_star_vs_delta_R_scan_US_hist->GetYaxis()->GetBinLowEdge(bin_min);
    float bin_edge_high = L0_L0bar_cos_theta_star_vs_delta_R_scan_US_hist->GetYaxis()->GetBinUpEdge(bin_max);

    //cout<<bin_edge_low<<endl;
    //cout<<bin_edge_high<<endl;

    float Delta_R_cut = bin_edge_high - (bin_edge_high-bin_edge_low)/2;




    //-----------------------------------------------------------------------------------------

    TCanvas *L0_L0bar_cosThetaProdPlane_eta_no_corr_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_eta_no_corr_can_%i", Delta_R_bin), Form("L0_L0bar_cosThetaProdPlane_eta_no_corr_can_%i", Delta_R_bin), 1200, 1000);

    L0_L0bar_cosThetaProdPlane_eta_no_corr_can->cd();

    TH1D *L0_L0bar_cosThetaProdPlane_eta_US_hist = L0_L0bar_cos_theta_star_vs_delta_R_scan_US_hist->ProjectionX( Form("proj_US_%i", Delta_R_bin), bin_min, bin_max);

    L0_L0bar_cosThetaProdPlane_eta_US_hist->GetXaxis()->SetTitle("cos(#theta*)");
    L0_L0bar_cosThetaProdPlane_eta_US_hist->GetXaxis()->CenterTitle();
    //L0_L0bar_cosThetaProdPlane_eta_US_hist->GetXaxis()->SetTextSizePixels(30);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0_L0bar_cosThetaProdPlane_eta_US_hist->GetYaxis()->CenterTitle();
    //L0_L0bar_cosThetaProdPlane_eta_US_hist->GetYaxis()->SetTextSizePixels(30);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->GetYaxis()->SetMaxDigits(3);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->SetMarkerSize(1.5);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->SetMarkerStyle(20);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->SetMarkerColor(kRed);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->SetLineColor(kRed);
    double nLLbar = L0_L0bar_cosThetaProdPlane_eta_US_hist->Integral();
    L0_L0bar_cosThetaProdPlane_eta_US_hist->Sumw2();
    //L0_L0bar_cosThetaProdPlane_eta_US_hist->Divide(L0_L0bar_cosThetaProdPlane_eta_eff);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->Scale(1./L0_L0bar_cosThetaProdPlane_eta_US_hist->GetXaxis()->GetBinWidth(1));
    //L0_L0bar_cosThetaProdPlane_eta_US_hist->Scale(1./L0_L0bar_cosThetaProdPlane_eta_US_hist->Integral());
    L0_L0bar_cosThetaProdPlane_eta_US_hist->SetMinimum(0);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->Draw("p e");

    //ME here just for plotting, used loser
    TH1D *L0_L0bar_cosThetaProdPlane_eta_ME_hist = L0_L0bar_cos_theta_star_vs_delta_R_scan_US_ME_hist->ProjectionX( Form("proj_US_ME_%i", Delta_R_bin), bin_min, bin_max);

    L0_L0bar_cosThetaProdPlane_eta_ME_hist->SetMarkerSize(1.5);
    L0_L0bar_cosThetaProdPlane_eta_ME_hist->SetMarkerStyle(24);
    L0_L0bar_cosThetaProdPlane_eta_ME_hist->SetMarkerColor(1);
    L0_L0bar_cosThetaProdPlane_eta_ME_hist->SetLineColor(1);
    L0_L0bar_cosThetaProdPlane_eta_ME_hist->Sumw2();
    L0_L0bar_cosThetaProdPlane_eta_ME_hist->Scale(nLLbar/L0_L0bar_cosThetaProdPlane_eta_ME_hist->Integral()); //scale ME to US
    L0_L0bar_cosThetaProdPlane_eta_ME_hist->Scale(1./L0_L0bar_cosThetaProdPlane_eta_ME_hist->GetXaxis()->GetBinWidth(1)); //US is scaled by bin width, ME needs to be scaled as well
    L0_L0bar_cosThetaProdPlane_eta_ME_hist->SetMinimum(0);
    L0_L0bar_cosThetaProdPlane_eta_ME_hist->Draw("p e same");

    TF1 *fitL0_L0bar_US_ThetaStar_no_corr_ME = new TF1("fitL0_L0bar_US_ThetaStar_no_corr_ME", "[0]*(1 + [1]*x)", -1, 1);
    fitL0_L0bar_US_ThetaStar_no_corr_ME->SetParameters(100, 0.5);

    L0_L0bar_cosThetaProdPlane_eta_ME_hist->Fit(fitL0_L0bar_US_ThetaStar_no_corr_ME, "s i 0 r");

    float P_L0_L0bar_no_corr_ME = fitL0_L0bar_US_ThetaStar_no_corr_ME->GetParameter(1)/(L0_alpha*L0bar_alpha);
    float P_L0_L0bar_no_corr_ME_err = fitL0_L0bar_US_ThetaStar_no_corr_ME->GetParError(1)/(L0_alpha*L0bar_alpha);


    TH1D *L0_L0bar_cosThetaProdPlane_eta_LS_hist = L0_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_hist->ProjectionX( Form("proj_LS_%i", Delta_R_bin), bin_min, bin_max);

    L0_L0bar_cosThetaProdPlane_eta_LS_hist->SetMarkerSize(1.5);
    L0_L0bar_cosThetaProdPlane_eta_LS_hist->SetMarkerStyle(21);
    L0_L0bar_cosThetaProdPlane_eta_LS_hist->SetMarkerColor(kBlue);
    double nLLbar_back = L0_L0bar_cosThetaProdPlane_eta_LS_hist->Integral();
    L0_L0bar_cosThetaProdPlane_eta_LS_hist->Sumw2();
    L0_L0bar_cosThetaProdPlane_eta_LS_hist->Scale(1./L0_L0bar_cosThetaProdPlane_eta_LS_hist->GetXaxis()->GetBinWidth(1));
    //L0_L0bar_cosThetaProdPlane_eta_LS_hist->Divide(L0_L0bar_cosThetaProdPlane_eta_eff);
    L0_L0bar_cosThetaProdPlane_eta_LS_hist->Draw("p e same");


    TH1D *L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist = L0_L0bar_cos_theta_star_vs_delta_R_scan_US_LS_ME_hist->ProjectionX( Form("proj_LS_ME_%i", Delta_R_bin), bin_min, bin_max);

    L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist->SetMarkerSize(1.5);
    L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist->SetMarkerStyle(25);
    L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist->SetMarkerColor(kMagenta+1);
    L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist->SetLineColor(kMagenta+1);
    L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist->Sumw2();
    L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist->Scale(nLLbar_back/L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist->Integral()); //scale ME_LS to background
    L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist->Scale(1./L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist->GetXaxis()->GetBinWidth(1)); //background is scaled by bin width, ME_LS needs to be scaled as well
    L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist->SetMinimum(0);
    L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist->Draw("p e same");


    TF1 *fitL0_L0bar_US_ThetaStar_no_corr = new TF1("fitL0_L0bar_US_ThetaStar_no_corr", "[0]*(1 + [1]*x)", -1, 1);
    fitL0_L0bar_US_ThetaStar_no_corr->SetParameters(100, 0.5);

    //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->Fit(fitL0_L0bar_US_ThetaStar_no_corr, "s i 0 r");

    float P_L0_L0bar_no_corr = fitL0_L0bar_US_ThetaStar_no_corr->GetParameter(1)/(L0_alpha*L0bar_alpha);
    float P_L0_L0bar_no_corr_err = fitL0_L0bar_US_ThetaStar_no_corr->GetParError(1)/(L0_alpha*L0bar_alpha);

    fitL0_L0bar_US_ThetaStar_no_corr->SetLineColor(1);
    //fitL0_L0bar_US_ThetaStar_no_corr->Draw("same");

    TLegend *L0_L0bar_leg = new TLegend(0.15, 0.45, 0.45, 0.69);
    L0_L0bar_leg->AddEntry(L0_L0bar_cosThetaProdPlane_eta_US_hist, "(US-US) p#pi");
    L0_L0bar_leg->AddEntry(L0_L0bar_cosThetaProdPlane_eta_ME_hist, "(US-US) ME");
    L0_L0bar_leg->AddEntry(L0_L0bar_cosThetaProdPlane_eta_LS_hist, "Combinatorial bckg.");
    L0_L0bar_leg->AddEntry(L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist, "Bckg. ME");
    //L0_L0bar_leg->AddEntry(fitL0_L0bar_US_ThetaStar_no_corr, "Linear fit to US");
    L0_L0bar_leg->SetBorderSize(0);
    L0_L0bar_leg->SetFillColorAlpha(0, 0.01);
    L0_L0bar_leg->Draw("same");

    TPaveText *L0_L0bar_text_no_corr = new TPaveText(0.5, 0.4, 0.85, 0.75, "NDC");
    L0_L0bar_text_no_corr->SetTextFont(42);
    //L0_L0bar_text_no_corr->AddText("STAR Internal");
    //L0_L0bar_text_no_corr->AddText("STAR preliminary");
    //((TText*)L0_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
    L0_L0bar_text_no_corr->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
    L0_L0bar_text_no_corr->AddText("Minimum bias, no correction");
    L0_L0bar_text_no_corr->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
    L0_L0bar_text_no_corr->AddText("|#it{y}| < 1");
    L0_L0bar_text_no_corr->AddText("0.5 < p_{T,1} < 5.0 GeV/#it{c}");
    L0_L0bar_text_no_corr->AddText("0.5 < p_{T,2} < 5.0 GeV/#it{c}");
    L0_L0bar_text_no_corr->AddText(Form("P = %.2f #pm %.2f", P_L0_L0bar_no_corr, fabs(P_L0_L0bar_no_corr_err)));
    L0_L0bar_text_no_corr->AddText(Form("P_{ME} = %.3f #pm %.3f", P_L0_L0bar_no_corr_ME, fabs(P_L0_L0bar_no_corr_ME_err)));
    L0_L0bar_text_no_corr->SetFillColorAlpha(0, 0.01);
    L0_L0bar_text_no_corr->Draw("same");

    L0_L0bar_cosThetaProdPlane_eta_no_corr_can->SaveAs(Form("./plots/Lambda/2D_analysis/L_correlations_delta_R_scan_bins/L0_L0bar_cosThetaProdPlane_eta_no_corr_less_Delta_R_%i.png", Delta_R_bin));

    //----------------------------------------------------------------------------------------------------

    TCanvas *L0_L0bar_cosThetaProdPlane_eta_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_eta_can_%i", Delta_R_bin), Form("L0_L0bar_cosThetaProdPlane_eta_can_%i", Delta_R_bin), 1200, 1000);

    L0_L0bar_cosThetaProdPlane_eta_can->cd();

    //ME histogram higher

    L0_L0bar_cosThetaProdPlane_eta_US_hist->Divide(L0_L0bar_cosThetaProdPlane_eta_ME_hist); // correct US using ME
    L0_L0bar_cosThetaProdPlane_eta_US_hist->Scale(nLLbar/L0_L0bar_cosThetaProdPlane_eta_US_hist->Integral()); //scale back to raw US
    L0_L0bar_cosThetaProdPlane_eta_US_hist->Scale(1./L0_L0bar_cosThetaProdPlane_eta_US_hist->GetXaxis()->GetBinWidth(1)); //bin width scaling
    L0_L0bar_cosThetaProdPlane_eta_US_hist->SetMinimum(0);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->Draw("p e");

    //L0_L0bar_cosThetaProdPlane_eta_ME_hist->Draw("same p e");

    TH1D *L0_L0bar_cosThetaProdPlane_eta_LS_hist_clone = (TH1D*)L0_L0bar_cosThetaProdPlane_eta_LS_hist->Clone("L0_L0bar_cosThetaProdPlane_eta_LS_hist_clone");
    L0_L0bar_cosThetaProdPlane_eta_LS_hist_clone->Divide(L0_L0bar_cosThetaProdPlane_eta_ME_hist); // correct US using ME

    L0_L0bar_cosThetaProdPlane_eta_LS_hist->Divide(L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist); //correct background using ME
    L0_L0bar_cosThetaProdPlane_eta_LS_hist->Scale(nLLbar_back/L0_L0bar_cosThetaProdPlane_eta_LS_hist->Integral()); //scale back to raw background
    L0_L0bar_cosThetaProdPlane_eta_LS_hist->Scale(1./L0_L0bar_cosThetaProdPlane_eta_LS_hist->GetXaxis()->GetBinWidth(1)); //bin width scaling
    L0_L0bar_cosThetaProdPlane_eta_LS_hist->SetMinimum(0);
    L0_L0bar_cosThetaProdPlane_eta_LS_hist->Draw("p e same");

    //fit dN/dcos(theta*)
    //signal + bacground
    TF1 *fitL0_L0bar_US_ThetaStar = new TF1("fitL0_L0bar_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
    fitL0_L0bar_US_ThetaStar->SetParameters(100, 0.5);

    //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->Fit(fitL0_L0bar_US_ThetaStar, "s i 0 r");

    float P_L0_L0bar = fitL0_L0bar_US_ThetaStar->GetParameter(1)/(L0_alpha*L0bar_alpha);
    float P_L0_L0bar_err = fitL0_L0bar_US_ThetaStar->GetParError(1)/(L0_alpha*L0bar_alpha);

    fitL0_L0bar_US_ThetaStar->SetLineColor(1);
    fitL0_L0bar_US_ThetaStar->Draw("same");

    //background
    TF1 *fitL0_L0bar_US_LS_ThetaStar = new TF1("fitL0_L0bar_US_LS_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
    fitL0_L0bar_US_LS_ThetaStar->SetParameters(100, 0.5);

    //fit_res_gaUS_LS_wrong_sign = L_inv_mass_US_LS->Fit(fitGaUS_LSsBack, "s i 0", "", 1.07, 1.4);
    L0_L0bar_cosThetaProdPlane_eta_LS_hist->Fit(fitL0_L0bar_US_LS_ThetaStar, "s i 0 r");

    float P_L0_L0bar_back = fitL0_L0bar_US_LS_ThetaStar->GetParameter(1)/(L0_alpha*L0bar_alpha);
    float P_L0_L0bar_back_err = fitL0_L0bar_US_LS_ThetaStar->GetParError(1)/(L0_alpha*L0bar_alpha);

    fitL0_L0bar_US_LS_ThetaStar->SetLineColor(1);
    fitL0_L0bar_US_LS_ThetaStar->Draw("same");


    TPaveText *L0_L0bar_text = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
    L0_L0bar_text->SetTextFont(42);
    //L0_L0bar_text->AddText("STAR Internal");
    //L0_L0bar_text->AddText("STAR preliminary");
    //((TText*)L0_L0bar_text->GetListOfLines()->Last())->SetTextColor(2);
    L0_L0bar_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
    L0_L0bar_text->AddText("Minimum bias");
    L0_L0bar_text->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
    L0_L0bar_text->AddText("|#it{y}| < 1");
    L0_L0bar_text->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
    L0_L0bar_text->AddText(Form("P_{tot} = %.2f #pm %.2f", P_L0_L0bar, fabs(P_L0_L0bar_err)));
    L0_L0bar_text->AddText(Form("P_{bckg} = %.2f #pm %.2f", P_L0_L0bar_back, fabs(P_L0_L0bar_back_err)));
    L0_L0bar_text->SetFillColorAlpha(0, 0.01);
    L0_L0bar_text->Draw("same");

    L0_L0bar_leg->Draw("same");

    L0_L0bar_cosThetaProdPlane_eta_can->SaveAs(Form("./plots/Lambda/2D_analysis/L_correlations_delta_R_scan_bins/L0_L0bar_cosThetaProdPlane_eta_less_Delta_R_%i.png", Delta_R_bin));

    //----------------------------------------------------------------------

    TCanvas *L0_L0bar_cosThetaProdPlane_eta_can_2 = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_eta_can_2_%i", Delta_R_bin), Form("L0_L0bar_cosThetaProdPlane_eta_can_2_%i", Delta_R_bin), 1200, 1000);

    L0_L0bar_cosThetaProdPlane_eta_can_2->cd();

    L0_L0bar_cosThetaProdPlane_eta_US_hist->Add(L0_L0bar_cosThetaProdPlane_eta_LS_hist, -1);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->SetMinimum(0);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->Write(Form("L0_L0bar_cosThetaProdPlane_Delta_R_%i", Delta_R_bin-1));
    L0_L0bar_cosThetaProdPlane_eta_US_hist->Draw("p e");


    TF1 *fitL0_L0bar_US_ThetaStar_2 = new TF1("fitL0_L0bar_US_ThetaStar_2", "[0]*(1 + [1]*x)", -1, 1);
    fitL0_L0bar_US_ThetaStar_2->SetParameters(100, 0.5);

    //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->Fit(fitL0_L0bar_US_ThetaStar_2, "s i 0 r");




    float P_L0_L0bar_2 = fitL0_L0bar_US_ThetaStar_2->GetParameter(1)/(L0_alpha*L0bar_alpha);
    float P_L0_L0bar_err_2 = fitL0_L0bar_US_ThetaStar_2->GetParError(1)/(L0_alpha*L0bar_alpha);

    fitL0_L0bar_US_ThetaStar_2->SetLineColor(1);
    fitL0_L0bar_US_ThetaStar_2->Draw("same");






    //----------------------------

    //background subtraction
    float nLLbar_fit = fitL0_L0bar_US_ThetaStar->GetParameter(0);
    float nLLbar_fit_err = fitL0_L0bar_US_ThetaStar->GetParError(0);

    float nLLbar_back_fit = fitL0_L0bar_US_LS_ThetaStar->GetParameter(0);
    float nLLbar_back_fit_err = fitL0_L0bar_US_LS_ThetaStar->GetParError(0);

    float P_L0_L0bar_from_fits = P_L0_L0bar + nLLbar_back_fit/(nLLbar_fit-nLLbar_back_fit)*(P_L0_L0bar - P_L0_L0bar_back);

    float P_L0_L0bar_from_fits_err = sqrt( P_L0_L0bar_err*P_L0_L0bar_err + nLLbar_back_fit/(nLLbar_fit-nLLbar_back_fit)*nLLbar_back_fit/(nLLbar_fit-nLLbar_back_fit)*P_L0_L0bar_err*P_L0_L0bar_err +
                                           nLLbar_back_fit/(nLLbar_fit-nLLbar_back_fit)*nLLbar_back_fit/(nLLbar_fit-nLLbar_back_fit)*P_L0_L0bar_back_err*P_L0_L0bar_back_err +
                                           nLLbar_back_fit*nLLbar_back_fit/(nLLbar_fit-nLLbar_back_fit)/(nLLbar_fit-nLLbar_back_fit)/(nLLbar_fit-nLLbar_back_fit)/(nLLbar_fit-nLLbar_back_fit)*(P_L0_L0bar - P_L0_L0bar_back)*(P_L0_L0bar - P_L0_L0bar_back)*nLLbar_fit_err*nLLbar_fit_err +
                                           (P_L0_L0bar - P_L0_L0bar_back)*(P_L0_L0bar - P_L0_L0bar_back)/(nLLbar_fit-nLLbar_back_fit)/(nLLbar_fit-nLLbar_back_fit)/(nLLbar_fit-nLLbar_back_fit)/(nLLbar_fit-nLLbar_back_fit)*nLLbar_back_fit_err*nLLbar_back_fit_err );


    //--------------------------------------------------------------------------------------------------

    //total

    //float SysErrTot_L0_L0bar = sys_err_tot_abs;


    //TPaveText *L0_L0bar_text_2 = new TPaveText(0.47, 0.15, 0.85, 0.53, "NDC"); //with polarization value
    TPaveText *L0_L0bar_text_2 = new TPaveText(0.47, 0.15, 0.8, 0.53, "NDC"); //without polarization value
    L0_L0bar_text_2->SetTextFont(42);
    //L0_L0bar_text_2->SetTextSize(15);
    //L0_L0bar_text_2->AddText("STAR");
    //L0_L0bar_text_2->AddText("STAR preliminary");
    //((TText*)L0_L0bar_text_2->GetListOfLines()->Last())->SetTextColor(2);
    L0_L0bar_text_2->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
    L0_L0bar_text_2->AddText("Minimum bias");
    L0_L0bar_text_2->AddText("#Lambda^{0}-#bar{#Lambda^{0}}");
    L0_L0bar_text_2->AddText("|#it{y}| < 1");
    L0_L0bar_text_2->AddText("0.5 < p_{T,1} < 5.0 GeV/#it{c}");
    L0_L0bar_text_2->AddText("0.5 < p_{T,2} < 5.0 GeV/#it{c}");
    L0_L0bar_text_2->AddText(Form("P_{#Lambda_{1}#Lambda_{2}} = %.3f #pm %.3f #pm %.3f", P_L0_L0bar_2, fabs(P_L0_L0bar_from_fits_err), fabs(sys_err_tot_abs*P_L0_L0bar_from_fits)));
    //L0_L0bar_text_2->AddText(Form("P_{topo} = %.2f", P_L0_L0bar_tight_topo_cuts));
    //L0_L0bar_text_2->AddText(Form("P_{pT} = %.2f", P_L0_L0bar_tight_pT_cuts));
    L0_L0bar_text_2->SetFillColorAlpha(0, 0.01);
    L0_L0bar_text_2->Draw("same");


    TLegend *L0_L0bar_2_leg = new TLegend(0.15, 0.3, 0.4, 0.49);
    //L0_L0bar_2_leg->SetTextSizePixels(15);
    L0_L0bar_2_leg->AddEntry(L0_L0bar_cosThetaProdPlane_eta_US_hist, "(US-US)-Bckg.");
    L0_L0bar_2_leg->AddEntry(fitL0_L0bar_US_ThetaStar_2, "Fit", "l");
    L0_L0bar_2_leg->SetBorderSize(0);
    L0_L0bar_2_leg->SetFillColorAlpha(0, 0.01);
    L0_L0bar_2_leg->Draw("same");

    L0_L0bar_cosThetaProdPlane_eta_can_2->SaveAs(Form("./plots/Lambda/2D_analysis/L_correlations_delta_R_scan_bins/L0_L0bar_cosThetaProdPlane_eta_subtract_less_Delta_R_%i.png", Delta_R_bin));

    //-----------------------


    PolarizationGraph_Delta_R->SetPoint(Delta_R_bin, Delta_R_cut, P_L0_L0bar_from_fits);
    PolarizationGraph_Delta_R->SetPointError(Delta_R_bin, (bin_edge_high-bin_edge_low)/2, fabs(P_L0_L0bar_from_fits_err));


    PolarizationGraph_Delta_R_sys_err->SetPoint(Delta_R_bin, Delta_R_cut, P_L0_L0bar_from_fits);
    PolarizationGraph_Delta_R_sys_err->SetPointError(Delta_R_bin, 0.02, fabs(sys_err_tot_abs));


    PolarizationGraph_Delta_R_sys_err_average->SetPoint(Delta_R_bin, Delta_R_cut, P_L0_L0bar_from_fits);
    PolarizationGraph_Delta_R_sys_err_average->SetPointError(Delta_R_bin, 0.02, fabs(sys_err_tot_abs));

    //____________________________________________________________________________________________________________________________________________________________________________________________________________


  }




  PolarizationGraph_Delta_R->RemovePoint(0);
  PolarizationGraph_Delta_R_sys_err->RemovePoint(0);
  PolarizationGraph_Delta_R_sys_err_average->RemovePoint(0);


  TCanvas *PolarizationGraph_Delta_R_can = new TCanvas("PolarizationGraph_Delta_R_can", "PolarizationGraph_Delta_R_can", 2000, 1200);
  PolarizationGraph_Delta_R_can->Divide(2,1);

  TH1F *DefaultHist = new TH1F("DefaultHist", "DefaultHist", 1, 0, 4);
  DefaultHist->GetXaxis()->SetTitle("#DeltaR");
  DefaultHist->GetXaxis()->CenterTitle();
  DefaultHist->GetYaxis()->SetTitle("P_{#Lambda_{1}#Lambda_{2}}");
  DefaultHist->GetYaxis()->CenterTitle();
  DefaultHist->GetYaxis()->SetRangeUser(-0.25, 0.5);
  DefaultHist->Draw();

  PolarizationGraph_Delta_R->SetMarkerStyle(20);
  PolarizationGraph_Delta_R->SetMarkerSize(2);
  PolarizationGraph_Delta_R->SetMarkerColor(kRed);
  PolarizationGraph_Delta_R->SetLineColor(kRed);
  PolarizationGraph_Delta_R->Write("PolarizationGraph_Delta_R_bins");
  PolarizationGraph_Delta_R->Draw("p e same");

  PolarizationGraph_Delta_R_sys_err->SetMarkerSize(2);
  PolarizationGraph_Delta_R_sys_err->SetFillColorAlpha(kRed, 0.25);
  PolarizationGraph_Delta_R_sys_err->SetLineColor(kRed);
  PolarizationGraph_Delta_R_sys_err->Write("PolarizationGraph_Delta_R_bins_sys_err");
  //PolarizationGraph_Delta_R_sys_err->Draw("2 same");


  PolarizationGraph_Delta_R_sys_err_average->SetMarkerSize(2);
  PolarizationGraph_Delta_R_sys_err_average->SetFillColorAlpha(kRed, 0.25);
  PolarizationGraph_Delta_R_sys_err_average->SetLineColor(kRed);
  PolarizationGraph_Delta_R_sys_err_average->Write("PolarizationGraph_Delta_R_bins_sys_err_average");
  PolarizationGraph_Delta_R_sys_err_average->Draw("2 same");


  TLine *ZeroLine_eta = new TLine(0,0.01,0,3.99);
  ZeroLine_eta->SetLineStyle(9);
  ZeroLine_eta->SetLineColor(1);
  //ZeroLine_eta->Draw("same");

  //TPaveText *Polarization_text = new TPaveText(0.12, 0.58, 0.49, 0.89, "NDC"); //with K0s
  TPaveText *Polarization_text = new TPaveText(0.55, 0.75, 0.89, 0.89, "NDC"); //without K0s
  Polarization_text->SetTextFont(42);
  //Polarization_text->AddText("STAR");
  //Polarization_text->AddText("STAR preliminary");
  //((TText*)Polarization_text->GetListOfLines()->Last())->SetTextColor(2);
  Polarization_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
  Polarization_text->AddText("Minimum bias");
  Polarization_text->AddText("|#it{y}| < 1");
  Polarization_text->AddText("0.5 < p_{T,1} < 5.0 GeV/#it{c}");
  Polarization_text->AddText("0.5 < p_{T,2} < 5.0 GeV/#it{c}");
  Polarization_text->SetFillColorAlpha(0, 0.01);
  Polarization_text->Draw("same");

  PolarizationGraph_Delta_R_can->SaveAs("./plots/Lambda/2D_analysis/L_polarization/L_polarization_delta_R_scan_bins.png");

  //__________________________________________________________________________________________________

  //L0_L0bar_Delta_R_US_hist->Add(L0_L0bar_Delta_R_LS_hist, -1);
  //L0_L0bar_Delta_R_US_hist->Write("L0_L0bar_Delta_R_hist");

  //__________________________________________________________________________________________________



  LLbarOutFile->Close();

  return true;
}
