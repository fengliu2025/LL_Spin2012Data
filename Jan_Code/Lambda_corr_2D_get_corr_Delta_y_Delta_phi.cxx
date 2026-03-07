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


//2024 values
const float L0_alpha = 0.747; //decay parameter of L0
const float L0_alpha_relat_err = 0.009/L0_alpha; //relative error of decay parameter

const float L0bar_alpha = -0.757; //decay paramteter of L0bar
const float L0bar_alpha_relat_err = 0.004/fabs(L0bar_alpha); //relative error of decay paramteter


pair<float,float> Get_P_LL_from_fit(TF1* fitL0_L0bar_US, TF1* fitL0_L0bar_US_LS, int pair_type)
{
  float alpha_square = 0;

  if( pair_type == 0 ) alpha_square = L0_alpha*L0bar_alpha;
  else if( pair_type == 1 ) alpha_square = L0_alpha*L0_alpha;
  else if( pair_type == 2 ) alpha_square = L0bar_alpha*L0bar_alpha;
  else
  {
    cout<<"Wrong Get_P_LL_from_fit function third parameter!"<<endl;

    pair<double,double> P_LL_pair(0, 0);

    return P_LL_pair;
  }

  float P_L0_L0bar_from_fits = fitL0_L0bar_US->GetParameter(1)/(alpha_square);
  float P_L0_L0bar_from_fits_err = fitL0_L0bar_US->GetParError(1)/(alpha_square);

  float nLLbar_fit = fitL0_L0bar_US->GetParameter(0);
  float nLLbar_fit_err = fitL0_L0bar_US->GetParError(0);


  float P_L0_L0bar_from_fits_back = fitL0_L0bar_US_LS->GetParameter(1)/(alpha_square);
  float P_L0_L0bar_from_fits_back_err = fitL0_L0bar_US_LS->GetParError(1)/(alpha_square);

  float nLLbar_back_fit = fitL0_L0bar_US_LS->GetParameter(0);
  float nLLbar_back_fit_err = fitL0_L0bar_US_LS->GetParError(0);


  float P_L0_L0bar_from_fits_signal = P_L0_L0bar_from_fits + nLLbar_back_fit/(nLLbar_fit-nLLbar_back_fit)*(P_L0_L0bar_from_fits - P_L0_L0bar_from_fits_back);

  float P_L0_L0bar_from_fits_signal_err = sqrt( P_L0_L0bar_from_fits_err*P_L0_L0bar_from_fits_err + nLLbar_back_fit/(nLLbar_fit-nLLbar_back_fit)*nLLbar_back_fit/(nLLbar_fit-nLLbar_back_fit)*P_L0_L0bar_from_fits_err*P_L0_L0bar_from_fits_err +
                                         nLLbar_back_fit/(nLLbar_fit-nLLbar_back_fit)*nLLbar_back_fit/(nLLbar_fit-nLLbar_back_fit)*P_L0_L0bar_from_fits_back_err*P_L0_L0bar_from_fits_back_err +
                                         nLLbar_back_fit*nLLbar_back_fit/(nLLbar_fit-nLLbar_back_fit)/(nLLbar_fit-nLLbar_back_fit)/(nLLbar_fit-nLLbar_back_fit)/(nLLbar_fit-nLLbar_back_fit)*(P_L0_L0bar_from_fits - P_L0_L0bar_from_fits_back)*(P_L0_L0bar_from_fits - P_L0_L0bar_from_fits_back)*nLLbar_fit_err*nLLbar_fit_err +
                                         (P_L0_L0bar_from_fits - P_L0_L0bar_from_fits_back)*(P_L0_L0bar_from_fits - P_L0_L0bar_from_fits_back)/(nLLbar_fit-nLLbar_back_fit)/(nLLbar_fit-nLLbar_back_fit)/(nLLbar_fit-nLLbar_back_fit)/(nLLbar_fit-nLLbar_back_fit)*nLLbar_back_fit_err*nLLbar_back_fit_err );



  pair<double,double> P_LL_pair(P_L0_L0bar_from_fits_signal, P_L0_L0bar_from_fits_signal_err);

  return P_LL_pair;

}


bool Lambda_corr_2D_get_corr_Delta_y_Delta_phi(const int cut_type = 0, const int energy = 510, const int year = 2017)
{
  //analyze stored Lambda pairs and save cos(theta*) histograms


  //_______________________________________________________________________________________________________________________________________________

  //systematic uncertainties
  //residual effect from PYTHIA closure test
  TFile *SysErrSlopeFile = new TFile("../Simulation/Closure_test/output/SysErr/SysErrSlope.root", "read");

  TH1F *SysErrSlope_delta_eta_delta_phi_hist[2];
  TH1F *ResidualPolarization_delta_eta_delta_phi_hist[2];

  for( unsigned int delta_eta_bin = 0; delta_eta_bin < 2; delta_eta_bin++ )
  {
    SysErrSlope_delta_eta_delta_phi_hist[delta_eta_bin] = (TH1F*)SysErrSlopeFile->Get(Form("SysErrSlope_delta_eta_delta_phi_hist_%i", delta_eta_bin));
    ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin] = (TH1F*)SysErrSlopeFile->Get(Form("ResidualPolarization_delta_eta_delta_phi_hist_%i", delta_eta_bin));
  }

  //output file with polarization graphs
  TFile *out_file = new TFile(Form("./output/Polarization/%i/Polarization_Delta_y_Delta_phi.root", year), "recreate");

  //systematic uncertainty histograms and values
  //alpha
  float sysErr_alpha_L0_L0bar = sqrt(L0_alpha_relat_err*L0_alpha_relat_err + L0bar_alpha_relat_err*L0bar_alpha_relat_err);
  float sysErr_alpha_L0_L0 = 2*L0_alpha_relat_err;
  float sysErr_alpha_L0bar_L0bar = 2*L0bar_alpha_relat_err;


  //DCA_PV sys. err. (calculated manually from DCA_PV variation in data)
  const float sysErr_DCA_abs = 0.004;

  //cuts variation
  //have to run this code with cuts_type = 1 and 2 first
  TFile *SysErrCutsTopo;
  TFile *SysErrCutsPt;
  TFile *SysErrCutsDCA;

  if(cut_type == 0 )
  {
    SysErrCutsTopo = new TFile(Form("./input/SysErr/%i/SysErr_tight_topo_cuts_Delta_y_Delta_phi_work.root", year), "read");

    SysErrCutsPt = new TFile(Form("./input/SysErr/%i/SysErr_tight_pT_cuts_Delta_y_Delta_phi_work.root", year), "read");

    SysErrCutsDCA = new TFile(Form("./input/SysErr/%i/SysErr_tight_DCA_cuts_Delta_y_Delta_phi_work.root", year), "read");
    
  }
  else
  {
    if( cut_type == 1 ) SysErrCutsTopo = new TFile(Form("./input/SysErr/%i/SysErr_tight_topo_cuts_Delta_y_Delta_phi_work.root", year), "recreate");
    if( cut_type == 2 ) SysErrCutsPt = new TFile(Form("./input/SysErr/%i/SysErr_tight_pT_cuts_Delta_y_Delta_phi_work.root", year), "recreate");
    if( cut_type == 3 ) SysErrCutsDCA = new TFile(Form("./input/SysErr/%i/SysErr_tight_DCA_cuts_Delta_y_Delta_phi_work.root", year), "recreate");
  }

  //tight cuts
  //topological cuts
  TF1 *fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US[2];
  TF1 *fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi_US[2];
  TF1 *fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US[2];

  TF1 *fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US_LS[2];
  TF1 *fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi_US_LS[2];
  TF1 *fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US_LS[2];

  TF1 *fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi[2];
  TF1 *fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi[2];
  TF1 *fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi[2];

  //daughter pT cuts
  TF1 *fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US[2];
  TF1 *fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi_US[2];
  TF1 *fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US[2];

  TF1 *fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US_LS[2];
  TF1 *fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi_US_LS[2];
  TF1 *fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US_LS[2];

  TF1 *fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi[2];
  TF1 *fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi[2];
  TF1 *fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi[2];

  //daughter DCA cuts
  TF1 *fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US[2];
  TF1 *fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi_US[2];
  TF1 *fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US[2];

  TF1 *fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US_LS[2];
  TF1 *fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi_US_LS[2];
  TF1 *fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US_LS[2];

  TF1 *fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi[2];
  TF1 *fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi[2];
  TF1 *fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi[2];

  for( unsigned int delta_eta_bin = 0; delta_eta_bin < 2; delta_eta_bin++)
  {
    if(cut_type == 0 )
    {
      fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin] = (TF1*)SysErrCutsTopo->Get(Form("fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US_%i", delta_eta_bin));
      fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin] = (TF1*)SysErrCutsTopo->Get(Form("fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi_US_%i", delta_eta_bin));
      fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin] = (TF1*)SysErrCutsTopo->Get(Form("fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US_%i", delta_eta_bin));

      fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin] = (TF1*)SysErrCutsTopo->Get(Form("fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US_LS_%i", delta_eta_bin));
      fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin] = (TF1*)SysErrCutsTopo->Get(Form("fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi_US_LS_%i", delta_eta_bin));
      fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin] = (TF1*)SysErrCutsTopo->Get(Form("fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US_LS_%i", delta_eta_bin));

      fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin] = (TF1*)SysErrCutsTopo->Get(Form("fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi_%i", delta_eta_bin));
      fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin] = (TF1*)SysErrCutsTopo->Get(Form("fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi_%i", delta_eta_bin));
      fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin] = (TF1*)SysErrCutsTopo->Get(Form("fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi_%i", delta_eta_bin));

      //-----------------------------------

      fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin] = (TF1*)SysErrCutsPt->Get(Form("fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US_%i", delta_eta_bin));
      fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin] = (TF1*)SysErrCutsPt->Get(Form("fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi_US_%i", delta_eta_bin));
      fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin] = (TF1*)SysErrCutsPt->Get(Form("fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US_%i", delta_eta_bin));

      fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin] = (TF1*)SysErrCutsPt->Get(Form("fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US_LS_%i", delta_eta_bin));
      fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin] = (TF1*)SysErrCutsPt->Get(Form("fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi_US_LS_%i", delta_eta_bin));
      fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin] = (TF1*)SysErrCutsPt->Get(Form("fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US_LS_%i", delta_eta_bin));

      fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin] = (TF1*)SysErrCutsPt->Get(Form("fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi_%i", delta_eta_bin));
      fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin] = (TF1*)SysErrCutsPt->Get(Form("fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi_%i", delta_eta_bin));
      fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin] = (TF1*)SysErrCutsPt->Get(Form("fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi_%i", delta_eta_bin));

      //-----------------------------------

      fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US[delta_eta_bin] = (TF1*)SysErrCutsDCA->Get(Form("fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US_%i", delta_eta_bin));
      fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi_US[delta_eta_bin] = (TF1*)SysErrCutsDCA->Get(Form("fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi_US_%i", delta_eta_bin));
      fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US[delta_eta_bin] = (TF1*)SysErrCutsDCA->Get(Form("fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US_%i", delta_eta_bin));

      fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin] = (TF1*)SysErrCutsDCA->Get(Form("fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US_LS_%i", delta_eta_bin));
      fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin] = (TF1*)SysErrCutsDCA->Get(Form("fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi_US_LS_%i", delta_eta_bin));
      fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin] = (TF1*)SysErrCutsDCA->Get(Form("fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US_LS_%i", delta_eta_bin));

      fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi[delta_eta_bin] = (TF1*)SysErrCutsDCA->Get(Form("fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_%i", delta_eta_bin));
      fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi[delta_eta_bin] = (TF1*)SysErrCutsDCA->Get(Form("fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi_%i", delta_eta_bin));
      fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi[delta_eta_bin] = (TF1*)SysErrCutsDCA->Get(Form("fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_%i", delta_eta_bin));

    }
    else
    {
      fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin] = new TF1(Form("fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin]->SetParameters(100, 0.5);

      fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin] = new TF1(Form("fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi_US_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin]->SetParameters(100, 0.5);

      fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin] = new TF1(Form("fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin]->SetParameters(100, 0.5);


      fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin] = new TF1(Form("fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US_LS_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin]->SetParameters(100, 0.5);

      fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin] = new TF1(Form("fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi_US_LS_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin]->SetParameters(100, 0.5);

      fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin] = new TF1(Form("fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US_LS_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin]->SetParameters(100, 0.5);


      fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin] = new TF1(Form("fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin]->SetParameters(100, 0.5);

      fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin] = new TF1(Form("fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin]->SetParameters(100, 0.5);

      fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin] = new TF1(Form("fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin]->SetParameters(100, 0.5);

      //-----------------------------------------------------------------------------------------

      fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin] = new TF1(Form("fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin]->SetParameters(100, 0.5);

      fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin] = new TF1(Form("fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi_US_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin]->SetParameters(100, 0.5);

      fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin] = new TF1(Form("fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin]->SetParameters(100, 0.5);


      fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin] = new TF1(Form("fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US_LS_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin]->SetParameters(100, 0.5);

      fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin] = new TF1(Form("fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi_US_LS_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin]->SetParameters(100, 0.5);

      fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin] = new TF1(Form("fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US_LS_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin]->SetParameters(100, 0.5);


      fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin] = new TF1(Form("fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin]->SetParameters(100, 0.5);

      fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin] = new TF1(Form("fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin]->SetParameters(100, 0.5);

      fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin] = new TF1(Form("fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin]->SetParameters(100, 0.5);

      //-----------------------------------------------------------------------------------------

      fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US[delta_eta_bin] = new TF1(Form("fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US[delta_eta_bin]->SetParameters(100, 0.5);

      fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi_US[delta_eta_bin] = new TF1(Form("fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi_US_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi_US[delta_eta_bin]->SetParameters(100, 0.5);

      fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US[delta_eta_bin] = new TF1(Form("fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US[delta_eta_bin]->SetParameters(100, 0.5);


      fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin] = new TF1(Form("fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US_LS_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin]->SetParameters(100, 0.5);

      fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin] = new TF1(Form("fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi_US_LS_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin]->SetParameters(100, 0.5);

      fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin] = new TF1(Form("fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US_LS_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin]->SetParameters(100, 0.5);


      fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi[delta_eta_bin] = new TF1(Form("fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi[delta_eta_bin]->SetParameters(100, 0.5);

      fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi[delta_eta_bin] = new TF1(Form("fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi[delta_eta_bin]->SetParameters(100, 0.5);

      fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi[delta_eta_bin] = new TF1(Form("fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_%i", delta_eta_bin), "[0]*(1 + [1]*x)", -1, 1);
      fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi[delta_eta_bin]->SetParameters(100, 0.5);
    }

  }

  //for weighted average of uncorrected errors (absolute errors)

  float SysErr_sum_from_fits_topo_cuts = 0;
  float SysErr_sum_from_fits_of_w_topo_cuts = 0;

  float SysErr_sum_from_fits_pT_cut = 0;
  float SysErr_sum_from_fits_of_w_pT_cut = 0;

  float SysErr_sum_from_fits_DCA_cut = 0;
  float SysErr_sum_from_fits_of_w_DCA_cut = 0;

  //----------------------------------

  float SysErr_sum_topo_cuts = 0;
  float SysErr_sum_of_w_topo_cuts = 0;

  float SysErr_sum_pT_cut = 0;
  float SysErr_sum_of_w_pT_cut = 0;

  float SysErr_sum_DCA_cut = 0;
  float SysErr_sum_of_w_DCA_cut = 0;

  float SysErr_sum_background = 0;
  float SysErr_sum_of_w_background = 0;

  float SysErr_sum_ME = 0;
  float SysErr_sum_of_w_ME = 0;

  //----------------------------------------------------------------------


  TFile *LLbarOutFile; //output file to store production plane histograms;

  if(cut_type == 0)
  {
    //Run12
    //work file with any recent updates
    LLbarOutFile = new TFile(Form("./output/ProdPlane/%i/ProdPlane_Lambda_ana_cuts.root", year), "read"); //analysis cuts
    
    //--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  }
  else if(cut_type == 1)
  {
    LLbarOutFile = new TFile(Form("./output/ProdPlane/%i/ProdPlane_Lambda_topo_tight_cuts.root", year), "read");

  }
  else if(cut_type == 2)
  {
    LLbarOutFile = new TFile(Form("./output/ProdPlane/%i/ProdPlane_Lambda_pT_tight_cuts.root", year), "read");

  }
  else if(cut_type == 3)
  {
    LLbarOutFile = new TFile(Form("./output/ProdPlane/%i/ProdPlane_Lambda_topo_tight_cuts.root", year), "read");

  }
  else
  {
    cout<<"Wrong cut type"<<endl;

    return false;
  }


  //data histograms
  TH2F *L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_hist = (TH2F*)LLbarOutFile->Get("L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_hist");

  TH2F *L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist = (TH2F*)LLbarOutFile->Get("L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist");

  //--------------------------------------------------------

  TH2F *L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_hist = (TH2F*)LLbarOutFile->Get("L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_hist");

  TH2F *L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist = (TH2F*)LLbarOutFile->Get("L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist");

  //--------------------------------------------------------

  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_hist = (TH2F*)LLbarOutFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_hist");

  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist = (TH2F*)LLbarOutFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist");

  //--------------------------------------------------------

  //mixed event
  TH2F *L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist = (TH2F*)LLbarOutFile->Get("L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist");

  TH2F *L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist = (TH2F*)LLbarOutFile->Get("L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist");

  //--------------------------------------------------------

  TH2F *L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist = (TH2F*)LLbarOutFile->Get("L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist");

  TH2F *L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist = (TH2F*)LLbarOutFile->Get("L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist");

  //--------------------------------------------------------

  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist = (TH2F*)LLbarOutFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist");

  TH2F *L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist = (TH2F*)LLbarOutFile->Get("L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist");

  //________________________________________________________________________________________



  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetLegendFont(42);

  //polarization graph for Delta y
  //first graph is for |Delta y| < 0.5, second is for 0.5 < |Delta y| < 2.0
  TGraphErrors *PolarizationGraph_delta_eta_delta_phi[2];
  TGraphErrors *PolarizationGraph_delta_eta_delta_phi_sys_err[2];
  TGraphErrors *PolarizationGraph_delta_eta_delta_phi_sys_err_average[2];


  TGraphErrors *PolarizationGraph_delta_eta_delta_phi_epsilon[2];
  TGraphErrors *PolarizationGraph_delta_eta_delta_phi_epsilon_sys_err[2];
  TGraphErrors *PolarizationGraph_delta_eta_delta_phi_epsilon_sys_err_average[2];


  TCanvas *PolarizationGraph_delta_eta_delta_phi_can = new TCanvas("PolarizationGraph_delta_eta_delta_phi_can", "PolarizationGraph_delta_eta_delta_phi_can", 2000, 1200);
  PolarizationGraph_delta_eta_delta_phi_can->Divide(2,1);

  TH1F *DefaultHist = new TH1F("DefaultHist", "DefaultHist", 10, -0.5, 0.5);
  DefaultHist->GetXaxis()->SetTitle("P_{#Lambda_{1}#Lambda_{2}}");
  DefaultHist->GetXaxis()->CenterTitle();
  DefaultHist->GetXaxis()->SetRangeUser(-0.42, 0.42);
  //DefaultHist->GetYaxis()->SetRangeUser(0.01, 4.99); //with K0s
  DefaultHist->GetYaxis()->SetRangeUser(0.01, 3.99); //without K0s
  //DefaultHist->GetYaxis()->SetNdivisions(505);
  //DefaultHist->GetYaxis()->SetNdivisions(06); //with K0s
  DefaultHist->GetYaxis()->SetNdivisions(06); //without K0s
  //DefaultHist->GetYaxis()->SetNdivisions(-5);
  DefaultHist->GetYaxis()->ChangeLabel(1,-1,-1,-1,-1,-1,"#Lambda#bar{#Lambda}");
  DefaultHist->GetYaxis()->ChangeLabel(2,-1,-1,-1,-1,-1,"#Lambda#Lambda");
  DefaultHist->GetYaxis()->ChangeLabel(3,-1,-1,-1,-1,-1,"#bar{#Lambda}#bar{#Lambda}");
  //DefaultHist->GetYaxis()->ChangeLabel(4,-1,-1,-1,-1,-1,"K_{s}^{0}K_{s}^{0}");

  TLine *ZeroLine_eta = new TLine(0,0.01,0,3.99);
  ZeroLine_eta->SetLineStyle(9);
  ZeroLine_eta->SetLineColor(1);
  //ZeroLine_eta->Draw("same");

  //---------------------------------------


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
  //Polarization_text->Draw("same");

  TPaveText *bin_text[2];

  TPaveText *Delta_eta_text_1 = new TPaveText(0.15, 0.16, 0.35, 0.25, "NDC");
  Delta_eta_text_1->SetTextFont(43);
  Delta_eta_text_1->SetTextSize(33);
  Delta_eta_text_1->SetFillColorAlpha(0, 0.01);
  Delta_eta_text_1->AddText("|#Delta#it{y}| < 0.5, |#Delta#phi| < #pi/3");

  TPaveText *Delta_eta_text_2 = new TPaveText(0.15, 0.11, 0.45, 0.2, "NDC");
  Delta_eta_text_2->SetTextFont(43);
  Delta_eta_text_2->SetTextSize(33);
  Delta_eta_text_2->SetFillColorAlpha(0, 0.01);
  Delta_eta_text_2->AddText("0.5 < |#Delta#it{y}| < 2.0 or |#Delta#phi| > #pi/3");


  out_file->cd();

  for( unsigned int delta_eta_bin = 1; delta_eta_bin < 3; delta_eta_bin++)
  {
    //create polarization graphs, defined earlier

    PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1] = new TGraphErrors(3);
    PolarizationGraph_delta_eta_delta_phi_sys_err[delta_eta_bin-1] = new TGraphErrors(3);
    PolarizationGraph_delta_eta_delta_phi_sys_err_average[delta_eta_bin-1] = new TGraphErrors(3);

    //-----------------------------------------------------------------------------------------

    TCanvas *L0_L0bar_cosThetaProdPlane_eta_no_corr_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_eta_no_corr_can_%i", delta_eta_bin), Form("L0_L0bar_cosThetaProdPlane_eta_no_corr_can_%i", delta_eta_bin), 1200, 1000);

    L0_L0bar_cosThetaProdPlane_eta_no_corr_can->cd();

    TH1D *L0_L0bar_cosThetaProdPlane_eta_US_hist = L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_hist->ProjectionX( Form("proj_US_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);



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



     //ME here just for plotting, used loser
    TH1D *L0_L0bar_cosThetaProdPlane_eta_ME_hist = L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist->ProjectionX( Form("proj_US_ME_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);

    L0_L0bar_cosThetaProdPlane_eta_ME_hist->SetMarkerSize(1.5);
    L0_L0bar_cosThetaProdPlane_eta_ME_hist->SetMarkerStyle(24);
    L0_L0bar_cosThetaProdPlane_eta_ME_hist->SetMarkerColor(1);
    L0_L0bar_cosThetaProdPlane_eta_ME_hist->SetLineColor(1);
    L0_L0bar_cosThetaProdPlane_eta_ME_hist->Sumw2();
    L0_L0bar_cosThetaProdPlane_eta_ME_hist->Scale(nLLbar/L0_L0bar_cosThetaProdPlane_eta_ME_hist->Integral()); //scale ME to US
    L0_L0bar_cosThetaProdPlane_eta_ME_hist->Scale(1./L0_L0bar_cosThetaProdPlane_eta_ME_hist->GetXaxis()->GetBinWidth(1)); //US is scaled by bin width, ME needs to be scaled as well
    L0_L0bar_cosThetaProdPlane_eta_ME_hist->SetMinimum(0);


    L0_L0bar_cosThetaProdPlane_eta_US_hist->GetYaxis()->SetRangeUser(0, L0_L0bar_cosThetaProdPlane_eta_ME_hist->GetMaximum()*1.05 );
    L0_L0bar_cosThetaProdPlane_eta_US_hist->Draw("p e");

    L0_L0bar_cosThetaProdPlane_eta_ME_hist->Draw("p e same");


    cout<<endl;
    //cout<<"L-Lbar ME/SE ratio = "<<L0_L0bar_cosThetaProdPlane_eta_ME_hist->GetBinContent(1)/L0_L0bar_cosThetaProdPlane_eta_ME_hist->GetBinError(1)/(L0_L0bar_cosThetaProdPlane_eta_US_hist->GetBinContent(1)/L0_L0bar_cosThetaProdPlane_eta_US_hist->GetBinError(1))<<endl;
    cout<<"L-Lbar ME/SE ratio = "<<L0_L0bar_cosThetaProdPlane_eta_US_hist->GetBinError(1)/L0_L0bar_cosThetaProdPlane_eta_ME_hist->GetBinError(1)<<endl;
    //cout<<"L-Lbar ME/SE ratio = "<<L0_L0bar_cosThetaProdPlane_eta_ME_hist->GetSumOfWeights()/L0_L0bar_cosThetaProdPlane_eta_US_hist->GetSumOfWeights()<<endl;
    cout<<endl;

    TF1 *fitL0_L0bar_US_ThetaStar_no_corr_ME = new TF1("fitL0_L0bar_US_ThetaStar_no_corr_ME", "[0]*(1 + [1]*x)", -1, 1);
    fitL0_L0bar_US_ThetaStar_no_corr_ME->SetParameters(100, 0.5);

    L0_L0bar_cosThetaProdPlane_eta_ME_hist->Fit(fitL0_L0bar_US_ThetaStar_no_corr_ME, "s i 0 r");

    float P_L0_L0bar_no_corr_ME = fitL0_L0bar_US_ThetaStar_no_corr_ME->GetParameter(1)/(L0_alpha*L0bar_alpha);
    float P_L0_L0bar_no_corr_ME_err = fitL0_L0bar_US_ThetaStar_no_corr_ME->GetParError(1)/(L0_alpha*L0bar_alpha);


    TH1D *L0_L0bar_cosThetaProdPlane_eta_LS_hist = L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist->ProjectionX( Form("proj_LS_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);

    L0_L0bar_cosThetaProdPlane_eta_LS_hist->SetMarkerSize(1.5);
    L0_L0bar_cosThetaProdPlane_eta_LS_hist->SetMarkerStyle(21);
    L0_L0bar_cosThetaProdPlane_eta_LS_hist->SetMarkerColor(kBlue);
    double nLLbar_back = L0_L0bar_cosThetaProdPlane_eta_LS_hist->Integral();
    L0_L0bar_cosThetaProdPlane_eta_LS_hist->Sumw2();
    L0_L0bar_cosThetaProdPlane_eta_LS_hist->Scale(1./L0_L0bar_cosThetaProdPlane_eta_LS_hist->GetXaxis()->GetBinWidth(1));
    //L0_L0bar_cosThetaProdPlane_eta_LS_hist->Divide(L0_L0bar_cosThetaProdPlane_eta_eff);
    L0_L0bar_cosThetaProdPlane_eta_LS_hist->Draw("p e same");


    TH1D *L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist = L0_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist->ProjectionX( Form("proj_LS_ME_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);

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


    //cout<<nLLbar<<endl;
    //cout<<nLLbar/fitL0_L0bar_US_ThetaStar_no_corr->GetParameter(0)<<endl;
    //cout<<endl;

    float P_L0_L0bar_no_corr = fitL0_L0bar_US_ThetaStar_no_corr->GetParameter(1)/(L0_alpha*L0bar_alpha);
    float P_L0_L0bar_no_corr_err = fitL0_L0bar_US_ThetaStar_no_corr->GetParError(1)/(L0_alpha*L0bar_alpha);

    fitL0_L0bar_US_ThetaStar_no_corr->SetLineColor(1);
    //fitL0_L0bar_US_ThetaStar_no_corr->Draw("same");

    TLegend *L0_L0bar_leg_no_corr = new TLegend(0.15, 0.3, 0.45, 0.55);
    L0_L0bar_leg_no_corr->AddEntry(L0_L0bar_cosThetaProdPlane_eta_US_hist, "(US-US) p#pi");
    L0_L0bar_leg_no_corr->AddEntry(L0_L0bar_cosThetaProdPlane_eta_ME_hist, "(US-US) ME");
    L0_L0bar_leg_no_corr->AddEntry(L0_L0bar_cosThetaProdPlane_eta_LS_hist, "Combinatorial bckg.");
    L0_L0bar_leg_no_corr->AddEntry(L0_L0bar_cosThetaProdPlane_eta_ME_LS_hist, "Bckg. ME");
    //L0_L0bar_leg_no_corr->AddEntry(fitL0_L0bar_US_ThetaStar_no_corr, "Linear fit to US");
    L0_L0bar_leg_no_corr->SetBorderSize(0);
    L0_L0bar_leg_no_corr->SetFillColorAlpha(0, 0.01);
    L0_L0bar_leg_no_corr->Draw("same");

    TPaveText *L0_L0bar_text_no_corr = new TPaveText(0.5, 0.3, 0.85, 0.65, "NDC");
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
    L0_L0bar_text_no_corr->AddText("DCA_{#Lambda} < 1 cm");
    L0_L0bar_text_no_corr->AddText(Form("P = %.2f #pm %.2f", P_L0_L0bar_no_corr, fabs(P_L0_L0bar_no_corr_err)));
    L0_L0bar_text_no_corr->AddText(Form("P_{ME} = %.3f #pm %.3f", P_L0_L0bar_no_corr_ME, fabs(P_L0_L0bar_no_corr_ME_err)));
    L0_L0bar_text_no_corr->SetFillColorAlpha(0, 0.01);
    L0_L0bar_text_no_corr->Draw("same");

    if(delta_eta_bin == 1) Delta_eta_text_1->Draw("same");
    if(delta_eta_bin == 2) Delta_eta_text_2->Draw("same");

    L0_L0bar_cosThetaProdPlane_eta_no_corr_can->SaveAs(Form("./plots/Lambda/L_correlations_delta_eta_delta_phi/L0_L0bar_cosThetaProdPlane_eta_no_corr_less_delta_eta_delta_phi_%i.png", delta_eta_bin));

    //----------------------------------------------------------------------------------------------------

    TCanvas *L0_L0bar_cosThetaProdPlane_eta_can = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_eta_can_%i", delta_eta_bin), Form("L0_L0bar_cosThetaProdPlane_eta_can_%i", delta_eta_bin), 1200, 1000);

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

    float chi_2_LLbar = L0_L0bar_cosThetaProdPlane_eta_US_hist->Chisquare(fitL0_L0bar_US_ThetaStar);

    cout<<endl;
    cout<<"Chi2(LLbar) = "<<chi_2_LLbar<<endl;
    cout<<"Chi2(LLbar)/NDF = "<<chi_2_LLbar/8.<<endl; //NDF should be 10 - 2 (Npoints - N free params)
    cout<<endl;

    float P_L0_L0bar = fitL0_L0bar_US_ThetaStar->GetParameter(1)/(L0_alpha*L0bar_alpha);
    float P_L0_L0bar_err = fitL0_L0bar_US_ThetaStar->GetParError(1)/(L0_alpha*L0bar_alpha);

    fitL0_L0bar_US_ThetaStar->SetLineColor(1);
    fitL0_L0bar_US_ThetaStar->Draw("same");

    //background
    TF1 *fitL0_L0bar_US_LS_ThetaStar = new TF1("fitL0_L0bar_US_LS_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
    fitL0_L0bar_US_LS_ThetaStar->SetParameters(100, 0.5);

    //fit_res_gaUS_LS_wrong_sign = L_inv_mass_US_LS->Fit(fitGaUS_LSsBack, "s i 0", "", 1.07, 1.4);
    L0_L0bar_cosThetaProdPlane_eta_LS_hist->Fit(fitL0_L0bar_US_LS_ThetaStar, "s i 0 r");

    float chi_2_LLbar_back = L0_L0bar_cosThetaProdPlane_eta_LS_hist->Chisquare(fitL0_L0bar_US_LS_ThetaStar);

    cout<<endl;
    cout<<"Chi2(LLbar back) = "<<chi_2_LLbar_back<<endl;
    cout<<"Chi2(LLbar back)/NDF = "<<chi_2_LLbar_back/8.<<endl; //NDF should be 10 - 2 (Npoints - N free params)
    cout<<endl;

    float P_L0_L0bar_back = fitL0_L0bar_US_LS_ThetaStar->GetParameter(1)/(L0_alpha*L0bar_alpha);
    float P_L0_L0bar_back_err = fitL0_L0bar_US_LS_ThetaStar->GetParError(1)/(L0_alpha*L0bar_alpha);

    fitL0_L0bar_US_LS_ThetaStar->SetLineColor(1);
    fitL0_L0bar_US_LS_ThetaStar->SetLineStyle(7);
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

    TLegend *L0_L0bar_leg = new TLegend(0.15, 0.3, 0.45, 0.55);
    L0_L0bar_leg->AddEntry(L0_L0bar_cosThetaProdPlane_eta_US_hist, "(US-US) p#pi");
    L0_L0bar_leg->AddEntry(fitL0_L0bar_US_ThetaStar, "(US-US) Fit");
    L0_L0bar_leg->AddEntry(L0_L0bar_cosThetaProdPlane_eta_LS_hist, "Combinatorial bckg.");
    L0_L0bar_leg->AddEntry(fitL0_L0bar_US_LS_ThetaStar, "Bckg. fit");
    //L0_L0bar_leg->AddEntry(fitL0_L0bar_US_ThetaStar_no_corr, "Linear fit to US");
    L0_L0bar_leg->SetBorderSize(0);
    L0_L0bar_leg->SetFillColorAlpha(0, 0.01);
    L0_L0bar_leg->Draw("same");

    if(delta_eta_bin == 1) Delta_eta_text_1->Draw("same");
    if(delta_eta_bin == 2) Delta_eta_text_2->Draw("same");

    L0_L0bar_cosThetaProdPlane_eta_can->SaveAs(Form("./plots/Lambda/L_correlations_delta_eta_delta_phi/L0_L0bar_cosThetaProdPlane_eta_less_delta_eta_delta_phi_%i.png", delta_eta_bin));

    //----------------------------------------------------------------------

    TCanvas *L0_L0bar_cosThetaProdPlane_eta_can_2 = new TCanvas(Form("L0_L0bar_cosThetaProdPlane_eta_can_2_%i", delta_eta_bin), Form("L0_L0bar_cosThetaProdPlane_eta_can_2_%i", delta_eta_bin), 1200, 1000);

    L0_L0bar_cosThetaProdPlane_eta_can_2->cd();

    L0_L0bar_cosThetaProdPlane_eta_US_hist->Add(L0_L0bar_cosThetaProdPlane_eta_LS_hist, -1);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->SetMinimum(0);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->Write(Form("L0_L0bar_cosThetaProdPlane_delta_y_delta_phi_%i", delta_eta_bin-1));
    L0_L0bar_cosThetaProdPlane_eta_US_hist->Draw("p e");


    TF1 *fitL0_L0bar_US_ThetaStar_2 = new TF1("fitL0_L0bar_US_ThetaStar_2", "[0]*(1 + [1]*x)", -1, 1);
    fitL0_L0bar_US_ThetaStar_2->SetParameters(100, 0.5);

    //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
    L0_L0bar_cosThetaProdPlane_eta_US_hist->Fit(fitL0_L0bar_US_ThetaStar_2, "s i 0 r");

    //store fit result for systematic errors
    //tight topo cuts
    if(cut_type == 1)
    {
      SysErrCutsTopo->cd();

      fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParameters(fitL0_L0bar_US_ThetaStar->GetParameter(0), fitL0_L0bar_US_ThetaStar->GetParameter(1));
      fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParError(0, fitL0_L0bar_US_ThetaStar->GetParError(0));
      fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParError(1, fitL0_L0bar_US_ThetaStar->GetParError(1));
      fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->Write();


      fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParameters(fitL0_L0bar_US_LS_ThetaStar->GetParameter(0), fitL0_L0bar_US_LS_ThetaStar->GetParameter(1));
      fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParError(0, fitL0_L0bar_US_LS_ThetaStar->GetParError(0));
      fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParError(1, fitL0_L0bar_US_LS_ThetaStar->GetParError(1));
      fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->Write();


      fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParameters(fitL0_L0bar_US_ThetaStar_2->GetParameter(0), fitL0_L0bar_US_ThetaStar_2->GetParameter(1));
      fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParError(0, fitL0_L0bar_US_ThetaStar_2->GetParError(0));
      fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParError(1, fitL0_L0bar_US_ThetaStar_2->GetParError(1));
      fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->Write();
    }

    //tight pT cuts
    if(cut_type == 2)
    {
      SysErrCutsPt->cd();

      fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParameters(fitL0_L0bar_US_ThetaStar->GetParameter(0), fitL0_L0bar_US_ThetaStar->GetParameter(1));
      fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParError(0, fitL0_L0bar_US_ThetaStar->GetParError(0));
      fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParError(1, fitL0_L0bar_US_ThetaStar->GetParError(1));
      fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->Write();


      fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParameters(fitL0_L0bar_US_LS_ThetaStar->GetParameter(0), fitL0_L0bar_US_LS_ThetaStar->GetParameter(1));
      fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParError(0, fitL0_L0bar_US_LS_ThetaStar->GetParError(0));
      fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParError(1, fitL0_L0bar_US_LS_ThetaStar->GetParError(1));
      fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->Write();

      fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParameters(fitL0_L0bar_US_ThetaStar_2->GetParameter(0), fitL0_L0bar_US_ThetaStar_2->GetParameter(1));
      fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParError(0, fitL0_L0bar_US_ThetaStar_2->GetParError(0));
      fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParError(1, fitL0_L0bar_US_ThetaStar_2->GetParError(1));
      fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->Write();
    }

    //tight pT cuts
    if(cut_type == 3)
    {
      SysErrCutsDCA->cd();

      fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParameters(fitL0_L0bar_US_ThetaStar->GetParameter(0), fitL0_L0bar_US_ThetaStar->GetParameter(1));
      fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParError(0, fitL0_L0bar_US_ThetaStar->GetParError(0));
      fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParError(1, fitL0_L0bar_US_ThetaStar->GetParError(1));
      fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->Write();


      fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParameters(fitL0_L0bar_US_LS_ThetaStar->GetParameter(0), fitL0_L0bar_US_LS_ThetaStar->GetParameter(1));
      fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParError(0, fitL0_L0bar_US_LS_ThetaStar->GetParError(0));
      fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParError(1, fitL0_L0bar_US_LS_ThetaStar->GetParError(1));
      fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->Write();

      fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParameters(fitL0_L0bar_US_ThetaStar_2->GetParameter(0), fitL0_L0bar_US_ThetaStar_2->GetParameter(1));
      fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParError(0, fitL0_L0bar_US_ThetaStar_2->GetParError(0));
      fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParError(1, fitL0_L0bar_US_ThetaStar_2->GetParError(1));
      fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->Write();
    }

    float P_L0_L0bar_2 = fitL0_L0bar_US_ThetaStar_2->GetParameter(1)/(L0_alpha*L0bar_alpha);
    float P_L0_L0bar_err_2 = fitL0_L0bar_US_ThetaStar_2->GetParError(1)/(L0_alpha*L0bar_alpha);

    fitL0_L0bar_US_ThetaStar_2->SetLineColor(1);
    fitL0_L0bar_US_ThetaStar_2->Draw("same");


    //calculate total systematic uncertainty for L-Lbar

    //slope
    //float SysErrSlope_L0_L0bar = SysErrSlope_delta_eta_delta_phi_hist[delta_eta_bin-1]->GetBinContent(1);
    float SysErrSlope_L0_L0bar = fabs(ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin-1]->GetBinContent(1))/fabs(P_L0_L0bar_2);

    SysErr_sum_ME += fabs(ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin-1]->GetBinContent(1))/fabs(ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin-1]->GetBinError(1));
    SysErr_sum_of_w_ME += 1./fabs(ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin-1]->GetBinError(1));


    //--------------------------------------------------------------

    //background subtraction

    pair<float,float> P_L0_L0bar_from_fits_pair = Get_P_LL_from_fit(fitL0_L0bar_US_ThetaStar,fitL0_L0bar_US_LS_ThetaStar,0);

    float P_L0_L0bar_from_fits = P_L0_L0bar_from_fits_pair.first;
    float P_L0_L0bar_from_fits_err = P_L0_L0bar_from_fits_pair.second;

    //--------------------------------------------------------------

    //statistical error correction for systematic error
    float SysErrBackground_L0_L0bar_corr = sqrt( fabs( P_L0_L0bar_err_2*P_L0_L0bar_err_2 - P_L0_L0bar_from_fits_err*P_L0_L0bar_from_fits_err ) );

    float SysErrBackground_L0_L0bar_work = ( fabs( P_L0_L0bar_2 - P_L0_L0bar_from_fits) - SysErrBackground_L0_L0bar_corr )/fabs(P_L0_L0bar_from_fits);

    float SysErrBackground_L0_L0bar = 0;

    if( SysErrBackground_L0_L0bar_work > 0 ) SysErrBackground_L0_L0bar = SysErrBackground_L0_L0bar_work; //store sys. err. only if it is larger than statistical fluctuations

    float SysErrBackground_L0_L0bar_no_corr = fabs( P_L0_L0bar_2 - P_L0_L0bar_from_fits)/fabs(P_L0_L0bar_from_fits);

    //SysErr_sum_background += SysErrBackground_L0_L0bar_no_corr*fabs(P_L0_L0bar_2)/SysErrBackground_L0_L0bar_corr;
    //SysErr_sum_background += SysErrBackground_L0_L0bar*fabs(P_L0_L0bar_2)/SysErrBackground_L0_L0bar_corr;
    //SysErr_sum_of_w_background += 1./SysErrBackground_L0_L0bar_corr;

    SysErr_sum_background += SysErrBackground_L0_L0bar*fabs(P_L0_L0bar_from_fits)/fabs(P_L0_L0bar_from_fits_err);
    SysErr_sum_of_w_background += 1./fabs(P_L0_L0bar_from_fits_err);

    //----------------------------------------------------------------------------------------------------------------------------

    //cuts variation
    float SysErr_tight_cuts_L0_L0bar_from_fits = 0;

    float SysErr_tight_topo_cuts_L0_L0bar_from_fits = 0;
    float SysErr_tight_pT_cuts_L0_L0bar_from_fits = 0;
    float SysErr_tight_DCA_cuts_L0_L0bar_from_fits = 0;


    float SysErr_tight_cuts_L0_L0bar_from_fits_no_corr = 0;

    float SysErr_tight_topo_cuts_L0_L0bar_from_fits_no_corr = 0;
    float SysErr_tight_pT_cuts_L0_L0bar_from_fits_no_corr = 0;
    float SysErr_tight_DCA_cuts_L0_L0bar_from_fits_no_corr = 0;

    //------------------------------------------------------

    float SysErr_tight_cuts_L0_L0bar = 0;

    float SysErr_tight_topo_cuts_L0_L0bar = 0;
    float SysErr_tight_pT_cuts_L0_L0bar = 0;
    float SysErr_tight_DCA_cuts_L0_L0bar = 0;


    float SysErr_tight_cuts_L0_L0bar_no_corr = 0;

    float SysErr_tight_topo_cuts_L0_L0bar_no_corr = 0;
    float SysErr_tight_pT_cuts_L0_L0bar_no_corr = 0;
    float SysErr_tight_DCA_cuts_L0_L0bar_no_corr = 0;


    if(cut_type == 0)
    {
      //independent fit of US-US and background

      pair<float,float> P_L0_L0bar_tight_topo_cuts_from_fits_signal_pair = Get_P_LL_from_fit(fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1], fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1], 0);

      float P_L0_L0bar_tight_topo_cuts_from_fits_signal = P_L0_L0bar_tight_topo_cuts_from_fits_signal_pair.first;
      float P_L0_L0bar_tight_topo_cuts_from_fits_signal_err = P_L0_L0bar_tight_topo_cuts_from_fits_signal_pair.second;


      //statistical error correction for systematic error
      float SysErr_tight_topo_cuts_L0_L0bar_from_fits_corr = sqrt( fabs( P_L0_L0bar_from_fits_err*P_L0_L0bar_from_fits_err - P_L0_L0bar_tight_topo_cuts_from_fits_signal_err*P_L0_L0bar_tight_topo_cuts_from_fits_signal_err ) );

      float SysErr_tight_topo_cuts_L0_L0bar_from_fits_work = ( fabs( P_L0_L0bar_from_fits - P_L0_L0bar_tight_topo_cuts_from_fits_signal) - SysErr_tight_topo_cuts_L0_L0bar_from_fits_corr )/fabs(P_L0_L0bar_from_fits);

      if( SysErr_tight_topo_cuts_L0_L0bar_from_fits_work > 0 ) SysErr_tight_topo_cuts_L0_L0bar_from_fits = SysErr_tight_topo_cuts_L0_L0bar_from_fits_work; //store sys. err. only if it is larger than statistical fluctuations

      SysErr_tight_topo_cuts_L0_L0bar_from_fits_no_corr = fabs( P_L0_L0bar_from_fits - P_L0_L0bar_tight_topo_cuts_from_fits_signal)/fabs(P_L0_L0bar_from_fits);


      //SysErr_sum_from_fits_topo_cuts += SysErr_tight_topo_cuts_L0_L0bar_from_fits*fabs(P_L0_L0bar_from_fits)/SysErr_tight_topo_cuts_L0_L0bar_from_fits_corr;
      //SysErr_sum_from_fits_of_w_topo_cuts += 1./SysErr_tight_topo_cuts_L0_L0bar_from_fits_corr;

      SysErr_sum_from_fits_topo_cuts += SysErr_tight_topo_cuts_L0_L0bar_from_fits*fabs(P_L0_L0bar_from_fits)/fabs(P_L0_L0bar_from_fits_err);
      SysErr_sum_from_fits_of_w_topo_cuts += 1./fabs(P_L0_L0bar_from_fits_err);

      //----------------------------------------------------------------

      //fit after background subtraction

      float P_L0_L0bar_tight_topo_cuts = fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->GetParameter(1)/(L0_alpha*L0bar_alpha);
      float P_L0_L0bar_tight_topo_cuts_err = fitL0_L0bar_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->GetParError(1)/(L0_alpha*L0bar_alpha);

      //statistical error correction for systematic error
      float SysErr_tight_topo_cuts_L0_L0bar_corr = sqrt( fabs( P_L0_L0bar_err_2*P_L0_L0bar_err_2 - P_L0_L0bar_tight_topo_cuts_err*P_L0_L0bar_tight_topo_cuts_err ) );

      float SysErr_tight_topo_cuts_L0_L0bar_work = ( fabs( P_L0_L0bar_2 - P_L0_L0bar_tight_topo_cuts) - SysErr_tight_topo_cuts_L0_L0bar_corr )/fabs(P_L0_L0bar_2);

      if( SysErr_tight_topo_cuts_L0_L0bar_work > 0 ) SysErr_tight_topo_cuts_L0_L0bar = SysErr_tight_topo_cuts_L0_L0bar_work; //store sys. err. only if it is larger than statistical fluctuations

      SysErr_tight_topo_cuts_L0_L0bar_no_corr = fabs( P_L0_L0bar_2 - P_L0_L0bar_tight_topo_cuts)/fabs(P_L0_L0bar_2);


      //SysErr_sum_topo_cuts += SysErr_tight_topo_cuts_L0_L0bar*fabs(P_L0_L0bar_2)/SysErr_tight_topo_cuts_L0_L0bar_corr;
      //SysErr_sum_of_w_topo_cuts += 1./SysErr_tight_topo_cuts_L0_L0bar_corr;

      SysErr_sum_topo_cuts += SysErr_tight_topo_cuts_L0_L0bar*fabs(P_L0_L0bar_2)/fabs(P_L0_L0bar_err_2);
      SysErr_sum_of_w_topo_cuts += 1./fabs(P_L0_L0bar_err_2);


      //___________________________________________________________________________________________________________________________________

      //independent fit of US-US and background
      pair<float,float> P_L0_L0bar_tight_pT_cuts_from_fits_signal_pair = Get_P_LL_from_fit(fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1], fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1], 0);

      float P_L0_L0bar_tight_pT_cuts_from_fits_signal = P_L0_L0bar_tight_pT_cuts_from_fits_signal_pair.first;
      float P_L0_L0bar_tight_pT_cuts_from_fits_signal_err = P_L0_L0bar_tight_pT_cuts_from_fits_signal_pair.second;


      //statistical error correction for systematic error
      float SysErr_tight_pT_cuts_L0_L0bar_from_fits_corr = sqrt( fabs( P_L0_L0bar_from_fits_err*P_L0_L0bar_from_fits_err - P_L0_L0bar_tight_pT_cuts_from_fits_signal_err*P_L0_L0bar_tight_pT_cuts_from_fits_signal_err ) );

      float SysErr_tight_pT_cuts_L0_L0bar_from_fits_work = ( fabs( P_L0_L0bar_from_fits - P_L0_L0bar_tight_pT_cuts_from_fits_signal) - SysErr_tight_pT_cuts_L0_L0bar_from_fits_corr )/fabs(P_L0_L0bar_from_fits);

      if( SysErr_tight_pT_cuts_L0_L0bar_from_fits_work > 0 ) SysErr_tight_pT_cuts_L0_L0bar_from_fits = SysErr_tight_pT_cuts_L0_L0bar_from_fits_work; //store sys. err. only if it is larger than statistical fluctuations

      SysErr_tight_pT_cuts_L0_L0bar_from_fits_no_corr = fabs( P_L0_L0bar_from_fits - P_L0_L0bar_tight_pT_cuts_from_fits_signal)/fabs(P_L0_L0bar_from_fits);


      //SysErr_sum_from_fits_pT_cut += SysErr_tight_pT_cuts_L0_L0bar_from_fits*fabs(P_L0_L0bar_from_fits)/SysErr_tight_pT_cuts_L0_L0bar_from_fits_corr;
      //SysErr_sum_from_fits_of_w_pT_cut += 1./SysErr_tight_pT_cuts_L0_L0bar_from_fits_corr;

      SysErr_sum_from_fits_pT_cut += SysErr_tight_pT_cuts_L0_L0bar_from_fits*fabs(P_L0_L0bar_from_fits)/fabs(P_L0_L0bar_from_fits_err);
      SysErr_sum_from_fits_of_w_pT_cut += 1./fabs(P_L0_L0bar_from_fits_err);

      //----------------------------------------------------------------


      float P_L0_L0bar_tight_pT_cuts = fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->GetParameter(1)/(L0_alpha*L0bar_alpha);
      float P_L0_L0bar_tight_pT_cuts_err = fitL0_L0bar_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->GetParError(1)/(L0_alpha*L0bar_alpha);

      //statistical error correction for systematic error
      float SysErr_tight_pT_cuts_L0_L0bar_corr = sqrt( fabs( P_L0_L0bar_err_2*P_L0_L0bar_err_2 - P_L0_L0bar_tight_pT_cuts_err*P_L0_L0bar_tight_pT_cuts_err ) );

      float SysErr_tight_pT_cuts_L0_L0bar_work = ( fabs( P_L0_L0bar_2 - P_L0_L0bar_tight_pT_cuts) - SysErr_tight_pT_cuts_L0_L0bar_corr )/fabs(P_L0_L0bar_2);

      if( SysErr_tight_pT_cuts_L0_L0bar_work > 0 ) SysErr_tight_pT_cuts_L0_L0bar = SysErr_tight_pT_cuts_L0_L0bar_work; //store sys. err. only if it is larger than statistical fluctuations

      SysErr_tight_pT_cuts_L0_L0bar_no_corr = fabs( P_L0_L0bar_2 - P_L0_L0bar_tight_pT_cuts)/fabs(P_L0_L0bar_2);


      //SysErr_sum_pT_cut += SysErr_tight_pT_cuts_L0_L0bar*fabs(P_L0_L0bar_2)/SysErr_tight_pT_cuts_L0_L0bar_corr;
      //SysErr_sum_of_w_pT_cut += 1./SysErr_tight_pT_cuts_L0_L0bar_corr;

      SysErr_sum_pT_cut += SysErr_tight_pT_cuts_L0_L0bar*fabs(P_L0_L0bar_2)/fabs(P_L0_L0bar_err_2);
      SysErr_sum_of_w_pT_cut += 1./fabs(P_L0_L0bar_err_2);

      //___________________________________________________________________________________________________________________________________

      //independent fit of US-US and background
      pair<float,float> P_L0_L0bar_tight_DCA_cuts_from_fits_signal_pair = Get_P_LL_from_fit(fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1], fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1], 0);

      float P_L0_L0bar_tight_DCA_cuts_from_fits_signal = P_L0_L0bar_tight_DCA_cuts_from_fits_signal_pair.first;
      float P_L0_L0bar_tight_DCA_cuts_from_fits_signal_err = P_L0_L0bar_tight_DCA_cuts_from_fits_signal_pair.second;


      //statistical error correction for systematic error
      float SysErr_tight_DCA_cuts_L0_L0bar_from_fits_corr = sqrt( fabs( P_L0_L0bar_from_fits_err*P_L0_L0bar_from_fits_err - P_L0_L0bar_tight_DCA_cuts_from_fits_signal_err*P_L0_L0bar_tight_DCA_cuts_from_fits_signal_err ) );

      float SysErr_tight_DCA_cuts_L0_L0bar_from_fits_work = ( fabs( P_L0_L0bar_from_fits - P_L0_L0bar_tight_DCA_cuts_from_fits_signal) - SysErr_tight_DCA_cuts_L0_L0bar_from_fits_corr )/fabs(P_L0_L0bar_from_fits);

      if( SysErr_tight_DCA_cuts_L0_L0bar_from_fits_work > 0 ) SysErr_tight_DCA_cuts_L0_L0bar_from_fits = SysErr_tight_DCA_cuts_L0_L0bar_from_fits_work; //store sys. err. only if it is larger than statistical fluctuations

      SysErr_tight_DCA_cuts_L0_L0bar_from_fits_no_corr = fabs( P_L0_L0bar_from_fits - P_L0_L0bar_tight_DCA_cuts_from_fits_signal)/fabs(P_L0_L0bar_from_fits);


      //SysErr_sum_from_fits_DCA_cut += SysErr_tight_DCA_cuts_L0_L0bar_from_fits*fabs(P_L0_L0bar_from_fits)/SysErr_tight_DCA_cuts_L0_L0bar_from_fits_corr;
      //SysErr_sum_from_fits_of_w_DCA_cut += 1./SysErr_tight_DCA_cuts_L0_L0bar_from_fits_corr;

      SysErr_sum_from_fits_DCA_cut += SysErr_tight_DCA_cuts_L0_L0bar_from_fits*fabs(P_L0_L0bar_from_fits)/fabs(P_L0_L0bar_from_fits_err);
      SysErr_sum_from_fits_of_w_DCA_cut += 1./fabs(P_L0_L0bar_from_fits_err);

      //----------------------------------------------------------------


      float P_L0_L0bar_tight_DCA_cuts = fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->GetParameter(1)/(L0_alpha*L0bar_alpha);
      float P_L0_L0bar_tight_DCA_cuts_err = fitL0_L0bar_tight_DCA_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->GetParError(1)/(L0_alpha*L0bar_alpha);

      //statistical error correction for systematic error
      float SysErr_tight_DCA_cuts_L0_L0bar_corr = sqrt( fabs( P_L0_L0bar_err_2*P_L0_L0bar_err_2 - P_L0_L0bar_tight_DCA_cuts_err*P_L0_L0bar_tight_DCA_cuts_err ) );

      float SysErr_tight_DCA_cuts_L0_L0bar_work = ( fabs( P_L0_L0bar_2 - P_L0_L0bar_tight_DCA_cuts) - SysErr_tight_DCA_cuts_L0_L0bar_corr )/fabs(P_L0_L0bar_2);

      if( SysErr_tight_DCA_cuts_L0_L0bar_work > 0 ) SysErr_tight_DCA_cuts_L0_L0bar = SysErr_tight_DCA_cuts_L0_L0bar_work; //store sys. err. only if it is larger than statistical fluctuations

      SysErr_tight_DCA_cuts_L0_L0bar_no_corr = fabs( P_L0_L0bar_2 - P_L0_L0bar_tight_DCA_cuts)/fabs(P_L0_L0bar_2);


      //SysErr_sum_DCA_cut += SysErr_tight_DCA_cuts_L0_L0bar*fabs(P_L0_L0bar_2)/SysErr_tight_DCA_cuts_L0_L0bar_corr;
      //SysErr_sum_of_w_DCA_cut += 1./SysErr_tight_DCA_cuts_L0_L0bar_corr;

      SysErr_sum_DCA_cut += SysErr_tight_DCA_cuts_L0_L0bar*fabs(P_L0_L0bar_2)/fabs(P_L0_L0bar_err_2);
      SysErr_sum_of_w_DCA_cut += 1./fabs(P_L0_L0bar_err_2);

      //___________________________________________________________________________________________________________________________________

      SysErr_tight_cuts_L0_L0bar = sqrt(SysErr_tight_topo_cuts_L0_L0bar*SysErr_tight_topo_cuts_L0_L0bar + SysErr_tight_pT_cuts_L0_L0bar*SysErr_tight_pT_cuts_L0_L0bar + SysErr_tight_DCA_cuts_L0_L0bar*SysErr_tight_DCA_cuts_L0_L0bar);

      SysErr_tight_cuts_L0_L0bar_no_corr = sqrt(SysErr_tight_topo_cuts_L0_L0bar_no_corr*SysErr_tight_topo_cuts_L0_L0bar_no_corr + SysErr_tight_pT_cuts_L0_L0bar_no_corr*SysErr_tight_pT_cuts_L0_L0bar_no_corr +
                                                SysErr_tight_DCA_cuts_L0_L0bar_no_corr*SysErr_tight_DCA_cuts_L0_L0bar_no_corr);

      //----------------------------

      SysErr_tight_cuts_L0_L0bar_from_fits = sqrt(SysErr_tight_topo_cuts_L0_L0bar_from_fits*SysErr_tight_topo_cuts_L0_L0bar_from_fits + SysErr_tight_pT_cuts_L0_L0bar_from_fits*SysErr_tight_pT_cuts_L0_L0bar_from_fits +
                                                  SysErr_tight_DCA_cuts_L0_L0bar_from_fits*SysErr_tight_DCA_cuts_L0_L0bar_from_fits);

      SysErr_tight_cuts_L0_L0bar_from_fits_no_corr = sqrt(SysErr_tight_topo_cuts_L0_L0bar_from_fits_no_corr*SysErr_tight_topo_cuts_L0_L0bar_from_fits_no_corr + SysErr_tight_pT_cuts_L0_L0bar_from_fits_no_corr*SysErr_tight_pT_cuts_L0_L0bar_from_fits_no_corr +
                                                          SysErr_tight_DCA_cuts_L0_L0bar_from_fits_no_corr*SysErr_tight_DCA_cuts_L0_L0bar_from_fits_no_corr );


    }

    //total
    //float SysErrTot_L0_L0bar = sqrt( sysErr_alpha_L0_L0bar*sysErr_alpha_L0_L0bar + SysErrSlope_L0_L0bar*SysErrSlope_L0_L0bar + SysErrBackground_L0_L0bar*SysErrBackground_L0_L0bar  );
    //float SysErrTot_L0_L0bar = sqrt( sysErr_alpha_L0_L0bar*sysErr_alpha_L0_L0bar + SysErrSlope_L0_L0bar*SysErrSlope_L0_L0bar + SysErrBackground_L0_L0bar*SysErrBackground_L0_L0bar + SysErr_tight_cuts_L0_L0bar*SysErr_tight_cuts_L0_L0bar );
    float SysErrTot_L0_L0bar = sqrt( sysErr_alpha_L0_L0bar*sysErr_alpha_L0_L0bar + SysErrSlope_L0_L0bar*SysErrSlope_L0_L0bar + SysErrBackground_L0_L0bar*SysErrBackground_L0_L0bar + SysErr_tight_cuts_L0_L0bar_from_fits*SysErr_tight_cuts_L0_L0bar_from_fits );

    //--------------------------------------------------------------------------------------------------------------

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
    L0_L0bar_text_2->AddText(Form("P_{#Lambda_{1}#Lambda_{2}} = %.3f #pm %.3f", P_L0_L0bar_2, fabs(P_L0_L0bar_err_2)));
    L0_L0bar_text_2->AddText(Form("P_{#Lambda_{1}#Lambda_{2},fit} = %.3f #pm %.3f", P_L0_L0bar_from_fits, fabs(P_L0_L0bar_from_fits_err)));
    //L0_L0bar_text_2->AddText(Form("P_{#Lambda_{1}#Lambda_{2}} = %.3f #pm %.3f #pm %.3f", P_L0_L0bar_2, fabs(P_L0_L0bar_err_2), fabs(SysErrTot_L0_L0bar*P_L0_L0bar_2)));
    //L0_L0bar_text_2->AddText(Form("P_{#Lambda_{1}#Lambda_{2}} = %.3f #pm %.3f #pm %.3f", P_L0_L0bar_from_fits, fabs(P_L0_L0bar_from_fits_err), fabs(P_L0_L0bar_from_fits*P_L0_L0bar_2)));
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

    if(delta_eta_bin == 1) Delta_eta_text_1->Draw("same");
    if(delta_eta_bin == 2) Delta_eta_text_2->Draw("same");

    L0_L0bar_cosThetaProdPlane_eta_can_2->SaveAs(Form("./plots/Lambda/L_correlations_delta_eta_delta_phi/L0_L0bar_cosThetaProdPlane_eta_subtract_less_delta_eta_delta_phi_%i.png", delta_eta_bin));

    //-----------------------

    //PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->SetPoint(1, P_L0_L0bar_2, 1);
    //PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->SetPointError(1, fabs(P_L0_L0bar_err_2), 0);

    PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->SetPoint(1, P_L0_L0bar_from_fits, 1);
    PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->SetPointError(1, fabs(P_L0_L0bar_from_fits_err), 0);


    //PolarizationGraph_delta_eta_delta_phi_sys_err[delta_eta_bin-1]->SetPoint(1, P_L0_L0bar_2, 1);
    //PolarizationGraph_delta_eta_delta_phi_sys_err[delta_eta_bin-1]->SetPointError(1, fabs(SysErrTot_L0_L0bar*P_L0_L0bar_2), 0.045);

    PolarizationGraph_delta_eta_delta_phi_sys_err[delta_eta_bin-1]->SetPoint(1, P_L0_L0bar_from_fits, 1);
    PolarizationGraph_delta_eta_delta_phi_sys_err[delta_eta_bin-1]->SetPointError(1, fabs(SysErrTot_L0_L0bar*P_L0_L0bar_from_fits), 0.045);


    //____________________________________________________________________________________________________________________________________________________________________________________________________________

    TCanvas *L0_L0_cosThetaProdPlane_eta_no_corr_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_eta_no_corr_can_%i", delta_eta_bin), Form("L0_L0_cosThetaProdPlane_eta_no_corr_can_%i", delta_eta_bin), 1200, 1000);

    L0_L0_cosThetaProdPlane_eta_no_corr_can->cd();

    TH1D *L0_L0_cosThetaProdPlane_eta_US_hist = L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_hist->ProjectionX( Form("proj_US_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);

    L0_L0_cosThetaProdPlane_eta_US_hist->GetXaxis()->SetTitle("cos(#theta*)");
    L0_L0_cosThetaProdPlane_eta_US_hist->GetXaxis()->CenterTitle();
    //L0_L0_cosThetaProdPlane_eta_US_hist->GetXaxis()->SetTextSizePixels(30);
    L0_L0_cosThetaProdPlane_eta_US_hist->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0_L0_cosThetaProdPlane_eta_US_hist->GetYaxis()->CenterTitle();
    //L0_L0_cosThetaProdPlane_eta_US_hist->GetYaxis()->SetTextSizePixels(30);
    L0_L0_cosThetaProdPlane_eta_US_hist->GetYaxis()->SetMaxDigits(3);
    L0_L0_cosThetaProdPlane_eta_US_hist->SetMarkerSize(1.5);
    L0_L0_cosThetaProdPlane_eta_US_hist->SetMarkerStyle(20);
    L0_L0_cosThetaProdPlane_eta_US_hist->SetMarkerColor(kRed);
    L0_L0_cosThetaProdPlane_eta_US_hist->SetLineColor(kRed);
    double nLL = L0_L0_cosThetaProdPlane_eta_US_hist->Integral();
    L0_L0_cosThetaProdPlane_eta_US_hist->Sumw2();
    //L0_L0_cosThetaProdPlane_eta_US_hist->Divide(L0_L0_cosThetaProdPlane_eta_eff);
    L0_L0_cosThetaProdPlane_eta_US_hist->Scale(1./L0_L0_cosThetaProdPlane_eta_US_hist->GetXaxis()->GetBinWidth(1));
    //L0_L0_cosThetaProdPlane_eta_US_hist->Scale(1./L0_L0_cosThetaProdPlane_eta_US_hist->Integral());
    L0_L0_cosThetaProdPlane_eta_US_hist->SetMinimum(0);
    L0_L0_cosThetaProdPlane_eta_US_hist->Draw("p e");

    //ME here just for plotting, used loser
    TH1D *L0_L0_cosThetaProdPlane_eta_ME_hist = L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist->ProjectionX( Form("proj_US_ME_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);



    L0_L0_cosThetaProdPlane_eta_ME_hist->SetMarkerSize(1.5);
    L0_L0_cosThetaProdPlane_eta_ME_hist->SetMarkerStyle(24);
    L0_L0_cosThetaProdPlane_eta_ME_hist->SetMarkerColor(1);
    L0_L0_cosThetaProdPlane_eta_ME_hist->SetLineColor(1);
    L0_L0_cosThetaProdPlane_eta_ME_hist->Sumw2();
    L0_L0_cosThetaProdPlane_eta_ME_hist->Scale(nLL/L0_L0_cosThetaProdPlane_eta_ME_hist->Integral()); //scale ME to US
    L0_L0_cosThetaProdPlane_eta_ME_hist->Scale(1./L0_L0_cosThetaProdPlane_eta_ME_hist->GetXaxis()->GetBinWidth(1)); //US is scaled by bin width, ME needs to be scaled as well
    L0_L0_cosThetaProdPlane_eta_ME_hist->SetMinimum(0);
    L0_L0_cosThetaProdPlane_eta_ME_hist->Draw("p e same");

    cout<<endl;
    //cout<<"L-L ME/SE ratio = "<<L0_L0_cosThetaProdPlane_eta_ME_hist->Integral()/nLL<<endl;
    cout<<"L-L ME/SE ratio = "<<L0_L0_cosThetaProdPlane_eta_US_hist->GetBinError(1)/L0_L0_cosThetaProdPlane_eta_ME_hist->GetBinError(1)<<endl;
    cout<<endl;

    TF1 *fitL0_L0_US_ThetaStar_no_corr_ME = new TF1("fitL0_L0_US_ThetaStar_no_corr_ME", "[0]*(1 + [1]*x)", -1, 1);
    fitL0_L0_US_ThetaStar_no_corr_ME->SetParameters(100, 0.5);

    L0_L0_cosThetaProdPlane_eta_ME_hist->Fit(fitL0_L0_US_ThetaStar_no_corr_ME, "s i 0 r");

    float P_L0_L0_no_corr_ME = fitL0_L0_US_ThetaStar_no_corr_ME->GetParameter(1)/(L0_alpha*L0_alpha);
    float P_L0_L0_no_corr_ME_err = fitL0_L0_US_ThetaStar_no_corr_ME->GetParError(1)/(L0_alpha*L0_alpha);


    TH1D *L0_L0_cosThetaProdPlane_eta_LS_hist = L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist->ProjectionX( Form("proj_LS_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);

    L0_L0_cosThetaProdPlane_eta_LS_hist->SetMarkerSize(1.5);
    L0_L0_cosThetaProdPlane_eta_LS_hist->SetMarkerStyle(21);
    L0_L0_cosThetaProdPlane_eta_LS_hist->SetMarkerColor(kBlue);
    double nLL_back = L0_L0_cosThetaProdPlane_eta_LS_hist->Integral();
    L0_L0_cosThetaProdPlane_eta_LS_hist->Sumw2();
    L0_L0_cosThetaProdPlane_eta_LS_hist->Scale(1./L0_L0_cosThetaProdPlane_eta_LS_hist->GetXaxis()->GetBinWidth(1));
    //L0_L0_cosThetaProdPlane_eta_LS_hist->Divide(L0_L0_cosThetaProdPlane_eta_eff);
    L0_L0_cosThetaProdPlane_eta_LS_hist->Draw("p e same");


    TH1D *L0_L0_cosThetaProdPlane_eta_ME_LS_hist = L0_L0_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist->ProjectionX( Form("proj_LS_ME_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);

    L0_L0_cosThetaProdPlane_eta_ME_LS_hist->SetMarkerSize(1.5);
    L0_L0_cosThetaProdPlane_eta_ME_LS_hist->SetMarkerStyle(25);
    L0_L0_cosThetaProdPlane_eta_ME_LS_hist->SetMarkerColor(kMagenta+1);
    L0_L0_cosThetaProdPlane_eta_ME_LS_hist->SetLineColor(kMagenta+1);
    L0_L0_cosThetaProdPlane_eta_ME_LS_hist->Sumw2();
    L0_L0_cosThetaProdPlane_eta_ME_LS_hist->Scale(nLL_back/L0_L0_cosThetaProdPlane_eta_ME_LS_hist->Integral()); //scale ME_LS to background
    L0_L0_cosThetaProdPlane_eta_ME_LS_hist->Scale(1./L0_L0_cosThetaProdPlane_eta_ME_LS_hist->GetXaxis()->GetBinWidth(1)); //background is scaled by bin width, ME_LS needs to be scaled as well
    L0_L0_cosThetaProdPlane_eta_ME_LS_hist->SetMinimum(0);
    L0_L0_cosThetaProdPlane_eta_ME_LS_hist->Draw("p e same");


    TF1 *fitL0_L0_US_ThetaStar_no_corr = new TF1("fitL0_L0_US_ThetaStar_no_corr", "[0]*(1 + [1]*x)", -1, 1);
    fitL0_L0_US_ThetaStar_no_corr->SetParameters(100, 0.5);

    //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
    L0_L0_cosThetaProdPlane_eta_US_hist->Fit(fitL0_L0_US_ThetaStar_no_corr, "s i 0 r");

    float P_L0_L0_no_corr = fitL0_L0_US_ThetaStar_no_corr->GetParameter(1)/(L0_alpha*L0_alpha);
    float P_L0_L0_no_corr_err = fitL0_L0_US_ThetaStar_no_corr->GetParError(1)/(L0_alpha*L0_alpha);

    fitL0_L0_US_ThetaStar_no_corr->SetLineColor(1);
    //fitL0_L0_US_ThetaStar_no_corr->Draw("same");

    TLegend *L0_L0_leg = new TLegend(0.15, 0.45, 0.45, 0.69);
    L0_L0_leg->AddEntry(L0_L0_cosThetaProdPlane_eta_US_hist, "(US-US) p#pi");
    L0_L0_leg->AddEntry(L0_L0_cosThetaProdPlane_eta_ME_hist, "(US-US) ME");
    L0_L0_leg->AddEntry(L0_L0_cosThetaProdPlane_eta_LS_hist, "Combinatorial bckg.");
    L0_L0_leg->AddEntry(L0_L0_cosThetaProdPlane_eta_ME_LS_hist, "Bckg. ME");
    //L0_L0_leg->AddEntry(fitL0_L0_US_ThetaStar_no_corr, "Linear fit to US");
    L0_L0_leg->SetBorderSize(0);
    L0_L0_leg->SetFillColorAlpha(0, 0.01);
    L0_L0_leg->Draw("same");

    TPaveText *L0_L0_text_no_corr = new TPaveText(0.5, 0.4, 0.85, 0.75, "NDC");
    L0_L0_text_no_corr->SetTextFont(42);
    //L0_L0_text_no_corr->AddText("STAR Internal");
    //L0_L0_text_no_corr->AddText("STAR preliminary");
    //((TText*)L0_L0_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
    L0_L0_text_no_corr->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
    L0_L0_text_no_corr->AddText("Minimum bias, no correction");
    L0_L0_text_no_corr->AddText("#Lambda^{0}-#Lambda^{0}");
    L0_L0_text_no_corr->AddText("|#it{y}| < 1");
    L0_L0_text_no_corr->AddText("0.5 < p_{T,1} < 5.0 GeV/#it{c}");
    L0_L0_text_no_corr->AddText("0.5 < p_{T,2} < 5.0 GeV/#it{c}");
    L0_L0_text_no_corr->AddText(Form("P = %.2f #pm %.2f", P_L0_L0_no_corr, fabs(P_L0_L0_no_corr_err)));
    L0_L0_text_no_corr->AddText(Form("P_{ME} = %.3f #pm %.3f", P_L0_L0_no_corr_ME, fabs(P_L0_L0_no_corr_ME_err)));
    L0_L0_text_no_corr->SetFillColorAlpha(0, 0.01);
    L0_L0_text_no_corr->Draw("same");

    if(delta_eta_bin == 1) Delta_eta_text_1->Draw("same");
    if(delta_eta_bin == 2) Delta_eta_text_2->Draw("same");

    L0_L0_cosThetaProdPlane_eta_no_corr_can->SaveAs(Form("./plots/Lambda/L_correlations_delta_eta_delta_phi/L0_L0_cosThetaProdPlane_eta_no_corr_less_delta_eta_delta_phi_%i.png", delta_eta_bin));

    //----------------------------------------------------------------------------------------------------

    TCanvas *L0_L0_cosThetaProdPlane_eta_can = new TCanvas(Form("L0_L0_cosThetaProdPlane_eta_can_%i", delta_eta_bin), Form("L0_L0_cosThetaProdPlane_eta_can_%i", delta_eta_bin), 1200, 1000);

    L0_L0_cosThetaProdPlane_eta_can->cd();

    //ME histogram higher

    L0_L0_cosThetaProdPlane_eta_US_hist->Divide(L0_L0_cosThetaProdPlane_eta_ME_hist); // correct US using ME
    L0_L0_cosThetaProdPlane_eta_US_hist->Scale(nLL/L0_L0_cosThetaProdPlane_eta_US_hist->Integral()); //scale back to raw US
    L0_L0_cosThetaProdPlane_eta_US_hist->Scale(1./L0_L0_cosThetaProdPlane_eta_US_hist->GetXaxis()->GetBinWidth(1)); //bin width scaling
    L0_L0_cosThetaProdPlane_eta_US_hist->SetMinimum(0);
    L0_L0_cosThetaProdPlane_eta_US_hist->Draw("p e");

    //L0_L0_cosThetaProdPlane_eta_ME_hist->Draw("same p e");

    L0_L0_cosThetaProdPlane_eta_LS_hist->Divide(L0_L0_cosThetaProdPlane_eta_ME_LS_hist); //correct background using ME
    L0_L0_cosThetaProdPlane_eta_LS_hist->Scale(nLL_back/L0_L0_cosThetaProdPlane_eta_LS_hist->Integral()); //scale back to raw background
    L0_L0_cosThetaProdPlane_eta_LS_hist->Scale(1./L0_L0_cosThetaProdPlane_eta_LS_hist->GetXaxis()->GetBinWidth(1)); //bin width scaling
    L0_L0_cosThetaProdPlane_eta_LS_hist->SetMinimum(0);
    L0_L0_cosThetaProdPlane_eta_LS_hist->Draw("p e same");

    //fit dN/dcos(theta*)
    //signal + bacground
    TF1 *fitL0_L0_US_ThetaStar = new TF1("fitL0_L0_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
    fitL0_L0_US_ThetaStar->SetParameters(100, 0.5);

    //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
    L0_L0_cosThetaProdPlane_eta_US_hist->Fit(fitL0_L0_US_ThetaStar, "s i 0 r");

    float chi_2_LL = L0_L0_cosThetaProdPlane_eta_US_hist->Chisquare(fitL0_L0_US_ThetaStar);

    cout<<endl;
    cout<<"Chi2(LL) = "<<chi_2_LL<<endl;
    cout<<"Chi2(LL)/NDF = "<<chi_2_LL/8.<<endl; //NDF should be 10 - 2 (Npoints - N free params)
    cout<<endl;

    float P_L0_L0 = fitL0_L0_US_ThetaStar->GetParameter(1)/(L0_alpha*L0_alpha);
    float P_L0_L0_err = fitL0_L0_US_ThetaStar->GetParError(1)/(L0_alpha*L0_alpha);

    fitL0_L0_US_ThetaStar->SetLineColor(1);
    fitL0_L0_US_ThetaStar->Draw("same");

    //background
    TF1 *fitL0_L0_US_LS_ThetaStar = new TF1("fitL0_L0_US_LS_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
    fitL0_L0_US_LS_ThetaStar->SetParameters(100, 0.5);

    //fit_res_gaUS_LS_wrong_sign = L_inv_mass_US_LS->Fit(fitGaUS_LSsBack, "s i 0", "", 1.07, 1.4);
    L0_L0_cosThetaProdPlane_eta_LS_hist->Fit(fitL0_L0_US_LS_ThetaStar, "s i 0 r");

    float chi_2_LL_back = L0_L0_cosThetaProdPlane_eta_LS_hist->Chisquare(fitL0_L0_US_LS_ThetaStar);

    cout<<endl;
    cout<<"Chi2(LL,back) = "<<chi_2_LL_back<<endl;
    cout<<"Chi2(LL,back)/NDF = "<<chi_2_LL_back/8.<<endl; //NDF should be 10 - 2 (Npoints - N free params)
    cout<<endl;

    float P_L0_L0_back = fitL0_L0_US_LS_ThetaStar->GetParameter(1)/(L0_alpha*L0_alpha);
    float P_L0_L0_back_err = fitL0_L0_US_LS_ThetaStar->GetParError(1)/(L0_alpha*L0_alpha);

    fitL0_L0_US_LS_ThetaStar->SetLineColor(1);
    fitL0_L0_US_LS_ThetaStar->Draw("same");


    TPaveText *L0_L0_text = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
    L0_L0_text->SetTextFont(42);
    //L0_L0_text->AddText("STAR Internal");
    //L0_L0_text->AddText("STAR preliminary");
    //((TText*)L0_L0_text->GetListOfLines()->Last())->SetTextColor(2);
    L0_L0_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
    L0_L0_text->AddText("Minimum bias");
    L0_L0_text->AddText("#Lambda^{0}-#Lambda^{0}");
    L0_L0_text->AddText("|#it{y}| < 1");
    L0_L0_text->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
    L0_L0_text->AddText(Form("P_{tot} = %.2f #pm %.2f", P_L0_L0, fabs(P_L0_L0_err)));
    L0_L0_text->AddText(Form("P_{bckg} = %.2f #pm %.2f", P_L0_L0_back, fabs(P_L0_L0_back_err)));
    L0_L0_text->SetFillColorAlpha(0, 0.01);
    L0_L0_text->Draw("same");

    L0_L0_leg->Draw("same");

    if(delta_eta_bin == 1) Delta_eta_text_1->Draw("same");
    if(delta_eta_bin == 2) Delta_eta_text_2->Draw("same");

    L0_L0_cosThetaProdPlane_eta_can->SaveAs(Form("./plots/Lambda/L_correlations_delta_eta_delta_phi/L0_L0_cosThetaProdPlane_eta_less_delta_eta_delta_phi_%i.png", delta_eta_bin));

    //----------------------------------------------------------------------

    TCanvas *L0_L0_cosThetaProdPlane_eta_can_2 = new TCanvas(Form("L0_L0_cosThetaProdPlane_eta_can_2_%i", delta_eta_bin), Form("L0_L0_cosThetaProdPlane_eta_can_2_%i", delta_eta_bin), 1200, 1000);

    L0_L0_cosThetaProdPlane_eta_can_2->cd();

    L0_L0_cosThetaProdPlane_eta_US_hist->Add(L0_L0_cosThetaProdPlane_eta_LS_hist, -1);
    L0_L0_cosThetaProdPlane_eta_US_hist->SetMinimum(0);
    L0_L0_cosThetaProdPlane_eta_US_hist->Write(Form("L0_L0_cosThetaProdPlane_delta_y_delta_phi_%i", delta_eta_bin-1));
    L0_L0_cosThetaProdPlane_eta_US_hist->Draw("p e");


    TF1 *fitL0_L0_US_ThetaStar_2 = new TF1("fitL0_L0_US_ThetaStar_2", "[0]*(1 + [1]*x)", -1, 1);
    fitL0_L0_US_ThetaStar_2->SetParameters(100, 0.5);

    //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
    L0_L0_cosThetaProdPlane_eta_US_hist->Fit(fitL0_L0_US_ThetaStar_2, "s i 0 r");


    //store fit result for systematic errors
    //tight topo cuts
    if(cut_type == 1)
    {
      SysErrCutsTopo->cd();

      fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParameters(fitL0_L0_US_ThetaStar->GetParameter(0), fitL0_L0_US_ThetaStar->GetParameter(1));
      fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParError(0, fitL0_L0_US_ThetaStar->GetParError(0));
      fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParError(1, fitL0_L0_US_ThetaStar->GetParError(1));
      fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->Write();


      fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParameters(fitL0_L0_US_LS_ThetaStar->GetParameter(0), fitL0_L0_US_LS_ThetaStar->GetParameter(1));
      fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParError(0, fitL0_L0_US_LS_ThetaStar->GetParError(0));
      fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParError(1, fitL0_L0_US_LS_ThetaStar->GetParError(1));
      fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->Write();

      fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParameters(fitL0_L0_US_ThetaStar_2->GetParameter(0), fitL0_L0_US_ThetaStar_2->GetParameter(1));
      fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParError(0, fitL0_L0_US_ThetaStar_2->GetParError(0));
      fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParError(1, fitL0_L0_US_ThetaStar_2->GetParError(1));
      fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->Write();
    }

    //tight pT cuts
    if(cut_type == 2)
    {
      SysErrCutsPt->cd();

      fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParameters(fitL0_L0_US_ThetaStar->GetParameter(0), fitL0_L0_US_ThetaStar->GetParameter(1));
      fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParError(0, fitL0_L0_US_ThetaStar->GetParError(0));
      fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParError(1, fitL0_L0_US_ThetaStar->GetParError(1));
      fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->Write();


      fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParameters(fitL0_L0_US_LS_ThetaStar->GetParameter(0), fitL0_L0_US_LS_ThetaStar->GetParameter(1));
      fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParError(0, fitL0_L0_US_LS_ThetaStar->GetParError(0));
      fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParError(1, fitL0_L0_US_LS_ThetaStar->GetParError(1));
      fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->Write();

      fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParameters(fitL0_L0_US_ThetaStar_2->GetParameter(0), fitL0_L0_US_ThetaStar_2->GetParameter(1));
      fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParError(0, fitL0_L0_US_ThetaStar_2->GetParError(0));
      fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParError(1, fitL0_L0_US_ThetaStar_2->GetParError(1));
      fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->Write();
    }

    //tight pT cuts
    if(cut_type == 3)
    {
      SysErrCutsDCA->cd();

      fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParameters(fitL0_L0_US_ThetaStar->GetParameter(0), fitL0_L0_US_ThetaStar->GetParameter(1));
      fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParError(0, fitL0_L0_US_ThetaStar->GetParError(0));
      fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParError(1, fitL0_L0_US_ThetaStar->GetParError(1));
      fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->Write();


      fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParameters(fitL0_L0_US_LS_ThetaStar->GetParameter(0), fitL0_L0_US_LS_ThetaStar->GetParameter(1));
      fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParError(0, fitL0_L0_US_LS_ThetaStar->GetParError(0));
      fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParError(1, fitL0_L0_US_LS_ThetaStar->GetParError(1));
      fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->Write();

      fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParameters(fitL0_L0_US_ThetaStar_2->GetParameter(0), fitL0_L0_US_ThetaStar_2->GetParameter(1));
      fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParError(0, fitL0_L0_US_ThetaStar_2->GetParError(0));
      fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParError(1, fitL0_L0_US_ThetaStar_2->GetParError(1));
      fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->Write();
    }

    float P_L0_L0_2 = fitL0_L0_US_ThetaStar_2->GetParameter(1)/(L0_alpha*L0_alpha);
    float P_L0_L0_err_2 = fitL0_L0_US_ThetaStar_2->GetParError(1)/(L0_alpha*L0_alpha);

    fitL0_L0_US_ThetaStar_2->SetLineColor(1);
    fitL0_L0_US_ThetaStar_2->Draw("same");

    //-------------------------------------------------------------------------------------------

    //calculate total systematic uncertainty for L-L

    //slope
    //float SysErrSlope_L0_L0 = SysErrSlope_delta_eta_delta_phi_hist[delta_eta_bin-1]->GetBinContent(2);
    float SysErrSlope_L0_L0 = fabs(ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin-1]->GetBinContent(2))/fabs(P_L0_L0_2);
    //float SysErrSlope_L0_L0 = 0;

    SysErr_sum_ME += fabs(ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin-1]->GetBinContent(2))/fabs(ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin-1]->GetBinError(2));
    SysErr_sum_of_w_ME += 1./fabs(ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin-1]->GetBinError(2));

    //----------------------------

    //background subtraction
    pair<float,float> P_L0_L0_from_fits_pair = Get_P_LL_from_fit(fitL0_L0_US_ThetaStar,fitL0_L0_US_LS_ThetaStar,1);

    float P_L0_L0_from_fits = P_L0_L0_from_fits_pair.first;
    float P_L0_L0_from_fits_err = P_L0_L0_from_fits_pair.second;

    //statistical error correction for systematic error
    float SysErrBackground_L0_L0_corr = sqrt( fabs( P_L0_L0_err_2*P_L0_L0_err_2 - P_L0_L0_from_fits_err*P_L0_L0_from_fits_err ) );

    float SysErrBackground_L0_L0_work = ( fabs( P_L0_L0_2 - P_L0_L0_from_fits) - SysErrBackground_L0_L0_corr )/fabs(P_L0_L0_from_fits);

    float SysErrBackground_L0_L0 = 0;

    if( SysErrBackground_L0_L0_work > 0 ) SysErrBackground_L0_L0 = SysErrBackground_L0_L0_work; //store sys. err. only if it is larger than statistical fluctuations

    float SysErrBackground_L0_L0_no_corr = fabs( P_L0_L0_2 - P_L0_L0_from_fits)/fabs(P_L0_L0_from_fits);


    //SysErr_sum_background += SysErrBackground_L0_L0_no_corr*fabs(P_L0_L0_2)/SysErrBackground_L0_L0_corr;
    //SysErr_sum_background += SysErrBackground_L0_L0*fabs(P_L0_L0_2)/SysErrBackground_L0_L0_corr;
    //SysErr_sum_of_w_background += 1./SysErrBackground_L0_L0_corr;

    SysErr_sum_background += SysErrBackground_L0_L0*fabs(P_L0_L0_from_fits)/fabs(P_L0_L0_from_fits_err);
    SysErr_sum_of_w_background += 1./fabs(P_L0_L0_from_fits_err);

    //SysErr_sum_background += SysErrBackground_L0_L0_no_corr*fabs(P_L0_L0_2)/sqrt(P_L0_L0_err_2*P_L0_L0_err_2 + P_L0_L0_from_fits_err*P_L0_L0_from_fits_err);
    //SysErr_sum_of_w_background += 1./sqrt(P_L0_L0_err_2*P_L0_L0_err_2 + P_L0_L0_from_fits_err*P_L0_L0_from_fits_err);

    //----------------------------

    //cuts variation
    float SysErr_tight_cuts_L0_L0_from_fits = 0;

    float SysErr_tight_topo_cuts_L0_L0_from_fits = 0;
    float SysErr_tight_pT_cuts_L0_L0_from_fits = 0;
    float SysErr_tight_DCA_cuts_L0_L0_from_fits = 0;


    float SysErr_tight_cuts_L0_L0_from_fits_no_corr = 0;

    float SysErr_tight_topo_cuts_L0_L0_from_fits_no_corr = 0;
    float SysErr_tight_pT_cuts_L0_L0_from_fits_no_corr = 0;
    float SysErr_tight_DCA_cuts_L0_L0_from_fits_no_corr = 0;

    //------------------------------------------------------

    float SysErr_tight_cuts_L0_L0 = 0;

    float SysErr_tight_topo_cuts_L0_L0 = 0;
    float SysErr_tight_pT_cuts_L0_L0 = 0;
    float SysErr_tight_DCA_cuts_L0_L0 = 0;


    float SysErr_tight_cuts_L0_L0_no_corr = 0;

    float SysErr_tight_topo_cuts_L0_L0_no_corr = 0;
    float SysErr_tight_pT_cuts_L0_L0_no_corr = 0;
    float SysErr_tight_DCA_cuts_L0_L0_no_corr = 0;

    if(cut_type == 0)
    {
      //independent fit of US-US and background
      pair<float,float> P_L0_L0_tight_topo_cuts_from_fits_signal_pair = Get_P_LL_from_fit(fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1], fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1], 1);

      float P_L0_L0_tight_topo_cuts_from_fits_signal = P_L0_L0_tight_topo_cuts_from_fits_signal_pair.first;
      float P_L0_L0_tight_topo_cuts_from_fits_signal_err = P_L0_L0_tight_topo_cuts_from_fits_signal_pair.second;

      //statistical error correction for systematic error
      float SysErr_tight_topo_cuts_L0_L0_from_fits_corr = sqrt( fabs( P_L0_L0_from_fits_err*P_L0_L0_from_fits_err - P_L0_L0_tight_topo_cuts_from_fits_signal_err*P_L0_L0_tight_topo_cuts_from_fits_signal_err ) );

      float SysErr_tight_topo_cuts_L0_L0_from_fits_work = ( fabs( P_L0_L0_from_fits - P_L0_L0_tight_topo_cuts_from_fits_signal) - SysErr_tight_topo_cuts_L0_L0_from_fits_corr )/fabs(P_L0_L0_from_fits);

      if( SysErr_tight_topo_cuts_L0_L0_from_fits_work > 0 ) SysErr_tight_topo_cuts_L0_L0_from_fits = SysErr_tight_topo_cuts_L0_L0_from_fits_work; //store sys. err. only if it is larger than statistical fluctuations

      SysErr_tight_topo_cuts_L0_L0_from_fits_no_corr = fabs( P_L0_L0_from_fits - P_L0_L0_tight_topo_cuts_from_fits_signal)/fabs(P_L0_L0_from_fits);


      //SysErr_sum_from_fits_topo_cuts += SysErr_tight_topo_cuts_L0_L0_from_fits*fabs(P_L0_L0_from_fits)/SysErr_tight_topo_cuts_L0_L0_from_fits_corr;
      //SysErr_sum_from_fits_of_w_topo_cuts += 1./SysErr_tight_topo_cuts_L0_L0_from_fits_corr;

      SysErr_sum_from_fits_topo_cuts += SysErr_tight_topo_cuts_L0_L0_from_fits*fabs(P_L0_L0_from_fits)/fabs(P_L0_L0_from_fits_err);
      SysErr_sum_from_fits_of_w_topo_cuts += 1./fabs(P_L0_L0_from_fits_err);

      //----------------------------------------------------------------

      //fit after background subtraction

      float P_L0_L0_tight_topo_cuts = fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->GetParameter(1)/(L0_alpha*L0_alpha);
      float P_L0_L0_tight_topo_cuts_err = fitL0_L0_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->GetParError(1)/(L0_alpha*L0_alpha);

      //statistical error correction for systematic error
      float SysErr_tight_topo_cuts_L0_L0_corr = sqrt( fabs( P_L0_L0_err_2*P_L0_L0_err_2 - P_L0_L0_tight_topo_cuts_err*P_L0_L0_tight_topo_cuts_err ) );

      float SysErr_tight_topo_cuts_L0_L0_work = ( fabs( P_L0_L0_2 - P_L0_L0_tight_topo_cuts) - SysErr_tight_topo_cuts_L0_L0_corr )/fabs(P_L0_L0_2);

      if( SysErr_tight_topo_cuts_L0_L0_work > 0 ) SysErr_tight_topo_cuts_L0_L0 = SysErr_tight_topo_cuts_L0_L0_work; //store sys. err. only if it is larger than statistical fluctuations

      SysErr_tight_topo_cuts_L0_L0_no_corr = fabs( P_L0_L0_2 - P_L0_L0_tight_topo_cuts)/fabs(P_L0_L0_2);


      SysErr_sum_topo_cuts += SysErr_tight_topo_cuts_L0_L0*fabs(P_L0_L0_2)/fabs(P_L0_L0_err_2);
      SysErr_sum_of_w_topo_cuts += 1./fabs(P_L0_L0_err_2);


      //___________________________________________________________________________________________________________________________________

      //independent fit of US-US and background
      pair<float,float> P_L0_L0_tight_pT_cuts_from_fits_signal_pair = Get_P_LL_from_fit(fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1], fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1], 1);

      float P_L0_L0_tight_pT_cuts_from_fits_signal = P_L0_L0_tight_pT_cuts_from_fits_signal_pair.first;
      float P_L0_L0_tight_pT_cuts_from_fits_signal_err = P_L0_L0_tight_pT_cuts_from_fits_signal_pair.second;

      //statistical error correction for systematic error
      float SysErr_tight_pT_cuts_L0_L0_from_fits_corr = sqrt( fabs( P_L0_L0_from_fits_err*P_L0_L0_from_fits_err - P_L0_L0_tight_pT_cuts_from_fits_signal_err*P_L0_L0_tight_pT_cuts_from_fits_signal_err ) );

      float SysErr_tight_pT_cuts_L0_L0_from_fits_work = ( fabs( P_L0_L0_from_fits - P_L0_L0_tight_pT_cuts_from_fits_signal) - SysErr_tight_pT_cuts_L0_L0_from_fits_corr )/fabs(P_L0_L0_from_fits);

      if( SysErr_tight_pT_cuts_L0_L0_from_fits_work > 0 ) SysErr_tight_pT_cuts_L0_L0_from_fits = SysErr_tight_pT_cuts_L0_L0_from_fits_work; //store sys. err. only if it is larger than statistical fluctuations

      SysErr_tight_pT_cuts_L0_L0_from_fits_no_corr = fabs( P_L0_L0_from_fits - P_L0_L0_tight_pT_cuts_from_fits_signal)/fabs(P_L0_L0_from_fits);


      SysErr_sum_from_fits_pT_cut += SysErr_tight_pT_cuts_L0_L0_from_fits*fabs(P_L0_L0_from_fits)/fabs(P_L0_L0_from_fits_err);
      SysErr_sum_from_fits_of_w_pT_cut += 1./fabs(P_L0_L0_from_fits_err);

      //----------------------------------------------------------------


      float P_L0_L0_tight_pT_cuts = fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->GetParameter(1)/(L0_alpha*L0_alpha);
      float P_L0_L0_tight_pT_cuts_err = fitL0_L0_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->GetParError(1)/(L0_alpha*L0_alpha);

      //statistical error correction for systematic error
      float SysErr_tight_pT_cuts_L0_L0_corr = sqrt( fabs( P_L0_L0_err_2*P_L0_L0_err_2 - P_L0_L0_tight_pT_cuts_err*P_L0_L0_tight_pT_cuts_err ) );

      float SysErr_tight_pT_cuts_L0_L0_work = ( fabs( P_L0_L0_2 - P_L0_L0_tight_pT_cuts) - SysErr_tight_pT_cuts_L0_L0_corr )/fabs(P_L0_L0_2);

      if( SysErr_tight_pT_cuts_L0_L0_work > 0 ) SysErr_tight_pT_cuts_L0_L0 = SysErr_tight_pT_cuts_L0_L0_work; //store sys. err. only if it is larger than statistical fluctuations

      SysErr_tight_pT_cuts_L0_L0_no_corr = fabs( P_L0_L0_2 - P_L0_L0_tight_pT_cuts)/fabs(P_L0_L0_2);


      SysErr_sum_pT_cut += SysErr_tight_pT_cuts_L0_L0*fabs(P_L0_L0_2)/fabs(P_L0_L0_err_2);
      SysErr_sum_of_w_pT_cut += 1./fabs(P_L0_L0_err_2);

      //___________________________________________________________________________________________________________________________________

      //independent fit of US-US and background
      pair<float,float> P_L0_L0_tight_DCA_cuts_from_fits_signal_pair = Get_P_LL_from_fit(fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1], fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1], 1);

      float P_L0_L0_tight_DCA_cuts_from_fits_signal = P_L0_L0_tight_DCA_cuts_from_fits_signal_pair.first;
      float P_L0_L0_tight_DCA_cuts_from_fits_signal_err = P_L0_L0_tight_DCA_cuts_from_fits_signal_pair.second;


      //statistical error correction for systematic error
      float SysErr_tight_DCA_cuts_L0_L0_from_fits_corr = sqrt( fabs( P_L0_L0_from_fits_err*P_L0_L0_from_fits_err - P_L0_L0_tight_DCA_cuts_from_fits_signal_err*P_L0_L0_tight_DCA_cuts_from_fits_signal_err ) );

      float SysErr_tight_DCA_cuts_L0_L0_from_fits_work = ( fabs( P_L0_L0_from_fits - P_L0_L0_tight_DCA_cuts_from_fits_signal) - SysErr_tight_DCA_cuts_L0_L0_from_fits_corr )/fabs(P_L0_L0_from_fits);

      if( SysErr_tight_DCA_cuts_L0_L0_from_fits_work > 0 ) SysErr_tight_DCA_cuts_L0_L0_from_fits = SysErr_tight_DCA_cuts_L0_L0_from_fits_work; //store sys. err. only if it is larger than statistical fluctuations

      SysErr_tight_DCA_cuts_L0_L0_from_fits_no_corr = fabs( P_L0_L0_from_fits - P_L0_L0_tight_DCA_cuts_from_fits_signal)/fabs(P_L0_L0_from_fits);


      //SysErr_sum_from_fits_DCA_cut += SysErr_tight_DCA_cuts_L0_L0_from_fits*fabs(P_L0_L0_from_fits)/SysErr_tight_DCA_cuts_L0_L0_from_fits_corr;
      //SysErr_sum_from_fits_of_w_DCA_cut += 1./SysErr_tight_DCA_cuts_L0_L0_from_fits_corr;

      SysErr_sum_from_fits_DCA_cut += SysErr_tight_DCA_cuts_L0_L0_from_fits*fabs(P_L0_L0_from_fits)/fabs(P_L0_L0_from_fits_err);
      SysErr_sum_from_fits_of_w_DCA_cut += 1./fabs(P_L0_L0_from_fits_err);

      //----------------------------------------------------------------


      float P_L0_L0_tight_DCA_cuts = fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->GetParameter(1)/(L0_alpha*L0_alpha);
      float P_L0_L0_tight_DCA_cuts_err = fitL0_L0_tight_DCA_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->GetParError(1)/(L0_alpha*L0_alpha);

      //statistical error correction for systematic error
      float SysErr_tight_DCA_cuts_L0_L0_corr = sqrt( fabs( P_L0_L0_err_2*P_L0_L0_err_2 - P_L0_L0_tight_DCA_cuts_err*P_L0_L0_tight_DCA_cuts_err ) );

      float SysErr_tight_DCA_cuts_L0_L0_work = ( fabs( P_L0_L0_2 - P_L0_L0_tight_DCA_cuts) - SysErr_tight_DCA_cuts_L0_L0_corr )/fabs(P_L0_L0_2);

      if( SysErr_tight_DCA_cuts_L0_L0_work > 0 ) SysErr_tight_DCA_cuts_L0_L0 = SysErr_tight_DCA_cuts_L0_L0_work; //store sys. err. only if it is larger than statistical fluctuations

      SysErr_tight_DCA_cuts_L0_L0_no_corr = fabs( P_L0_L0_2 - P_L0_L0_tight_DCA_cuts)/fabs(P_L0_L0_2);


      //SysErr_sum_DCA_cut += SysErr_tight_DCA_cuts_L0_L0*fabs(P_L0_L0_2)/SysErr_tight_DCA_cuts_L0_L0_corr;
      //SysErr_sum_of_w_DCA_cut += 1./SysErr_tight_DCA_cuts_L0_L0_corr;

      SysErr_sum_DCA_cut += SysErr_tight_DCA_cuts_L0_L0*fabs(P_L0_L0_2)/fabs(P_L0_L0_err_2);
      SysErr_sum_of_w_DCA_cut += 1./fabs(P_L0_L0_err_2);

      //___________________________________________________________________________________________________________________________________

      SysErr_tight_cuts_L0_L0 = sqrt(SysErr_tight_topo_cuts_L0_L0*SysErr_tight_topo_cuts_L0_L0 + SysErr_tight_pT_cuts_L0_L0*SysErr_tight_pT_cuts_L0_L0 + SysErr_tight_DCA_cuts_L0_L0*SysErr_tight_DCA_cuts_L0_L0);

      SysErr_tight_cuts_L0_L0_no_corr = sqrt(SysErr_tight_topo_cuts_L0_L0_no_corr*SysErr_tight_topo_cuts_L0_L0_no_corr + SysErr_tight_pT_cuts_L0_L0_no_corr*SysErr_tight_pT_cuts_L0_L0_no_corr +
                                                SysErr_tight_DCA_cuts_L0_L0_no_corr*SysErr_tight_DCA_cuts_L0_L0_no_corr);

      //----------------------------

      SysErr_tight_cuts_L0_L0_from_fits = sqrt(SysErr_tight_topo_cuts_L0_L0_from_fits*SysErr_tight_topo_cuts_L0_L0_from_fits + SysErr_tight_pT_cuts_L0_L0_from_fits*SysErr_tight_pT_cuts_L0_L0_from_fits +
                                                  SysErr_tight_DCA_cuts_L0_L0_from_fits*SysErr_tight_DCA_cuts_L0_L0_from_fits);

      SysErr_tight_cuts_L0_L0_from_fits_no_corr = sqrt(SysErr_tight_topo_cuts_L0_L0_from_fits_no_corr*SysErr_tight_topo_cuts_L0_L0_from_fits_no_corr + SysErr_tight_pT_cuts_L0_L0_from_fits_no_corr*SysErr_tight_pT_cuts_L0_L0_from_fits_no_corr +
                                                          SysErr_tight_DCA_cuts_L0_L0_from_fits_no_corr*SysErr_tight_DCA_cuts_L0_L0_from_fits_no_corr );


    }


    //total
    //float SysErrTot_L0_L0 = sqrt( sysErr_alpha_L0_L0*sysErr_alpha_L0_L0 + SysErrSlope_L0_L0*SysErrSlope_L0_L0 + SysErrBackground_L0_L0*SysErrBackground_L0_L0  );
    //float SysErrTot_L0_L0 = sqrt( sysErr_alpha_L0_L0*sysErr_alpha_L0_L0 + SysErrSlope_L0_L0*SysErrSlope_L0_L0 + SysErrBackground_L0_L0*SysErrBackground_L0_L0 + SysErr_tight_cuts_L0_L0*SysErr_tight_cuts_L0_L0  );
    float SysErrTot_L0_L0 = sqrt( sysErr_alpha_L0_L0*sysErr_alpha_L0_L0 + SysErrSlope_L0_L0*SysErrSlope_L0_L0 + SysErrBackground_L0_L0*SysErrBackground_L0_L0 + SysErr_tight_cuts_L0_L0_from_fits*SysErr_tight_cuts_L0_L0_from_fits  );

    //TPaveText *L0_L0_text_2 = new TPaveText(0.47, 0.15, 0.85, 0.53, "NDC"); //with polarization value
    TPaveText *L0_L0_text_2 = new TPaveText(0.47, 0.15, 0.8, 0.53, "NDC"); //without polarization value
    L0_L0_text_2->SetTextFont(42);
    //L0_L0_text_2->SetTextSize(15);
    //L0_L0_text_2->AddText("STAR");
    //L0_L0_text_2->AddText("STAR preliminary");
    //((TText*)L0_L0_text_2->GetListOfLines()->Last())->SetTextColor(2);
    L0_L0_text_2->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
    L0_L0_text_2->AddText("Minimum bias");
    L0_L0_text_2->AddText("#Lambda^{0}-#Lambda^{0}");
    L0_L0_text_2->AddText("|#it{y}| < 1");
    L0_L0_text_2->AddText("0.5 < p_{T,1} < 5.0 GeV/#it{c}");
    L0_L0_text_2->AddText("0.5 < p_{T,2} < 5.0 GeV/#it{c}");
    L0_L0_text_2->AddText(Form("P_{#Lambda_{1}#Lambda_{2}} = %.3f #pm %.3f", P_L0_L0_2, fabs(P_L0_L0_err_2)));
    L0_L0_text_2->AddText(Form("P_{#Lambda_{1}#Lambda_{2}} = %.3f #pm %.3f", P_L0_L0_from_fits, fabs(P_L0_L0_from_fits_err)));
    //L0_L0_text_2->AddText(Form("P_{topo} = %.2f", P_L0_L0_tight_topo_cuts));
    //L0_L0_text_2->AddText(Form("P_{pT} = %.2f", P_L0_L0_tight_pT_cuts));
    L0_L0_text_2->SetFillColorAlpha(0, 0.01);
    L0_L0_text_2->Draw("same");


    TLegend *L0_L0_2_leg = new TLegend(0.15, 0.3, 0.4, 0.49);
    //L0_L0_2_leg->SetTextSizePixels(15);
    L0_L0_2_leg->AddEntry(L0_L0_cosThetaProdPlane_eta_US_hist, "(US-US)-Bckg.");
    L0_L0_2_leg->AddEntry(fitL0_L0_US_ThetaStar_2, "Fit", "l");
    L0_L0_2_leg->SetBorderSize(0);
    L0_L0_2_leg->SetFillColorAlpha(0, 0.01);
    L0_L0_2_leg->Draw("same");

    if(delta_eta_bin == 1) Delta_eta_text_1->Draw("same");
    if(delta_eta_bin == 2) Delta_eta_text_2->Draw("same");

    L0_L0_cosThetaProdPlane_eta_can_2->SaveAs(Form("./plots/Lambda/L_correlations_delta_eta_delta_phi/L0_L0_cosThetaProdPlane_eta_subtract_less_delta_eta_delta_phi_%i.png", delta_eta_bin));

    //----------------------------------


    //PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->SetPoint(2, P_L0_L0_2, 2);
    //PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->SetPointError(2, fabs(P_L0_L0_err_2), 0);

    PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->SetPoint(2, P_L0_L0_from_fits, 2);
    PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->SetPointError(2, fabs(P_L0_L0_from_fits_err), 0);


    //PolarizationGraph_delta_eta_delta_phi_sys_err[delta_eta_bin-1]->SetPoint(2, P_L0_L0_2, 2);
    //PolarizationGraph_delta_eta_delta_phi_sys_err[delta_eta_bin-1]->SetPointError(2, fabs(SysErrTot_L0_L0*P_L0_L0_2), 0.045);

    PolarizationGraph_delta_eta_delta_phi_sys_err[delta_eta_bin-1]->SetPoint(2, P_L0_L0_from_fits, 2);
    PolarizationGraph_delta_eta_delta_phi_sys_err[delta_eta_bin-1]->SetPointError(2, fabs(SysErrTot_L0_L0*P_L0_L0_from_fits), 0.045);


    //____________________________________________________________________________________________________________________________________________________________________________________________________________

    TCanvas *L0bar_L0bar_cosThetaProdPlane_eta_no_corr_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_eta_no_corr_can_%i", delta_eta_bin), Form("L0bar_L0bar_cosThetaProdPlane_eta_no_corr_can_%i", delta_eta_bin), 1200, 1000);

    L0bar_L0bar_cosThetaProdPlane_eta_no_corr_can->cd();

    TH1D *L0bar_L0bar_cosThetaProdPlane_eta_US_hist = L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_hist->ProjectionX( Form("proj_US_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);

    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->GetXaxis()->SetTitle("cos(#theta*)");
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->GetXaxis()->CenterTitle();
    //L0bar_L0bar_cosThetaProdPlane_eta_US_hist->GetXaxis()->SetTextSizePixels(30);
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->GetYaxis()->SetTitle("d#it{N}/d cos(#theta*)");
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->GetYaxis()->CenterTitle();
    //L0bar_L0bar_cosThetaProdPlane_eta_US_hist->GetYaxis()->SetTextSizePixels(30);
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->GetYaxis()->SetMaxDigits(3);
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->SetMarkerSize(1.5);
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->SetMarkerStyle(20);
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->SetMarkerColor(kRed);
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->SetLineColor(kRed);
    double nLbarLbar = L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Integral();
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Sumw2();
    //L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Divide(L0bar_L0bar_cosThetaProdPlane_eta_eff);
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Scale(1./L0bar_L0bar_cosThetaProdPlane_eta_US_hist->GetXaxis()->GetBinWidth(1));
    //L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Scale(1./L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Integral());
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->SetMinimum(0);
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Draw("p e");

    //ME here just for plotting, used loser
    TH1D *L0bar_L0bar_cosThetaProdPlane_eta_ME_hist = L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_ME_hist->ProjectionX( Form("proj_US_ME_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);



    L0bar_L0bar_cosThetaProdPlane_eta_ME_hist->SetMarkerSize(1.5);
    L0bar_L0bar_cosThetaProdPlane_eta_ME_hist->SetMarkerStyle(24);
    L0bar_L0bar_cosThetaProdPlane_eta_ME_hist->SetMarkerColor(1);
    L0bar_L0bar_cosThetaProdPlane_eta_ME_hist->SetLineColor(1);
    L0bar_L0bar_cosThetaProdPlane_eta_ME_hist->Sumw2();
    L0bar_L0bar_cosThetaProdPlane_eta_ME_hist->Scale(nLbarLbar/L0bar_L0bar_cosThetaProdPlane_eta_ME_hist->Integral()); //scale ME to US
    L0bar_L0bar_cosThetaProdPlane_eta_ME_hist->Scale(1./L0bar_L0bar_cosThetaProdPlane_eta_ME_hist->GetXaxis()->GetBinWidth(1)); //US is scaled by bin width, ME needs to be scaled as well
    L0bar_L0bar_cosThetaProdPlane_eta_ME_hist->SetMinimum(0);
    L0bar_L0bar_cosThetaProdPlane_eta_ME_hist->Draw("p e same");

    cout<<endl;
    //cout<<"Lbar-Lbar ME/SE ratio = "<<L0bar_L0bar_cosThetaProdPlane_eta_ME_hist->Integral()/nLbarLbar<<endl;
    cout<<"Lbar-Lbar ME/SE ratio = "<<L0bar_L0bar_cosThetaProdPlane_eta_US_hist->GetBinError(1)/L0bar_L0bar_cosThetaProdPlane_eta_ME_hist->GetBinError(1)<<endl;
    cout<<endl;

    TF1 *fitL0bar_L0bar_US_ThetaStar_no_corr_ME = new TF1("fitL0bar_L0bar_US_ThetaStar_no_corr_ME", "[0]*(1 + [1]*x)", -1, 1);
    fitL0bar_L0bar_US_ThetaStar_no_corr_ME->SetParameters(100, 0.5);

    L0bar_L0bar_cosThetaProdPlane_eta_ME_hist->Fit(fitL0bar_L0bar_US_ThetaStar_no_corr_ME, "s i 0 r");

    float P_L0bar_L0bar_no_corr_ME = fitL0bar_L0bar_US_ThetaStar_no_corr_ME->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
    float P_L0bar_L0bar_no_corr_ME_err = fitL0bar_L0bar_US_ThetaStar_no_corr_ME->GetParError(1)/(L0bar_alpha*L0bar_alpha);


    TH1D *L0bar_L0bar_cosThetaProdPlane_eta_LS_hist = L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_hist->ProjectionX( Form("proj_LS_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);

    L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->SetMarkerSize(1.5);
    L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->SetMarkerStyle(21);
    L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->SetMarkerColor(kBlue);
    double nLbarLbar_back = L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->Integral();
    L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->Sumw2();
    L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->Scale(1./L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->GetXaxis()->GetBinWidth(1));
    //L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->Divide(L0bar_L0bar_cosThetaProdPlane_eta_eff);
    L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->Draw("p e same");


    TH1D *L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist = L0bar_L0bar_cos_theta_star_vs_delta_eta_delta_phi_US_LS_ME_hist->ProjectionX( Form("proj_LS_ME_%i", delta_eta_bin), delta_eta_bin, delta_eta_bin);

    L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist->SetMarkerSize(1.5);
    L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist->SetMarkerStyle(25);
    L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist->SetMarkerColor(kMagenta+1);
    L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist->SetLineColor(kMagenta+1);
    L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist->Sumw2();
    L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist->Scale(nLbarLbar_back/L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist->Integral()); //scale ME_LS to background
    L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist->Scale(1./L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist->GetXaxis()->GetBinWidth(1)); //background is scaled by bin width, ME_LS needs to be scaled as well
    L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist->SetMinimum(0);
    L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist->Draw("p e same");


    TF1 *fitL0bar_L0bar_US_ThetaStar_no_corr = new TF1("fitL0bar_L0bar_US_ThetaStar_no_corr", "[0]*(1 + [1]*x)", -1, 1);
    fitL0bar_L0bar_US_ThetaStar_no_corr->SetParameters(100, 0.5);

    //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Fit(fitL0bar_L0bar_US_ThetaStar_no_corr, "s i 0 r");

    float P_L0bar_L0bar_no_corr = fitL0bar_L0bar_US_ThetaStar_no_corr->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
    float P_L0bar_L0bar_no_corr_err = fitL0bar_L0bar_US_ThetaStar_no_corr->GetParError(1)/(L0bar_alpha*L0bar_alpha);

    fitL0bar_L0bar_US_ThetaStar_no_corr->SetLineColor(1);
    //fitL0bar_L0bar_US_ThetaStar_no_corr->Draw("same");

    TLegend *L0bar_L0bar_leg = new TLegend(0.15, 0.45, 0.45, 0.69);
    L0bar_L0bar_leg->AddEntry(L0bar_L0bar_cosThetaProdPlane_eta_US_hist, "(US-US) p#pi");
    L0bar_L0bar_leg->AddEntry(L0bar_L0bar_cosThetaProdPlane_eta_ME_hist, "(US-US) ME");
    L0bar_L0bar_leg->AddEntry(L0bar_L0bar_cosThetaProdPlane_eta_LS_hist, "Combinatorial bckg.");
    L0bar_L0bar_leg->AddEntry(L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist, "Bckg. ME");
    //L0bar_L0bar_leg->AddEntry(fitL0bar_L0bar_US_ThetaStar_no_corr, "Linear fit to US");
    L0bar_L0bar_leg->SetBorderSize(0);
    L0bar_L0bar_leg->SetFillColorAlpha(0, 0.01);
    L0bar_L0bar_leg->Draw("same");

    TPaveText *L0bar_L0bar_text_no_corr = new TPaveText(0.5, 0.4, 0.85, 0.75, "NDC");
    L0bar_L0bar_text_no_corr->SetTextFont(42);
    //L0bar_L0bar_text_no_corr->AddText("STAR Internal");
    //L0bar_L0bar_text_no_corr->AddText("STAR preliminary");
    //((TText*)L0bar_L0bar_text_no_corr->GetListOfLines()->Last())->SetTextColor(2);
    L0bar_L0bar_text_no_corr->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
    L0bar_L0bar_text_no_corr->AddText("Minimum bias, no correction");
    L0bar_L0bar_text_no_corr->AddText("#bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
    L0bar_L0bar_text_no_corr->AddText("|#it{y}| < 1");
    L0bar_L0bar_text_no_corr->AddText("0.5 < p_{T,1} < 5.0 GeV/#it{c}");
    L0bar_L0bar_text_no_corr->AddText("0.5 < p_{T,2} < 5.0 GeV/#it{c}");
    L0bar_L0bar_text_no_corr->AddText(Form("P = %.2f #pm %.2f", P_L0bar_L0bar_no_corr, fabs(P_L0bar_L0bar_no_corr_err)));
    L0bar_L0bar_text_no_corr->AddText(Form("P_{ME} = %.3f #pm %.3f", P_L0bar_L0bar_no_corr_ME, fabs(P_L0bar_L0bar_no_corr_ME_err)));
    L0bar_L0bar_text_no_corr->SetFillColorAlpha(0, 0.01);
    L0bar_L0bar_text_no_corr->Draw("same");

    if(delta_eta_bin == 1) Delta_eta_text_1->Draw("same");
    if(delta_eta_bin == 2) Delta_eta_text_2->Draw("same");

    L0bar_L0bar_cosThetaProdPlane_eta_no_corr_can->SaveAs(Form("./plots/Lambda/L_correlations_delta_eta_delta_phi/L0bar_L0bar_cosThetaProdPlane_eta_no_corr_less_delta_eta_delta_phi_%i.png", delta_eta_bin));

    //----------------------------------------------------------------------------------------------------

    TCanvas *L0bar_L0bar_cosThetaProdPlane_eta_can = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_eta_can_%i", delta_eta_bin), Form("L0bar_L0bar_cosThetaProdPlane_eta_can_%i", delta_eta_bin), 1200, 1000);

    L0bar_L0bar_cosThetaProdPlane_eta_can->cd();

    //ME histogram higher

    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Divide(L0bar_L0bar_cosThetaProdPlane_eta_ME_hist); // correct US using ME
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Scale(nLbarLbar/L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Integral()); //scale back to raw US
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Scale(1./L0bar_L0bar_cosThetaProdPlane_eta_US_hist->GetXaxis()->GetBinWidth(1)); //bin width scaling
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->SetMinimum(0);
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Draw("p e");

    //L0bar_L0bar_cosThetaProdPlane_eta_ME_hist->Draw("same p e");

    L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->Divide(L0bar_L0bar_cosThetaProdPlane_eta_ME_LS_hist); //correct background using ME
    L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->Scale(nLbarLbar_back/L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->Integral()); //scale back to raw background
    L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->Scale(1./L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->GetXaxis()->GetBinWidth(1)); //bin width scaling
    L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->SetMinimum(0);
    L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->Draw("p e same");

    //fit dN/dcos(theta*)
    //signal + bacground
    TF1 *fitL0bar_L0bar_US_ThetaStar = new TF1("fitL0bar_L0bar_US_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
    fitL0bar_L0bar_US_ThetaStar->SetParameters(100, 0.5);

    //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Fit(fitL0bar_L0bar_US_ThetaStar, "s i 0 r");

    float chi_2_LbarLbar = L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Chisquare(fitL0bar_L0bar_US_ThetaStar);

    cout<<endl;
    cout<<"Chi2(LbarLbar) = "<<chi_2_LbarLbar<<endl;
    cout<<"Chi2(LbarLbar)/NDF = "<<chi_2_LbarLbar/8.<<endl; //NDF should be 10 - 2 (Npoints - N free params)
    cout<<endl;

    float P_L0bar_L0bar = fitL0bar_L0bar_US_ThetaStar->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
    float P_L0bar_L0bar_err = fitL0bar_L0bar_US_ThetaStar->GetParError(1)/(L0bar_alpha*L0bar_alpha);

    fitL0bar_L0bar_US_ThetaStar->SetLineColor(1);
    fitL0bar_L0bar_US_ThetaStar->Draw("same");

    //background
    TF1 *fitL0bar_L0bar_US_LS_ThetaStar = new TF1("fitL0bar_L0bar_US_LS_ThetaStar", "[0]*(1 + [1]*x)", -1, 1);
    fitL0bar_L0bar_US_LS_ThetaStar->SetParameters(100, 0.5);

    //fit_res_gaUS_LS_wrong_sign = L_inv_mass_US_LS->Fit(fitGaUS_LSsBack, "s i 0", "", 1.07, 1.4);
    L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->Fit(fitL0bar_L0bar_US_LS_ThetaStar, "s i 0 r");

    float chi_2_LbarLbar_back = L0bar_L0bar_cosThetaProdPlane_eta_LS_hist->Chisquare(fitL0bar_L0bar_US_LS_ThetaStar);

    cout<<endl;
    cout<<"Chi2(LbarLbar,back) = "<<chi_2_LbarLbar_back<<endl;
    cout<<"Chi2(LbarLbar,back)/NDF = "<<chi_2_LbarLbar_back/8.<<endl; //NDF should be 10 - 2 (Npoints - N free params)
    cout<<endl;

    float P_L0bar_L0bar_back = fitL0bar_L0bar_US_LS_ThetaStar->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
    float P_L0bar_L0bar_back_err = fitL0bar_L0bar_US_LS_ThetaStar->GetParError(1)/(L0bar_alpha*L0bar_alpha);

    fitL0bar_L0bar_US_LS_ThetaStar->SetLineColor(1);
    fitL0bar_L0bar_US_LS_ThetaStar->Draw("same");


    TPaveText *L0bar_L0bar_text = new TPaveText(0.6, 0.4, 0.85, 0.69, "NDC");
    L0bar_L0bar_text->SetTextFont(42);
    //L0bar_L0bar_text->AddText("STAR Internal");
    //L0bar_L0bar_text->AddText("STAR preliminary");
    //((TText*)L0bar_L0bar_text->GetListOfLines()->Last())->SetTextColor(2);
    L0bar_L0bar_text->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
    L0bar_L0bar_text->AddText("Minimum bias");
    L0bar_L0bar_text->AddText("#bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
    L0bar_L0bar_text->AddText("|#it{y}| < 1");
    L0bar_L0bar_text->AddText("0.5 < p_{T} < 5.0 GeV/#it{c}");
    L0bar_L0bar_text->AddText(Form("P_{tot} = %.2f #pm %.2f", P_L0bar_L0bar, fabs(P_L0bar_L0bar_err)));
    L0bar_L0bar_text->AddText(Form("P_{bckg} = %.2f #pm %.2f", P_L0bar_L0bar_back, fabs(P_L0bar_L0bar_back_err)));
    L0bar_L0bar_text->SetFillColorAlpha(0, 0.01);
    L0bar_L0bar_text->Draw("same");

    L0bar_L0bar_leg->Draw("same");

    if(delta_eta_bin == 1) Delta_eta_text_1->Draw("same");
    if(delta_eta_bin == 2) Delta_eta_text_2->Draw("same");

    L0bar_L0bar_cosThetaProdPlane_eta_can->SaveAs(Form("./plots/Lambda/L_correlations_delta_eta_delta_phi/L0bar_L0bar_cosThetaProdPlane_eta_less_delta_eta_delta_phi_%i.png", delta_eta_bin));

    //----------------------------------------------------------------------

    TCanvas *L0bar_L0bar_cosThetaProdPlane_eta_can_2 = new TCanvas(Form("L0bar_L0bar_cosThetaProdPlane_eta_can_2_%i", delta_eta_bin), Form("L0bar_L0bar_cosThetaProdPlane_eta_can_2_%i", delta_eta_bin), 1200, 1000);

    L0bar_L0bar_cosThetaProdPlane_eta_can_2->cd();

    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Add(L0bar_L0bar_cosThetaProdPlane_eta_LS_hist, -1);
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->SetMinimum(0);
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Write(Form("L0bar_L0bar_cosThetaProdPlane_delta_y_delta_phi_%i", delta_eta_bin-1));
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Draw("p e");


    TF1 *fitL0bar_L0bar_US_ThetaStar_2 = new TF1("fitL0bar_L0bar_US_ThetaStar_2", "[0]*(1 + [1]*x)", -1, 1);
    fitL0bar_L0bar_US_ThetaStar_2->SetParameters(100, 0.5);

    //fit_res_gaus_wrong_sign = L_inv_mass_US->Fit(fitGaussBack, "s i 0", "", 1.07, 1.4);
    L0bar_L0bar_cosThetaProdPlane_eta_US_hist->Fit(fitL0bar_L0bar_US_ThetaStar_2, "s i 0 r");

    //------------------------------------------------------------

    //store fit result for systematic errors
    //tight topo cuts
    if(cut_type == 1)
    {
      SysErrCutsTopo->cd();

      fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParameters(fitL0bar_L0bar_US_ThetaStar->GetParameter(0), fitL0bar_L0bar_US_ThetaStar->GetParameter(1));
      fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParError(0, fitL0bar_L0bar_US_ThetaStar->GetParError(0));
      fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParError(1, fitL0bar_L0bar_US_ThetaStar->GetParError(1));
      fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->Write();


      fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParameters(fitL0bar_L0bar_US_LS_ThetaStar->GetParameter(0), fitL0bar_L0bar_US_LS_ThetaStar->GetParameter(1));
      fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParError(0, fitL0bar_L0bar_US_LS_ThetaStar->GetParError(0));
      fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParError(1, fitL0bar_L0bar_US_LS_ThetaStar->GetParError(1));
      fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->Write();


      fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParameters(fitL0bar_L0bar_US_ThetaStar_2->GetParameter(0), fitL0bar_L0bar_US_ThetaStar_2->GetParameter(1));
      fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParError(0, fitL0bar_L0bar_US_ThetaStar_2->GetParError(0));
      fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParError(1, fitL0bar_L0bar_US_ThetaStar_2->GetParError(1));
      fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->Write();
    }

    //tight pT cuts
    if(cut_type == 2)
    {
      SysErrCutsPt->cd();

      fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParameters(fitL0bar_L0bar_US_ThetaStar->GetParameter(0), fitL0bar_L0bar_US_ThetaStar->GetParameter(1));
      fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParError(0, fitL0bar_L0bar_US_ThetaStar->GetParError(0));
      fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParError(1, fitL0bar_L0bar_US_ThetaStar->GetParError(1));
      fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->Write();


      fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParameters(fitL0bar_L0bar_US_LS_ThetaStar->GetParameter(0), fitL0bar_L0bar_US_LS_ThetaStar->GetParameter(1));
      fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParError(0, fitL0bar_L0bar_US_LS_ThetaStar->GetParError(0));
      fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParError(1, fitL0bar_L0bar_US_LS_ThetaStar->GetParError(1));
      fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->Write();


      fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParameters(fitL0bar_L0bar_US_ThetaStar_2->GetParameter(0), fitL0bar_L0bar_US_ThetaStar_2->GetParameter(1));
      fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParError(0, fitL0bar_L0bar_US_ThetaStar_2->GetParError(0));
      fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParError(1, fitL0bar_L0bar_US_ThetaStar_2->GetParError(1));
      fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->Write();
    }

    //tight pT cuts
    if(cut_type == 3)
    {
      SysErrCutsDCA->cd();

      fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParameters(fitL0bar_L0bar_US_ThetaStar->GetParameter(0), fitL0bar_L0bar_US_ThetaStar->GetParameter(1));
      fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParError(0, fitL0bar_L0bar_US_ThetaStar->GetParError(0));
      fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->SetParError(1, fitL0bar_L0bar_US_ThetaStar->GetParError(1));
      fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1]->Write();


      fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParameters(fitL0bar_L0bar_US_LS_ThetaStar->GetParameter(0), fitL0bar_L0bar_US_LS_ThetaStar->GetParameter(1));
      fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParError(0, fitL0bar_L0bar_US_LS_ThetaStar->GetParError(0));
      fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->SetParError(1, fitL0bar_L0bar_US_LS_ThetaStar->GetParError(1));
      fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1]->Write();

      fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParameters(fitL0bar_L0bar_US_ThetaStar_2->GetParameter(0), fitL0bar_L0bar_US_ThetaStar_2->GetParameter(1));
      fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParError(0, fitL0bar_L0bar_US_ThetaStar_2->GetParError(0));
      fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->SetParError(1, fitL0bar_L0bar_US_ThetaStar_2->GetParError(1));
      fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->Write();
    }

    float P_L0bar_L0bar_2 = fitL0bar_L0bar_US_ThetaStar_2->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
    float P_L0bar_L0bar_err_2 = fitL0bar_L0bar_US_ThetaStar_2->GetParError(1)/(L0bar_alpha*L0bar_alpha);

    fitL0bar_L0bar_US_ThetaStar_2->SetLineColor(1);
    fitL0bar_L0bar_US_ThetaStar_2->Draw("same");

    //--------------------------------------------------------------------

    //calculate total systematic uncertainty for L-Lbar

    //slope
    //float SysErrSlope_L0bar_L0bar = SysErrSlope_delta_eta_delta_phi_hist[delta_eta_bin-1]->GetBinContent(3);
    float SysErrSlope_L0bar_L0bar = fabs(ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin-1]->GetBinContent(3))/fabs(P_L0bar_L0bar_2);

    SysErr_sum_ME += fabs(ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin-1]->GetBinContent(3))/fabs(ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin-1]->GetBinError(3));
    SysErr_sum_of_w_ME += 1./fabs(ResidualPolarization_delta_eta_delta_phi_hist[delta_eta_bin-1]->GetBinError(3));

    //----------------------------

    //background subtraction
    pair<float,float> P_L0bar_L0bar_from_fits_pair = Get_P_LL_from_fit(fitL0bar_L0bar_US_ThetaStar,fitL0bar_L0bar_US_LS_ThetaStar,2);

    float P_L0bar_L0bar_from_fits = P_L0bar_L0bar_from_fits_pair.first;
    float P_L0bar_L0bar_from_fits_err = P_L0bar_L0bar_from_fits_pair.second;

    //statistical error correction for systematic error
    float SysErrBackground_L0bar_L0bar_corr = sqrt( fabs( P_L0bar_L0bar_err_2*P_L0bar_L0bar_err_2 - P_L0bar_L0bar_from_fits_err*P_L0bar_L0bar_from_fits_err ) );

    float SysErrBackground_L0bar_L0bar_work = ( fabs( P_L0bar_L0bar_2 - P_L0bar_L0bar_from_fits) - SysErrBackground_L0bar_L0bar_corr )/fabs(P_L0bar_L0bar_from_fits_err);

    float SysErrBackground_L0bar_L0bar = 0;

    if( SysErrBackground_L0bar_L0bar_work > 0 ) SysErrBackground_L0bar_L0bar = SysErrBackground_L0bar_L0bar_work; //store sys. err. only if it is larger than statistical fluctuations

    float SysErrBackground_L0bar_L0bar_no_corr = fabs( P_L0bar_L0bar_2 - P_L0bar_L0bar_from_fits)/fabs(P_L0bar_L0bar_from_fits); //no stat. err. corr.


    //SysErr_sum_background += SysErrBackground_L0bar_L0bar_no_corr*fabs(P_L0bar_L0bar_2)/SysErrBackground_L0bar_L0bar_corr;
    //SysErr_sum_background += SysErrBackground_L0bar_L0bar*fabs(P_L0bar_L0bar_2)/SysErrBackground_L0bar_L0bar_corr;
    //SysErr_sum_of_w_background += 1./SysErrBackground_L0bar_L0bar_corr;

    SysErr_sum_background += SysErrBackground_L0bar_L0bar*fabs(P_L0bar_L0bar_from_fits)/fabs(P_L0bar_L0bar_from_fits_err);
    SysErr_sum_of_w_background += 1./fabs(P_L0bar_L0bar_from_fits_err);

    //SysErr_sum_background += SysErrBackground_L0bar_L0bar_no_corr*fabs(P_L0bar_L0bar_2)/sqrt(P_L0bar_L0bar_err_2*P_L0bar_L0bar_err_2 + P_L0bar_L0bar_from_fits_err*P_L0bar_L0bar_from_fits_err);
    //SysErr_sum_of_w_background += 1./sqrt(P_L0bar_L0bar_err_2*P_L0bar_L0bar_err_2 + P_L0bar_L0bar_from_fits_err*P_L0bar_L0bar_from_fits_err);

    //-------------------------------------------------------------

    //cuts variation
    float SysErr_tight_cuts_L0bar_L0bar_from_fits = 0;

    float SysErr_tight_topo_cuts_L0bar_L0bar_from_fits = 0;
    float SysErr_tight_pT_cuts_L0bar_L0bar_from_fits = 0;
    float SysErr_tight_DCA_cuts_L0bar_L0bar_from_fits = 0;


    float SysErr_tight_cuts_L0bar_L0bar_from_fits_no_corr = 0;

    float SysErr_tight_topo_cuts_L0bar_L0bar_from_fits_no_corr = 0;
    float SysErr_tight_pT_cuts_L0bar_L0bar_from_fits_no_corr = 0;
    float SysErr_tight_DCA_cuts_L0bar_L0bar_from_fits_no_corr = 0;

    //------------------------------------------------------

    float SysErr_tight_cuts_L0bar_L0bar = 0;

    float SysErr_tight_topo_cuts_L0bar_L0bar = 0;
    float SysErr_tight_pT_cuts_L0bar_L0bar = 0;
    float SysErr_tight_DCA_cuts_L0bar_L0bar = 0;


    float SysErr_tight_cuts_L0bar_L0bar_no_corr = 0;

    float SysErr_tight_topo_cuts_L0bar_L0bar_no_corr = 0;
    float SysErr_tight_pT_cuts_L0bar_L0bar_no_corr = 0;
    float SysErr_tight_DCA_cuts_L0bar_L0bar_no_corr = 0;

    if(cut_type == 0)
    {
      //independent fit of US-US and background
      pair<float,float> P_L0bar_L0bar_tight_topo_cuts_from_fits_signal_pair = Get_P_LL_from_fit(fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1], fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1], 2);

      float P_L0bar_L0bar_tight_topo_cuts_from_fits_signal = P_L0bar_L0bar_tight_topo_cuts_from_fits_signal_pair.first;
      float P_L0bar_L0bar_tight_topo_cuts_from_fits_signal_err = P_L0bar_L0bar_tight_topo_cuts_from_fits_signal_pair.second;

      //statistical error correction for systematic error
      float SysErr_tight_topo_cuts_L0bar_L0bar_from_fits_corr = sqrt( fabs( P_L0bar_L0bar_from_fits_err*P_L0bar_L0bar_from_fits_err - P_L0bar_L0bar_tight_topo_cuts_from_fits_signal_err*P_L0bar_L0bar_tight_topo_cuts_from_fits_signal_err ) );

      float SysErr_tight_topo_cuts_L0bar_L0bar_from_fits_work = ( fabs( P_L0bar_L0bar_from_fits - P_L0bar_L0bar_tight_topo_cuts_from_fits_signal) - SysErr_tight_topo_cuts_L0bar_L0bar_from_fits_corr )/fabs(P_L0bar_L0bar_from_fits);

      if( SysErr_tight_topo_cuts_L0bar_L0bar_from_fits_work > 0 ) SysErr_tight_topo_cuts_L0bar_L0bar_from_fits = SysErr_tight_topo_cuts_L0bar_L0bar_from_fits_work; //store sys. err. only if it is larger than statistical fluctuations

      SysErr_tight_topo_cuts_L0bar_L0bar_from_fits_no_corr = fabs( P_L0bar_L0bar_from_fits - P_L0bar_L0bar_tight_topo_cuts_from_fits_signal)/fabs(P_L0bar_L0bar_from_fits);


      SysErr_sum_from_fits_topo_cuts += SysErr_tight_topo_cuts_L0bar_L0bar_from_fits*fabs(P_L0bar_L0bar_from_fits)/fabs(P_L0bar_L0bar_from_fits_err);
      SysErr_sum_from_fits_of_w_topo_cuts += 1./fabs(P_L0bar_L0bar_from_fits_err);

      //----------------------------------------------------------------

      //fit after background subtraction

      float P_L0bar_L0bar_tight_topo_cuts = fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
      float P_L0bar_L0bar_tight_topo_cuts_err = fitL0bar_L0bar_tight_topo_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->GetParError(1)/(L0bar_alpha*L0bar_alpha);

      //statistical error correction for systematic error
      float SysErr_tight_topo_cuts_L0bar_L0bar_corr = sqrt( fabs( P_L0bar_L0bar_err_2*P_L0bar_L0bar_err_2 - P_L0bar_L0bar_tight_topo_cuts_err*P_L0bar_L0bar_tight_topo_cuts_err ) );

      float SysErr_tight_topo_cuts_L0bar_L0bar_work = ( fabs( P_L0bar_L0bar_2 - P_L0bar_L0bar_tight_topo_cuts) - SysErr_tight_topo_cuts_L0bar_L0bar_corr )/fabs(P_L0bar_L0bar_2);

      if( SysErr_tight_topo_cuts_L0bar_L0bar_work > 0 ) SysErr_tight_topo_cuts_L0bar_L0bar = SysErr_tight_topo_cuts_L0bar_L0bar_work; //store sys. err. only if it is larger than statistical fluctuations

      SysErr_tight_topo_cuts_L0bar_L0bar_no_corr = fabs( P_L0bar_L0bar_2 - P_L0bar_L0bar_tight_topo_cuts)/fabs(P_L0bar_L0bar_2);


      SysErr_sum_topo_cuts += SysErr_tight_topo_cuts_L0bar_L0bar*fabs(P_L0bar_L0bar_2)/fabs(P_L0bar_L0bar_err_2);
      SysErr_sum_of_w_topo_cuts += 1./fabs(P_L0bar_L0bar_err_2);


      //___________________________________________________________________________________________________________________________________

      //independent fit of US-US and background
      pair<float,float> P_L0bar_L0bar_tight_pT_cuts_from_fits_signal_pair = Get_P_LL_from_fit(fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1], fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1], 2);

      float P_L0bar_L0bar_tight_pT_cuts_from_fits_signal = P_L0bar_L0bar_tight_pT_cuts_from_fits_signal_pair.first;
      float P_L0bar_L0bar_tight_pT_cuts_from_fits_signal_err = P_L0bar_L0bar_tight_pT_cuts_from_fits_signal_pair.second;

      //statistical error correction for systematic error
      float SysErr_tight_pT_cuts_L0bar_L0bar_from_fits_corr = sqrt( fabs( P_L0bar_L0bar_from_fits_err*P_L0bar_L0bar_from_fits_err - P_L0bar_L0bar_tight_pT_cuts_from_fits_signal_err*P_L0bar_L0bar_tight_pT_cuts_from_fits_signal_err ) );

      float SysErr_tight_pT_cuts_L0bar_L0bar_from_fits_work = ( fabs( P_L0bar_L0bar_from_fits - P_L0bar_L0bar_tight_pT_cuts_from_fits_signal) - SysErr_tight_pT_cuts_L0bar_L0bar_from_fits_corr )/fabs(P_L0bar_L0bar_from_fits);

      if( SysErr_tight_pT_cuts_L0bar_L0bar_from_fits_work > 0 ) SysErr_tight_pT_cuts_L0bar_L0bar_from_fits = SysErr_tight_pT_cuts_L0bar_L0bar_from_fits_work; //store sys. err. only if it is larger than statistical fluctuations

      SysErr_tight_pT_cuts_L0bar_L0bar_from_fits_no_corr = fabs( P_L0bar_L0bar_from_fits - P_L0bar_L0bar_tight_pT_cuts_from_fits_signal)/fabs(P_L0bar_L0bar_from_fits);


      SysErr_sum_from_fits_pT_cut += SysErr_tight_pT_cuts_L0bar_L0bar_from_fits*fabs(P_L0bar_L0bar_from_fits)/fabs(P_L0bar_L0bar_from_fits_err);
      SysErr_sum_from_fits_of_w_pT_cut += 1./fabs(P_L0bar_L0bar_from_fits_err);

      //----------------------------------------------------------------


      float P_L0bar_L0bar_tight_pT_cuts = fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
      float P_L0bar_L0bar_tight_pT_cuts_err = fitL0bar_L0bar_tight_pT_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->GetParError(1)/(L0bar_alpha*L0bar_alpha);

      //statistical error correction for systematic error
      float SysErr_tight_pT_cuts_L0bar_L0bar_corr = sqrt( fabs( P_L0bar_L0bar_err_2*P_L0bar_L0bar_err_2 - P_L0bar_L0bar_tight_pT_cuts_err*P_L0bar_L0bar_tight_pT_cuts_err ) );

      float SysErr_tight_pT_cuts_L0bar_L0bar_work = ( fabs( P_L0bar_L0bar_2 - P_L0bar_L0bar_tight_pT_cuts) - SysErr_tight_pT_cuts_L0bar_L0bar_corr )/fabs(P_L0bar_L0bar_2);

      if( SysErr_tight_pT_cuts_L0bar_L0bar_work > 0 ) SysErr_tight_pT_cuts_L0bar_L0bar = SysErr_tight_pT_cuts_L0bar_L0bar_work; //store sys. err. only if it is larger than statistical fluctuations

      SysErr_tight_pT_cuts_L0bar_L0bar_no_corr = fabs( P_L0bar_L0bar_2 - P_L0bar_L0bar_tight_pT_cuts)/fabs(P_L0bar_L0bar_2);


      SysErr_sum_pT_cut += SysErr_tight_pT_cuts_L0bar_L0bar*fabs(P_L0bar_L0bar_2)/fabs(P_L0bar_L0bar_err_2);
      SysErr_sum_of_w_pT_cut += 1./fabs(P_L0bar_L0bar_err_2);

      //___________________________________________________________________________________________________________________________________

      //independent fit of US-US and background
      pair<float,float> P_L0bar_L0bar_tight_DCA_cuts_from_fits_signal_pair = Get_P_LL_from_fit(fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US[delta_eta_bin-1], fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi_US_LS[delta_eta_bin-1], 2);

      float P_L0bar_L0bar_tight_DCA_cuts_from_fits_signal = P_L0bar_L0bar_tight_DCA_cuts_from_fits_signal_pair.first;
      float P_L0bar_L0bar_tight_DCA_cuts_from_fits_signal_err = P_L0bar_L0bar_tight_DCA_cuts_from_fits_signal_pair.second;


      //statistical error correction for systematic error
      float SysErr_tight_DCA_cuts_L0bar_L0bar_from_fits_corr = sqrt( fabs( P_L0bar_L0bar_from_fits_err*P_L0bar_L0bar_from_fits_err - P_L0bar_L0bar_tight_DCA_cuts_from_fits_signal_err*P_L0bar_L0bar_tight_DCA_cuts_from_fits_signal_err ) );

      float SysErr_tight_DCA_cuts_L0bar_L0bar_from_fits_work = ( fabs( P_L0bar_L0bar_from_fits - P_L0bar_L0bar_tight_DCA_cuts_from_fits_signal) - SysErr_tight_DCA_cuts_L0bar_L0bar_from_fits_corr )/fabs(P_L0bar_L0bar_from_fits);

      if( SysErr_tight_DCA_cuts_L0bar_L0bar_from_fits_work > 0 ) SysErr_tight_DCA_cuts_L0bar_L0bar_from_fits = SysErr_tight_DCA_cuts_L0bar_L0bar_from_fits_work; //store sys. err. only if it is larger than statistical fluctuations

      SysErr_tight_DCA_cuts_L0bar_L0bar_from_fits_no_corr = fabs( P_L0bar_L0bar_from_fits - P_L0bar_L0bar_tight_DCA_cuts_from_fits_signal)/fabs(P_L0bar_L0bar_from_fits);


      //SysErr_sum_from_fits_DCA_cut += SysErr_tight_DCA_cuts_L0bar_L0bar_from_fits*fabs(P_L0bar_L0bar_from_fits)/SysErr_tight_DCA_cuts_L0bar_L0bar_from_fits_corr;
      //SysErr_sum_from_fits_of_w_DCA_cut += 1./SysErr_tight_DCA_cuts_L0bar_L0bar_from_fits_corr;

      SysErr_sum_from_fits_DCA_cut += SysErr_tight_DCA_cuts_L0bar_L0bar_from_fits*fabs(P_L0bar_L0bar_from_fits)/fabs(P_L0bar_L0bar_from_fits_err);
      SysErr_sum_from_fits_of_w_DCA_cut += 1./fabs(P_L0bar_L0bar_from_fits_err);

      //----------------------------------------------------------------


      float P_L0bar_L0bar_tight_DCA_cuts = fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->GetParameter(1)/(L0bar_alpha*L0bar_alpha);
      float P_L0bar_L0bar_tight_DCA_cuts_err = fitL0bar_L0bar_tight_DCA_cuts_Delta_y_Delta_phi[delta_eta_bin-1]->GetParError(1)/(L0bar_alpha*L0bar_alpha);

      //statistical error correction for systematic error
      float SysErr_tight_DCA_cuts_L0bar_L0bar_corr = sqrt( fabs( P_L0bar_L0bar_err_2*P_L0bar_L0bar_err_2 - P_L0bar_L0bar_tight_DCA_cuts_err*P_L0bar_L0bar_tight_DCA_cuts_err ) );

      float SysErr_tight_DCA_cuts_L0bar_L0bar_work = ( fabs( P_L0bar_L0bar_2 - P_L0bar_L0bar_tight_DCA_cuts) - SysErr_tight_DCA_cuts_L0bar_L0bar_corr )/fabs(P_L0bar_L0bar_2);

      if( SysErr_tight_DCA_cuts_L0bar_L0bar_work > 0 ) SysErr_tight_DCA_cuts_L0bar_L0bar = SysErr_tight_DCA_cuts_L0bar_L0bar_work; //store sys. err. only if it is larger than statistical fluctuations

      SysErr_tight_DCA_cuts_L0bar_L0bar_no_corr = fabs( P_L0bar_L0bar_2 - P_L0bar_L0bar_tight_DCA_cuts)/fabs(P_L0bar_L0bar_2);


      //SysErr_sum_DCA_cut += SysErr_tight_DCA_cuts_L0bar_L0bar*fabs(P_L0bar_L0bar_2)/SysErr_tight_DCA_cuts_L0bar_L0bar_corr;
      //SysErr_sum_of_w_DCA_cut += 1./SysErr_tight_DCA_cuts_L0bar_L0bar_corr;

      SysErr_sum_DCA_cut += SysErr_tight_DCA_cuts_L0bar_L0bar*fabs(P_L0bar_L0bar_2)/fabs(P_L0bar_L0bar_err_2);
      SysErr_sum_of_w_DCA_cut += 1./fabs(P_L0bar_L0bar_err_2);

      //___________________________________________________________________________________________________________________________________

      SysErr_tight_cuts_L0bar_L0bar = sqrt(SysErr_tight_topo_cuts_L0bar_L0bar*SysErr_tight_topo_cuts_L0bar_L0bar + SysErr_tight_pT_cuts_L0bar_L0bar*SysErr_tight_pT_cuts_L0bar_L0bar + SysErr_tight_DCA_cuts_L0bar_L0bar*SysErr_tight_DCA_cuts_L0bar_L0bar);

      SysErr_tight_cuts_L0bar_L0bar_no_corr = sqrt(SysErr_tight_topo_cuts_L0bar_L0bar_no_corr*SysErr_tight_topo_cuts_L0bar_L0bar_no_corr + SysErr_tight_pT_cuts_L0bar_L0bar_no_corr*SysErr_tight_pT_cuts_L0bar_L0bar_no_corr +
                                                SysErr_tight_DCA_cuts_L0bar_L0bar_no_corr*SysErr_tight_DCA_cuts_L0bar_L0bar_no_corr);

      //----------------------------

      SysErr_tight_cuts_L0bar_L0bar_from_fits = sqrt(SysErr_tight_topo_cuts_L0bar_L0bar_from_fits*SysErr_tight_topo_cuts_L0bar_L0bar_from_fits + SysErr_tight_pT_cuts_L0bar_L0bar_from_fits*SysErr_tight_pT_cuts_L0bar_L0bar_from_fits +
                                                  SysErr_tight_DCA_cuts_L0bar_L0bar_from_fits*SysErr_tight_DCA_cuts_L0bar_L0bar_from_fits);

      SysErr_tight_cuts_L0bar_L0bar_from_fits_no_corr = sqrt(SysErr_tight_topo_cuts_L0bar_L0bar_from_fits_no_corr*SysErr_tight_topo_cuts_L0bar_L0bar_from_fits_no_corr + SysErr_tight_pT_cuts_L0bar_L0bar_from_fits_no_corr*SysErr_tight_pT_cuts_L0bar_L0bar_from_fits_no_corr +
                                                          SysErr_tight_DCA_cuts_L0bar_L0bar_from_fits_no_corr*SysErr_tight_DCA_cuts_L0bar_L0bar_from_fits_no_corr );

    }

    //total
    //float SysErrTot_L0bar_L0bar = sqrt( sysErr_alpha_L0bar_L0bar*sysErr_alpha_L0bar_L0bar + SysErrSlope_L0bar_L0bar*SysErrSlope_L0bar_L0bar + SysErrBackground_L0bar_L0bar*SysErrBackground_L0bar_L0bar  );
    //float SysErrTot_L0bar_L0bar = sqrt( sysErr_alpha_L0bar_L0bar*sysErr_alpha_L0bar_L0bar + SysErrSlope_L0bar_L0bar*SysErrSlope_L0bar_L0bar + SysErrBackground_L0bar_L0bar*SysErrBackground_L0bar_L0bar + SysErr_tight_cuts_L0bar_L0bar*SysErr_tight_cuts_L0bar_L0bar );
    float SysErrTot_L0bar_L0bar = sqrt( sysErr_alpha_L0bar_L0bar*sysErr_alpha_L0bar_L0bar + SysErrSlope_L0bar_L0bar*SysErrSlope_L0bar_L0bar + SysErrBackground_L0bar_L0bar*SysErrBackground_L0bar_L0bar + SysErr_tight_cuts_L0bar_L0bar_from_fits*SysErr_tight_cuts_L0bar_L0bar_from_fits );

    //TPaveText *L0bar_L0bar_text_2 = new TPaveText(0.47, 0.15, 0.85, 0.53, "NDC"); //with polarization value
    TPaveText *L0bar_L0bar_text_2 = new TPaveText(0.47, 0.15, 0.8, 0.53, "NDC"); //without polarization value
    L0bar_L0bar_text_2->SetTextFont(42);
    //L0bar_L0bar_text_2->SetTextSize(15);
    //L0bar_L0bar_text_2->AddText("STAR");
    //L0bar_L0bar_text_2->AddText("STAR preliminary");
    //((TText*)L0bar_L0bar_text_2->GetListOfLines()->Last())->SetTextColor(2);
    L0bar_L0bar_text_2->AddText(Form("%i p+p #sqrt{s} = %i GeV", year, energy));
    L0bar_L0bar_text_2->AddText("Minimum bias");
    L0bar_L0bar_text_2->AddText("#bar{#Lambda^{0}}-#bar{#Lambda^{0}}");
    L0bar_L0bar_text_2->AddText("|#it{y}| < 1");
    L0bar_L0bar_text_2->AddText("0.5 < p_{T,1} < 5.0 GeV/#it{c}");
    L0bar_L0bar_text_2->AddText("0.5 < p_{T,2} < 5.0 GeV/#it{c}");
    L0bar_L0bar_text_2->AddText(Form("P_{#Lambda_{1}#Lambda_{2}} = %.3f #pm %.3f", P_L0bar_L0bar_2, fabs(P_L0bar_L0bar_err_2)));
    L0bar_L0bar_text_2->AddText(Form("P_{#Lambda_{1}#Lambda_{2}} = %.3f #pm %.3f", P_L0bar_L0bar_from_fits, fabs(P_L0bar_L0bar_from_fits_err)));
    //L0bar_L0bar_text_2->AddText(Form("P_{topo} = %.2f", P_L0bar_L0bar_tight_topo_cuts));
    //L0bar_L0bar_text_2->AddText(Form("P_{pT} = %.2f", P_L0bar_L0bar_tight_pT_cuts));
    L0bar_L0bar_text_2->SetFillColorAlpha(0, 0.01);
    L0bar_L0bar_text_2->Draw("same");


    TLegend *L0bar_L0bar_2_leg = new TLegend(0.15, 0.3, 0.4, 0.49);
    //L0bar_L0bar_2_leg->SetTextSizePixels(15);
    L0bar_L0bar_2_leg->AddEntry(L0bar_L0bar_cosThetaProdPlane_eta_US_hist, "(US-US)-Bckg.");
    L0bar_L0bar_2_leg->AddEntry(fitL0bar_L0bar_US_ThetaStar_2, "Fit", "l");
    L0bar_L0bar_2_leg->SetBorderSize(0);
    L0bar_L0bar_2_leg->SetFillColorAlpha(0, 0.01);
    L0bar_L0bar_2_leg->Draw("same");

    if(delta_eta_bin == 1) Delta_eta_text_1->Draw("same");
    if(delta_eta_bin == 2) Delta_eta_text_2->Draw("same");

    L0bar_L0bar_cosThetaProdPlane_eta_can_2->SaveAs(Form("./plots/Lambda/L_correlations_delta_eta_delta_phi/L0bar_L0bar_cosThetaProdPlane_eta_subtract_less_delta_eta_delta_phi_%i.png", delta_eta_bin));


    //---------------------------------------------



    //PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->SetPoint(3, P_L0bar_L0bar_2, 3);
    //PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->SetPointError(3, fabs(P_L0bar_L0bar_err_2), 0);

    PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->SetPoint(3, P_L0bar_L0bar_from_fits, 3);
    PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->SetPointError(3, fabs(P_L0bar_L0bar_from_fits_err), 0);


    //PolarizationGraph_delta_eta_delta_phi_sys_err[delta_eta_bin-1]->SetPoint(3, P_L0bar_L0bar_2, 3);
    //PolarizationGraph_delta_eta_delta_phi_sys_err[delta_eta_bin-1]->SetPointError(3, fabs(SysErrTot_L0bar_L0bar*P_L0bar_L0bar_2), 0.045);

    PolarizationGraph_delta_eta_delta_phi_sys_err[delta_eta_bin-1]->SetPoint(3, P_L0bar_L0bar_from_fits, 3);
    PolarizationGraph_delta_eta_delta_phi_sys_err[delta_eta_bin-1]->SetPointError(3, fabs(SysErrTot_L0bar_L0bar*P_L0bar_L0bar_from_fits), 0.045);


    //____________________________________________________________________________________________________________________________________________________________________________________________________________

    //save polarization graphs
    out_file->cd();


    PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->Write(Form("PolarizationGraph_delta_y_delta_phi_%i", delta_eta_bin-1));
    PolarizationGraph_delta_eta_delta_phi_sys_err[delta_eta_bin-1]->Write(Form("PolarizationGraph_delta_y_delta_phi_sys_err_%i", delta_eta_bin-1));


    //____________________________________________________________________________________________________________________________________________________________________________________________________________



    //----------------------------------L pair stats---------------------------------
    cout<<endl;
    cout<<"Bin "<<delta_eta_bin<<endl;
    cout<<endl;
    cout<<"N L-Lbar pairs from hist: "<<nLLbar<<endl;
    cout<<"N L-Lbar background pairs from hist: "<<nLLbar_back<<endl;
    cout<<"Signal/Background "<<(nLLbar - nLLbar_back)/nLLbar_back<<endl;
    cout<<endl;
    cout<<"N L-L pairs from hist: "<<nLL<<endl;
    cout<<"N L-L background pairs from hist: "<<nLL_back<<endl;
    cout<<"Signal/Background "<<(nLL - nLL_back)/nLL_back<<endl;
    cout<<endl;
    cout<<"N Lbar-Lbar pairs from hist: "<<nLbarLbar<<endl;
    cout<<"N Lbar-Lbar background pairs from hist: "<<nLbarLbar_back<<endl;
    cout<<"Signal/Background "<<(nLbarLbar - nLbarLbar_back)/nLbarLbar_back<<endl;
    cout<<endl;
    cout<<endl;

    cout<<"P LLbar = "<<P_L0_L0bar_2<<" +- "<<P_L0_L0bar_err_2<<" +- "<<SysErrTot_L0_L0bar*P_L0_L0bar_2<<endl;
    cout<<"P LL = "<<P_L0_L0_2<<" +- "<<P_L0_L0_err_2<<" +- "<<SysErrTot_L0_L0*P_L0_L0_2<<endl;
    cout<<"P LbarLbar = "<<P_L0bar_L0bar_2<<" +- "<<P_L0bar_L0bar_err_2<<" +- "<<SysErrTot_L0bar_L0bar*P_L0bar_L0bar_2<<endl;
    cout<<endl;

    cout<<"From fits:"<<endl;
    cout<<"P LLbar = "<<P_L0_L0bar_from_fits<<" +- "<<P_L0_L0bar_from_fits_err<<" +- "<<SysErrTot_L0_L0bar*P_L0_L0bar_2<<endl;
    cout<<"P LL = "<<P_L0_L0_from_fits<<" +- "<<P_L0_L0_from_fits_err<<" +- "<<SysErrTot_L0_L0*P_L0_L0_2<<endl;
    cout<<"P LbarLbar = "<<P_L0bar_L0bar_from_fits<<" +- "<<P_L0bar_L0bar_from_fits_err<<" +- "<<SysErrTot_L0bar_L0bar*P_L0bar_L0bar_2<<endl;
    cout<<endl;

    cout<<"Systematic uncetainties by source:"<<endl;
    cout<<"pair | s_topo | s_pT | s_DCA | s_bckg | s_ME | s_alpha"<<endl;
    cout<<"LLbar: "<<P_L0_L0bar_2*SysErr_tight_topo_cuts_L0_L0bar<<" | "<<P_L0_L0bar_2*SysErr_tight_pT_cuts_L0_L0bar<<" | "<<P_L0_L0bar_2*SysErr_tight_DCA_cuts_L0_L0bar<<" | "
    <<P_L0_L0bar_2*SysErrBackground_L0_L0bar<<" | "<<P_L0_L0bar_2*SysErrSlope_L0_L0bar<<" | "<<P_L0_L0bar_2*sysErr_alpha_L0_L0bar<<endl;

    cout<<"LL: "<<P_L0_L0_2*SysErr_tight_topo_cuts_L0_L0<<" | "<<P_L0_L0_2*SysErr_tight_pT_cuts_L0_L0<<" | "<<P_L0_L0_2*SysErr_tight_DCA_cuts_L0_L0<<" | "
    <<P_L0_L0_2*SysErrBackground_L0_L0<<" | "<<P_L0_L0_2*SysErrSlope_L0_L0<<" | "<<P_L0_L0_2*sysErr_alpha_L0_L0<<endl;

    cout<<"LbarLbar: "<<P_L0bar_L0bar_2*SysErr_tight_topo_cuts_L0bar_L0bar<<" | "<<P_L0bar_L0bar_2*SysErr_tight_pT_cuts_L0bar_L0bar<<" | "<<P_L0bar_L0bar_2*SysErr_tight_DCA_cuts_L0bar_L0bar<<" | "
    <<P_L0bar_L0bar_2*SysErrBackground_L0bar_L0bar<<" | "<<P_L0bar_L0bar_2*SysErrSlope_L0bar_L0bar<<" | "<<P_L0bar_L0bar_2*sysErr_alpha_L0bar_L0bar<<endl;
    cout<<endl;

    cout<<"Systematic uncetainties by source from fits:"<<endl;
    cout<<"pair | s_topo | s_pT | s_DCA | s_bckg | s_ME | s_alpha"<<endl;
    cout<<"LLbar: "<<P_L0_L0bar_from_fits*SysErr_tight_topo_cuts_L0_L0bar_from_fits<<" | "<<P_L0_L0bar_from_fits*SysErr_tight_pT_cuts_L0_L0bar_from_fits<<" | "<<P_L0_L0bar_from_fits*SysErr_tight_DCA_cuts_L0_L0bar_from_fits<<" | "
    <<P_L0_L0bar_from_fits*SysErrBackground_L0_L0bar<<" | "<<P_L0_L0bar_from_fits*SysErrSlope_L0_L0bar<<" | "<<P_L0_L0bar_from_fits*sysErr_alpha_L0_L0bar<<endl;

    cout<<"LL: "<<P_L0_L0_from_fits*SysErr_tight_topo_cuts_L0_L0_from_fits<<" | "<<P_L0_L0_from_fits*SysErr_tight_pT_cuts_L0_L0_from_fits<<" | "<<P_L0_L0_from_fits*SysErr_tight_DCA_cuts_L0_L0_from_fits<<" | "
    <<P_L0_L0_from_fits*SysErrBackground_L0_L0<<" | "<<P_L0_L0_from_fits*SysErrSlope_L0_L0<<" | "<<P_L0_L0_from_fits*sysErr_alpha_L0_L0<<endl;

    cout<<"LbarLbar: "<<P_L0bar_L0bar_from_fits*SysErr_tight_topo_cuts_L0bar_L0bar_from_fits<<" | "<<P_L0bar_L0bar_from_fits*SysErr_tight_pT_cuts_L0bar_L0bar_from_fits<<" | "<<P_L0bar_L0bar_from_fits*SysErr_tight_DCA_cuts_L0bar_L0bar_from_fits<<" | "
    <<P_L0bar_L0bar_from_fits*SysErrBackground_L0bar_L0bar<<" | "<<P_L0bar_L0bar_from_fits*SysErrSlope_L0bar_L0bar<<" | "<<P_L0bar_L0bar_from_fits*sysErr_alpha_L0bar_L0bar<<endl;
    cout<<endl;


  }

  cout<<"Relative sys. err. of alpha:"<<endl;
  cout<<"LLbar: "<<sysErr_alpha_L0_L0bar<<endl;
  cout<<"LL: "<<sysErr_alpha_L0_L0<<endl;
  cout<<"LbarLbar: "<<sysErr_alpha_L0bar_L0bar<<endl;
  cout<<endl;


  float SysErr_topo_weights = SysErr_sum_topo_cuts/SysErr_sum_of_w_topo_cuts;
  float SysErr_pT_weights = SysErr_sum_pT_cut/SysErr_sum_of_w_pT_cut;
  float SysErr_DCA_weights = SysErr_sum_DCA_cut/SysErr_sum_of_w_DCA_cut;
  float SysErr_bckg_weights = SysErr_sum_background/SysErr_sum_of_w_background;
  float SysErr_ME_weights = SysErr_sum_ME/SysErr_sum_of_w_ME;

  cout<<"Systematic uncetainty weighted averages by source:"<<endl;
  cout<<"s_topo | s_pT | s_DCA | s_bckg | s_ME"<<endl;
  cout<<SysErr_topo_weights<<" | "<<SysErr_pT_weights<<" | "<<SysErr_DCA_weights<<" | "<<SysErr_bckg_weights<<" | "<<SysErr_ME_weights<<endl;
  cout<<endl;

  float SysErrTot_weights_work = sqrt(SysErr_topo_weights*SysErr_topo_weights + SysErr_pT_weights*SysErr_pT_weights + SysErr_DCA_weights*SysErr_DCA_weights + SysErr_bckg_weights*SysErr_bckg_weights + SysErr_ME_weights*SysErr_ME_weights);

  cout<<"Total uncertainty from weighter average:"<<endl;
  cout<<SysErrTot_weights_work<<endl;
  cout<<endl;

  //-------------------------------------

  float SysErr_topo_weights_from_fits = SysErr_sum_from_fits_topo_cuts/SysErr_sum_from_fits_of_w_topo_cuts;
  float SysErr_pT_weights_from_fits = SysErr_sum_from_fits_pT_cut/SysErr_sum_from_fits_of_w_pT_cut;
  float SysErr_DCA_weights_from_fits = SysErr_sum_from_fits_DCA_cut/SysErr_sum_from_fits_of_w_DCA_cut;

  cout<<"Systematic uncetainty weighted averages by source:"<<endl;
  cout<<"s_topo | s_pT | s_DCA | s_bckg | s_ME"<<endl;
  cout<<SysErr_topo_weights_from_fits<<" | "<<SysErr_pT_weights_from_fits<<" | "<<SysErr_DCA_weights_from_fits<<" | "<<SysErr_bckg_weights<<" | "<<SysErr_ME_weights<<endl;
  cout<<endl;

  float SysErrTot_weights_from_fits_work = sqrt(SysErr_topo_weights_from_fits*SysErr_topo_weights_from_fits + SysErr_pT_weights_from_fits*SysErr_pT_weights_from_fits + SysErr_DCA_weights_from_fits*SysErr_DCA_weights_from_fits +
                                                SysErr_bckg_weights*SysErr_bckg_weights + SysErr_ME_weights*SysErr_ME_weights);

  cout<<"Total uncertainty from fits from weighter average:"<<endl;
  cout<<SysErrTot_weights_from_fits_work<<endl;
  cout<<endl;

  out_file->cd();

  for( unsigned int delta_eta_bin = 1; delta_eta_bin < 3; delta_eta_bin++)
  {

    //L-Lbar

    //float L0L0bar_SysErrTot_weights_work_relat = SysErrTot_weights_work/PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->GetPointX(1); //point 1 is LLbar
    float L0L0bar_SysErrTot_weights_work_relat = SysErrTot_weights_from_fits_work/PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->GetPointX(1); //point 1 is LLbar

    float L0L0bar_SysErr_DCA_relat = sysErr_DCA_abs/PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->GetPointX(1);
    //float L0L0bar_SysErr_embed_relat = sysErr_embed_abs/PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->GetPointX(1);

    float L0L0bar_SysErrTot_weights = sqrt( L0L0bar_SysErrTot_weights_work_relat*L0L0bar_SysErrTot_weights_work_relat + sysErr_alpha_L0_L0bar*sysErr_alpha_L0_L0bar + L0L0bar_SysErr_DCA_relat*L0L0bar_SysErr_DCA_relat );



    PolarizationGraph_delta_eta_delta_phi_sys_err_average[delta_eta_bin-1]->SetPoint(1, PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->GetPointX(1), 1);
    PolarizationGraph_delta_eta_delta_phi_sys_err_average[delta_eta_bin-1]->SetPointError(1, L0L0bar_SysErrTot_weights*PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->GetPointX(1), 0.045);


    //L-L

    //float L0L0_SysErrTot_weights_work_relat = SysErrTot_weights_work/PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->GetPointX(2); //point 1 is LLbar
    float L0L0_SysErrTot_weights_work_relat = SysErrTot_weights_from_fits_work/PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->GetPointX(2); //point 1 is LLbar

    float L0L0_SysErr_DCA_relat = sysErr_DCA_abs/PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->GetPointX(2);
    //float L0L0_SysErr_embed_relat = sysErr_embed_abs/PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->GetPointX(2);

    float L0L0_SysErrTot_weights = sqrt( L0L0_SysErrTot_weights_work_relat*L0L0_SysErrTot_weights_work_relat + sysErr_alpha_L0_L0*sysErr_alpha_L0_L0 + L0L0_SysErr_DCA_relat*L0L0_SysErr_DCA_relat );


    PolarizationGraph_delta_eta_delta_phi_sys_err_average[delta_eta_bin-1]->SetPoint(2, PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->GetPointX(2), 2);
    PolarizationGraph_delta_eta_delta_phi_sys_err_average[delta_eta_bin-1]->SetPointError(2, L0L0_SysErrTot_weights*PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->GetPointX(2), 0.045);

    //Lbar-Lbar

    //float L0barL0bar_SysErrTot_weights_work_relat = SysErrTot_weights_work/PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->GetPointX(3); //point 1 is LLbar
    float L0barL0bar_SysErrTot_weights_work_relat = SysErrTot_weights_from_fits_work/PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->GetPointX(3); //point 1 is LLbar

    float L0barL0bar_SysErr_DCA_relat = sysErr_DCA_abs/PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->GetPointX(3);
    //float L0barL0bar_SysErr_embed_relat = sysErr_embed_abs/PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->GetPointX(3);

    float L0barL0bar_SysErrTot_weights = sqrt( L0barL0bar_SysErrTot_weights_work_relat*L0barL0bar_SysErrTot_weights_work_relat + sysErr_alpha_L0bar_L0bar*sysErr_alpha_L0bar_L0bar + L0barL0bar_SysErr_DCA_relat*L0barL0bar_SysErr_DCA_relat );


    PolarizationGraph_delta_eta_delta_phi_sys_err_average[delta_eta_bin-1]->SetPoint(3, PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->GetPointX(3), 3);
    PolarizationGraph_delta_eta_delta_phi_sys_err_average[delta_eta_bin-1]->SetPointError(3, L0barL0bar_SysErrTot_weights*PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->GetPointX(3), 0.045);


    PolarizationGraph_delta_eta_delta_phi_sys_err_average[delta_eta_bin-1]->Write(Form("PolarizationGraph_delta_y_delta_phi_sys_err_average_%i", delta_eta_bin-1));

    //----------------------------------------------------------------------------------------

    //plot polarization graphs in bins
    PolarizationGraph_delta_eta_delta_phi_can->cd(delta_eta_bin);

    //gPad->SetLeftMargin(0.17);
    gPad->SetRightMargin(0.05);

    DefaultHist->Draw();

    PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->SetMarkerStyle(20);
    PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->SetMarkerSize(2);
    PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->SetMarkerColor(kRed);
    PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->SetLineColor(kRed);
    PolarizationGraph_delta_eta_delta_phi[delta_eta_bin-1]->Draw("p e same");

    PolarizationGraph_delta_eta_delta_phi_sys_err[delta_eta_bin-1]->SetMarkerSize(2);
    PolarizationGraph_delta_eta_delta_phi_sys_err[delta_eta_bin-1]->SetFillColorAlpha(kRed, 0.25);
    PolarizationGraph_delta_eta_delta_phi_sys_err[delta_eta_bin-1]->SetLineColor(kRed);
    //PolarizationGraph_delta_eta_delta_phi_sys_err[delta_eta_bin-1]->Draw("2 same");

    PolarizationGraph_delta_eta_delta_phi_sys_err_average[delta_eta_bin-1]->SetMarkerSize(2);
    PolarizationGraph_delta_eta_delta_phi_sys_err_average[delta_eta_bin-1]->SetFillColorAlpha(kRed, 0.25);
    PolarizationGraph_delta_eta_delta_phi_sys_err_average[delta_eta_bin-1]->SetLineColor(kRed);
    if(year != 2016) PolarizationGraph_delta_eta_delta_phi_sys_err_average[delta_eta_bin-1]->Draw("2 same");

    ZeroLine_eta->Draw("same");

    bin_text[delta_eta_bin-1] = new TPaveText(0.12, 0.85, 0.49, 0.89, "NDC");
    bin_text[delta_eta_bin-1]->SetTextFont(42);
    if(delta_eta_bin == 1) bin_text[delta_eta_bin-1]->AddText("|#Delta#it{y}| < 0.5, |#Delta#phi| < #pi/3");
    if(delta_eta_bin == 2) bin_text[delta_eta_bin-1]->AddText("0.5 < |#Delta#it{y}| < 2.0 or |#Delta#phi| > #pi/3");
    bin_text[delta_eta_bin-1]->SetFillColorAlpha(0, 0.01);
    bin_text[delta_eta_bin-1]->Draw("same");

    if(delta_eta_bin == 1 ) Polarization_text->Draw("same");

  }


  PolarizationGraph_delta_eta_delta_phi_can->SaveAs("./plots/Lambda/L_polarization/L_polarization_delta_eta_delta_phi.png");

  //__________________________________________________________________________________________________




  LLbarOutFile->Close();
  out_file->Close();

  return true;
}
