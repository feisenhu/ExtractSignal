// Task for Peak extraction
//
#include "TLatex.h"
#include "TMath.h"
#include "TGaxis.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
// #include "TRatioPlot.h"
#include "TObjArray.h"
#include "TSystem.h"
#include "TLegend.h"
#include "ExtractSignal.h"
#include <iostream>

Bool_t DoDebug                = kFALSE;

Bool_t DoRebin                = kTRUE;    // Do Rebinning of TH2 Histograms

// Use Either PirmSec or SecSec, not both!
Bool_t PrimSecPairing         = kTRUE;
Bool_t SecSecPairing          = !PrimSecPairing;

Bool_t PrimarySoBRatio        = kFALSE;    // Plot S/B Ratio for PrimaryPairs and find best region to apply primary mass cut

Bool_t ExtPairSig             = kFALSE;    // Extract Pair Signal
Bool_t ExtFourPairSig         = kTRUE;    // Extract Four Pair Signal
Bool_t ExtFourPairPionSig     = kFALSE;
Bool_t ExtFourPairEtaSig      = kTRUE;
  Bool_t DoPtProj               = kFALSE;    // Do Pt Projection
  Bool_t DoMassProj             = kTRUE;    // Do Mass Projection
    Bool_t DoFit                  = kTRUE;    // Fitting Backround
    Bool_t DoTemplateFit          = kTRUE;    // true=fit data with a*MCSignalHist+b*MCBackgrounfHist , false = standard TF1 gauss exp  and pol2 fits ; Do Fit of true reconstructed Singal and reconstructed Background
    // Bool_t DoSignalOverBackground = kTRUE;    // Do Signal over Background Ratio (based on DoMassProj)

  Bool_t ExtGen                 = kFALSE;    // Extract Generated Signal
  Bool_t ExtGenSmeared          = kFALSE;    // Extract Generated Smeared Signal
  Bool_t ExtRec                 = kTRUE;    // Extract Reconstructed Signal

Bool_t ExtDielectronPairSpectra = kFALSE;  // Extract the TH2 mass vs pt histogram from the reconstructed dielectron pairs

Bool_t ExtPreFilterSpectra    = kFALSE;     // Extract the Prefilter histograms before and after the prefilters are applied

Bool_t ExtMassCut             = kFALSE;    // Extract Mass Cut Histogram
Bool_t DoCutEff               = kFALSE;    // Plotting Cut Efficiency

Bool_t DoCompGenSmearRec      = kFALSE;    // Comparison between GeneratedSmeared and Reconstructed in Pt, Eta and Phi on single electron level
Bool_t DoLossBySecondaryCuts  = kFALSE;    // Plotting the loss of secondaries, for variation of cuts
Bool_t DoPlotPIDComponents    = kFALSE;    // Plotting each particle contibution as a projection to the number of sigmas for each detector

double MassCutPrimary   = 0.547862;
double MassCutSecondary = 0.01;

// String is used to keep order in output path
// TString TrainNumber = "local_test";
// TString TrainNumber = "merged_309_LHC18h1_310_LHC17d1_lowField_311_LHC17h3_WithOutPreFilter_NewMissMatchSignals";
// TString TrainNumber = "merged_312_LHC18h1_313_LHC17d1_LF_314_LHC17h3_widerPreFilter";
// TString TrainNumber = "merged_315_LHC18h1_316_LHC17d1_LF_317_LHC17h3_without_Prefilter_primprecuts_pionmasscut";
// TString TrainNumber = "merged_318_LHC18h1_319_LHC17d1_LF_320_LHC17h3_Prefilter_primprecuts_pionmasscut";
// TString TrainNumber = "merged_321_LHC18h1_322_LHC17d1_LF_323_LHC17h3_SecSecPreFilter";
// TString TrainNumber = "merged_325_17d1_LF_326_18h1_327_17h3_CorrectedV0selection";
// TString TrainNumber = "merged_332_LHC17h3_333_LHC18h1_334LHC17d1_LF_RemoveTracksInPreFilter_DataLikeSignals";
// TString TrainNumber = "merged_335_LHC17h3_336_LHC18h1_337_LHC17d1_LF_PrimSecPrefilter_0.1-0.165";
// TString TrainNumber = "merged_338_LHC17h3_339_LHC18h1_340_LHC17d1_LF_separated_PIDHistos_oldPreFiltercuts";
// TString TrainNumber = "merged_341_LHC17h3_342_LHC18h1_343_LHC17d1_LF_AllTrackPairsPt75Eta0.8Cut";
// TString TrainNumber = "merged_347_LHC17d1_LF_348_LHC18h1_351_LHC17h3_openMassCuts";
// TString TrainNumber = "merged_352_LHC17d1_LF_353_LHC18h1_354_LHC17h3_TrainWithSecSecFourPairing";
// TString TrainNumber = "merged_364_LHC17d1_LF_365_LHC17h3_366_LHC18h1_GammaGamma_withoutPreFilter_NoMassCut";
// TString TrainNumber = "merged_376_LHC17d1_LF_377_LHC17h3_378_LHC18h1_Dalitz_withoutPreFilter_NoMassCut";
// TString TrainNumber = "merged_403_LHC18h1_child3_404_LHC17h3_405_LHC17d1_LF_407_LHC18h1_child1+2_Dalitz_withPreFilter_withMasscut0-0.35_lowerSplitLevel";
// TString TrainNumber = "merged_408_LHC18H1_409_LHC17h3_410_LHC17d1_LF_Dalitz_GammaGamma_noMassCut_withPrefilter";
// TString TrainNumber = "merged_412_LHC18h1_413_LHC17h3_414_LHC17d1_LF_OnlyRec_DalitzGamma_withPrefilter_0.1-0.2MassCut";
TString TrainNumber = "merged_419_LHC18h1_420_LHC17h3_421_LHC17d1_LF_OnlyRec_DalitzGammaGamma_widerSecSecPrefilter_0-0.35MassCut";


void ExtractSignal(){

  // std::vector<Double_t> vec_rebin_pt = {0.00,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,1.00,1.10,1.20,1.30,1.40,1.50,1.60,1.80,2.00,2.40,2.80,3.20,4.0}; // pt Intervalls for rebinning
  // pT rebin vector with rebin factors
  std::vector<Int_t> vec_rebin_pt = {2}; // pt Intervalls for rebinning

  // mass rebin vector with rebin factors dependent on pT-projecvec_bin_pttion intervalls
  // std::vector<Int_t> vec_rebin_mass = {1,  1,  1,  1,  1,  1,  1, 1, 1};
  // std::vector<Int_t> vec_rebin_mass = {1,  1,  1,  1,  1,  1,  1};
  // std::vector<Int_t> vec_rebin_mass = {5,  5,  5,  5,  5,  5,  5, 5, 5};
  // std::vector<Int_t> vec_rebin_mass = {5,  5,  5,  5,  5,  5,  5};
  // std::vector<Int_t> vec_rebin_mass = {5,  5,  5,  5,  5};
  // std::vector<Int_t> vec_rebin_mass = {10, 10, 10, 10, 10, 10, 10, 10, 10};
  // std::vector<Int_t> vec_rebin_mass = {10, 10, 10, 10, 10, 10, 10};
  std::vector<Int_t> vec_rebin_mass = {20, 20, 20, 20, 20, 20, 20};
  // std::vector<Int_t> vec_rebin_mass = {20, 20, 20, 20, 20};
  // std::vector<Int_t> vec_rebin_mass = {50, 50, 50, 50, 50, 50, 50};
  // std::vector<Int_t> vec_rebin_mass = {50, 50, 50, 50, 50};
  // std::vector<Int_t> vec_rebin_mass = {75, 75, 75, 75, 75, 75, 75};
  // std::vector<Int_t> vec_rebin_mass = {75, 75, 75, 75, 75};

  // std::vector<Double_t> vec_rebin_mass = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1,                                                                                                                                             /*10 MeV steps*/
                                          // 0.11, 0.112, 0.114, 0.116, 0.118, 0.12, 0.122, 0.124, 0.126, 0.128, 0.13, 0.132, 0.134, 0.135, 0.136, 0.138, 0.14, 0.142, 0.144, 0.146, 0.148, 0.15, 0.152, 0.154, 0.156, 0.158, 0.16,                       /*2  MeV steps*/
                                          // 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52,  /*10 MeV steps*/
                                          // 0.522, 0.524, 0.526, 0.528, 0.53, 0.532, 0.534, 0.536, 0.538, 0.54, 0.542, 0.544, 0.546, 0.548, 0.549, 0.55, 0.552, 0.554, 0.556, 0.558, 0.56, 0.562, 0.564, 0.566, 0.568, 0.57, 0.572, 0.574, 0.576, 0.578, /*2 MeV steps*/
                                          // 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.95, 1.00}; // mass intervalls for rebinning

 // 10 MeV Binning
 // std::vector<Double_t> vec_rebin_mass = {0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49,
                                         // 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00}; // mass intervalls for rebinning
 // 20 MeV Binning + Pion/Eta mass region 10 MeV Bins
 // std::vector<Double_t> vec_rebin_mass = {0.00, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.15, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32, 0.34, 0.36, 0.38, 0.40, 0.42, 0.44, 0.46, 0.48, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.60, 0.62, 0.64, 0.66, 0.68, 0.70, 0.72, 0.74, 0.76, 0.78, 0.80, 0.82, 0.84, 0.86, 0.87, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98, 1.00}; // mass intervalls for rebinning
 // 40 MeV Binning
 // std::vector<Double_t> vec_rebin_mass = {0.00, 0.04, 0.08, 0.12, 0.16, 0.20, 0.24, 0.28, 0.32, 0.36, 0.40, 0.44, 0.48, 0.52, 0.56, 0.60, 0.64, 0.68, 0.72, 0.76, 0.80, 0.84, 0.88, 0.92, 0.96, 1.00}; // mass intervalls for rebinning

  // Intervalls for projection in pt slices (if number of proj is changed, eddit vec_rebin_mass)
  // std::vector<Double_t> vec_proj_bin_pt = {0.0, 0.2, 0.3, 0.4, 0.5, 0.8, 1.0, 2.0, 4.0};
  std::vector<Double_t> vec_proj_bin_pt = {0.0, 0.3, 0.5, 0.8, 1.0, 2.0, 4.0};
  // std::vector<Double_t> vec_proj_bin_pt = {0.0, 0.2, 0.3, 0.4, 0.5};
  std::vector<Double_t> vec_proj_bin_mass = {0.0, 1.0}; // Intervalls for projection in mass slices

                                                                                if(DoDebug) Printf("%d", __LINE__);

  // Open file to read in the required histograms
  // TFile *fFile = new TFile("/data/feisenhut/1_EtaReconstruction/AnalysisResults.root"); // This file should exist
  // TFile *fFile = new TFile("/u/feisenhut/Documents/Uni/Masterarbeit/Sicherung/LegoTrainsOutputs/DQ_pp_MC_AOD/merged_309_LHC18h1_310_LHC17d1_lowField_311_LHC17h3_WithOutPreFilter_NewMissMatchSignals/AnalysisResults.root");
  // TFile *fFile = new TFile("/u/feisenhut/Documents/Uni/Masterarbeit/Sicherung/LegoTrainsOutputs/DQ_pp_MC_AOD/merged_312_LHC18h1_313_LHC17d1_LF_314_LHC17h3_widerPreFilter/AnalysisResults.root");
  // TFile *fFile = new TFile("/u/feisenhut/Documents/Uni/Masterarbeit/Sicherung/LegoTrainsOutputs/DQ_pp_MC_AOD/merged_315_LHC18h1_316_LHC17d1_LF_317_LHC17h3_without_Prefilter_primprecuts_pionmasscut/AnalysisResults.root");
  // TFile *fFile = new TFile("/u/feisenhut/Documents/Uni/Masterarbeit/Sicherung/LegoTrainsOutputs/DQ_pp_MC_AOD/merged_318_LHC18h1_319_LHC17d1_LF_320_LHC17h3_Prefilter_primprecuts_pionmasscut/AnalysisResults.root");
  // TFile *fFile = new TFile("/u/feisenhut/Documents/Uni/Masterarbeit/Sicherung/LegoTrainsOutputs/DQ_pp_MC_AOD/merged_321_LHC18h1_322_LHC17d1_LF_323_LHC17h3_SecSecPreFilter/AnalysisResults.root");
  // TFile *fFile = new TFile("/u/feisenhut/Documents/Uni/Masterarbeit/Sicherung/LegoTrainsOutputs/DQ_pp_MC_AOD/merged_325_17d1_LF_326_18h1_327_17h3_CorrectedV0selection/AnalysisResults.root");
  // TFile *fFile = new TFile("/u/feisenhut/Documents/Uni/Masterarbeit/Sicherung/LegoTrainsOutputs/DQ_pp_MC_AOD/merged_332_LHC17h3_333_LHC18h1_334LHC17d1_LF_RemoveTracksInPreFilter_DataLikeSignals/AnalysisResults.root");
  // TFile *fFile = new TFile("/u/feisenhut/Documents/Uni/Masterarbeit/Sicherung/LegoTrainsOutputs/DQ_pp_MC_AOD/merged_335_LHC17h3_336_LHC18h1_337_LHC17d1_LF_PrimSecPrefilter_0.1-0.165/AnalysisResults.root");
  // TFile *fFile = new TFile("/u/feisenhut/Documents/Uni/Masterarbeit/Sicherung/LegoTrainsOutputs/DQ_pp_MC_AOD/merged_338_LHC17h3_339_LHC18h1_340_LHC17d1_LF_separated_PIDHistos_oldPreFiltercuts/AnalysisResults.root");
  // TFile *fFile = new TFile("/u/feisenhut/Documents/Uni/Masterarbeit/Sicherung/LegoTrainsOutputs/DQ_pp_MC_AOD/merged_341_LHC17h3_342_LHC18h1_343_LHC17d1_LF_AllTrackPairsPt75Eta0.8Cut/AnalysisResults.root");
  // TFile *fFile = new TFile("/u/feisenhut/Documents/Uni/Masterarbeit/Sicherung/LegoTrainsOutputs/DQ_pp_MC_AOD/merged_347_LHC17d1_LF_348_LHC18h1_351_LHC17h3_openMassCuts/AnalysisResults.root");
  // TFile *fFile = new TFile("/u/feisenhut/Documents/Uni/Masterarbeit/Sicherung/LegoTrainsOutputs/DQ_pp_MC_AOD/merged_352_LHC17d1_LF_353_LHC18h1_354_LHC17h3_TrainWithSecSecFourPairing/AnalysisResults.root");
  // TFile *fFile = new TFile("/u/feisenhut/Documents/Uni/Masterarbeit/Sicherung/LegoTrainsOutputs/DQ_pp_MC_AOD/merged_364_LHC17d1_LF_365_LHC17h3_366_LHC18h1_GammaGamma_withoutPreFilter_NoMassCut/AnalysisResults.root");
  // TFile *fFile = new TFile("/u/feisenhut/Documents/Uni/Masterarbeit/Sicherung/LegoTrainsOutputs/DQ_pp_MC_AOD/merged_376_LHC17d1_LF_377_LHC17h3_378_LHC18h1_Dalitz_withoutPreFilter_NoMassCut/AnalysisResults.root");
  // TFile *fFile = new TFile("/u/feisenhut/Documents/Uni/Masterarbeit/Sicherung/LegoTrainsOutputs/DQ_pp_MC_AOD/merged_403_LHC18h1_child3_404_LHC17h3_405_LHC17d1_LF_407_LHC18h1_child1+2_Dalitz_withPreFilter_withMasscut0-0.35_lowerSplitLevel/AnalysisResults.root");
  // TFile *fFile = new TFile("/u/feisenhut/Documents/Uni/Masterarbeit/Sicherung/LegoTrainsOutputs/DQ_pp_MC_AOD/merged_408_LHC18H1_409_LHC17h3_410_LHC17d1_LF_Dalitz_GammaGamma_noMassCut_withPrefilter/AnalysisResults.root");
  // TFile *fFile = new TFile(Form("/u/feisenhut/Documents/Uni/Masterarbeit/Sicherung/LegoTrainsOutputs/DQ_pp_MC_AOD/%s/AnalysisResults.root",TrainNumber.Data()));
  TFile *fFile = new TFile(Form("~/Dokumente/Uni/feisenhut/LegoTrainsOutputs/DQ_pp_MC_AOD/%s/AnalysisResults.root",TrainNumber.Data()));

                                                                                if(DoDebug) Printf("%d", __LINE__);


  // Names of the reconstructed primary cuts
  // TString names_Prim_Cuts=("noPID;JPID_sum1_pt75_sec_wImpParXY;JPID_sum2_pt75_sec_wImpParXYZ;JPID_sum3_pt75_sec_wImpParXYZ_NclsTPC;JPID_sum4_pt75_sec_wImpParXYZ_NclsTPC_TPCchi2Cl;JPID_sum5_pt75_sec_wImpParXYZ_NclsTPC_TPCchi2Cl_SharedCls;JPID_sum_pt75");
  // TString names_Prim_Cuts=("noPID;JPID_sum1_pt75_sec_OnlyImpParXY;JPID_sum2_pt75_sec_OnlyImpParZ;JPID_sum3_pt75_sec_OnlyNclsTPC;JPID_sum4_pt75_sec_OnlyTPCchi2Cl;JPID_sum5_pt75_sec_OnlySharedCls;JPID_sum6_pt75_sec_kV0;JPID_sum_pt75");
  // TString names_Prim_Cuts=("JPID_sum_pt75;JPID_sum2_pt75_sec_OnlyCosOpenAngle;JPID_sum3_pt75_sec_OnlyChi2NDF;JPID_sum4_pt75_sec_OnlyLegDist;JPID_sum5_pt75_sec_OnlyR;JPID_sum6_pt75_sec_OnlyPsiPair;JPID_sum7_pt75_sec_OnlyM;JPID_sum8_pt75_sec_OnlyArmPt;JPID_sum9_pt75_sec_OnlyArmAlpha;JPID_sum1_pt75_sec_kV0");
  // TString names_Prim_Cuts=("JPID_sum_pt75;JPID_sum1_pt75_sec_wCos;JPID_sum2_pt75_sec_wCosChi2;JPID_sum3_pt75_sec_wCosChi2LegDist;JPID_sum4_pt75_sec_wCosChi2LegDistR;JPID_sum5_pt75_sec_wCosChi2LegDistRPsiPair;JPID_sum6_pt75_sec_wCosChi2LegDistRPsiPairM;JPID_sum7_pt75_sec_wCosChi2LegDistRPsiPairMArmPt;JPID_sum8_pt75_sec_wCosChi2LegDistRPsiPairMArmPtAlpha");
  // TString names_Prim_Cuts=("JPID_sum_pt75;JPID_sum2_pt75_sec_OnlyCosOpenAngle;JPID_sum7_pt75_sec_OnlyM;JPID_sum8_pt75_sec_OnlyArmPt;JPID_sum9_pt75_sec_OnlyArmAlpha;JPID_sum1_pt75_sec_kV0");
  // TString names_Prim_Cuts=("JPID_sum_pt75;JPID_sum6_pt75_sec_wCosChi2LegDistRPsiPairM;JPID_sum7_pt75_sec_wCosChi2LegDistRPsiPairMArmPt;JPID_sum8_pt75_sec_wCosChi2LegDistRPsiPairMArmPtAlpha");
  TString names_Prim_Cuts=("JPID_sum_pt75;JPID_sum1_pt75_sec_kV0");





  // Names of the reconstructed secondary cuts
  // TString names_Sec_Cuts=("noPID;TESTCuts_wImpParXY;TESTCuts_wImpParXYZ;TESTCuts_wImpParXYZ_NclsTPC;TESTCuts_wImpParXYZ_NclsTPC_TPCchi2Cl;TESTCuts_wImpParXYZ_NclsTPC_TPCchi2Cl_SharedCls;JPID_sum_pt75_secondary");
  // TString names_Sec_Cuts=("noPID;TESTCutsOnlyImpParXY;TESTCutsOnlyImpParZ;TESTCutsOnlyNclsTPC;TESTCutsOnlyTPCchi2Cl;TESTCutsOnlySharedCls;TESTCutskV0;JPID_sum_pt75_secondary");

  // TString names_Sec_Cuts=("noPID;kV0OnlyCosOpenAngle;kV0OnlyChi2NDF;kV0OnlyLegDist;kV0OnlyR;kV0OnlyPsiPair;kV0OnlyM;kV0OnlyArmPt;kV0OnlyArmAlpha;kV0");
  // TString names_Sec_Cuts=("noPID;kV0wCosOpenAngle;kV0wChi2NDF;kV0wCosChi2LegDist;wCosChi2LegDistR;wCosChi2LegDistRPsiPair;kV0wCosChi2LegDistRPsiPairM;kV0wCosChi2LegDistRPsiPairMArmPt;kV0wCosChi2LegDistRPsiPairMArmPtAlpha");
  // TString names_Sec_Cuts=("noPID;kV0OnlyCosOpenAngle;kV0OnlyM;kV0OnlyArmPt;kV0OnlyArmAlpha;kV0");
  // TString names_Sec_Cuts=("noPID;kV0wPairM;kV0wPairMArmPt;kV0wPairMArmPtAlpha");
  // TString names_Sec_Cuts=("noPID;kV0");                 // use for old Root files e.g. merged_309
  TString names_Sec_Cuts=("noPID;pairkV0");
  TString names_Sec_Cuts_SecSecPairing=("noPID_V0_standard;track_V0_standard");

                                                                                if(DoDebug) Printf("%d", __LINE__);

  // Read in the different Lists from root file
  TList *EfficiencyHistList       = (TList*)fFile->Get("efficiency0");
  TList *SingleList               = (TList*)EfficiencyHistList->FindObject("SingleElectrons");
  TList *PairList                 = (TList*)EfficiencyHistList->FindObject("Pairs");
  // TList *FourPairListPrimSec             = (TList*)EfficiencyHistList->FindObject("4 el. Pairs");                                                // Before merged_352
  TList * FourPairListPrimSec     ; if(PrimSecPairing) FourPairListPrimSec        = (TList*)EfficiencyHistList->FindObject("4Pairs_PrimSec");       // Since merged_352
  TList * FourPairListSecSec      ; if(SecSecPairing)  FourPairListSecSec         = (TList*)EfficiencyHistList->FindObject("4Pairs_SecSec");        // Since merged_352
  TList *SupportList              ; SupportList                 = (TList*)EfficiencyHistList->FindObject("Support");
  TList *SinglePrimGenList        ; SinglePrimGenList           = (TList*)SingleList->FindObject("Generated_Primary");
  TList *SingleSecGenList         ; SingleSecGenList            = (TList*)SingleList->FindObject("Generated_Secondary");
  TList *SinglePrimGenSmearedList ; SinglePrimGenSmearedList    = (TList*)SingleList->FindObject("GeneratedSmeared_Primary");
  TList *SingleSecGenSmearedList  ; SingleSecGenSmearedList     = (TList*)SingleList->FindObject("GeneratedSmeared_Secondary");
  TList *PairPrimGenList          ; PairPrimGenList             = (TList*)PairList->FindObject("Generated_Primary");
  TList *PairPrimGenSmearedList   ; PairPrimGenSmearedList      = (TList*)PairList->FindObject("GeneratedSmeared_Primary");
  TList *PairSecGenList           ; PairSecGenList              = (TList*)PairList->FindObject("Generated_Secondary");
  TList *PairSecGenSmearedList    ; PairSecGenSmearedList       = (TList*)PairList->FindObject("GeneratedSmeared_Secondary");
  TList *FourPairGenList          ; if(PrimSecPairing) FourPairGenList        = (TList*)FourPairListPrimSec->FindObject("Generated");
  TList *FourPairGenSmearedList   ; if(PrimSecPairing) FourPairGenSmearedList = (TList*)FourPairListPrimSec->FindObject("GeneratedSmeared");
  TList *SupportSecondaryTrackCutsList_wo_Cuts  ; SupportSecondaryTrackCutsList_wo_Cuts = (TList*)SupportList->FindObject("SecondaryTrackCuts: noPID");
  TList *SupportSecondaryTrackCutsList_w_Cuts   ; SupportSecondaryTrackCutsList_w_Cuts  = (TList*)SupportList->FindObject("SecondaryTrackCuts: pairkV0");


  TObjArray *arrRecPrimaryFolderNames =  names_Prim_Cuts.Tokenize(";");
  TObjArray *arrRecSecondaryFolderNames =  names_Sec_Cuts.Tokenize(";");

  // TObjArray *arrRecSecSecPairingFolderNames =  names_Sec_Cuts_SecSecPairing.Tokenize(";");
  TObjArray *arrRecSecSecPairingFolderNames =  names_Sec_Cuts.Tokenize(";");
  const unsigned int nRecPrimCuts=arrRecPrimaryFolderNames->GetEntriesFast();
  const unsigned int nRecSecCuts=arrRecSecondaryFolderNames->GetEntriesFast();

  std::vector <TList*> SingleRecPrimList;
  std::vector <TList*> SingleRecSecList;
  std::vector <TList*> PairRecPrimList;
  std::vector <TList*> PairRecSecList;
  std::vector <TList*> FourPairRecList;
  std::vector <TList*> FourPairRecListSecSec;

                                                                                if(DoDebug) Printf("%d", __LINE__);

  // Create Single Reconstructed Lists  (Primary)
  for (unsigned int iRecCuts = 0; iRecCuts < nRecPrimCuts; iRecCuts++) {
    TString recPrimaryList_name = arrRecPrimaryFolderNames->At(iRecCuts)->GetName();
    recPrimaryList_name += "_Primary";
    TList* TempSingleRecPrimList = (TList*)SingleList->FindObject(recPrimaryList_name.Data());
    TempSingleRecPrimList->SetName(Form("Single%sList",recPrimaryList_name.Data()));
    TempSingleRecPrimList->SetOwner();
    SingleRecPrimList.push_back(TempSingleRecPrimList);
  }
                                                                                if(DoDebug) Printf("%d", __LINE__);
  // Create Single Reconstructed Lists  (Secondary)
  for (unsigned int iRecCuts = 0; iRecCuts < nRecSecCuts; iRecCuts++) {
    TString recSecondaryList_name = arrRecSecondaryFolderNames->At(iRecCuts)->GetName();
    recSecondaryList_name += "_Secondary";
    TList* TempSingleRecSecList = (TList*)SingleList->FindObject(recSecondaryList_name.Data());
    TempSingleRecSecList->SetName(Form("Single%sList",recSecondaryList_name.Data()));
    TempSingleRecSecList->SetOwner();
    SingleRecSecList.push_back(TempSingleRecSecList);
  }
                                                                                if(DoDebug) Printf("%d", __LINE__);
  // Create Pair Reconstructed Lists  (Primary)
  for (unsigned int iRecCuts = 0; iRecCuts < nRecPrimCuts; iRecCuts++) {
    TString recPrimaryList_name = arrRecPrimaryFolderNames->At(iRecCuts)->GetName();
    recPrimaryList_name += "_Primary";
    TList* TempPairRecPrimList = (TList*)PairList->FindObject(recPrimaryList_name.Data());
    TempPairRecPrimList->SetName(Form("Pair%sList",recPrimaryList_name.Data()));
    TempPairRecPrimList->SetOwner();
    PairRecPrimList.push_back(TempPairRecPrimList);
  }
                                                                                if(DoDebug) Printf("%d", __LINE__);
  // Create Pair Reconstructed Lists  (Secondary)
  for (unsigned int iRecCuts = 0; iRecCuts < nRecSecCuts; iRecCuts++) {
    TString recSecondaryList_name = arrRecSecondaryFolderNames->At(iRecCuts)->GetName();
    recSecondaryList_name += "_Secondary";
    TList* TempPairRecSecList = (TList*)PairList->FindObject(recSecondaryList_name.Data());
    TempPairRecSecList->SetName(Form("Pair%sList",recSecondaryList_name.Data()));
    TempPairRecSecList->SetOwner();
    PairRecSecList.push_back(TempPairRecSecList);
  }
                                                                                if(DoDebug) Printf("%d", __LINE__);
  // Create Four Reconstructed Lists
  if(PrimSecPairing) {
    for (unsigned int iRecCuts = 0; iRecCuts < nRecPrimCuts; iRecCuts++) {
      TString recList_name = arrRecPrimaryFolderNames->At(iRecCuts)->GetName();
      TList* TempFourPairRecList = (TList*)FourPairListPrimSec->FindObject(recList_name.Data());
      TempFourPairRecList->SetName(Form("FourPair%sList",recList_name.Data()));
      TempFourPairRecList->SetOwner();
      FourPairRecList.push_back(TempFourPairRecList);
    }
  }
                                                                                if(DoDebug) Printf("%d", __LINE__);
  // Create Four SecSec Reconstructed Lists
  if(SecSecPairing) {
    for (unsigned int iRecCuts = 0; iRecCuts < nRecSecCuts; iRecCuts++) {
      TString recListSecSec_name = arrRecSecSecPairingFolderNames->At(iRecCuts)->GetName();
      TList* TempFourPairRecListSecSec = (TList*)FourPairListSecSec->FindObject(recListSecSec_name.Data());
      TempFourPairRecListSecSec->SetName(Form("FourPair%sListSecSec",recListSecSec_name.Data()));
      TempFourPairRecListSecSec->SetOwner();
      FourPairRecListSecSec.push_back(TempFourPairRecListSecSec);
    }
  }

                                                                                if(DoDebug) Printf("%d", __LINE__);

  // Names of single primary electron histograms
  TString VecGenSmearedPrimSingleHistoNames = "Ngen_Neg_eleFinalState";
  TString VecRecPrimSingleHistoNames = "Nrec_Neg_eleFinalState";

  // Names of single secondary electron histograms
  TString VecGenSmearedSecSingleHistoNames = "Ngen_Neg_eleSecondaryFromPhoton";
  TString VecRecSecSingleHistoNames = "Nrec_Neg_eleSecondaryFromPhoton";

  // Names of secondary pair histograms
  // TString VecRecSecPairHistoNames = "Nrec_pair_sameMother_photon_secondary";
  TString VecRecSecPairHistoNames = "Nrec_pair_NotSameMother_secondary";

  // Names of two pair histograms
  // TString VecPairHistoNames = "Ngen_pair_sameMother_photon_finalstate;Ngen_pair_sameMother_photon_secondary;Ngen_pair_sameMother_photon_secondaryfromMaterial;Ngen_pair_sameMother_photon_secondaryfromWD;Ngen_pair_UndefinedMother_photon_secondary;Ngen_pair_DifferentMother_photon_secondary;Ngen_pair_sameMother_photon_secondary_eta;Ngen_sameMother_eta;Ngen_pair_conversion;Ngen_pair_random;Ngen_pair_sameMother_pion";
  // \/\/\/\/ use for old Root files e.g. merged_309 \/\/\/\/
  // TString VecPrimPairHistoNames    = "Ngen_pair_sameMother_finalstate;Ngen_pair_DifferentMother_finalstate;Ngen_pair_UndefinedMother_finalstate;Ngen_pair_sameMother_photon_finalstate;Ngen_pair_sameMother_eta_finalstate;Ngen_pair_DifferentMother_eta_finalstate;Ngen_pair_UndefinedMother_eta_finalstate";
  // TString VecPrimPairHistoNames    = "Ngen_pair_sameMother_finalstate;Ngen_pair_DifferentMother_finalstate;Ngen_pair_UndefinedMother_finalstate;Ngen_pair_sameMother_photon_finalstate;Ngen_pair_sameMother_eta_finalstate;Ngen_pair_DifferentMother_eta_finalstate;Ngen_pair_UndefinedMother_eta_finalstate;Ngen_pair_sameMother_pion_finalstate;Ngen_pair_DifferentMother_pion_finalstate;Ngen_pair_UndefinedMother_pion_finalstate";
  // TString VecSecPairHistoNames     = "Ngen_pair_sameMother_photon_secondary;Ngen_pair_sameMother_photon_secondary_pion;Ngen_pair_sameMother_photon_secondary_eta;Ngen_pair_conversion_secondary;Ngen_pair_random_secondary;Ngen_pair_NotSameMother_secondary";
  // in future root files (older than merged_335)
  TString VecPrimPairHistoNames    = "Ngen_anyPair_anyPart_UndefinedMother_finalstate;Ngen_elePair_sameMother_finalstate;Ngen_elePair_DifferentMother_finalstate;Ngen_elePair_UndefinedMother_finalstate;Ngen_elePair_sameMother_photon_finalstate;Ngen_elePair_sameMother_eta_finalstate;Ngen_elePair_DifferentMother_eta_finalstate;Ngen_elePair_UndefinedMother_eta_finalstate;Ngen_elePair_sameMother_pion_finalstate;Ngen_elePair_DifferentMother_pion_finalstate;Ngen_elePair_UndefinedMother_pion_finalstate";
  TString VecSecPairHistoNames     = "Ngen_anyPair_anyPart_random_secondary;Ngen_elePair_sameMother_photon_secondary;Ngen_elePair_sameMother_photon_secondary_pion;Ngen_elePair_sameMother_photon_secondary_eta;Ngen_elePair_conversion_secondary;Ngen_elePair_random_secondary;Ngen_elePair_NotSameMother_secondary";

  // Names of reconstructed two pair histograms
    // \/\/\/\/ use for old Root files e.g. merged_309 \/\/\/\/
  // TString VecPrimRecPairHistoNames = "Nrec_pair_sameMother_finalstate;Nrec_pair_DifferentMother_finalstate;Nrec_pair_UndefinedMother_finalstate;Nrec_pair_sameMother_photon_finalstate;Nrec_pair_sameMother_eta_finalstate;Nrec_pair_DifferentMother_eta_finalstate;Nrec_pair_UndefinedMother_eta_finalstate";
  // TString VecPrimRecPairHistoNames = "Nrec_pair_sameMother_finalstate;Nrec_pair_DifferentMother_finalstate;Nrec_pair_UndefinedMother_finalstate;Nrec_pair_sameMother_photon_finalstate;Nrec_pair_sameMother_eta_finalstate;Nrec_pair_DifferentMother_eta_finalstate;Nrec_pair_UndefinedMother_eta_finalstate;Nrec_pair_sameMother_pion_finalstate;Nrec_pair_DifferentMother_pion_finalstate;Nrec_pair_UndefinedMother_pion_finalstate";
  // TString VecSecRecPairHistoNames  = "Nrec_pair_sameMother_photon_secondary;Nrec_pair_sameMother_photon_secondary_pion;Nrec_pair_sameMother_photon_secondary_eta;Nrec_pair_conversion_secondary;Nrec_pair_random_secondary;Nrec_pair_NotSameMother_secondary";
  // in future root files (older than merged_335)
  TString VecPrimRecPairHistoNames = "Nrec_anyPair_anyPart_UndefinedMother_finalstate;Nrec_elePair_sameMother_finalstate;Nrec_elePair_DifferentMother_finalstate;Nrec_elePair_UndefinedMother_finalstate;Nrec_elePair_sameMother_photon_finalstate;Nrec_elePair_sameMother_eta_finalstate;Nrec_elePair_DifferentMother_eta_finalstate;Nrec_elePair_UndefinedMother_eta_finalstate;Nrec_elePair_sameMother_pion_finalstate;Nrec_elePair_DifferentMother_pion_finalstate;Nrec_elePair_UndefinedMother_pion_finalstate";
  TString VecSecRecPairHistoNames  = "Nrec_anyPair_anyPart_random_secondary;Nrec_elePair_sameMother_photon_secondary;Nrec_elePair_sameMother_photon_secondary_pion;Nrec_elePair_sameMother_photon_secondary_eta;Nrec_elePair_conversion_secondary;Nrec_elePair_random_secondary;Nrec_elePair_NotSameMother_secondary";

  // Names of four pair histograms
  // TString VecGenFourPairHistoNames = "Ngen_ULSFourElePair1_FinalState";
  // TString VecGenFourPairHistoNames = "Ngen_FourElePair1_FinalState";
  // TString VecGenFourPairHistoNames = "Ngen_FourElePair1_FinalState;Ngen_FourElePair1_FinalState_Dalitz;Ngen_FourAnyPartPair1_FinalState";
  TString VecGenFourPairHistoNames = "Ngen_FourElePair1_FinalState;Ngen_FourElePair1_sameEta_Dalitz;Ngen_FourAnyPartPair1_FinalState_UndefinedMother;Ngen_FourElePair1_samePion_Dalitz";

  // Names of reconstructed four pair histograms
  // TString VecRecFourPairHistoNamesPrimSec = "Nrec_ULSFourElePair1_FinalState";
  // TString VecRecFourPairHistoNamesPrimSec = "Nrec_FourElePair1_FinalState";
  // TString VecRecFourPairHistoNamesPrimSec = "Nrec_FourElePair1_FinalState;Nrec_FourElePair1_FinalState_Dalitz;Nrec_FourAnyPartPair1_FinalState";
  TString VecRecFourPairHistoNamesPrimSec = "Nrec_FourElePair1_FinalState;Nrec_FourElePair1_sameEta_Dalitz;Nrec_FourAnyPartPair1_FinalState_UndefinedMother;Nrec_FourElePair1_samePion_Dalitz";
  TString VecRecFourPairHistoNamesSecSec = "Nrec_FourElePair1_Secondary_samePhoton;Nrec_FourElePair1_Secondary_samePhoton_sameEta;Nrec_FourAnyPartPair1_Secondary_UndefinedMother;Nrec_FourElePair1_Secondary_samePhoton_samePion";

  // Name of TH2 histograms where we do a mass projection and plot same and different mother pairs in one histogram
  // First sameMother_Name , then differntMother_Name  (fistSignal: sameMother + differnetMother, secondSignal:  ...)  DifferentMother Hisots shouldn't be from a defined MotherSignal
  // primary Hist names
  TString VecPrimGenSmearedMassCutHistoNames = "Ngen_pair_sameMother_finalstate;Ngen_pair_DifferentMother_finalstate;Ngen_pair_sameMother_pion_finalstate;Ngen_pair_DifferentMother_finalstate;Ngen_pair_sameMother_eta_finalstate;Ngen_pair_DifferentMother_finalstate";
  TString VecPrimRecMassCutHistoNames        = "Nrec_pair_sameMother_finalstate;Nrec_pair_DifferentMother_finalstate;Nrec_pair_sameMother_pion_finalstate;Nrec_pair_DifferentMother_finalstate;Nrec_pair_sameMother_eta_finalstate;Nrec_pair_DifferentMother_finalstate";
  // secondary Hist names
  TString VecSecGenSmearedMassCutHistoNames = "Ngen_pair_sameMother_photon_secondary;Ngen_pair_DifferentMother_photon_secondary";
  TString VecSecRecMassCutHistoNames        = "Nrec_pair_sameMother_photon_secondary;Nrec_pair_DifferentMother_photon_secondary";


  // Titles of mass cut histograms  + dummy "; ;"
  TString VecPrimMassCutHistoTitles = "Pair from same and different Mother, Primary; Pair from same Pion and different Mother, Primary; ;Pair from same Eta and different Mother, Primary; ";
  TString VecSecMassCutHistoTitles = "Pair from same Mother and different Mother, Secondary; ";

  // List names of MCSignalsList that should be plotted in PlotPIDComponents
  TString ListNamesSupportMCSignalsSecondary = "SecondaryTrackCuts: anyPair_anyPart_random_secondary";

  // Names of Prefilter histograms before and after Prefilter
  TString VecPrefilterHistNames = "HistBeforePrimSecPrefilter;HistBeforeSecSecPrefilter;HistAfterPrimSecPrefilter;HistAfterSecSecPrefilter";


  // Output Check what Settings are selected
  std::cout << "ExtPairSig " << ExtPairSig <<" ExtFourPairSig " << ExtFourPairSig << std::endl;
  std::cout << "ExtGen " << ExtGen <<" ExtGenSmeared "<< ExtGenSmeared << std::endl;
  std::cout << "DoPtProj " << DoPtProj <<" DoMassProj "<< DoMassProj << std::endl;

  // Get Number of Events
  TH1D* hNEvents =  (TH1D*) EfficiencyHistList->FindObject("events");
  nEvents = hNEvents->GetEntries();
                                                                                if(DoDebug) std::cout << "Number of Events = " << nEvents << std::endl;



                                                                                if(DoDebug) Printf("%d", __LINE__);


  if (ExtPairSig == kTRUE) {
    DoExtractPairSignal(VecPrimPairHistoNames, VecPrimRecPairHistoNames, PairPrimGenList, PairPrimGenSmearedList, PairRecPrimList[1], ExtGen, ExtGenSmeared, ExtRec, vec_rebin_pt, vec_rebin_mass, vec_proj_bin_pt, vec_proj_bin_mass, nEvents);
    DoExtractPairSignal(VecSecPairHistoNames,  VecSecRecPairHistoNames,  PairSecGenList,  PairSecGenSmearedList,  PairRecSecList[1],  ExtGen, ExtGenSmeared, ExtRec, vec_rebin_pt, vec_rebin_mass, vec_proj_bin_pt, vec_proj_bin_mass, nEvents);
  }

  if (PrimarySoBRatio == kTRUE) {
    PlotPrimarySoBRatio(VecPrimRecPairHistoNames, PairRecPrimList, vec_rebin_pt, vec_rebin_mass, "Pair");
  }


                                                                                if(DoDebug) Printf("%d", __LINE__);

  if (ExtFourPairSig == kTRUE) {
    std::vector <TH2D> VecFourGenPairHistos;
    std::vector <TH2D> VecFourGenSmearedPairHistos;
    std::vector <TH2D> VecFourRecPairHistos;
    // generated and Generated Smeared part
    TObjArray* arrGenFourPairHistoNames=VecGenFourPairHistoNames.Tokenize(";");
    TObjArray* arrRecFourPairHistoNames;
    if (PrimSecPairing == kTRUE)  arrRecFourPairHistoNames=VecRecFourPairHistoNamesPrimSec.Tokenize(";");
    if (SecSecPairing == kTRUE)   arrRecFourPairHistoNames=VecRecFourPairHistoNamesSecSec.Tokenize(";");
    const Int_t nGenFourPairTH2hist=arrGenFourPairHistoNames->GetEntriesFast();
                                                                                std::cout << __LINE__ << " Cout " << nGenFourPairTH2hist << std::endl;
    for (Int_t i = 0; i < nGenFourPairTH2hist; i++) {
      if (ExtGen == kTRUE|| ExtGenSmeared == kTRUE) {
        TString temp1 = arrGenFourPairHistoNames->At(i)->GetName();
        TH2D hTempGen4PairHist         = *(dynamic_cast<TH2D*>(FourPairGenList->FindObject(temp1)));
        TH2D hTempGenSmear4PairHist    = *(dynamic_cast<TH2D*>(FourPairGenSmearedList->FindObject(temp1)));
        hTempGen4PairHist     .Sumw2();
        hTempGenSmear4PairHist.Sumw2();
        // if (DoRebin == kTRUE) {
        //   if (!Rebin2DHistogram(hTempGen4PairHist     , vec_rebin_mass, vec_rebin_pt)) continue; // Rebin2DHistogram sould only fail if new binning has less than 2 bins
        //   if (!Rebin2DHistogram(hTempGenSmear4PairHist, vec_rebin_mass, vec_rebin_pt)) continue; // Rebin2DHistogram sould only fail if new binning has less than 2 bins
        // }
        VecFourGenPairHistos.push_back(hTempGen4PairHist);
        VecFourGenSmearedPairHistos.push_back(hTempGenSmear4PairHist);
      }
                                                                                if(DoDebug) Printf("%d", __LINE__);
      if (ExtRec == kTRUE) {
        TString temp2 = arrRecFourPairHistoNames->At(i)->GetName();
        TH2D hTempRec4PairHist;
        if (PrimSecPairing == kTRUE) hTempRec4PairHist = *(dynamic_cast<TH2D*>(FourPairRecList[1]->FindObject(temp2))); // Select histograms of second CutSettings
        if (SecSecPairing == kTRUE)  hTempRec4PairHist = *(dynamic_cast<TH2D*>(FourPairRecListSecSec[1]->FindObject(temp2))); // Select histograms of second CutSettings
        // hTempRec4PairHist     .Sumw2();
        // if (DoRebin == kTRUE) {
        //   if (!Rebin2DHistogram(hTempRec4PairHist     , vec_rebin_mass, vec_rebin_pt)) continue; // Rebin2DHistogram sould only fail if new binning has less than 2 bins
        // }
        VecFourRecPairHistos.push_back(hTempRec4PairHist);
      }
    }
                                                                                if(DoDebug) Printf("%d", __LINE__);
    if (ExtGen        == kTRUE) {DrawProjection(VecFourGenPairHistos       , vec_proj_bin_pt, vec_proj_bin_mass, vec_rebin_pt, vec_rebin_mass, ExtFourPairPionSig, ExtFourPairEtaSig, "FourPair/Generated"        , nEvents);}
    if (ExtGenSmeared == kTRUE) {DrawProjection(VecFourGenSmearedPairHistos, vec_proj_bin_pt, vec_proj_bin_mass, vec_rebin_pt, vec_rebin_mass, ExtFourPairPionSig, ExtFourPairEtaSig, "FourPair/GeneratedSmeared" , nEvents);}
                                                                                if(DoDebug) Printf("%d", __LINE__);
    if (ExtRec == kTRUE && PrimSecPairing) {DrawProjection(VecFourRecPairHistos       , vec_proj_bin_pt, vec_proj_bin_mass, vec_rebin_pt, vec_rebin_mass, ExtFourPairPionSig, ExtFourPairEtaSig, Form("FourPair/%s",FourPairRecList[1]->GetName())       , nEvents);}
    if (ExtRec == kTRUE && SecSecPairing)  {DrawProjection(VecFourRecPairHistos       , vec_proj_bin_pt, vec_proj_bin_mass, vec_rebin_pt, vec_rebin_mass, ExtFourPairPionSig, ExtFourPairEtaSig, Form("FourPair/%s",FourPairRecListSecSec[1]->GetName()) , nEvents);}


  }

                                                                                if(DoDebug) Printf("%d", __LINE__);

  if (ExtMassCut == kTRUE) {
    if(ExtGenSmeared == kTRUE){
      std::cout <<"Begin merge Prim GenSmeared Histograms" << std::endl;
      MergeMassCutHistos(PairPrimGenSmearedList,VecPrimGenSmearedMassCutHistoNames,VecPrimMassCutHistoTitles,"Pair/GeneratedSmeared", MassCutPrimary,!DoCutEff);
      std::cout <<"Begin merge Sec GenSmeared Histograms" << std::endl;
      MergeMassCutHistos(PairSecGenSmearedList,VecSecGenSmearedMassCutHistoNames,VecSecMassCutHistoTitles,"Pair/GeneratedSmeared", MassCutSecondary,DoCutEff);
    }
    if(ExtRec == kTRUE){
      std::cout <<"Begin merge Prim Rec Histograms" << std::endl;
      MergeMassCutHistos(PairRecPrimList[6],VecPrimRecMassCutHistoNames,VecPrimMassCutHistoTitles,"Pair/JPID_sum_pt75_Primary", MassCutPrimary,!DoCutEff);
      std::cout <<"Begin merge Sec Rec Histograms" << std::endl;
      MergeMassCutHistos(PairRecSecList[6],VecSecRecMassCutHistoNames,VecSecMassCutHistoTitles,"Pair/JPID_sum_pt75_secondary_Secondary", MassCutSecondary,DoCutEff);
    }
  }

                                                                                if(DoDebug) Printf("%d", __LINE__);

  // Comparison between secondary GeneratedSmeared and Reconstructed in Pt, Eta and Phi on single electron level
  if(DoCompGenSmearRec){
    // CompareGenSmearRec(VecGenSmearedSecSingleHistoNames,VecRecSecSingleHistoNames,SingleSecGenSmearedList,SingleRecSecList[0]);
  }

                                                                                if(DoDebug) Printf("%d", __LINE__);

  if(DoLossBySecondaryCuts){
                                                                                if(DoDebug) Printf("before LossBySecondaryCuts\n");

 // LossBySecondaryCuts("SingleElectrons",VecRecSecSingleHistoNames, SingleRecSecList, nRecPrimCuts, "Single", "_SecondaryList");
    LossBySecondaryCuts("Pair"           ,VecRecSecPairHistoNames  , PairRecSecList  , nRecSecCuts , "Pair"  , "_SecondaryList");
  }

  if (DoPlotPIDComponents) {
    PlotPIDComponents (SupportSecondaryTrackCutsList_w_Cuts, ListNamesSupportMCSignalsSecondary);
  }

  if(ExtDielectronPairSpectra) {
                                                                                if(DoDebug) Printf("%d", __LINE__);
    TCanvas *cSignalCanvas           = new TCanvas("cSignalCanvas","",800,800);
    SetCanvasStandardSettings(cSignalCanvas);
    cSignalCanvas->SetLogz(1);

    TString histName;
    TString DocumentPath;
    TObjArray* arrPrimSecHistNames;
    TObjArray* arrSecSecHistNames;
    TH2D* dielectronPairHist;
    arrPrimSecHistNames = VecRecFourPairHistoNamesPrimSec.Tokenize(";");
    arrSecSecHistNames  = VecRecFourPairHistoNamesSecSec.Tokenize(";");
    if(PrimSecPairing){
      histName = arrPrimSecHistNames->At(2)->GetName();
      DocumentPath = Form("Plots/%s/%s/DielectronPairSpectra/%s/",TrainNumber.Data(),"PrimSecPairing",FourPairRecList.at(1)->GetName());
      dielectronPairHist = (TH2D*)FourPairRecList.at(1)->FindObject(histName);
    }
    if(SecSecPairing){
      histName = arrSecSecHistNames->At(2)->GetName();
      DocumentPath = Form("Plots/%s/%s/DielectronPairSpectra/%s/",TrainNumber.Data(),"SecSecPairing",FourPairRecListSecSec.at(1)->GetName());
      dielectronPairHist = (TH2D*)FourPairRecListSecSec.at(1)->FindObject(histName);
    }
    gSystem->Exec(Form("mkdir -p %s",DocumentPath.Data()));
    cSignalCanvas->cd();
    gPad->SetLeftMargin(0.1);
    gPad->SetRightMargin(0.13);
    dielectronPairHist->Draw("ColZ");
    cSignalCanvas -> SaveAs(Form("%s%s.pdf",DocumentPath.Data(),histName.Data()));
  }

  if(ExtPreFilterSpectra) {
    TCanvas *cSignalCanvas           = new TCanvas("cSignalCanvas","",800,800);
    SetCanvasStandardSettings(cSignalCanvas);
    cSignalCanvas->SetLogz(1);

    TString DocumentPath = Form("Plots/%s/PrefilterHistos/",TrainNumber.Data());
    gSystem->Exec(Form("mkdir -p %s",DocumentPath.Data()));

    TH2D* hPrefilterHist;
    TObjArray* arrPrefilterNames = VecPrefilterHistNames.Tokenize(";");
    const Int_t nPrefilterHists = arrPrefilterNames->GetEntriesFast();
    for (Int_t i = 0; i < nPrefilterHists; i++) {
      hPrefilterHist = (TH2D*)FourPairListPrimSec->FindObject(arrPrefilterNames->At(i)->GetName());
      cSignalCanvas->cd();
      hPrefilterHist->Draw("ColZ");
      cSignalCanvas -> SaveAs(Form("%s%s.pdf",DocumentPath.Data(),arrPrefilterNames->At(i)->GetName()));
    }

  }


}


void DoExtractPairSignal(TString VecPairHistoNames, TString VecRecPairHistoNames, TList* PairGenList ,TList* PairGenSmearedList, TList* PairRecList, Bool_t ExtGen, Bool_t ExtGenSmeared, Bool_t ExtRec, std::vector<Int_t> vec_rebin_pt, std::vector<Int_t> vec_rebin_mass, std::vector<Double_t> vec_proj_bin_pt, std::vector<Double_t> vec_proj_bin_mass, Double_t nEvents){
  std::vector <TH2D> VecGenPairHistos;
  std::vector <TH2D> VecGenSmearedPairHistos;
  std::vector <TH2D> VecRecPairHistos;
  TObjArray* arrPairHistoNames=VecPairHistoNames.Tokenize(";");
  TObjArray* arrRecPairHistoNames=VecRecPairHistoNames.Tokenize(";");
  const Int_t nPairTH2hist=arrPairHistoNames->GetEntriesFast();
  std::cout << nPairTH2hist << std::endl;
  for (Int_t i = 0; i < nPairTH2hist; i++) {
    TString temp1 = arrPairHistoNames->At(i)->GetName();
    TString temp2 = arrRecPairHistoNames->At(i)->GetName();
    std::cout << temp1 << " " << temp2 << std::endl;

    TH2D hTempGenPairHist          = *(dynamic_cast<TH2D*>(PairGenList       ->FindObject(temp1)));
    TH2D hTempGenSmearPairHist     = *(dynamic_cast<TH2D*>(PairGenSmearedList->FindObject(temp1)));
    TH2D hTempRecPairHist          = *(dynamic_cast<TH2D*>(PairRecList       ->FindObject(temp2)));
    hTempGenPairHist     .Sumw2();
    hTempGenSmearPairHist.Sumw2();
    hTempRecPairHist     .Sumw2();

    // if (DoRebin == kTRUE){
    //   if(!Rebin2DHistogram(hTempGenPairHist      , vec_rebin_mass, vec_rebin_pt)) continue;  // Rebin2DHistogram sould only fail if new binning has less than 2 bins
    //   if(!Rebin2DHistogram(hTempGenSmearPairHist , vec_rebin_mass, vec_rebin_pt)) continue;  // Rebin2DHistogram sould only fail if new binning has less than 2 bins
    //   if(!Rebin2DHistogram(hTempRecPairHist      , vec_rebin_mass, vec_rebin_pt)) continue;  // Rebin2DHistogram sould only fail if new binning has less than 2 bins
    // }

    VecGenPairHistos.push_back(hTempGenPairHist);               // Fill vectors with selected histograms
    VecGenSmearedPairHistos.push_back(hTempGenSmearPairHist);   // Fill vectors with selected histograms
    VecRecPairHistos.push_back(hTempRecPairHist);               // Fill vectors with selected histograms
  }

  if (ExtGen        == kTRUE) {DrawProjection(VecGenPairHistos,        vec_proj_bin_pt, vec_proj_bin_mass, vec_rebin_pt, vec_rebin_mass, ExtFourPairPionSig, ExtFourPairEtaSig, "Pair/Generated"                      , nEvents);}
  if (ExtGenSmeared == kTRUE) {DrawProjection(VecGenSmearedPairHistos, vec_proj_bin_pt, vec_proj_bin_mass, vec_rebin_pt, vec_rebin_mass, ExtFourPairPionSig, ExtFourPairEtaSig, "Pair/GeneratedSmeared"               , nEvents);}
  if (ExtRec        == kTRUE) {DrawProjection(VecRecPairHistos,        vec_proj_bin_pt, vec_proj_bin_mass, vec_rebin_pt, vec_rebin_mass, ExtFourPairPionSig, ExtFourPairEtaSig, Form("Pair/%s",PairRecList->GetName()), nEvents);}
}


//  function to draw projection of x and y axis (mass and pt projections)
void DrawProjection(std::vector<TH2D> VecHistos, std::vector<Double_t> vec_proj_bin_pt, std::vector<Double_t> vec_proj_bin_mass, std::vector<Int_t> vec_rebin_pt, std::vector<Int_t> vec_rebin_mass, Bool_t ExtFourPairPionSig, Bool_t ExtFourPairEtaSig, TString PairCase, Double_t nEvents) {
                                                                                if(DoDebug) Printf("%d", __LINE__);
  TCanvas *cSignalCanvas           = new TCanvas("cSignalCanvas","",800,800);
  TCanvas *cProjPtPionCanvas       = new TCanvas("cProjPtPionCanvas","",1600,1600);
  TCanvas *cProjPtEtaCanvas        = new TCanvas("cProjPtEtaCanvas" ,"",1600,1600);
  TCanvas *cSignalPionCanvas       = new TCanvas("cSignalPionCanvas","",1600,1600);
  TCanvas *cSignalEtaCanvas        = new TCanvas("cSignalEtaCanvas" ,"",1600,1600);
  TCanvas *cSignificancePionCanvas = new TCanvas("cSignificancePionCanvas" ,"",1600,1600);
  TCanvas *cSignificanceEtaCanvas  = new TCanvas("cSignificanceEtaCanvas" ,"",1600,1600);
  TCanvas *cSoBPionCanvas          = new TCanvas("cSigOverBackPionCanvas" ,"",1600,1600);
  TCanvas *cSoBEtaCanvas           = new TCanvas("cSigOverBackEtaCanvas" ,"",1600,1600);
  // TCanvas *cSignalCanvasLog        = new TCanvas("cSignalCanvasLog","",800,800);
  SetCanvasStandardSettings(cSignalCanvas);
  SetCanvasStandardSettings(cProjPtPionCanvas);
  SetCanvasStandardSettings(cProjPtEtaCanvas);
  SetCanvasStandardSettings(cSignalPionCanvas);
  SetCanvasStandardSettings(cSignalEtaCanvas);
  SetCanvasStandardSettings(cSoBPionCanvas);
  SetCanvasStandardSettings(cSoBEtaCanvas);
  SetCanvasStandardSettings(cSignificancePionCanvas);
  SetCanvasStandardSettings(cSignificanceEtaCanvas);
  // SetCanvasStandardSettings(cSignalCanvasLog);
  cProjPtPionCanvas->DivideSquare(vec_proj_bin_pt.size()-1);
  cProjPtEtaCanvas->DivideSquare(vec_proj_bin_pt.size()-1);
  cSignalPionCanvas->DivideSquare(vec_proj_bin_pt.size()-1);
  cSignalEtaCanvas->DivideSquare(vec_proj_bin_pt.size()-1);
  cSoBPionCanvas->DivideSquare(vec_proj_bin_pt.size()-1);
  cSoBEtaCanvas->DivideSquare(vec_proj_bin_pt.size()-1);
  cSignificancePionCanvas->DivideSquare(vec_proj_bin_pt.size()-1);
  cSignificanceEtaCanvas->DivideSquare(vec_proj_bin_pt.size()-1);
  // cSignalCanvasLog->SetLogy();
  cSignalCanvas->SetLogy();
                                                                                if(DoDebug) Printf("%d", __LINE__);
  TH1D* projX;
  TH1D* projY;

  Double_t startTextX         = .48;
  Double_t startTextY         = .84;
  Double_t textHeight         = 0.01;
  Double_t textSize           = 25;

  Double_t lowerPionMassRange = 0.0;
  Double_t upperPionMassRange = 0.3;
  Double_t lowerEtaMassRange = 0.3;
  Double_t upperEtaMassRange = 0.8;
                                                                                if(DoDebug) Printf("%d", __LINE__);
  TString PairingType;
  if (PrimSecPairing == kTRUE)  PairingType = "PrimSecPairing";
  if (SecSecPairing  == kTRUE)  PairingType = "SecSecPairing";

  TString path_Pol_or_Template;
  if (DoTemplateFit == kTRUE)  path_Pol_or_Template = "TemplateFit";
  if (DoTemplateFit == kFALSE) path_Pol_or_Template = "PolFit";

  for (size_t iHist = 0; iHist < VecHistos.size(); iHist++) {              // loop over all Histograms
    TString DocumentPathProjX   = Form("Plots/%s/%s/%s/%s/%s/Mass_Projectoion/",TrainNumber.Data(),PairingType.Data(),PairCase.Data(),path_Pol_or_Template.Data(),VecHistos.at(iHist).GetName());  // define path to mass projection
    TString DocumentPathProjY   = Form("Plots/%s/%s/%s/%s/%s/Pt_Projectoion/"  ,TrainNumber.Data(),PairingType.Data(),PairCase.Data(),path_Pol_or_Template.Data(),VecHistos.at(iHist).GetName());    // define path to pt projections
    TString DocumentPathSignals = Form("Plots/%s/%s/%s/%s/%s/Signals/"         ,TrainNumber.Data(),PairingType.Data(),PairCase.Data(),path_Pol_or_Template.Data(),VecHistos.at(iHist).GetName());  // define path to mass projection
    gSystem->Exec(Form("mkdir -p %s",DocumentPathProjY.Data()));
    gSystem->Exec(Form("mkdir -p %s",DocumentPathProjX.Data()));
    gSystem->Exec(Form("mkdir -p %s",DocumentPathSignals.Data()));

    // Define hist to draw S/B in pt-slices
    std::vector<Double_t> Vec_SigBack_Pion;
    std::vector<Double_t> Vec_SigBack_Eta;
                                                                                if(DoDebug) Printf("%d", __LINE__);

    // TH1D* hSignalOverBackgroundPtPion = new TH1D("Signal over Background in pt slices for Pions", "Signal over Backgroun in pt slices for Pions ;Pt [GeV] ; Ratio S/B", vec_proj_bin_pt.size()-1, vec_proj_bin_pt.data());
    // TH1D* hSignalOverBackgroundPtEta = new TH1D("Signal over Background in pt slices for Eta", "Signal over Backgroun in pt slices for Eta; Pt [GeV] ; Ratio S/B", vec_proj_bin_pt.size()-1, vec_proj_bin_pt.data());
    // SetHistoStandardSettingsRed(hSignalOverBackgroundPtPion);
    // SetHistoStandardSettingsBlue(hSignalOverBackgroundPtEta);
    for (size_t j = 0; j < vec_proj_bin_pt.size()-1;  j++) {      // loop over all pt projection intervalls
      for (size_t k = 0; k < vec_proj_bin_mass.size()-1; k++) {     // loop over all mass projection intervalls
        TLatex *pT_Intervall    = new TLatex(startTextX, startTextY    , Form("%g < #font[12]{p}_{T,eeee} < %g  #frac{GeV}{#font[12]{c}}",vec_proj_bin_pt.at(j) ,vec_proj_bin_pt.at(j+1)));
        TLatex *mass_Intervall  = new TLatex(startTextX, startTextY-.05, Form("#font[12]{m}_{eeee}  Intervall = %g - %g  #frac{GeV}{#font[12]{c^{2}}}",vec_proj_bin_mass.at(k),vec_proj_bin_mass.at(k+1)));
        SetTextSettings(pT_Intervall,textSize);
        SetTextSettings(mass_Intervall,textSize);
        Int_t startbinX = VecHistos.at(iHist).GetXaxis()->FindBin(vec_proj_bin_mass[k]);    // select start bin of mass projection
        Int_t endbinX   = VecHistos.at(iHist).GetXaxis()->FindBin(vec_proj_bin_mass[k+1]);  // select end bin of mass projection
        Int_t startbinY = VecHistos.at(iHist).GetYaxis()->FindBin(vec_proj_bin_pt[j]);      // select start bin of pt projection
        Int_t endbinY   = VecHistos.at(iHist).GetYaxis()->FindBin(vec_proj_bin_pt[j+1]);    // select end bin of pt projection
        projX = (TH1D*) VecHistos.at(iHist).ProjectionX(Form("%s_ProjMass%g:%g_pt%g:%g",VecHistos.at(iHist).GetName(),vec_proj_bin_mass.at(k),vec_proj_bin_mass.at(k+1),vec_proj_bin_pt.at(j),vec_proj_bin_pt.at(j+1)),startbinY,endbinY)->Clone();  // Mass projection Histogram
        projY = (TH1D*) VecHistos.at(iHist).ProjectionY(Form("%s_ProjPt%g:%g_mass%g:%g",VecHistos.at(iHist).GetName(),vec_proj_bin_pt.at(j),vec_proj_bin_pt.at(j+1),vec_proj_bin_mass.at(k),vec_proj_bin_mass.at(k+1)),startbinX,endbinX)->Clone();  // Pt projection Histogram
        projX->RebinX(vec_rebin_mass.at(j));
        // projX->RebinY();
        // projX->Sumw2();
        // projY->Sumw2();

        // Normalise Histograms to Number of Events
        TF1* fTF1dummy = new TF1("dummy","1",-100,100);
        projX->Divide(fTF1dummy,nEvents);
        projY->Divide(fTF1dummy,nEvents);

        SetHistoStandardSettings(projX);
        SetHistoStandardSettings(projY);

        // projX->Scale(1,"width");
        // projY->Scale(1,"width");

                                                                                if(DoDebug) Printf("%d", __LINE__);

        //////////////////////////////////////////////////////////////////////////////
        // MC, TRUE SIGNAL AND BACKGROUNG: USED TO FIT DATA ; data = a*SignalHist + b*BackgroundHist
        //////////////////////////////////////////////////////////////////////////////
        TH1D* fMCprojSB_pion;
        TH1D* fMCprojSB_eta;
        // fMCprojSB_eta = (TH1D*) VecHistos.at(0).ProjectionX(Form("%s_ProjMass%g:%g_pt%g:%g",VecHistos.at(0).GetName(),vec_proj_bin_mass.at(k),vec_proj_bin_mass.at(k+1),vec_proj_bin_pt.at(j),vec_proj_bin_pt.at(j+1)) ,startbinY,endbinY)->Clone();  // Mass projection Histogram
        fMCprojSB_pion = (TH1D*) VecHistos.at(iHist).ProjectionX(Form("%s_ProjMass%g:%g_pt%g:%g",VecHistos.at(iHist).GetName(),vec_proj_bin_mass.at(k),vec_proj_bin_mass.at(k+1),vec_proj_bin_pt.at(j),vec_proj_bin_pt.at(j+1)) ,startbinY,endbinY)->Clone();  // Mass projection Histogram
        fMCprojSB_eta  = (TH1D*) VecHistos.at(iHist).ProjectionX(Form("%s_ProjMass%g:%g_pt%g:%g",VecHistos.at(iHist).GetName(),vec_proj_bin_mass.at(k),vec_proj_bin_mass.at(k+1),vec_proj_bin_pt.at(j),vec_proj_bin_pt.at(j+1)) ,startbinY,endbinY)->Clone();  // Mass projection Histogram
        fMCprojS_pion  = (TH1D*) VecHistos.at(3).ProjectionX(Form("%s_ProjMass%g:%g_pt%g:%g",VecHistos.at(3).GetName(),vec_proj_bin_mass.at(k),vec_proj_bin_mass.at(k+1),vec_proj_bin_pt.at(j),vec_proj_bin_pt.at(j+1)),startbinY,endbinY)->Clone();  // Mass projection Histogram
        fMCprojS_eta   = (TH1D*) VecHistos.at(1).ProjectionX(Form("%s_ProjMass%g:%g_pt%g:%g",VecHistos.at(1).GetName(),vec_proj_bin_mass.at(k),vec_proj_bin_mass.at(k+1),vec_proj_bin_pt.at(j),vec_proj_bin_pt.at(j+1)),startbinY,endbinY)->Clone();  // Mass projection Histogram

        // fMCprojSB_pion->Sumw2();
        // fMCprojSB_eta->Sumw2();
        // fMCprojS_pion->Sumw2();
        // fMCprojS_eta->Sumw2();

        // Normalise to Number of events
        fMCprojSB_pion  -> Divide(fTF1dummy,nEvents);
        fMCprojSB_eta   -> Divide(fTF1dummy,nEvents);
        fMCprojS_pion   -> Divide(fTF1dummy,nEvents);
        fMCprojS_eta    -> Divide(fTF1dummy,nEvents);

        fMCprojSB_pion -> Add(fMCprojS_pion,-1);
        fMCprojB_pion = fMCprojSB_pion;
        fMCprojSB_eta -> Add(fMCprojS_eta,-1);
        fMCprojB_eta = fMCprojSB_eta;
        SetHistoStandardSettingsBlue(fMCprojS_pion);
        SetHistoStandardSettingsBlue(fMCprojS_eta);
        SetHistoStandardSettingsGreen(fMCprojB_pion);
        SetHistoStandardSettingsGreen(fMCprojB_eta);
        fMCprojS_pion ->RebinX(vec_rebin_mass.at(j));
        fMCprojS_eta  ->RebinX(vec_rebin_mass.at(j));
        fMCprojB_pion ->RebinX(vec_rebin_mass.at(j));
        fMCprojB_eta  ->RebinX(vec_rebin_mass.at(j));

        TF1* fMCHistFit_pion = new TF1(Form("template1_ProjMass%g:%g_pt%g:%g",vec_proj_bin_mass.at(k),vec_proj_bin_mass.at(k+1),vec_proj_bin_pt.at(j),vec_proj_bin_pt.at(j+1)),template1,0.,1.,2);
        TF1* fMCHistFit_eta  = new TF1(Form("template2_ProjMass%g:%g_pt%g:%g",vec_proj_bin_mass.at(k),vec_proj_bin_mass.at(k+1),vec_proj_bin_pt.at(j),vec_proj_bin_pt.at(j+1)),template2,0.,1.,2);
        SetFitSettingsEta(fMCHistFit_pion);
        SetFitSettingsEta(fMCHistFit_eta);

        auto legendInfo = new TLegend(0.52,0.22,0.87,0.42);
        TLegendEntry *entry4=legendInfo->AddEntry("collisionSystem"   ,"pp, #sqrt{s} = 13 TeV, |#eta_{e}| < 0.8","");
        TLegendEntry *entry1=legendInfo->AddEntry("SinglePtPrim","#font[12]{p}_{T,e}^{prim} > 0.075 GeV/#font[12]{c}","");
        TLegendEntry *entry2=legendInfo->AddEntry("SinglePtSec" ,"#font[12]{p}_{T,e}^{sec} > 0.02 GeV/#font[12]{c}","");
        // TLegendEntry *entry2=legendInfo->AddEntry("PairPt"  ,"#font[12]{p}_{T,ee}^{prim} > 0.075 GeV/#font[12]{c} , #font[12]{p}_{T,ee}^{sec} > 0.02 GeV/#font[12]{c}","");
        TLegendEntry *entry3=legendInfo->AddEntry("FourPairPt"   ,Form("%g < #font[12]{p}_{T,eeee} < %g  GeV/#font[12]{c}",vec_proj_bin_pt.at(j) ,vec_proj_bin_pt.at(j+1)),"");
        legendInfo->SetBorderSize(0);
        legendInfo->SetFillColorAlpha(0, 0.0);
        legendInfo->SetTextSize(0.04);



        // Draw Mass projections
        if (DoMassProj == kTRUE) {

          // ------------------------
          // Plot in pion peak region
          if(ExtFourPairPionSig == kTRUE){
            TPad* fProjPad = new TPad("pad1","pad1",0.0,0.0,1.0,1.0,0);
            TPad* fSigPad  = new TPad("pad1","pad1",0.0,0.0,1.0,1.0,0);
            TLatex *massPion_Intervall  = new TLatex(startTextX, startTextY-.05, Form("#font[12]{m}_{ee}  Intervall = %g - %g  #frac{GeV}{#font[12]{c^{2}}}",lowerPionMassRange,upperPionMassRange));
            TH1D* projXpion = (TH1D*) projX->Clone(Form("%s_ProjMass%g:%g_pt%g:%g",VecHistos.at(iHist).GetName(),lowerPionMassRange,upperPionMassRange,vec_proj_bin_pt.at(j),vec_proj_bin_pt.at(j+1)));
            TH1D* fSignalPion = (TH1D*) projXpion->Clone(Form("%s_SignalPion%g:%g_pt%g:%g",VecHistos.at(iHist).GetName(),lowerPionMassRange,upperPionMassRange,vec_proj_bin_pt.at(j),vec_proj_bin_pt.at(j+1)));
            SetTextSettings(massPion_Intervall,textSize);
            SetMultipleHistoStandardSettings(projXpion);
            fSignalPion->SetTitle(Form("Signal in Pt Intervall %g-%g #frac{GeV}{#font[12]{c}}",vec_proj_bin_pt.at(j),vec_proj_bin_pt.at(j+1)));
            fProjPad->cd();
            // // projXpion->SetAxisRange(0, (projXpion->GetMaximum()+projXpion->GetBinError(projXpion->GetMaximumBin()))*1.2, "Y");
            projXpion->GetXaxis()->SetRangeUser(lowerPionMassRange,upperPionMassRange);
            projXpion->GetYaxis()->SetRangeUser(0.,(projXpion->GetMaximum()+projXpion->GetBinError(projXpion->GetMaximumBin()))*1.25);
            projXpion->GetYaxis()->SetTitle("Counts");
            projXpion->GetYaxis()->SetTitleOffset(3.0);
            projXpion->Draw("E1");

            auto legend = new TLegend(0.65,0.72,0.9,0.89);

            // TF1 *fPionFit = new TF1("fPionFit","pol2+gaus(3)",lowerPionMassRange,upperPionMassRange);
            TF1 *fPionFit = new TF1("fPionFit","(x<[1])*([0]*(TMath::Exp(-0.5*((x-[1])/[2])^2)+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*TMath::Exp(-0.5*((x-[1])/[2])^2))+pol2(4)+expo(7)",lowerPionMassRange,upperPionMassRange);
            // TF1* fSignalFit = new TF1("fSignalPion","gaus",lowerPionMassRange,upperPionMassRange);
            TF1* fSignalFit = new TF1("fSignalPion","(x<[1])*([0]*(TMath::Exp(-0.5*((x-[1])/[2])^2)+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*TMath::Exp(-0.5*((x-[1])/[2])^2))",lowerPionMassRange,upperPionMassRange);
            // TF1* fBackgroundFit = new TF1("fBackgroundPion","pol2",lowerPionMassRange,upperPionMassRange);
            TF1* fBackgroundFitPion = new TF1("fBackgroundFitPion",BackgroundFunctionPion,lowerPionMassRange,upperPionMassRange,5);
            SetFitSettingsPion(fPionFit);
            SetFitSettingsPionSignal(fSignalFit);
            SetFitSettingsPionBackground(fBackgroundFitPion);
            if(DoFit == kTRUE){
              if(DoTemplateFit == kFALSE) {
                fBackgroundFitPion->SetParameters(1,1,1,1,1);
                reject = kTRUE;
                                                                                                // if(DoDebug) Printf("%d", __LINE__);
                projXpion->Fit("fBackgroundFitPion" ,"QRMNE0","",lowerPionMassRange,upperPionMassRange);
                                                                                                // if(DoDebug) Printf("%d", __LINE__);
                projXpion->Fit("fBackgroundFitPion" ,"QRMNE0","",lowerPionMassRange,upperPionMassRange);
                                                                                                // if(DoDebug) Printf("%d", __LINE__);
                projXpion->Fit("fBackgroundFitPion" ,"QRMNE0","",lowerPionMassRange,upperPionMassRange);
                                                                                                // if(DoDebug) Printf("%d", __LINE__);
                reject = kFALSE;

                TLegendEntry *entry1=legend->AddEntry("Data","MC Data","P");
                TLegendEntry *entry2=legend->AddEntry("B","B Fit","l");
                TLegendEntry *entry3=legend->AddEntry("S+B","S+B Fit","l");
                legend->SetBorderSize(0);
                legend->SetFillColorAlpha(0, 0.0);
                legend->SetTextSize(0.04);
                entry1->SetMarkerColor(kBlack);       entry1->SetMarkerStyle(5);    entry1->SetLineColor(kBlack);
                entry2->SetLineColor(kGreen+2);       entry2->SetLineStyle(2);
                entry3->SetLineColor(kBlue);          entry3->SetLineStyle(2);

                fPionFit->SetParameters(projXpion->GetMaximum()/2,0.135,0.025,0.4,fBackgroundFitPion->GetParameter(0),fBackgroundFitPion->GetParameter(1),fBackgroundFitPion->GetParameter(2));
                // fPionFit->SetParameters(0,0.018,0.01,projXpion->GetMaximum()/2,0.135,0.025,0.4);
                // fPionFit->SetParLimits(4,0.130,0.140);
                // fPionFit->SetParLimits(5,0.0,0.02);
                fPionFit->SetParLimits(1,0.13,0.14);
                fPionFit->SetParLimits(2,0.0,0.1);
                fPionFit->SetParLimits(3,0,2);
                fPionFit->FixParameter(4,fBackgroundFitPion->GetParameter(0));
                fPionFit->FixParameter(5,fBackgroundFitPion->GetParameter(1));
                fPionFit->FixParameter(6,fBackgroundFitPion->GetParameter(2));
                fPionFit->FixParameter(7,fBackgroundFitPion->GetParameter(3));
                fPionFit->FixParameter(8,fBackgroundFitPion->GetParameter(4));
                                                                                                // if(DoDebug) Printf("%d", __LINE__);
                projXpion->Fit("fPionFit","QRMNE0");
                                                                                                // if(DoDebug) Printf("%d", __LINE__);
                projXpion->Fit("fPionFit","QRMNE0");
                                                                                                // if(DoDebug) Printf("%d", __LINE__);
                projXpion->Fit("fPionFit","QRMNE0");
                                                                                                // if(DoDebug) Printf("%d", __LINE__);

                fSignalFit->SetParameters(fPionFit->GetParameter(0),fPionFit->GetParameter(1),fPionFit->GetParameter(2),fPionFit->GetParameter(3));
                fPionFit->Draw("same");
                fSignalFit    -> Draw("same");
                fBackgroundFitPion-> Draw("same");
              }

              if (DoTemplateFit == kTRUE) {
                TLegendEntry *entry1=legend->AddEntry("Data","MC Data","P");
                TLegendEntry *entry2=legend->AddEntry("S","Signal"    ,"P");
                TLegendEntry *entry3=legend->AddEntry("B","Background","P");
                // TLegendEntry *entry4=legend->AddEntry("Fit","S+B Template Fit","l");
                TLegendEntry *entry5=legend->AddEntry("Fit","B Fit Pol2+Expo","l");
                legend->SetBorderSize(0);
                legend->SetFillColorAlpha(0, 0.0);
                legend->SetTextSize(0.04);
                entry1->SetMarkerColor(kBlack);       entry1->SetMarkerStyle(5);    entry1->SetLineColor(kBlack);
                entry2->SetMarkerColor(kBlue);        entry2->SetMarkerStyle(2);    entry2->SetLineColor(kBlue);
                entry3->SetMarkerColor(kGreen+2);     entry3->SetMarkerStyle(4);    entry3->SetLineColor(kGreen+2);
                // entry4->SetLineColor(kBlue);          entry4->SetLineStyle(2);
                entry5->SetLineColor(kGreen+2);       entry5->SetLineStyle(2);

                fMCHistFit_pion->SetParameter(0,0.9);
                fMCHistFit_pion->SetParameter(1,0.9);
                                                                                                // if(DoDebug) Printf("%d", __LINE__);
                projXpion->Fit(fMCHistFit_pion,"QIN");
                                                                                                // if(DoDebug) Printf("%d", __LINE__);

                fBackgroundFitPion->SetParameters(1,1,1,1,1);
                                                                                                // if(DoDebug) Printf("%d", __LINE__);
                fMCprojB_pion->Fit("fBackgroundFitPion" ,"QRMNE0","",lowerPionMassRange,upperPionMassRange);
                                                                                                // if(DoDebug) Printf("%d", __LINE__);
                fMCprojB_pion->Fit("fBackgroundFitPion" ,"QRMNE0","",lowerPionMassRange,upperPionMassRange);
                                                                                                // if(DoDebug) Printf("%d", __LINE__);
                fMCprojB_pion->Fit("fBackgroundFitPion" ,"QRMNE0","",lowerPionMassRange,upperPionMassRange);
                                                                                                // if(DoDebug) Printf("%d", __LINE__);
                fBackgroundFitPion-> Draw("same");

                fMCprojB_pion   -> Draw("same E1");
                fMCprojS_pion   -> Draw("same E1");
                projXpion   -> Draw("same P");
                // fMCHistFit_pion -> DrawClone("same");
              }
            }
            // pT_Intervall  -> Draw("same");
            legendInfo    -> Draw("same");
            legend        -> Draw("same");
            cSignalCanvas -> cd();
            fProjPad      -> Draw();
            cSignalCanvas -> SaveAs(Form("%s%s_Pion_rebin%i.pdf",DocumentPathProjX.Data(), projXpion->GetName(),vec_rebin_mass.at(0)));

            cProjPtPionCanvas->cd(j+1);
            fProjPad->Draw();


            // Substract background
            if (DoTemplateFit == kFALSE) {
              for (int iBin = 1; iBin < projXpion->GetNbinsX(); iBin++) {
                Double_t binCenter = projXpion->GetBinCenter(iBin);
                Double_t fitValue  = fBackgroundFitPion->Eval(binCenter);
                Double_t newBinContent = projXpion->GetBinContent(iBin) - fitValue;
                fSignalPion->SetBinContent(iBin, newBinContent);
              }
              SetMultipleHistoStandardSettings(fSignalPion);

              TAxis *xAxis = fSignalPion->GetXaxis();
              // Double_t lowerIntEdge = fSignalFit->GetParameter(1)-3*fSignalFit->GetParameter(2);
              Double_t lowerIntEdge = 0.1;
              Double_t upperIntEdge = fSignalFit->GetParameter(1)+3*fSignalFit->GetParameter(2);
              Int_t binLowerSigmaEdge   = xAxis->FindBin(lowerIntEdge);
              Int_t binUpperSigmaEdge   = xAxis->FindBin(upperIntEdge);
              Double_t fSignalValue     = fSignalPion->Integral(binLowerSigmaEdge,binUpperSigmaEdge);
              // Double_t fBackgroundValue = fBackgroundFitPion->Integral(fSignalFit->GetParameter(1)-3*fSignalFit->GetParameter(2),fSignalFit->GetParameter(1)+3*fSignalFit->GetParameter(2));
              Double_t fBackgroundValue = 0;
              for (int iSigmaBins = 0; iSigmaBins < (binUpperSigmaEdge-binLowerSigmaEdge); iSigmaBins++) {
                fBackgroundValue = fBackgroundValue + fBackgroundFitPion->Eval(projXpion->GetBinCenter(binLowerSigmaEdge+iSigmaBins));
              }
              // std::cout << fBackgroundValue << " " << fBackgroundValue1 << std::endl;
              Double_t fSignalOverBackground = fSignalValue/fBackgroundValue;
              Double_t fSignificance = fSignalValue/TMath::Sqrt(fSignalValue + fBackgroundValue);
              TLatex *fTextSigBack       = new TLatex(startTextX+0.15, startTextY     , Form("#frac{S}{B} = %g", fSignalOverBackground));
              TLatex *fTextSignificance  = new TLatex(startTextX+0.15, startTextY-0.1 , Form("#frac{S}{#sqrt{S+B}} = %g", fSignificance));
              SetTextSettings(fTextSigBack,textSize);
              SetTextSettings(fTextSignificance,textSize);

              fSigPad->cd();
              fSignalPion->GetXaxis()->SetRangeUser(lowerPionMassRange,upperPionMassRange);  // Range of x-Axis should be the same of the fit region
              TLine* fLowerSigmaLine = new TLine(lowerIntEdge,(fSignalPion->GetMinimum()-fSignalPion->GetBinError(fSignalPion->GetMinimumBin()))*1.05,lowerIntEdge,(fSignalPion->GetMaximum()+fSignalPion->GetBinError(fSignalPion->GetMaximumBin()))*1.05);
              TLine* fUpperSigmaLine = new TLine(upperIntEdge,(fSignalPion->GetMinimum()-fSignalPion->GetBinError(fSignalPion->GetMinimumBin()))*1.05,upperIntEdge,(fSignalPion->GetMaximum()+fSignalPion->GetBinError(fSignalPion->GetMaximumBin()))*1.05);
              SetSigmaLineSettings(fLowerSigmaLine);
              SetSigmaLineSettings(fUpperSigmaLine);
              fSignalPion->Draw();
              fSignalFit->Draw("same");
              fTextSigBack->Draw("same");
              fTextSignificance->Draw("same");
              fLowerSigmaLine->Draw("same");
              fUpperSigmaLine->Draw("same");
              cSignalCanvas->cd();
              fSigPad->Draw();
              cSignalCanvas->SaveAs(Form("%s%s.pdf",DocumentPathSignals.Data(), fSignalPion->GetName()));
              cSignalPionCanvas->cd(j+1);
              fSigPad->Draw();
              Vec_SigBack_Pion.push_back(fSignalOverBackground);
            }


            if (DoTemplateFit == kTRUE) {
              // //----------------------------
              // // S/B & S/sqrt(S+2B) integrated in an intervall
              // Double_t lowerIntEdge = fMCprojS_pion->GetXaxis()->FindBin(0.1);
              // Double_t upperIntEdge = fMCprojS_pion->GetXaxis()->FindBin(0.16);
              // Double_t fSignalValue     = fMCHistFit_pion->GetParameter(0)*(fMCprojS_pion->Integral(lowerIntEdge,upperIntEdge));
              // Double_t fBackgroundValue = fMCHistFit_pion->GetParameter(1)*(fMCprojB_pion->Integral(lowerIntEdge,upperIntEdge));
              //
              // Double_t fSignalOverBackground = fSignalValue/fBackgroundValue;
              // Double_t fSignificance = fSignalValue/TMath::Sqrt(fSignalValue + fBackgroundValue);
              // if (fSignalOverBackground > 100) fSignalOverBackground = 0;
              // TLatex *fTextSigBack       = new TLatex(startTextX+0.15, startTextY     , Form("#frac{S}{B} = %g", fSignalOverBackground));
              // TLatex *fTextSignificance  = new TLatex(startTextX+0.15, startTextY-0.1 , Form("#frac{S}{#sqrt{S+B}} = %g", fSignificance));
              // SetTextSettings(fTextSigBack,textSize);
              // SetTextSettings(fTextSignificance,textSize);
              //
              // fSigPad->cd();
              // fMCprojS_pion->GetXaxis()->SetRangeUser(lowerPionMassRange,upperPionMassRange);  // Range of x-Axis should be the same of the fit region
              // TLine* fLowerIntLine = new TLine(lowerIntEdge,(fMCprojS_pion->GetMinimum()-fMCprojS_pion->GetBinError(fMCprojS_pion->GetMinimumBin()))*1.05,upperIntEdge,(fMCprojS_pion->GetMaximum()+fMCprojS_pion->GetBinError(fMCprojS_pion->GetMaximumBin()))*1.05);
              // TLine* fUpperIntLine = new TLine(lowerIntEdge,(fMCprojS_pion->GetMinimum()-fMCprojS_pion->GetBinError(fMCprojS_pion->GetMinimumBin()))*1.05,upperIntEdge,(fMCprojS_pion->GetMaximum()+fMCprojS_pion->GetBinError(fMCprojS_pion->GetMaximumBin()))*1.05);
              // SetSigmaLineSettings(fLowerIntLine);
              // SetSigmaLineSettings(fUpperIntLine);
              // TF1 *fDummy = new TF1 ("fDummy", "1", -100., 100.);
              // fMCprojS_pion->Multiply(fDummy,fMCHistFit_pion->GetParameter(0));
              // fMCprojS_pion->Draw();
              // fTextSigBack->Draw("same");
              // fTextSignificance->Draw("same");
              // fLowerIntLine->Draw("same");
              // fUpperIntLine->Draw("same");
              // cSignalCanvas->cd();
              // fSigPad->Draw();
              // cSignalCanvas -> SaveAs(Form("%s%s_rebin%i.pdf",DocumentPathSignals.Data(), fMCprojS_pion->GetName(), vec_rebin_mass.at(1)));
              // // cSignalPionCanvas->cd(j+1);
              // // fSigPad->Draw();
              // Vec_SigBack_Pion.push_back(fSignalOverBackground);

              // S/B & S/sqrt(S+B) bin wise
              Int_t startBin = fMCprojS_pion->GetXaxis()->FindBin(0.0);
              Int_t endBin   = fMCprojS_pion->GetXaxis()->FindBin(0.3);
              Double_t fBinSignalValue;
              Double_t fBinBackgroundValue;
              Double_t fBinSignalOverBackground;
              Double_t fBinSignificance;
              TH1D* hSignalOverBackground = new TH1D("SignalOverBackgroundPion", "   ; m_{eeee} (GeV/#font[12]{c^{2}}); Ratio #frac{S}{B}"                , endBin-startBin+1,fMCprojS_pion->GetXaxis()->GetBinLowEdge(startBin),fMCprojS_pion->GetXaxis()->GetBinUpEdge(endBin));
              TH1D* hSignificance         = new TH1D("SignificancePion"        , "   ; m_{eeee} (GeV/#font[12]{c^{2}}); Significance #frac{S}{#sqrt{S+B}}", endBin-startBin+1,fMCprojS_pion->GetXaxis()->GetBinLowEdge(startBin),fMCprojS_pion->GetXaxis()->GetBinUpEdge(endBin));
              hSignalOverBackground -> GetYaxis()->SetTitleOffset(1.2);
              hSignificance         -> GetYaxis()->SetTitleOffset(1.1);

              for (Int_t iBin = 0; iBin < endBin-startBin+1; iBin++) {
                fBinSignalValue     =  fMCHistFit_pion->GetParameter(0)*(fMCprojS_pion->GetBinContent(iBin+startBin));
                fBinBackgroundValue =  fMCHistFit_pion->GetParameter(1)*(fMCprojB_pion->GetBinContent(iBin+startBin));
                if (fBinBackgroundValue != 0)      fBinSignalOverBackground = fBinSignalValue/fBinBackgroundValue;
                else if (fBinBackgroundValue == 0) fBinSignalOverBackground = 0;
                if (fBinBackgroundValue != 0 && fBinSignalValue != 0 && fBinBackgroundValue > 0)      fBinSignificance = fBinSignalValue/TMath::Sqrt(fBinSignalValue + fBinBackgroundValue);
                else fBinSignificance = 0;
                hSignalOverBackground ->SetBinContent(iBin,fBinSignalOverBackground);
                hSignificance         ->SetBinContent(iBin,fBinSignificance);
              }
              cSignalCanvas->cd();
              hSignalOverBackground->Draw();
              // pT_Intervall->Draw("same");
              legendInfo    -> Draw("same");
              hSignalOverBackground -> SaveAs(Form("%s%s_pt%g:%g_rebin%i.root",DocumentPathSignals.Data(), hSignalOverBackground->GetName(), vec_proj_bin_pt.at(j),vec_proj_bin_pt.at(j+1), vec_rebin_mass.at(0)));
              // cSignalCanvas -> SaveAs(Form("%s%s_pt%g:%g_rebin%i.pdf",DocumentPathSignals.Data(), hSignalOverBackground->GetName(), vec_proj_bin_pt.at(j),vec_proj_bin_pt.at(j+1), vec_rebin_mass.at(0)));
              hSignificance->Draw();
              // pT_Intervall->Draw("same");
              legendInfo    -> Draw("same");
              hSignificance -> SaveAs(Form("%s%s_pt%g:%g_rebin%i.root",DocumentPathSignals.Data(), hSignificance->GetName(), vec_proj_bin_pt.at(j),vec_proj_bin_pt.at(j+1), vec_rebin_mass.at(0)));
              // cSignalCanvas -> SaveAs(Form("%s%s_pt%g:%g_rebin%i.pdf",DocumentPathSignals.Data(), hSignificance->GetName(), vec_proj_bin_pt.at(j),vec_proj_bin_pt.at(j+1), vec_rebin_mass.at(0)));
              cSoBPionCanvas->cd(j+1);
              hSignalOverBackground->Draw();
              // pT_Intervall->Draw("same");
              legendInfo    -> Draw("same");
              cSignificancePionCanvas->cd(j+1);
              hSignificance->Draw();
              // pT_Intervall->Draw("same");
              legendInfo    -> Draw("same");

            }
          }

                                                                                if(DoDebug) Printf("%d", __LINE__);

          // ------------------------
          // Plot in eta peak region
          if(ExtFourPairEtaSig == kTRUE){
            TPad* fProjPad = new TPad("pad1","pad1",0.0,0.0,1.0,1.0,0);
            TPad* fSigPad  = new TPad("pad1","pad1",0.0,0.0,1.0,1.0,0);
            TLatex *massEta_Intervall  = new TLatex(startTextX, startTextY-.05, Form("#font[12]{m}_{eeee}  Intervall = %g - %g  #frac{GeV}{#font[12]{c^{2}}}",lowerEtaMassRange,upperEtaMassRange));
            TH1D* projXeta  = (TH1D*) projX->Clone(Form("%s_ProjMass%g:%g_pt%g:%g",VecHistos.at(iHist).GetName(),lowerEtaMassRange,upperEtaMassRange,vec_proj_bin_pt.at(j),vec_proj_bin_pt.at(j+1)));
            TH1D* fSignalEta = (TH1D*) projXeta->Clone(Form("%s_SignalEta%g:%g_pt%g:%g",VecHistos.at(iHist).GetName(),lowerEtaMassRange,upperEtaMassRange,vec_proj_bin_pt.at(j),vec_proj_bin_pt.at(j+1)));
            SetTextSettings(massEta_Intervall,textSize);
            SetMultipleHistoStandardSettings(projXeta);
            fSignalEta->SetTitle(Form("Signal in Pt Intervall %g-%g #frac{GeV}{#font[12]{c}}",vec_proj_bin_pt.at(j),vec_proj_bin_pt.at(j+1)));
            fProjPad->cd();
            // // projXeta->SetAxisRange(0, (projXeta->GetMaximum()+projXeta->GetBinError(projXeta->GetMaximumBin()))*1.2, "Y");
            projXeta->GetXaxis()->SetRangeUser(lowerEtaMassRange,upperEtaMassRange);
            projXeta->GetYaxis()->SetRangeUser(0.,(projXeta->GetMaximum()+projXeta->GetBinError(projXeta->GetMaximumBin()))*1.25);
            projXeta->GetYaxis()->SetTitle("Counts");
            projXeta->GetYaxis()->SetTitleOffset(3.5);
            projXeta->Draw("E1");

            // auto legend = new TLegend(0.11,0.72,0.3,0.89);
            auto legend = new TLegend(0.65,0.72,0.9,0.89);

            // TF1 *fEtaFit = new TF1("fEtaFit","pol2+gaus(3)",lowerEtaMassRange,upperEtaMassRange);
            TF1 *fEtaFit = new TF1("fEtaFit","(x<[1])*([0]*(TMath::Exp(-0.5*((x-[1])/[2])^2)+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*TMath::Exp(-0.5*((x-[1])/[2])^2))+pol2(4)+expo(7)",lowerEtaMassRange,upperEtaMassRange);
            // TF1* fSignalFit = new TF1("fSignalEta","gaus",lowerEtaMassRange,upperEtaMassRange);
            TF1* fSignalFit = new TF1("fSignalEta","(x<[1])*([0]*(TMath::Exp(-0.5*((x-[1])/[2])^2)+TMath::Exp((x-[1])/[3])*(1.-TMath::Exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*TMath::Exp(-0.5*((x-[1])/[2])^2))",lowerEtaMassRange,upperEtaMassRange);
            // TF1* fBackgroundFit = new TF1("fBackgroundFit",BackgroundFunctionEta,lowerEtaMassRange,upperEtaMassRange,3);
            TF1* fBackgroundFit = new TF1("fBackgroundFit",BackgroundFunctionEta,lowerEtaMassRange,upperEtaMassRange,5);
            // TF1* fBackgroundFit = new TF1("fBackgroundFit","pol2",lowerEtaMassRange,upperEtaMassRange);
            SetFitSettingsEta(fEtaFit);
            SetFitSettingsEta(fSignalFit);
            SetFitSettingsEtaBackground(fBackgroundFit);
            if(DoFit == kTRUE){
              if (DoTemplateFit == kFALSE) {
                fBackgroundFit->SetParameters(1,1,1,1,1);
                reject = kTRUE;
                                                                                                                // if(DoDebug) Printf("%d", __LINE__);
                projXeta->Fit("fBackgroundFit" ,"QRMNE0","",lowerEtaMassRange,upperEtaMassRange);
                                                                                                                // if(DoDebug) Printf("%d", __LINE__);
                projXeta->Fit("fBackgroundFit" ,"QRMNE0","",lowerEtaMassRange,upperEtaMassRange);
                                                                                                                // if(DoDebug) Printf("%d", __LINE__);
                projXeta->Fit("fBackgroundFit" ,"QRMNE0","",lowerEtaMassRange,upperEtaMassRange);
                                                                                                                // if(DoDebug) Printf("%d", __LINE__);
                reject = kFALSE;

                TLegendEntry *entry1=legend->AddEntry("Data","MC Data","P");
                TLegendEntry *entry2=legend->AddEntry("B","B Fit","l");
                TLegendEntry *entry3=legend->AddEntry("S+B","S+B Fit","l");
                legend->SetBorderSize(0);
                legend->SetFillColorAlpha(0, 0.0);
                legend->SetTextSize(0.04);
                entry1->SetMarkerColor(kBlack);       entry1->SetMarkerStyle(5);    entry1->SetLineColor(kBlack);
                entry2->SetLineColor(kGreen+2);       entry2->SetLineStyle(2);
                entry3->SetLineColor(kBlue);          entry3->SetLineStyle(2);

                // fEtaFit->SetParameters(0,0.018,0.01,projXeta->GetMaximum()/2,0.547,0.01,0.4);
                fEtaFit->SetParameters(projXeta->GetMaximum()/2,0.547,0.01,0.4,fBackgroundFit->GetParameter(0),fBackgroundFit->GetParameter(1),fBackgroundFit->GetParameter(2));
                // // fEtaFit->SetParLimits(4,0.53,0.565);  // (*)
                // // fEtaFit->SetParLimits(5,0.0,0.02);    // (*)
                fEtaFit->SetParLimits(1,0.54,0.56);
                fEtaFit->SetParLimits(2,0.0,0.1);
                fEtaFit->SetParLimits(3,0,2);
                fEtaFit->FixParameter(4,fBackgroundFit->GetParameter(0));
                fEtaFit->FixParameter(5,fBackgroundFit->GetParameter(1));
                fEtaFit->FixParameter(6,fBackgroundFit->GetParameter(2));
                fEtaFit->FixParameter(7,fBackgroundFit->GetParameter(3));
                fEtaFit->FixParameter(8,fBackgroundFit->GetParameter(4));
                                                                                                // if(DoDebug) Printf("%d", __LINE__);
                projXeta->Fit("fEtaFit","QRMNE0");
                                                                                                // if(DoDebug) Printf("%d", __LINE__);
                projXeta->Fit("fEtaFit","QRMNE0");
                                                                                                // if(DoDebug) Printf("%d", __LINE__);
                projXeta->Fit("fEtaFit","QRMNE0");
                                                                                                // if(DoDebug) Printf("%d", __LINE__);

                fSignalFit->SetParameters(fEtaFit->GetParameter(0),fEtaFit->GetParameter(1),fEtaFit->GetParameter(2),fEtaFit->GetParameter(3));
                fEtaFit->Draw("same");
                fSignalFit    -> Draw("same");
                fBackgroundFit-> Draw("same");
              }

              if (DoTemplateFit == kTRUE) {
                TLegendEntry *entry1=legend->AddEntry("Data","MC Data","P");
                TLegendEntry *entry2=legend->AddEntry("S","Signal"    ,"P");
                TLegendEntry *entry3=legend->AddEntry("B","Background","P");
                // TLegendEntry *entry4=legend->AddEntry("Fit","S+B Template Fit","l");
                TLegendEntry *entry5=legend->AddEntry("Fit","B Fit Pol2+Expo","l");
                legend->SetBorderSize(0);
                legend->SetFillColorAlpha(0, 0.0);
                legend->SetTextSize(0.04);
                entry1->SetMarkerColor(kBlack);       entry1->SetMarkerStyle(5);    entry1->SetLineColor(kBlack);
                entry2->SetMarkerColor(kBlue);        entry2->SetMarkerStyle(2);    entry2->SetLineColor(kBlue);
                entry3->SetMarkerColor(kGreen+2);     entry3->SetMarkerStyle(4);    entry3->SetLineColor(kGreen+2);
                // entry4->SetLineColor(kBlue);          entry4->SetLineStyle(2);
                entry5->SetLineColor(kGreen+2);       entry5->SetLineStyle(2);

                fMCHistFit_eta->SetParameter(0,0.9);
                fMCHistFit_eta->SetParameter(1,0.9);
                                                                                                // if(DoDebug) Printf("%d", __LINE__);
                projXeta->Fit(fMCHistFit_eta,"QIN");
                                                                                                // if(DoDebug) Printf("%d", __LINE__);

                fBackgroundFit->SetParameters(1,1,1,1,1);
                // fBackgroundFit->SetParameters(1/nEvents,1,1,1,1);
                // fBackgroundFit->SetParameters(1/nEvents,1/nEvents,1/nEvents,1,1);
                // fBackgroundFit->SetParameters(1/nEvents,1/nEvents,1/nEvents,1/nEvents,1/nEvents);
                // fBackgroundFit->SetParameters(-0.7,0.05,-0.002,-0.35,-0.075);
                                                                                                // if(DoDebug) Printf("%d", __LINE__);
                fMCprojB_eta->Fit("fBackgroundFit" ,"QRMNE0","",lowerEtaMassRange,upperEtaMassRange);
                                                                                                // if(DoDebug) Printf("%d", __LINE__);
                fMCprojB_eta->Fit("fBackgroundFit" ,"QRMNE0","",lowerEtaMassRange,upperEtaMassRange);
                                                                                                // if(DoDebug) Printf("%d", __LINE__);
                fMCprojB_eta->Fit("fBackgroundFit" ,"QRMNE0","",lowerEtaMassRange,upperEtaMassRange);
                                                                                                // if(DoDebug) Printf("%d", __LINE__);
                fBackgroundFit-> Draw("same");

                fMCprojB_eta   -> Draw("same E1");
                fMCprojS_eta   -> Draw("same E1");
                projXeta   -> Draw("same P");
                // fMCHistFit_eta -> DrawClone("same");
              }

            }
            // pT_Intervall  -> Draw("same");
            legendInfo    -> Draw("same");
            legend        -> Draw("same");
            cSignalCanvas -> cd();
            fProjPad      -> Draw();
            cSignalCanvas -> SaveAs(Form("%s%s_Eta_rebin%i.pdf",DocumentPathProjX.Data(), projXeta->GetName(),vec_rebin_mass.at(0)));

            cProjPtEtaCanvas->cd(j+1);
            fProjPad->Draw();


            // Substract background
            if (DoTemplateFit == kFALSE) {
              for (int iBin = 1; iBin <= projXeta->GetNbinsX(); iBin++) {
                Double_t binCenter = projXeta->GetBinCenter(iBin);
                Double_t fitValue  = fBackgroundFit->Eval(binCenter);
                Double_t newBinContent = projXeta->GetBinContent(iBin) - fitValue;
                fSignalEta->SetBinContent(iBin, newBinContent);
              }
              SetMultipleHistoStandardSettings(fSignalEta);

              TAxis *xAxis = fSignalEta->GetXaxis();
              // Double_t lowerIntEdge = fSignalFit->GetParameter(1)-3*fSignalFit->GetParameter(2);
              Double_t lowerIntEdge = 0.45;
              Double_t upperIntEdge = fSignalFit->GetParameter(1)+3*fSignalFit->GetParameter(2);
              Int_t binLowerSigmaEdge   = xAxis->FindBin(lowerIntEdge);
              Int_t binUpperSigmaEdge   = xAxis->FindBin(upperIntEdge);
              Double_t fSignalValue     = fSignalEta->Integral(binLowerSigmaEdge,binUpperSigmaEdge);
              // Double_t fBackgroundValue = fBackgroundFit->Integral(fSignalFit->GetParameter(1)-3*fSignalFit->GetParameter(2),fSignalFit->GetParameter(1)+3*fSignalFit->GetParameter(2));
              Double_t fBackgroundValue = 0;
              for (int iSigmaBins = 0; iSigmaBins < (binUpperSigmaEdge-binLowerSigmaEdge); iSigmaBins++) {
                fBackgroundValue = fBackgroundValue + fBackgroundFit->Eval(projXeta->GetBinCenter(binLowerSigmaEdge+iSigmaBins));
              }
              Double_t fSignalOverBackground = fSignalValue/fBackgroundValue;
              Double_t fSignificance = fSignalValue/TMath::Sqrt(fSignalValue + fBackgroundValue);
              TLatex *fTextSigBack       = new TLatex(startTextX+0.15, startTextY     , Form("#frac{S}{B} = %g", fSignalOverBackground));
              TLatex *fTextSignificance  = new TLatex(startTextX+0.15, startTextY-0.1 , Form("#frac{S}{#sqrt{S+B}} = %g", fSignificance));
              SetTextSettings(fTextSigBack,textSize);
              SetTextSettings(fTextSignificance,textSize);

              fSigPad->cd();
              fSignalEta->GetXaxis()->SetRangeUser(lowerEtaMassRange,upperEtaMassRange);  // Range of x-Axis should be the same of the fit region
              // fSignalEta->GetYaxis()->SetRangeUser((fSignalEta->GetMinimum()-fSignalEta->GetBinError(fSignalEta->GetMinimumBin()))*1.2,(fSignalEta->GetMaximum()+fSignalEta->GetBinError(fSignalEta->GetMaximumBin()))*1.2);
              TLine* fLowerSigmaLine = new TLine(lowerIntEdge,(fSignalEta->GetMinimum()-fSignalEta->GetBinError(fSignalEta->GetMinimumBin()))*1.05,lowerIntEdge,(fSignalEta->GetMaximum()+fSignalEta->GetBinError(fSignalEta->GetMaximumBin()))*1.05);
              TLine* fUpperSigmaLine = new TLine(upperIntEdge,(fSignalEta->GetMinimum()-fSignalEta->GetBinError(fSignalEta->GetMinimumBin()))*1.05,upperIntEdge,(fSignalEta->GetMaximum()+fSignalEta->GetBinError(fSignalEta->GetMaximumBin()))*1.05);
              SetSigmaLineSettings(fLowerSigmaLine);
              SetSigmaLineSettings(fUpperSigmaLine);
              fSignalEta->Draw();
              fSignalFit->Draw("same");
              fTextSigBack->Draw("same");
              fTextSignificance->Draw("same");
              fLowerSigmaLine->Draw("same");
              fUpperSigmaLine->Draw("same");
              // pT_Intervall->Draw("same");
              legendInfo    -> Draw("same");
              cSignalCanvas->cd();
              fSigPad->Draw();
              cSignalCanvas -> SaveAs(Form("%s%s.pdf",DocumentPathSignals.Data(), fSignalEta->GetName()));
              cSignalEtaCanvas->cd(j+1);
              fSigPad->Draw();
              Vec_SigBack_Eta.push_back(fSignalOverBackground);
            }

            if (DoTemplateFit == kTRUE) {
              //----------------------------
              // // S/B & S/sqrt(S+2B) integrated in an intervall
              // Double_t lowerIntEdge = fMCprojS_eta->GetXaxis()->FindBin(0.45);
              // Double_t upperIntEdge = fMCprojS_eta->GetXaxis()->FindBin(0.61);
              // Double_t fSignalValue     = fMCHistFit_eta->GetParameter(0)*(fMCprojS_eta->Integral(lowerIntEdge,upperIntEdge));
              // Double_t fBackgroundValue = fMCHistFit_eta->GetParameter(1)*(fMCprojB_eta->Integral(lowerIntEdge,upperIntEdge));
              //
              // Double_t fSignalOverBackground = fSignalValue/fBackgroundValue;
              // Double_t fSignificance = fSignalValue/TMath::Sqrt(fSignalValue + fBackgroundValue);
              // if (fSignalOverBackground > 100) fSignalOverBackground = 0;
              // TLatex *fTextSigBack    = new TLatex(startTextX+0.15, startTextY    , Form("#frac{S}{B} = %g", fSignalOverBackground));
              // TLatex *fTextSignificance  = new TLatex(startTextX+0.15, startTextY-0.1 , Form("#frac{S}{#sqrt{S+B}} = %g", fSignificance));
              // SetTextSettings(fTextSigBack,textSize);
              // SetTextSettings(fTextSignificance,textSize);
              //
              // fSigPad->cd();
              // fMCprojS_eta->GetXaxis()->SetRangeUser(lowerEtaMassRange,upperEtaMassRange);  // Range of x-Axis should be the same of the fit region
              // TLine* fLowerIntLine = new TLine(lowerIntEdge,(fMCprojS_eta->GetMinimum()-fMCprojS_eta->GetBinError(fMCprojS_eta->GetMinimumBin()))*1.05,upperIntEdge,(fMCprojS_eta->GetMaximum()+fMCprojS_eta->GetBinError(fMCprojS_eta->GetMaximumBin()))*1.05);
              // TLine* fUpperIntLine = new TLine(lowerIntEdge,(fMCprojS_eta->GetMinimum()-fMCprojS_eta->GetBinError(fMCprojS_eta->GetMinimumBin()))*1.05,upperIntEdge,(fMCprojS_eta->GetMaximum()+fMCprojS_eta->GetBinError(fMCprojS_eta->GetMaximumBin()))*1.05);
              // SetSigmaLineSettings(fLowerIntLine);
              // SetSigmaLineSettings(fUpperIntLine);
              // TF1 *fDummy = new TF1 ("fDummy", "1", -100., 100.);
              // fMCprojS_eta->Multiply(fDummy,fMCHistFit_eta->GetParameter(0));
              // fMCprojS_eta->Draw();
              // fTextSigBack->Draw("same");
              // fTextSignificance->Draw("same");
              // fLowerIntLine->Draw("same");
              // fUpperIntLine->Draw("same");
              // pT_Intervall->Draw("same");
              // cSignalCanvas->cd();
              // fSigPad->Draw();
              // cSignalCanvas -> SaveAs(Form("%s%s_rebin%i.pdf",DocumentPathSignals.Data(), fMCprojS_eta->GetName(), vec_rebin_mass.at(1)));
              // // cSignalEtaCanvas->cd(j+1);
              // // fSigPad->Draw();
              // Vec_SigBack_Eta.push_back(fSignalOverBackground);

              // S/B & S/sqrt(S+B) bin wise
              Int_t startBin = fMCprojS_eta->GetXaxis()->FindBin(0.3);
              Int_t endBin   = fMCprojS_eta->GetXaxis()->FindBin(0.8);
              Double_t fBinSignalValue;
              Double_t fBinBackgroundValue;
              Double_t fBinSignalOverBackground;
              Double_t fBinSignificance;
              TH1D* hSignalOverBackground = new TH1D("SignalOverBackgroundEta", "  ; m_{eeee} #frac{GeV}{#font[12]{c^{2}}}; Ratio #frac{S}{B}"                 , endBin-startBin+1,fMCprojS_eta->GetXaxis()->GetBinLowEdge(startBin),fMCprojS_eta->GetXaxis()->GetBinUpEdge(endBin));
              TH1D* hSignificance         = new TH1D("SignificanceEta"        , "  ; m_{eeee} #frac{GeV}{#font[12]{c^{2}}}; Significance #frac{S}{#sqrt{S+B}}" , endBin-startBin+1,fMCprojS_eta->GetXaxis()->GetBinLowEdge(startBin),fMCprojS_eta->GetXaxis()->GetBinUpEdge(endBin));
              hSignalOverBackground -> GetYaxis()->SetTitleOffset(1.2);
              hSignificance         -> GetYaxis()->SetTitleOffset(1.1);


              for (Int_t iBin = 0; iBin < endBin-startBin+1; iBin++) {
                fBinSignalValue     =  fMCHistFit_eta->GetParameter(0)*(fMCprojS_eta->GetBinContent(iBin+startBin));
                fBinBackgroundValue =  fMCHistFit_eta->GetParameter(1)*(fMCprojB_eta->GetBinContent(iBin+startBin));
                if (fBinBackgroundValue != 0)      fBinSignalOverBackground = fBinSignalValue/fBinBackgroundValue;
                else if (fBinBackgroundValue == 0) fBinSignalOverBackground = 0;
                if (fBinBackgroundValue != 0 && fBinSignalValue != 0 && fBinBackgroundValue > 0)      fBinSignificance = fBinSignalValue/TMath::Sqrt(fBinSignalValue + fBinBackgroundValue);
                else fBinSignificance = 0;
                hSignalOverBackground ->SetBinContent(iBin,fBinSignalOverBackground);
                hSignificance         ->SetBinContent(iBin,fBinSignificance);
              }
              cSignalCanvas->cd();
              hSignalOverBackground->Draw();
              // pT_Intervall->Draw("same");
              legendInfo    -> Draw("same");
              hSignalOverBackground -> SaveAs(Form("%s%s_pt%g:%g_rebin%i.root",DocumentPathSignals.Data(), hSignalOverBackground->GetName(), vec_proj_bin_pt.at(j),vec_proj_bin_pt.at(j+1), vec_rebin_mass.at(0)));
              // cSignalCanvas -> SaveAs(Form("%s%s_pt%g:%g_rebin%i.pdf",DocumentPathSignals.Data(), hSignalOverBackground->GetName(), vec_proj_bin_pt.at(j),vec_proj_bin_pt.at(j+1), vec_rebin_mass.at(0)));
              hSignificance->Draw();
              // pT_Intervall->Draw("same");
              legendInfo    -> Draw("same");
              hSignificance -> SaveAs(Form("%s%s_pt%g:%g_rebin%i.root",DocumentPathSignals.Data(), hSignificance->GetName(), vec_proj_bin_pt.at(j),vec_proj_bin_pt.at(j+1), vec_rebin_mass.at(0)));
              // cSignalCanvas -> SaveAs(Form("%s%s_pt%g:%g_rebin%i.pdf",DocumentPathSignals.Data(), hSignificance->GetName(), vec_proj_bin_pt.at(j),vec_proj_bin_pt.at(j+1), vec_rebin_mass.at(0)));
              cSoBEtaCanvas->cd(j+1);
              hSignalOverBackground->Draw();
              // pT_Intervall->Draw("same");
              legendInfo    -> Draw("same");
              cSignificanceEtaCanvas->cd(j+1);
              hSignificance->Draw();
              // pT_Intervall->Draw("same");
              legendInfo    -> Draw("same");

            }
          }
                                                                                if(DoDebug) Printf("%d", __LINE__);

          // // ------------------------
          // // Plot with pion and eta peak
          // cSignalCanvas->cd();
          // projX -> Draw();
          // if(DoFit == kTRUE){
          //   // TF1 *fFullFit = new TF1("fFullFit","pol2+gaus(3)+gaus(6)",0.0,1.0);
          //   // TF1 *fPionFit = new TF1("fPionFit","pol2+gaus(3)",0.04,0.3);
          //   // TF1 *fEtaFit = new TF1("fEtaFit","pol2+gaus(3)",0.3,0.8);
          //   TF1 *fPionFit = new TF1("fPionFit","pol2+(x<[4])*([3]*(TMath::Exp(-0.5*((x-[4])/[5])^2)+TMath::Exp((x-[4])/[6])*(1.-TMath::Exp(-0.5*((x-[4])/[5])^2))))+(x>=[4])*([3]*TMath::Exp(-0.5*((x-[4])/[5])^2))",0.04,0.3);
          //   TF1 *fEtaFit = new TF1("fEtaFit","pol2+(x<[4])*([3]*(TMath::Exp(-0.5*((x-[4])/[5])^2)+TMath::Exp((x-[4])/[6])*(1.-TMath::Exp(-0.5*((x-[4])/[5])^2))))+(x>=[4])*([3]*TMath::Exp(-0.5*((x-[4])/[5])^2))",0.3,0.8);
          //
          //
          //   SetFitSettingsPion(fPionFit);
          //   SetFitSettingsEta(fEtaFit);
          //
          //   // fFullFit -> SetParameters(0,0.018,0.01,projX->GetMaximum()/2,0.135,0.02,projX->GetMaximum()/2,0.549,0.025);
          //   // fPionFit -> SetParameters(0,0.018,0.01,projX->GetMaximum()/2,0.135,0.025);
          //   // fEtaFit  -> SetParameters(0,0.018,0.01,projX->GetMaximum()/2,0.549,0.025);
          //   fPionFit->SetParameters(0,0.018,0.01,projX->GetMaximum()/2,0.135,0.025,0.4);
          //   fEtaFit ->SetParameters(0,0.018,0.01,projX->GetMaximum()/2,0.547,0.01 ,0.4);
          //   // fFullFit -> SetParLimits(4,0.1,0.15);
          //   // fFullFit -> SetParLimits(7,0.54,0.56);
          //   // fFullFit -> SetParLimits(8,0.001,0.01);
          //   fPionFit->SetParLimits(4,0.13,0.14);  // (*)
          //   fEtaFit->SetParLimits(4,0.53,0.565);  // (*)
          //   fPionFit->SetParLimits(5,0.0,0.02);   // (*)
          //   fEtaFit->SetParLimits(5,0.0,0.02);    // (*)
          //
          //   // projX    -> Fit("fFullFit","QRMNE0");
          //   projX    -> Fit("fPionFit","QRMNE0");
          //   projX    -> Fit("fPionFit","QRMNE0");
          //   projX    -> Fit("fPionFit","QRMNE0");
          //   projX    -> Fit("fEtaFit","QRMNE0");
          //   projX    -> Fit("fEtaFit","QRMNE0");
          //   projX    -> Fit("fEtaFit","QRMNE0");
          //   // fFullFit -> Draw("same");
          //   fPionFit -> Draw("same");
          //   fEtaFit  -> Draw("same");
          // }
          // pT_Intervall  -> Draw("same");
          // // mass_Intervall -> Draw("same");
          // cSignalCanvas -> SaveAs(Form("%s%s.pdf",DocumentPathProjX.Data(), projX->GetName()));



        } // end of DoMassProj
                                                                                if(DoDebug) Printf("%d", __LINE__);

        // Draw Pt projections
        if (DoPtProj == kTRUE) {
          projY -> Draw();
          // pT_Intervall  -> Draw("same");
          legendInfo    -> Draw("same");
          mass_Intervall -> Draw("same");
          cSignalCanvas -> SaveAs(Form("%s%s.pdf",DocumentPathProjY.Data(), projY->GetName()));
        }
                                                                                // if(DoDebug) Printf("%d", __LINE__);
      } // end loop over all mass projection intervalls
                                                                                // if(DoDebug) Printf("%d", __LINE__);
    } // end loop over all pt projection intervalls
                                                                                // if(DoDebug) Printf("%d", __LINE__);


    if (DoMassProj == kTRUE && ExtFourPairPionSig == kTRUE) {
      cProjPtPionCanvas       -> SaveAs(Form("%s%s_Pion_rebin%i_%ipT-Intervalls.pdf",DocumentPathProjX.Data()  , "AllPtProjections"        , vec_rebin_mass.at(0), (unsigned int) vec_proj_bin_pt.size()-1));
      cSoBPionCanvas          -> SaveAs(Form("%s%s_Pion_rebin%i_%ipT-Intervalls.pdf",DocumentPathSignals.Data(), "AllSignalOverBackground" , vec_rebin_mass.at(0), (unsigned int) vec_proj_bin_pt.size()-1));
      cSignificancePionCanvas -> SaveAs(Form("%s%s_Pion_rebin%i_%ipT-Intervalls.pdf",DocumentPathSignals.Data(), "AllSignificances"        , vec_rebin_mass.at(0), (unsigned int) vec_proj_bin_pt.size()-1));
    }
    if (DoMassProj == kTRUE && ExtFourPairEtaSig == kTRUE) {
      cProjPtEtaCanvas        -> SaveAs(Form("%s%s_Eta_rebin%i_%ipT-Intervalls.pdf",DocumentPathProjX.Data()  , "AllPtProjections"         , vec_rebin_mass.at(0), (unsigned int) vec_proj_bin_pt.size()-1));
      cSoBEtaCanvas           -> SaveAs(Form("%s%s_Eta_rebin%i_%ipT-Intervalls.pdf",DocumentPathSignals.Data(), "AllSignalOverBackground"  , vec_rebin_mass.at(0), (unsigned int) vec_proj_bin_pt.size()-1));
      cSignificanceEtaCanvas  -> SaveAs(Form("%s%s_Eta_rebin%i_%ipT-Intervalls.pdf",DocumentPathSignals.Data(), "AllSignificances"         , vec_rebin_mass.at(0), (unsigned int) vec_proj_bin_pt.size()-1));
    }

                                                                                if(DoDebug) Printf("%d", __LINE__);

    // if (DoSignalOverBackground == kTRUE) {
    //   // Draw Signal over Background ratio for the pion and eta signal
    //   cSignalCanvas->cd();
    //   for (unsigned int iBin = 1; iBin < vec_proj_bin_pt.size(); iBin++) {
    //     if(ExtFourPairPionSig == kTRUE)  hSignalOverBackgroundPtPion->SetBinContent(iBin, Vec_SigBack_Pion.at(iBin-1));
    //     if(ExtFourPairEtaSig == kTRUE) hSignalOverBackgroundPtEta ->SetBinContent(iBin, Vec_SigBack_Eta.at(iBin-1));
    //   }
    //
    //   auto legend = new TLegend(0.9,0.85,0.75,0.9);
    //   TLegendEntry *entry1=legend->AddEntry("S/B Eta","Eta Signal","l");        // Legend with Line, Eta  S/B
    //   TLegendEntry *entry2=legend->AddEntry("S/B Pion","Pion Signal","l");      // Legend with Line, Pion S/B
    //   entry1->SetLineColor(kBlue); // blue
    //   entry2->SetLineColor(kRed); // red
    //
    //                                                                               if(DoDebug) Printf("%d", __LINE__);
    //   hSignalOverBackgroundPtEta->Draw("");
    //   legend->Draw("same");
    //   cSignalCanvas -> SaveAs(Form("%s%s.pdf",DocumentPathSignals.Data(), "SoB_Ratio_Eta"));
    //   hSignalOverBackgroundPtPion->Draw("");
    //   legend->Draw("same");
    //   cSignalCanvas -> SaveAs(Form("%s%s.pdf",DocumentPathSignals.Data(), "SoB_Ratio_Pion"));
    //
    //                                                                               if(DoDebug) Printf("%d", __LINE__);
    //   // // Setting Y axis range
    //   // double fMinHist;
    //   // double fMaxHist;
    //   // if (hSignalOverBackgroundPtPion->GetMaximum() >= hSignalOverBackgroundPtEta->GetMaximum()) {fMaxHist = hSignalOverBackgroundPtPion->GetMaximum()*1.2;}
    //   // else {fMaxHist = hSignalOverBackgroundPtEta->GetMaximum()*1.2;}
    //   // if (hSignalOverBackgroundPtPion->GetMinimum() <= hSignalOverBackgroundPtEta->GetMinimum()) {fMinHist = hSignalOverBackgroundPtPion->GetMinimum()*1,2;}
    //   // else {fMinHist = hSignalOverBackgroundPtEta->GetMinimum()*1,2;}
    //   // hSignalOverBackgroundPtPion->GetYaxis()->SetRangeUser(fMinHist,fMaxHist);
    //
    //   hSignalOverBackgroundPtPion->Draw("");
    //   hSignalOverBackgroundPtEta->Draw("same");
    //   legend->Draw("same");
    //   cSignalCanvas -> SaveAs(Form("%s%s.pdf",DocumentPathSignals.Data(), "SoB_Ratio_Pion_Eta"));
    //                                                                               if(DoDebug) Printf("%d", __LINE__);
    // } // end of DoSignalOverBackground
  } // end loop over all Histograms
  delete cSignalCanvas;
  delete cProjPtPionCanvas;
  delete cProjPtEtaCanvas;
  delete cSignalPionCanvas;
  delete cSignalEtaCanvas;
  delete cSignificancePionCanvas;
  delete cSignificanceEtaCanvas;
  delete cSoBPionCanvas;
  delete cSoBEtaCanvas;
} // end of DrawProjection void








void PlotPrimarySoBRatio(TString VecRecPairHistoNames, std::vector<TList*> PairRecPrimList, std::vector<Int_t> vec_rebin_pt, std::vector<Int_t> vec_rebin_mass, TString PairCase) {
                                                                                if(DoDebug)  std::cout << "Start plotting of S/B ratio for primary pairs" << std::endl;
  TCanvas *cSoBCanvas = new TCanvas("cSoBCanvas","",800,800);
  SetCanvasStandardSettings(cSoBCanvas);
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  cSoBCanvas->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
  pad2->SetTopMargin(0);
  pad2->Draw();
  pad2->cd();

  pad1->SetTopMargin(pad1->GetTopMargin()*0.5);
  // pad1->SetBottomMargin(pad1->GetBottomMargin()*4.0);
  pad2->SetTopMargin(pad2->GetTopMargin()*2.0);
  pad2->SetBottomMargin(pad2->GetBottomMargin()*2.0);
  // Printf("Pad1 Top Margin: %d", pad1->GetTopMargin());
  // Printf("Pad1 Bottom Margin: %d", pad1->GetBottomMargin());
  // Printf("Pad2 Top Margin: %d", pad2->GetTopMargin());
  // Printf("Pad2 Bottom Margin: %d", pad2->GetBottomMargin());

  for (size_t iRecPairList = 0; iRecPairList < PairRecPrimList.size(); iRecPairList++) {
    TList* PairRecList = PairRecPrimList.at(iRecPairList);
    TObjArray* arrPrimaryPairHistoNames=VecRecPairHistoNames.Tokenize(";");
    const Int_t nPairTH2hist=arrPrimaryPairHistoNames->GetEntriesFast();
    TString DocumentPath = Form("Plots/%s/%s/%s/%s/",TrainNumber.Data(),PairCase.Data(),"SoB",PairRecList->GetName());
    gSystem->Exec(Form("mkdir -p %s",DocumentPath.Data()));

    // TH2 Histograms
    TH2D* hFullULSHist;
    TH2D* hEtaSignal;
    TH2D* hPionSignal;
    TH2D* hULSminusSignals;
    // TH1 Projections
    TH1D* hProjFullULSHist;
    TH1D* hProjEtaSignal;
    TH1D* hProjPionSignal;
    TH1D* hProjULSminusSignals;


    for (Int_t i = 0; i < nPairTH2hist; i++) {
      TString temp = arrPrimaryPairHistoNames->At(i)->GetName();
      // selecting specific Histograms
      if (temp == "Nrec_anyPair_anyPart_UndefinedMother_finalstate") hFullULSHist       = (dynamic_cast<TH2D*>(PairRecList->FindObject(temp)));
      if (temp == "Nrec_anyPair_anyPart_UndefinedMother_finalstate") hULSminusSignals   = (dynamic_cast<TH2D*>(PairRecList->FindObject(temp)));
      if (temp == "Nrec_elePair_sameMother_eta_finalstate" )         hEtaSignal         = (dynamic_cast<TH2D*>(PairRecList->FindObject(temp)));
      if (temp == "Nrec_elePair_sameMother_pion_finalstate")         hPionSignal        = (dynamic_cast<TH2D*>(PairRecList->FindObject(temp)));
    }

    Double_t fProjStartBin = hFullULSHist->GetYaxis()->FindBin(0.);
    Double_t fProjEndBin   = hFullULSHist->GetYaxis()->FindBin(4.);

    hProjFullULSHist     = hFullULSHist     -> ProjectionX("ProjXFullULSHist",     fProjStartBin, fProjEndBin);
    hProjEtaSignal       = hEtaSignal       -> ProjectionX("ProjXEtaSignal",       fProjStartBin, fProjEndBin);
    hProjPionSignal      = hPionSignal      -> ProjectionX("ProjXPionSignal",      fProjStartBin, fProjEndBin);
    hProjULSminusSignals = hULSminusSignals -> ProjectionX("ProjXULSminusSignals", fProjStartBin, fProjEndBin);

    hProjULSminusSignals->Add(hProjEtaSignal,-1);
    // hProjULSminusSignals->Add(hProjPionSignal,-1);

    SetHistoStandardSettings(hProjFullULSHist);
    SetHistoStandardSettingsBlue(hProjEtaSignal);
    SetHistoStandardSettingsRed(hProjPionSignal);
    SetHistoStandardSettingsGreen(hProjULSminusSignals);

    hProjEtaSignal->GetXaxis()->SetTitleOffset(hProjEtaSignal->GetXaxis()->GetTitleOffset()*2.0);

    hProjFullULSHist     -> RebinX(vec_rebin_mass.at(0));
    hProjEtaSignal       -> RebinX(vec_rebin_mass.at(0));
    hProjPionSignal      -> RebinX(vec_rebin_mass.at(0));
    hProjULSminusSignals -> RebinX(vec_rebin_mass.at(0));

    hProjPionSignal->GetYaxis()->SetTitle("Counts");
    hProjEtaSignal->GetYaxis()->SetTitle("Counts");
    hProjFullULSHist->GetYaxis()->SetTitle("Counts");

    // Set Legend + Legend settings
    auto legend = new TLegend(0.35,0.75,0.85,0.92);
    TLegendEntry *entry1=legend->AddEntry("name","Sum","p");        // Legend with Line, same Mother
    TLegendEntry *entry2=legend->AddEntry("name","#pi^{0}-Signal, primary e^{+}e^{-} pair from #pi^{0}-Dalitz decay","p");   // Legend with Line, different Mother
    TLegendEntry *entry3=legend->AddEntry("name","#eta-Signal, primary e^{+}e^{-} pair from #eta-Dalitz decay","p");   // Legend with Line, different Mother
    TLegendEntry *entry4=legend->AddEntry("name","Background, (Sum - #eta-Signal)","p");   // Legend with Line, different Mother
    entry1->SetMarkerColor(kBlack);     entry1->SetMarkerStyle(20);
    entry2->SetMarkerColor(kRed);       entry2->SetMarkerStyle(27);
    entry3->SetMarkerColor(kBlue);      entry3->SetMarkerStyle(2);
    entry4->SetMarkerColor(kGreen+2);   entry4->SetMarkerStyle(4);
    legend->SetBorderSize(0);
    legend->SetFillColorAlpha(0, 0.0);
    legend->SetTextSize(0.03);

    auto legendInfo = new TLegend(0.095,0.75,0.42,0.92);
    TLegendEntry *entry5=legendInfo->AddEntry("collisionSystem"   ,"pp, #sqrt{s} = 13 TeV, |#eta_{e}| < 0.8","");
    TLegendEntry *entry6=legendInfo->AddEntry("SinglePtPrim","#font[12]{p}_{T,e}^{prim} > 0.075 GeV/#font[12]{c}","");
    TLegendEntry *entry7=legendInfo->AddEntry("SinglePtSec","#font[12]{p}_{T,ee}^{prim} > 0.075 GeV/#font[12]{c}","");
    legendInfo->SetBorderSize(0);
    legendInfo->SetFillColorAlpha(0, 0.0);
    legendInfo->SetTextSize(0.03);

    pad1->cd();
    pad1->SetLogy();
    hProjPionSignal->Draw("E1");
    hProjULSminusSignals->Draw("same E1");
    hProjEtaSignal->DrawClone("same E1");
    hProjFullULSHist->Draw("same E1");
    legend->Draw("same");
    legendInfo->Draw("same");

    pad2->cd();
    hProjEtaSignal->Divide(hProjULSminusSignals);
    // hProjEtaSignal->SetMarkerStyle(2);
    hProjEtaSignal->GetYaxis()->SetTitle("Ratio S/B");
    hProjEtaSignal->Draw("ep E1");

    cSoBCanvas->SaveAs(Form("%sPrimaryPairSoB.pdf",DocumentPath.Data()));
  }
  delete cSoBCanvas;
}








void MergeMassCutHistos(TList* PairList, TString VecMassCutHistoNames, TString VecMassCutHistoTitles, TString PairCase, double MassCut, Bool_t DoCutEff){
  // Create path to save output
  TString DocumentPath = Form("Plots/%s/%s/MassCuts/",TrainNumber.Data(),PairCase.Data());
  gSystem->Exec(Form("mkdir -p %s",DocumentPath.Data()));

  //Initialise Canvas + setting canvas settings
  TCanvas *cSignalCanvas = new TCanvas("cSignalCanvas","",800,800);
  TCanvas *cLogYCanvas = new TCanvas("cLogYCanvas","",800,800);
  SetCanvasStandardSettings(cSignalCanvas);
  SetCanvasStandardSettings(cLogYCanvas);

  TObjArray* arrMassCutHistoNames=VecMassCutHistoNames.Tokenize(";");
  TObjArray* arrHistTitles = VecMassCutHistoTitles.Tokenize(";");
  // Number of Histograms that become compared
  const Int_t nMassCutTH2hist=arrMassCutHistoNames->GetEntriesFast();
  std::cout << nMassCutTH2hist << std::endl;

  for (Int_t i = 0; i < nMassCutTH2hist; i+=2) {  // i+=2 caused by the order of names in arrMassCutHistoNames. (fistSignal: sameMother + differnetMother, secondSignal:  ...)
    TString temp1 = arrMassCutHistoNames->At(i)->GetName();    // using sameMotherHisto
    TString temp2 = arrMassCutHistoNames->At(i+1)->GetName();  // using DifferentMotherHisto
    std::cout << temp1 <<" " << temp2 << std::endl;
    // select both Histos from PairList out of File
    TH2D hTempSameMotherHist     = *(dynamic_cast<TH2D*>(PairList->FindObject(temp1)));
    TH2D hTempDiffMotherHist     = *(dynamic_cast<TH2D*>(PairList->FindObject(temp2)));
    hTempSameMotherHist.Sumw2();
    hTempDiffMotherHist.Sumw2();

    // Doing the projection onto x Axis
    TH1D* hXProj_SameMotherHist  = hTempSameMotherHist.ProjectionX();
    TH1D* hXProj_DiffMotherHist  = hTempDiffMotherHist.ProjectionX();
    hXProj_SameMotherHist->Sumw2();
    hXProj_DiffMotherHist->Sumw2();

    SetHistoStandardSettings(hXProj_SameMotherHist);
    SetHistoStandardSettingsRed(hXProj_DiffMotherHist);

    // Set new range on x-Axis for the projection Hist
    hXProj_SameMotherHist->SetAxisRange(0,0.8,"x");
    hXProj_DiffMotherHist->SetAxisRange(0,0.8,"x");

    TString name = arrMassCutHistoNames->At(i)->GetName();
    // TString status = "";
    // if (name.Contains("finalstate")) {status = "Primaries";}
    // if (name.Contains("secondary")) {status = "Secondaries";}

    // Set Legend + Legend settings
    auto legend = new TLegend(0.9,0.85,0.75,0.9);
    TLegendEntry *entry1=legend->AddEntry("name","same Mother","l");        // Legend with Line, same Mother
    TLegendEntry *entry2=legend->AddEntry("name","different Mother","l");   // Legend with Line, different Mother
    entry1->SetLineColor(4); // blue
    entry2->SetLineColor(2); // red

    // Giving the Projection Histogram a title
    hXProj_SameMotherHist->SetTitle(Form("%s", arrHistTitles->At(i)->GetName()));

    // Draw a line at MassCut
    double cutPosition = MassCut;
    double fMinHist;
    double fMaxHist;


    // Setting Y axis range and use same range for a vertical cut line
    if (hXProj_SameMotherHist->GetMaximum() >= hXProj_DiffMotherHist->GetMaximum()) {fMaxHist = hXProj_SameMotherHist->GetMaximum()*1.2;}
    else {fMaxHist = hXProj_DiffMotherHist->GetMaximum()*1.2;}
    if (hXProj_SameMotherHist->GetMinimum() <= hXProj_DiffMotherHist->GetMinimum()) {fMinHist = hXProj_SameMotherHist->GetMinimum();}
    else {fMinHist = hXProj_DiffMotherHist->GetMinimum();}
    hXProj_SameMotherHist->GetYaxis()->SetRangeUser(fMinHist,fMaxHist);
    TLine cutLine = TLine(cutPosition, fMinHist, cutPosition, fMaxHist);
    cutLine.SetLineColor(kGreen+2);
    cutLine.SetLineWidth(2);
    // Select linear canvas and Draw both Projections into same plot
    cSignalCanvas->cd();
    hXProj_SameMotherHist->Draw("");
    hXProj_DiffMotherHist->Draw("same");
    cutLine.Draw();
    legend->Draw();
    cSignalCanvas->SaveAs(Form("%sSameDiffernetMother%s.pdf",DocumentPath.Data(),name.Data()));
    // cSignalCanvas->SaveAs(Form("%sSameDiffernetMother%s.png",DocumentPath.Data(),name.Data()));


    // Setting Y axis range and use same range for a vertical cut line
    if (hXProj_SameMotherHist->GetMaximum() >= hXProj_DiffMotherHist->GetMaximum()) {fMaxHist = hXProj_SameMotherHist->GetMaximum()*10;}
    else {fMaxHist = hXProj_DiffMotherHist->GetMaximum()*10;}
    if (hXProj_SameMotherHist->GetMinimum() <= hXProj_DiffMotherHist->GetMinimum()) {fMinHist = hXProj_SameMotherHist->GetMinimum()+0.1;}
    else {fMinHist = hXProj_DiffMotherHist->GetMinimum()+0.1;}
    hXProj_SameMotherHist->GetYaxis()->SetRangeUser(fMinHist,fMaxHist);
    cutLine = TLine(cutPosition, fMinHist, cutPosition, fMaxHist);
    cutLine.SetLineColor(kGreen+2);
    cutLine.SetLineWidth(2);
    // Select log y canvas and Draw both Projections into same plot
    cLogYCanvas->cd();
    hXProj_SameMotherHist->Draw("");
    hXProj_DiffMotherHist->Draw("same");
    cutLine.Draw();
    legend->Draw();
    cLogYCanvas->SetLogy();
    cLogYCanvas->SaveAs(Form("%sLogYSameDiffernetMother%s.pdf",DocumentPath.Data(),name.Data()));
    // cLogYCanvas->SaveAs(Form("%sLogYSameDiffernetMother%s.png",DocumentPath.Data(),name.Data()));

    if (DoCutEff == kTRUE) {
      // Initialise Histogram
      TH1D* hSameMotherCutEfficiency     = new TH1D("Same Mother Cut Efficiency","Cut Efficiency ;Mass Cut (GeV/c^{2});Efficiency",10,0.,0.1);//,AliDielectronVarManager::kPt);
      TH1D* hDiffMotherCutEfficiency     = new TH1D("Differnt Mother Cut Efficiency","Cut Efficiency ;Mass Cut (GeV/c^{2});Efficiency",10,0.,0.1);//,AliDielectronVarManager::kPt);
      TH1D* hDiffToSameMotherCutEfficiency     = new TH1D("Differnt Mother to Same Mother Cut Efficiency","Cut Efficiency ;Mass Cut [GeV/c^{2}];Efficiency",10,0.,0.1);//,AliDielectronVarManager::kPt);
      TH1D* hContaminationDiffPairs      = new TH1D("Contamination of Pairs from different Mother Pairs","Contamination of Pairs from different Mother Pairs ;Mass Cut [GeV/c^{2}];Contamination",10,0.,0.1);//,AliDielectronVarManager::kPt);
      hSameMotherCutEfficiency->Sumw2();
      hDiffMotherCutEfficiency->Sumw2();
      hDiffToSameMotherCutEfficiency->Sumw2();
      hContaminationDiffPairs->Sumw2();
      SetHistoStandardSettings(hSameMotherCutEfficiency);
      SetHistoStandardSettingsRed(hDiffMotherCutEfficiency);
      SetHistoStandardSettings(hDiffToSameMotherCutEfficiency);
      SetHistoStandardSettings(hContaminationDiffPairs);
      hDiffToSameMotherCutEfficiency->SetLineColor(kOrange-3);
      hSameMotherCutEfficiency->GetYaxis()->SetRangeUser(0.0,1.1);
      hContaminationDiffPairs ->GetYaxis()->SetRangeUser(0.0,1.1);

      std::vector<Double_t> NMassCuts = {0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1};
      for (unsigned int iCut = 1; iCut < NMassCuts.size(); iCut++) {
        Double_t NTotalSameMother     = hXProj_SameMotherHist->Integral();       // Total number of pairs in same Mother Histogram
        Double_t NTotalDiffMother     = hXProj_DiffMotherHist->Integral();       // Total number of pairs in different Mother Histogram
        Double_t NInsideCutSameMother = hXProj_SameMotherHist->Integral(hXProj_SameMotherHist->GetXaxis()->FindBin(NMassCuts.at(0)),hXProj_SameMotherHist->GetXaxis()->FindBin(NMassCuts.at(iCut)));       // Total number of pairs inside the cut with same Mother
        Double_t NInsideCutDiffMother = hXProj_DiffMotherHist->Integral(hXProj_DiffMotherHist->GetXaxis()->FindBin(NMassCuts.at(0)),hXProj_DiffMotherHist->GetXaxis()->FindBin(NMassCuts.at(iCut)));       // Total number of pairs inside the cut with different Mother

        Double_t sameMotherEff = NInsideCutSameMother/NTotalSameMother;
        // Double_t diffMotherEff = NInsideCutDiffMother/NTotalDiffMother;
        Double_t diffMotherEff = NInsideCutDiffMother/NTotalDiffMother;
        Double_t diff_sameMotherEff = NInsideCutDiffMother/NTotalSameMother;
        Double_t Contamination = NInsideCutDiffMother/(NInsideCutSameMother+NInsideCutDiffMother);

        std::cout << " NTotalSameMother     = " << NTotalSameMother <<     ", NTotalDiffMother     = " << NTotalDiffMother << std::endl;
        std::cout << " NInsideCutSameMother = " << NInsideCutSameMother << ", NInsideCutDiffMother = " << NInsideCutDiffMother << std::endl;
        std::cout << " sameMotherEff        = " << sameMotherEff <<        ", diffMotherEff        = " << diffMotherEff << std::endl << std::endl;


        hSameMotherCutEfficiency->SetBinContent(iCut,sameMotherEff);
        hDiffMotherCutEfficiency->SetBinContent(iCut,diffMotherEff);
        hDiffToSameMotherCutEfficiency->SetBinContent(iCut,diff_sameMotherEff);
        hContaminationDiffPairs ->SetBinContent(iCut,Contamination);
        // std::cout << hSameMotherCutEfficiency->GetBinContent(iCut) << std::endl;
        // std::cout << hDiffMotherCutEfficiency->GetBinContent(iCut) << std::endl;
      }

      // Set Legend + Legend settings
      auto legend = new TLegend(0.4,0.4,0.9,0.6);
      TLegendEntry *entry1=legend->AddEntry("name","#frac{same Mother Pairs inside cut}{Total number same Mother pairs}","l");        // Legend with Line, same Mother
      TLegendEntry *entry2=legend->AddEntry("name","#frac{diff Mother Pairs inside cut}{Total number diff Mother pairs}","l");   // Legend with Line, different Mother
      TLegendEntry *entry3=legend->AddEntry("name","#frac{diff Mother Pairs inside cut}{Total number same Mother pairs}","l");   // Legend with Line, different Mother
      entry1->SetLineColor(kBlue); // blue
      entry2->SetLineColor(kRed); // red
      entry3->SetLineColor(kOrange-3);

      cSignalCanvas->cd();
      hSameMotherCutEfficiency->Draw("c");
      hDiffMotherCutEfficiency->Draw("csame");
      hDiffToSameMotherCutEfficiency->Draw("csame");
      legend->Draw();
      cSignalCanvas->SaveAs(Form("%sCutEfficiency.pdf",DocumentPath.Data()));

      // Set Legend + Legend settings
      auto legendCont = new TLegend(0.15,0.85,0.45,0.9);
      TLegendEntry *entryCont=legendCont->AddEntry("name","#frac{diff Mother Pairs inside cut}{Total number of Mother pairs inside Cut}","l");        // Legend with Line, same Mother
      entryCont->SetLineColor(4); // blue

      cSignalCanvas->cd();
      hContaminationDiffPairs->Draw("c");
      legendCont->Draw();
      cSignalCanvas->SaveAs(Form("%sContaminationOfDiffMotherPairs.pdf",DocumentPath.Data()));
    }
  }
}

//_______________________________________________________________________________________________
Bool_t Rebin2DHistogram(TH2D& hIn, std::vector<Double_t> fBinsRebin_Mee, std::vector<Double_t> fBinsRebin_Ptee)
{
  /// Rebin a 2D histogram based on two bin vectors with variable binning.
  /// If one of the binning vectors in this class is not sufficiently filled (<2 entries), there are 2 options:
  /// (1) if the vector has size 1 and the stored number is an integer and an exact divisor of the original number
  ///     of bins, this number is taken as rebin factor.
  /// (2) if (1) is not fulfilled, the original binning of the respective axis is kept.
  ///
  TString sname(hIn.GetName());  // Otherwise potential memory leak


  if (fBinsRebin_Mee.size()<2) {
    // Int_t n_bins_temp = hIn.GetNbinsX();
    // Double_t bins_temp[n_bins_temp+1];
    // //    cout << " sizeof(bins_temp)/sizeof(*bins_temp) = " << sizeof(bins_temp)/sizeof(*bins_temp) << endl;
    // //    cout << " hIn.GetXaxis()->GetBinLowEdge(n_bins_temp)   = " << hIn.GetXaxis()->GetBinLowEdge(n_bins_temp) << endl;
    // //    cout << " hIn.GetXaxis()->GetBinLowEdge(n_bins_temp+1) = " << hIn.GetXaxis()->GetBinLowEdge(n_bins_temp+1) << endl;
    // hIn.GetXaxis()->GetLowEdge(bins_temp); // Fills the binning of the original histogram into the array. // problems to fill the last bin! - why!?!
    // //    cout << " bins_temp[1] = " << bins_temp[1] << " bins_temp[99] = " << bins_temp[99] << " bins_temp[100] = " << bins_temp[100] << endl;
    // bins_temp[n_bins_temp] = hIn.GetXaxis()->GetBinLowEdge(n_bins_temp+1); // fill last bin edge manually...
    // //    cout << " bins_temp[1] = " << bins_temp[1] << " bins_temp[99] = " << bins_temp[99] << " bins_temp[100] = " << bins_temp[100] << endl;
    // Int_t n_rebin=1;
    // if (   fBinsRebin_Mee.size()==1 && fBinsRebin_Mee[0]!=0.
    //     && fBinsRebin_Mee[0]==Int_t(fBinsRebin_Mee[0])
    //     && n_bins_temp%Int_t(fBinsRebin_Mee[0])==0 ) {
    //   n_rebin = fBinsRebin_Mee[0];
    //   // Info(Form("LmRebinner: rebin X axis by factor %i.", n_rebin));
    // }
    // // else Error(Form("LmRebinner: rebinning X axis with factor %f failed", fBinsRebin_Mee[0]));
    // SetRebinningX(n_bins_temp, bins_temp, n_rebin,fBinsRebin_Mee);

    std::cout <<  " x-Axis NBins of Hist < 2 " << std::endl;
    std::cout <<  " ####################### END PROGRAM ####################### " << std::endl;
    return 0;
  }
  if (fBinsRebin_Ptee.size()<2) {
    // Int_t n_bins_temp = hIn.GetNbinsY();
    // Double_t bins_temp[n_bins_temp+1];
    // hIn.GetYaxis()->GetLowEdge(bins_temp); // Fills the binning of the original histogram into the array. // problems to fill the last bin! - why!?!
    // bins_temp[n_bins_temp] = hIn.GetYaxis()->GetBinLowEdge(n_bins_temp+1); // fill last bin edge manually...
    // Int_t n_rebin=1;
    // if (   fBinsRebin_Ptee.size()==1 && fBinsRebin_Ptee[0]!=0.
    //     && fBinsRebin_Ptee[0]==Int_t(fBinsRebin_Ptee[0])
    //     && n_bins_temp%Int_t(fBinsRebin_Ptee[0])==0 ) {
    //   n_rebin = fBinsRebin_Ptee[0];
    //   // Info(Form("LmRebinner: rebin Y axis by factor %i.", n_rebin));
    // }
    // // else Error(Form("LmRebinner: rebinning Y axis with factor %f failed", fBinsRebin_Ptee[0]));
    // SetRebinningY(n_bins_temp, bins_temp, n_rebin,fBinsRebin_Ptee);

    std::cout <<  " y-Axis NBins of Hist < 2 " << std::endl;
    std::cout <<  " ####################### END PROGRAM ####################### " << std::endl;
    return 0;
  }

  TH2D hOut("hOut", hIn.GetTitle(), fBinsRebin_Mee.size()-1, fBinsRebin_Mee.data(), fBinsRebin_Ptee.size()-1, fBinsRebin_Ptee.data());
  hOut.Sumw2();

  hOut.GetXaxis()->SetTitle(hIn.GetXaxis()->GetTitle());
  hOut.GetYaxis()->SetTitle(hIn.GetYaxis()->GetTitle());
  hOut.GetZaxis()->SetTitle(hIn.GetZaxis()->GetTitle());
  Double_t bin_content_hIn     =  0;
  Double_t bin_error_hIn       =  0;
  Double_t bin_center_mee_hIn  = -1;
  Double_t bin_center_ptee_hIn = -1;

  Double_t bin_content_hOut    = -1;
  Double_t bin_error_hOut_2    = -1;
  Int_t    bin_mee_hOut        = -1;
  Int_t    bin_ptee_hOut       = -1;

  // The error in this loop is quadratically added. This means that after this loop you have to sqrt the error bin-by-bin
  // to get the real uncertainty
  for (Int_t i_mee = 0; i_mee <= hIn.GetNbinsX()+1; ++i_mee){
    for (Int_t j_ptee = 0; j_ptee <= hIn.GetNbinsY()+1; ++j_ptee){
      bin_content_hIn = hIn.GetBinContent(i_mee, j_ptee); // input bin content
      bin_error_hIn   = hIn.GetBinError(i_mee, j_ptee);   // input bin error

      bin_center_mee_hIn  = hIn.GetXaxis()->GetBinCenter(i_mee);
      bin_center_ptee_hIn = hIn.GetYaxis()->GetBinCenter(j_ptee);

      bin_mee_hOut    = hOut.GetXaxis()->FindBin(bin_center_mee_hIn);
      bin_ptee_hOut   = hOut.GetYaxis()->FindBin(bin_center_ptee_hIn);

      bin_content_hOut = hOut.GetBinContent(bin_mee_hOut, bin_ptee_hOut);
      bin_error_hOut_2 = hOut.GetBinError(bin_mee_hOut, bin_ptee_hOut);

      // set hOut content
      hOut.SetBinContent(bin_mee_hOut, bin_ptee_hOut, bin_content_hOut + bin_content_hIn);
      hOut.SetBinError(bin_mee_hOut, bin_ptee_hOut, bin_error_hOut_2 + bin_error_hIn * bin_error_hIn);
    }
  }

  // sqrt of errors bin-by-bin
  for (Int_t i_mee = 0; i_mee <= hOut.GetNbinsX()+1; ++i_mee){
    for (Int_t j_ptee = 0; j_ptee <= hOut.GetNbinsY()+1; ++j_ptee){
      hOut.SetBinError(i_mee, j_ptee, TMath::Sqrt(hOut.GetBinError(i_mee, j_ptee)));
    }
  }

  hIn = hOut;
  hIn.SetName(sname.Data());


  return 1;
}

//_______________________________________________________________________________________________
void SetRebinningX(Int_t n_bins, const Double_t* bins, Int_t n_rebin, std::vector<Double_t> fBinsRebin)
{
  fBinsRebin.clear(); // to be able to overwrite
  for (int i=0; i<n_bins+1; i+=n_rebin) { fBinsRebin.push_back(bins[i]); }
  // SetEnableRebinning();
}

//_______________________________________________________________________________________________
void SetRebinningY(Int_t n_bins, const Double_t* bins, Int_t n_rebin, std::vector<Double_t> fBinsRebin)
{
  fBinsRebin.clear(); // to be able to overwrite
  for (int i=0; i<n_bins+1; i+=n_rebin) { fBinsRebin.push_back(bins[i]); }
  // SetEnableRebinning();
}

// void CompareGenSmearRec(TString VecGenSmearedSingleHistoNames, TString VecRecSingleHistoNames, TList* SingleSecGenSmearedList, TList* SingleSecRecList){
//   TString DocumentPathProj = Form("Plots/%s/%s/PtEtaPhiProj/",TrainNumber.Data(),"SingleElectrons");    // define path to pt projections
//   gSystem->Exec(Form("mkdir -p %s",DocumentPathProj.Data()));
//
//   TCanvas *cSignalCanvas = new TCanvas("cSignalCanvas","",800,800);
//   TCanvas *cLogYCanvas = new TCanvas("cLogYCanvas","",800,800);
//   SetCanvasStandardSettings(cSignalCanvas);
//   SetCanvasStandardSettings(cLogYCanvas);
//   cLogYCanvas->SetLogy();
//
//
//   TObjArray* arrGenSmearHistoNames=VecGenSmearedSingleHistoNames.Tokenize(";");
//   TObjArray* arrRecHistoNames=VecRecSingleHistoNames.Tokenize(";");
//   const Int_t nSingleTH3hist=arrGenSmearHistoNames->GetEntriesFast();
//   for (Int_t i = 0; i < nSingleTH3hist; i++) {
//     TString temp1 = arrGenSmearHistoNames->At(i)->GetName();
//     TString temp2 = arrRecHistoNames->At(i)->GetName();
//
//     // TH3D hTempGenSingleHist          = *(dynamic_cast<TH3D*>(SingleSecGenList       ->FindObject(temp1)));
//     TH3D hTempGenSmearSingleHist     = *(dynamic_cast<TH3D*>(SingleSecGenSmearedList->FindObject(temp1)));
//     TH3D hTempRecSingleHist          = *(dynamic_cast<TH3D*>(SingleSecRecList ->FindObject(temp2)));
//     // hTempGenSingleHist     .Sumw2();
//     hTempGenSmearSingleHist.Sumw2();
//     hTempRecSingleHist     .Sumw2();
//
//     TH1D* projPtGenSmear  = hTempGenSmearSingleHist.ProjectionX();  // Pt  projection Histogram
//     TH1D* projPtRec       = hTempRecSingleHist     .ProjectionX();  // Pt  projection Histogram
//     TH1D* projEtaGenSmear = hTempGenSmearSingleHist.ProjectionY();  // Eta projection Histogram
//     TH1D* projEtaRec      = hTempRecSingleHist     .ProjectionY();  // Eta projection Histogram
//     TH1D* projPhiGenSmear = hTempGenSmearSingleHist.ProjectionZ();  // Phi projection Histogram
//     TH1D* projPhiRec      = hTempRecSingleHist     .ProjectionZ();  // Phi projection Histogram
//     projPtGenSmear ->Sumw2();
//     projPtRec      ->Sumw2();
//     projEtaGenSmear->Sumw2();
//     projEtaRec     ->Sumw2();
//     projPhiGenSmear->Sumw2();
//     projPhiRec     ->Sumw2();
//     SetHistoStandardSettings(projPtGenSmear);
//     SetHistoStandardSettings(projEtaGenSmear);
//     SetHistoStandardSettings(projPhiGenSmear);
//     SetHistoStandardSettingsRed(projPtRec);
//     SetHistoStandardSettingsRed(projEtaRec);
//     SetHistoStandardSettingsRed(projPhiRec);
//     projPtGenSmear ->Scale(1,"width");
//     projEtaGenSmear->Scale(1,"width");
//     projPhiGenSmear->Scale(1,"width");
//     projPtRec      ->Scale(1,"width");
//     projEtaRec     ->Scale(1,"width");
//     projPhiRec     ->Scale(1,"width");
//
//     projPtGenSmear->SetTitle("Pt Projection");
//     projEtaGenSmear->SetTitle("Eta Projection");
//     projPhiGenSmear->SetTitle("Phi Projection");
//     projPtRec->SetTitle("Pt Projection");
//     projEtaRec->SetTitle("Eta Projection");
//     projPhiRec->SetTitle("Phi Projection");
//
//     // Double_t norm = 1;
//     // projPtGenSmear ->Scale(norm/projPtGenSmear ->Integral()/*, "width"*/);
//     // projPtRec      ->Scale(norm/projPtRec      ->Integral()/*, "width"*/);
//     // projEtaGenSmear->Scale(norm/projEtaGenSmear->Integral()/*, "width"*/);
//     // projEtaRec     ->Scale(norm/projEtaRec     ->Integral()/*, "width"*/);
//     // projPhiGenSmear->Scale(norm/projPhiGenSmear->Integral()/*, "width"*/);
//     // projPhiRec     ->Scale(norm/projPhiRec     ->Integral()/*, "width"*/);
//
//     auto ratioPt  = new TRatioPlot(projPtRec , projPtGenSmear );
//     auto ratioEta = new TRatioPlot(projEtaRec, projEtaGenSmear);
//     auto ratioPhi = new TRatioPlot(projPhiRec, projPhiGenSmear);
//
//     auto legend = new TLegend(0.9 ,0.85,0.65,0.9);
//     TLegendEntry *entry1=legend->AddEntry("name","generated smeared secondaries");        // Legend with Line, same Mother
//     TLegendEntry *entry2=legend->AddEntry("name","reconstructed secondaries");   // Legend with Line, different Mother
//     entry1->SetLineColor(kBlue); // blue
//     entry2->SetLineColor(kRed);  // red
//
//     cSignalCanvas->cd();
//     projEtaGenSmear -> Draw();
//     projEtaRec -> Draw("csame");
//     ratioEta -> Draw();
//     legend->Draw();
//     cSignalCanvas -> SaveAs(Form("%s%s.pdf",DocumentPathProj.Data(), "Eta_ComparisonGenSmearRec"));
//
//     projPhiGenSmear -> Draw();
//     projPhiRec -> Draw("same");
//     ratioPhi -> Draw();
//     legend->Draw();
//     cSignalCanvas -> SaveAs(Form("%s%s.pdf",DocumentPathProj.Data(), "Phi_ComparisonGenSmearRec"));
//
//     cLogYCanvas->cd();
//     projPtGenSmear -> Draw();
//     projPtRec -> Draw("same");
//     ratioPt -> Draw();
//     legend->Draw();
//     cLogYCanvas -> SaveAs(Form("%s%s.pdf",DocumentPathProj.Data(), "Pt_ComparisonGenSmearRec"));
//   }
// }

void LossBySecondaryCuts(TString PairCase, TString RecSecCaseNames, std::vector <TList*> HistCaseList, const unsigned int NRecCuts,  TString Replace1, TString Replace2){
                                                                                // if (DoDebug) printf("HistCaseList size: %d \n",HistCaseList->size());

  // LossBySecondaryCuts()
  TString DocumentPath = Form("Plots/%s/%s/",TrainNumber.Data(),PairCase.Data());
  gSystem->Exec(Form("mkdir -p %s",DocumentPath.Data()));
  TH1D* hCutEntiesHist = new TH1D("CutEntries","Entries of different cuts ;Cuts;#Entries",NRecCuts+1,0,0);  // Histogram get filled bin by bin with number of entries of different cuts
  TObjArray* arrHistName = RecSecCaseNames.Tokenize(":");
  unsigned int nHists=arrHistName->GetEntriesFast();
  Int_t nEntries;

  if (PairCase == "SingleElectrons") {
    TH3D hTempHist;
    for (size_t iCut = 0; iCut < NRecCuts; iCut++) {         // iteration for all cuts
      for (unsigned int iHist = 0; iHist < nHists; iHist++) {   // iteration for all histograms in list
        TString temp = arrHistName->At(iHist)->GetName();
        hTempHist = *(dynamic_cast<TH3D*>(HistCaseList[iCut]->FindObject(temp)));
        nEntries = hTempHist.GetEntries();
        TString binName = (((TString)HistCaseList[iCut]->GetName()).ReplaceAll(Form("%s",Replace1.Data()),"")).ReplaceAll(Form("%s",Replace2.Data()),"");
        hCutEntiesHist->Fill(binName,nEntries);
      }
    }
  }
  if (PairCase == "Pair") {
    TH2D hTempHist;
    for (size_t iCut = 0; iCut < NRecCuts; iCut++) {         // iteration for all cuts
      for (unsigned int iHist = 0; iHist < nHists; iHist++) {   // iteration for all histograms in list
        TString temp = arrHistName->At(iHist)->GetName();
        hTempHist = *(dynamic_cast<TH2D*>(HistCaseList[iCut]->FindObject(temp)));
        nEntries = hTempHist.GetEntries();
        TString binName = (((TString)HistCaseList[iCut]->GetName()).ReplaceAll(Form("%s",Replace1.Data()),"")).ReplaceAll(Form("%s",Replace2.Data()),"");
        hCutEntiesHist->Fill(binName,nEntries);
      }
    }
  }
  TCanvas *cSignalCanvas = new TCanvas("cSignalCanvas","",800,800);
  SetCanvasStandardSettings(cSignalCanvas);
  hCutEntiesHist -> Draw("HIST");
  cSignalCanvas -> SaveAs(Form("%s%s.pdf",DocumentPath.Data(), hCutEntiesHist->GetName()));
}


void PlotPIDComponents (TList* iSupportList, TString MCSingalList) {
  // TCanvas *cSignalCanvas = new TCanvas("cSignalCanvas","",800,800);
  // SetCanvasStandardSettings(cSignalCanvas);
  TObjArray* arrMCSignalLists=MCSingalList.Tokenize(";");
  const Int_t nMCSignalLists=arrMCSignalLists->GetEntriesFast();
  for ( Int_t iMCSigList = 0; iMCSigList < nMCSignalLists; iMCSigList++) {
    TString DocumentPathProjX   = Form("Plots/%s/%s/",TrainNumber.Data(), "Support"/*, arrMCSignalLists->At(iMCSigList)->GetName()*/);  // define path to mass projection
    gSystem->Exec(Form("mkdir -p %s",DocumentPathProjX.Data()));

    TList* iMCSignalList = (TList*)iSupportList->FindObject(arrMCSignalLists->At(iMCSigList)->GetName());
    TList* ITSnSigmaList = (TList*)iMCSignalList->FindObject("ITS_nSigmas");
    TList* TPCnSigmaList = (TList*)iMCSignalList->FindObject("TPC_nSigmas");
    TList* TOFnSigmaList = (TList*)iMCSignalList->FindObject("TOF_nSigmas");
    std::vector<TList*> vecDetectorLists;
    std::vector<TString> vecnSigmaPart;
    // vecDetectorLists.push_back(ITSnSigmaList);
    vecDetectorLists.push_back(TPCnSigmaList);
    // vecDetectorLists.push_back(TOFnSigmaList);
    vecnSigmaPart.push_back((TString)"elec");
    vecnSigmaPart.push_back((TString)"muon");
    vecnSigmaPart.push_back((TString)"pion");
    vecnSigmaPart.push_back((TString)"kaon");
    vecnSigmaPart.push_back((TString)"prot");
    for (unsigned int iDetector = 0; iDetector < vecDetectorLists.size(); iDetector++) {                                // Loop over the 3 detector folders
      for (Int_t iSelectedPIDSigma = 0; iSelectedPIDSigma < 5; iSelectedPIDSigma++) {         // Loop over the 5 nSigma Sections one for each: ele, muon, pion, kaon, proton
        TCanvas *cSignalCanvas = new TCanvas("cSignalCanvas","",800,800);
        TCanvas *cSignalCanvasLog = new TCanvas("cSignalCanvasLog","",800,800);
        SetCanvasStandardSettings(cSignalCanvas);
        SetCanvasStandardSettings(cSignalCanvasLog);
        cSignalCanvasLog->SetLogy();
        cSignalCanvas->cd();

        TH2D hTempnSigmaHistAll     = *(dynamic_cast<TH2D*>(vecDetectorLists.at(iDetector)->At(iSelectedPIDSigma*6+0)));
        TH2D hTempnSigmaHistEle     = *(dynamic_cast<TH2D*>(vecDetectorLists.at(iDetector)->At(iSelectedPIDSigma*6+1)));
        TH2D hTempnSigmaHistMuon    = *(dynamic_cast<TH2D*>(vecDetectorLists.at(iDetector)->At(iSelectedPIDSigma*6+2)));
        TH2D hTempnSigmaHistPion    = *(dynamic_cast<TH2D*>(vecDetectorLists.at(iDetector)->At(iSelectedPIDSigma*6+3)));
        TH2D hTempnSigmaHistKaon    = *(dynamic_cast<TH2D*>(vecDetectorLists.at(iDetector)->At(iSelectedPIDSigma*6+4)));
        TH2D hTempnSigmaHistProton  = *(dynamic_cast<TH2D*>(vecDetectorLists.at(iDetector)->At(iSelectedPIDSigma*6+5)));
        TH1D* hProjYTempnSigmaHistAll    = hTempnSigmaHistAll   .ProjectionY(     "nSigma_projection"                        );
        TH1D* hProjYTempnSigmaHistEle    = hTempnSigmaHistEle   .ProjectionY(Form("nSigma_projection_%s",vecnSigmaPart.at(0).Data()));
        TH1D* hProjYTempnSigmaHistMuon   = hTempnSigmaHistMuon  .ProjectionY(Form("nSigma_projection_%s",vecnSigmaPart.at(1).Data()));
        TH1D* hProjYTempnSigmaHistPion   = hTempnSigmaHistPion  .ProjectionY(Form("nSigma_projection_%s",vecnSigmaPart.at(2).Data()));
        TH1D* hProjYTempnSigmaHistKaon   = hTempnSigmaHistKaon  .ProjectionY(Form("nSigma_projection_%s",vecnSigmaPart.at(3).Data()));
        TH1D* hProjYTempnSigmaHistProton = hTempnSigmaHistProton.ProjectionY(Form("nSigma_projection_%s",vecnSigmaPart.at(4).Data()));
        SetHistoStandardSettings(hProjYTempnSigmaHistAll);
        SetHistoStandardSettingsBlue(hProjYTempnSigmaHistEle);
        SetHistoStandardSettingsRed(hProjYTempnSigmaHistPion);
        SetHistoStandardSettingsGreen(hProjYTempnSigmaHistMuon);
        SetHistoStandardSettingsOrange(hProjYTempnSigmaHistKaon);
        SetHistoStandardSettingsPink(hProjYTempnSigmaHistProton);

        auto legend = new TLegend(0.6,0.75,0.9,0.9);
        TLegendEntry *entry1=legend->AddEntry("name","all contributions","l");        // Legend with Line, same Mother
        TLegendEntry *entry2=legend->AddEntry("name","electrons","l");   // Legend with Line, different Mother
        TLegendEntry *entry3=legend->AddEntry("name","muons","l");   // Legend with Line, different Mother
        TLegendEntry *entry4=legend->AddEntry("name","pions","l");   // Legend with Line, different Mother
        TLegendEntry *entry5=legend->AddEntry("name","kaons","l");   // Legend with Line, different Mother
        TLegendEntry *entry6=legend->AddEntry("name","protons","l");   // Legend with Line, different Mother
        entry1->SetLineColor(kBlack);
        entry2->SetLineColor(kBlue);
        entry3->SetLineColor(kGreen+2);
        entry4->SetLineColor(kRed);
        entry5->SetLineColor(kOrange+2);
        entry6->SetLineColor(kPink+2);
        // double yMaxValue;
        // if(iSelectedPIDSigma == 0) yMaxValue = hProjYTempnSigmaHistAll->GetMaximum()*1.25;
        // hProjYTempnSigmaHistAll->GetYaxis()->SetRangeUser(0,yMaxValue);
        hProjYTempnSigmaHistAll->Draw("");
        hProjYTempnSigmaHistEle->Draw("same");
        hProjYTempnSigmaHistMuon->Draw("same");
        hProjYTempnSigmaHistPion->Draw("same");
        hProjYTempnSigmaHistKaon->Draw("same");
        hProjYTempnSigmaHistProton->Draw("same");
        legend->Draw("same");
        cSignalCanvas -> SaveAs(Form("%s%s%s_%s.pdf",DocumentPathProjX.Data(), vecDetectorLists.at(iDetector)->GetName(), hProjYTempnSigmaHistAll->GetName(), vecnSigmaPart.at(iSelectedPIDSigma).Data()));
        cSignalCanvasLog->cd();
        double yMaxValueLog;
        if(iSelectedPIDSigma == 0) yMaxValueLog = hProjYTempnSigmaHistAll->GetMaximum()*10;
        hProjYTempnSigmaHistAll->GetYaxis()->SetRangeUser(0.1,yMaxValueLog);
        hProjYTempnSigmaHistAll->Draw("");
        hProjYTempnSigmaHistEle->Draw("same");
        hProjYTempnSigmaHistMuon->Draw("same");
        hProjYTempnSigmaHistPion->Draw("same");
        hProjYTempnSigmaHistKaon->Draw("same");
        hProjYTempnSigmaHistProton->Draw("same");
        legend->Draw("same");
        cSignalCanvasLog -> SaveAs(Form("%s%s%s_log_%s.pdf",DocumentPathProjX.Data(), vecDetectorLists.at(iDetector)->GetName(), hProjYTempnSigmaHistAll->GetName(), vecnSigmaPart.at(iSelectedPIDSigma).Data()));
      }
    }

  }
}


// void ExtractDielectronPairSpectra(TString listName){
//   FourPairRecList->FindObject()
// }
