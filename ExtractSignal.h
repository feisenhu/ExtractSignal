#include "TLatex.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include <iostream>

void DrawProjection(std::vector<TH2D> VecPairHistos, std::vector<Double_t> vec_bin_pt, std::vector<Double_t> vec_bin_mass, std::vector<Int_t> vec_rebin_pt, std::vector<Int_t> vec_rebin_mass, Bool_t ExtFourPairPionSig, Bool_t ExtFourPairEtaSig, TString PairCase, Double_t nEvents);
void PlotPrimarySoBRatio(TString VecRecPairHistoNames, std::vector<TList*> PairRecPrimList, std::vector<Int_t> vec_rebin_pt, std::vector<Int_t> vec_rebin_mass, TString PairCase);
void MergeMassCutHistos(TList* PairList, TString VecMassCutHistoNames, TString VecMassCutHistoTitles, TString PairCase, double MassCut, Bool_t DoCutEff);
void DoExtractPairSignal(TString VecPairHistoNames, TString VecRecPairHistoNames, TList* PairGenList ,TList* PairGenSmearedList, TList* PairRecList, Bool_t ExtGen, Bool_t ExtGenSmeared, Bool_t ExtRec, std::vector<Int_t> vec_rebin_pt, std::vector<Int_t> vec_rebin_mass, std::vector<Double_t> vec_proj_bin_pt, std::vector<Double_t> vec_proj_bin_mass, Double_t nEvents);
void SetRebinningX(Int_t n_bins, const Double_t* bins, Int_t n_rebin, std::vector<Double_t> fBinsRebin);
void SetRebinningY(Int_t n_bins, const Double_t* bins, Int_t n_rebin, std::vector<Double_t> fBinsRebin);
// void CompareGenSmearRec(TString VecGenSmearedSingleHistoNames, TString VecRecSingleHistoNames, TList* SingleSecGenSmearedList, TList* SingleSecRecNoCutsList);
void LossBySecondaryCuts(TString PairCase, TString RecSecCaseNames, std::vector <TList*> HistCaseList, const unsigned int NRecCuts,  TString Replace1, TString Replace2);
void PlotPIDComponents (TList* iSupportList, TString MCSingalList);
Bool_t Rebin2DHistogram(TH2D& hIn, std::vector<Double_t> fBinsRebin_Mee, std::vector<Double_t> fBinsRebin_Ptee);


TH1D* fMCprojS_pion;
TH1D* fMCprojB_pion;
TH1D* fMCprojS_eta;
TH1D* fMCprojB_eta;

Double_t nEvents;
Bool_t reject;

Double_t BackgroundFunctionPion (Double_t *x, Double_t *par) {
  if (reject && x[0] > 0.08 && x[0] < 0.15) {
    TF1::RejectPoint();
    return 0;
  }
  // return par[0] + par[1]*x[0];
  // return par[0]+par[1]*x[0]+par[2]*x[0]*x[0];
  return (par[0]+par[1]*x[0]+par[2]*x[0]*x[0] +  TMath::Exp(par[3]+par[4]*x[0]))/nEvents;
}

Double_t BackgroundFunctionEta (Double_t *x, Double_t *par) {
  if (reject && x[0] > 0.45 && x[0] < 0.6) {
    TF1::RejectPoint();
    return 0;
  }
  // return par[0] + par[1]*x[0];
  // return par[0]+par[1]*x[0]+par[2]*x[0]*x[0];
  return (par[0]+par[1]*x[0]+par[2]*x[0]*x[0] +  TMath::Exp(par[3]+par[4]*x[0]))/nEvents;
}

Double_t template1(Double_t *x, Double_t *par)
{
  Double_t y1 = fMCprojS_pion->GetBinContent(fMCprojS_pion->GetXaxis()->FindFixBin(x[0]));
  Double_t y2 = fMCprojB_pion->GetBinContent(fMCprojB_pion->GetXaxis()->FindFixBin(x[0]));
  return par[0]*y1 + par[1]*y2;
}

Double_t template2(Double_t *x, Double_t *par)
{
  Double_t y1 = fMCprojS_eta->GetBinContent(fMCprojS_eta->GetXaxis()->FindFixBin(x[0]));
  Double_t y2 = fMCprojB_eta->GetBinContent(fMCprojB_eta->GetXaxis()->FindFixBin(x[0]));
  return par[0]*y1 + par[1]*y2;
}


void SetTextSettings(TLatex* text, Double_t textSize){
  text->SetNDC();
  text->SetTextColor(1);
  text->SetTextSize(textSize);
  text->SetTextFont(43);

}

void SetCanvasStandardSettings(TCanvas *cCanv){
  //
  TGaxis::SetMaxDigits(3);
  gStyle->SetOptStat(0); // <- das hier macht dies box rechts oben weg
  cCanv->SetTopMargin(0.1);
  cCanv->SetBottomMargin(0.15);
  cCanv->SetRightMargin(0.1);
  cCanv->SetLeftMargin(0.15);
  cCanv->SetTickx();
  cCanv->SetTicky();
  cCanv->SetLogy(0);
  cCanv->SetLogx(0);
}

void SetHistoStandardSettings(TH1* histo, Double_t XOffset = 1.2, Double_t YOffset = 1.3){
  histo->GetXaxis()->SetTitleOffset(XOffset);
  histo->GetYaxis()->SetTitleOffset(YOffset);
  histo->GetXaxis()->SetTitleSize(25);
  histo->GetYaxis()->SetTitleSize(25);
  histo->GetXaxis()->SetLabelSize(20);
  histo->GetYaxis()->SetLabelSize(20);
  histo->GetXaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetTitleFont(43);
  histo->GetXaxis()->SetTitleFont(43);

  histo->SetMarkerStyle(20);
  histo->SetMarkerSize(1.0);
  histo->SetLineWidth(1.0);
  histo->SetLineColor(kBlack);
  histo->SetMarkerColor(kBlack);
}

void SetHistoStandardSettingsRed(TH1* histo, Double_t XOffset = 1.2, Double_t YOffset = 1.3){
  histo->GetXaxis()->SetTitleOffset(XOffset);
  histo->GetYaxis()->SetTitleOffset(YOffset);
  histo->GetXaxis()->SetTitleSize(25);
  histo->GetYaxis()->SetTitleSize(25);
  histo->GetXaxis()->SetLabelSize(20);
  histo->GetYaxis()->SetLabelSize(20);
  histo->GetXaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetTitleFont(43);
  histo->GetXaxis()->SetTitleFont(43);

  histo->SetMarkerStyle(27);
  histo->SetMarkerSize(1.0);
  histo->SetLineWidth(1.0);
  histo->SetLineColor(kRed);
  histo->SetMarkerColor(kRed);
}

void SetHistoStandardSettingsBlue(TH1* histo, Double_t XOffset = 1.2, Double_t YOffset = 1.3){
  histo->GetXaxis()->SetTitleOffset(XOffset);
  histo->GetYaxis()->SetTitleOffset(YOffset);
  histo->GetXaxis()->SetTitleSize(25);
  histo->GetYaxis()->SetTitleSize(25);
  histo->GetXaxis()->SetLabelSize(20);
  histo->GetYaxis()->SetLabelSize(20);
  histo->GetXaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetTitleFont(43);
  histo->GetXaxis()->SetTitleFont(43);

  histo->SetMarkerStyle(34);
  histo->SetMarkerSize(1.0);
  histo->SetLineWidth(1.0);
  histo->SetLineColor(kBlue);
  histo->SetMarkerColor(kBlue);
}

void SetHistoStandardSettingsGreen(TH1* histo, Double_t XOffset = 1.2, Double_t YOffset = 1.3){
  histo->GetXaxis()->SetTitleOffset(XOffset);
  histo->GetYaxis()->SetTitleOffset(YOffset);
  histo->GetXaxis()->SetTitleSize(25);
  histo->GetYaxis()->SetTitleSize(25);
  histo->GetXaxis()->SetLabelSize(20);
  histo->GetYaxis()->SetLabelSize(20);
  histo->GetXaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetTitleFont(43);
  histo->GetXaxis()->SetTitleFont(43);

  histo->SetMarkerStyle(4);
  histo->SetMarkerSize(1.0);
  histo->SetLineWidth(1.0);
  histo->SetLineColor(kGreen+2);
  histo->SetMarkerColor(kGreen+2);
}

void SetHistoStandardSettingsOrange(TH1* histo, Double_t XOffset = 1.2, Double_t YOffset = 1.3){
  histo->GetXaxis()->SetTitleOffset(XOffset);
  histo->GetYaxis()->SetTitleOffset(YOffset);
  histo->GetXaxis()->SetTitleSize(25);
  histo->GetYaxis()->SetTitleSize(25);
  histo->GetXaxis()->SetLabelSize(20);
  histo->GetYaxis()->SetLabelSize(20);
  histo->GetXaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetTitleFont(43);
  histo->GetXaxis()->SetTitleFont(43);

  histo->SetMarkerStyle(20);
  histo->SetMarkerSize(0.5);
  histo->SetLineWidth(1);
  histo->SetLineColor(kOrange+2);
  histo->SetMarkerColor(kOrange+2);
}

void SetHistoStandardSettingsPink(TH1* histo, Double_t XOffset = 1.2, Double_t YOffset = 1.3){
  histo->GetXaxis()->SetTitleOffset(XOffset);
  histo->GetYaxis()->SetTitleOffset(YOffset);
  histo->GetXaxis()->SetTitleSize(25);
  histo->GetYaxis()->SetTitleSize(25);
  histo->GetXaxis()->SetLabelSize(20);
  histo->GetYaxis()->SetLabelSize(20);
  histo->GetXaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetLabelFont(43);
  histo->GetYaxis()->SetTitleFont(43);
  histo->GetXaxis()->SetTitleFont(43);

  histo->SetMarkerStyle(20);
  histo->SetMarkerSize(0.5);
  histo->SetLineWidth(1);
  histo->SetLineColor(kPink+2);
  histo->SetMarkerColor(kPink+2);
}


void SetMultipleHistoStandardSettings(TH1* histo, Double_t XOffset = 2.5, Double_t YOffset = 1.3){
  histo->GetXaxis()->SetTitleOffset(XOffset);
  histo->GetYaxis()->SetTitleOffset(YOffset);
  histo->GetXaxis()->SetTitleSize(20);
  histo->GetYaxis()->SetTitleSize(20);
  histo->GetXaxis()->SetLabelSize(20);
  histo->GetYaxis()->SetLabelSize(20);
  histo->GetXaxis()->SetLabelFont(23);
  histo->GetYaxis()->SetLabelFont(23);
  histo->GetYaxis()->SetTitleFont(23);
  histo->GetXaxis()->SetTitleFont(23);

  histo->SetMarkerStyle(5);
  histo->SetMarkerSize(1.);
  histo->SetLineWidth(1.);
  histo->SetLineColor(kBlack);
  histo->SetMarkerColor(kBlack);
}

void SetFitSettingsPion(TF1* fFit){
  fFit->SetLineColor(kRed);
  fFit->SetLineWidth(2);
  fFit->SetLineStyle(1);
}

void SetFitSettingsPionSignal(TF1* fFit){
  // fFit->SetLineColor(kRed);
  fFit->SetLineWidth(2);
  fFit->SetLineStyle(1);
  fFit->SetLineColor(kRed);
}

void SetFitSettingsPionBackground(TF1* fFit){
  // fFit->SetLineColor(kRed);
  fFit->SetLineWidth(2);
  fFit->SetLineStyle(2);
  fFit->SetLineColor(kGreen+2);
}

void SetFitSettingsEta(TF1* fFit){
  fFit->SetLineColor(kBlue);
  fFit->SetLineWidth(3);
  fFit->SetLineStyle(2);
}

void SetFitSettingsEtaSignal(TF1* fFit){
  // fFit->SetLineColor(kRed);
  fFit->SetLineWidth(3);
  fFit->SetLineStyle(1);
  fFit->SetLineColor(kViolet);
}

void SetFitSettingsEtaBackground(TF1* fFit){
  fFit->SetLineWidth(2);
  fFit->SetLineStyle(2);
  fFit->SetLineColor(kGreen+2);
}

void SetSigmaLineSettings(TLine* line){
  line->SetLineColor(kRed);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
}
