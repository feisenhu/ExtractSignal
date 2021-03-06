// #include "TLatex.h"
// #include "TMath.h"
#include "TGaxis.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TFile.h"
// #include "TF1.h"
#include "TH1.h"
// #include "TH2.h"
// #include "TH3.h"
// // #include "TRatioPlot.h"
#include "TObjArray.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TStyle.h"
// #include "TString.h"
#include <iostream>

// Program to plot two histograms, which are separatly stored as root files, into the same histogram
// The root files should be named the same.

Bool_t DoDebug        = kFALSE;

// PDF or PNG
Bool_t output_PDF     = kTRUE;
Bool_t output_PNG     = kFALSE;

//Plot Signal over Background or Significance
Bool_t doSoB          = kFALSE;
Bool_t doSignificance = !doSoB;

// Bool to plot Dalitz and GammaGamma hisograms from 2 different  Trains (toal: 4 histos)
Bool_t CompDiffTrains = kTRUE;

// Chnage legend to first entry with Prefilter, second without Prefilter
// if Bool true: compare w/wo Prefilter, if false: compare Dalitz diff prim pair mass cuts
Bool_t compWandWOPreFilter = kTRUE;

// TString rebinCase = "rebin50";
TString rebinCase = "rebin20";

TString fFileHistNames;

// TString DataFolder = "local_test";
// TString DataFolder = "merged_352_LHC17d1_LF_353_LHC18h1_354_LHC17h3_TrainWithSecSecFourPairing";
// TString DataFolder = "merged_403_LHC18h1_child3_404_LHC17h3_405_LHC17d1_LF_407_LHC18h1_child1+2_Dalitz_withPreFilter_withMasscut0-0.35_lowerSplitLevel";
// TString DataFolder = "merged_408_LHC18H1_409_LHC17h3_410_LHC17d1_LF_Dalitz_GammaGamma_noMassCut_withPrefilter";
// TString DataFolder = "merged_412_LHC18h1_413_LHC17h3_414_LHC17d1_LF_OnlyRec_DalitzGamma_withPrefilter_0.1-0.2MassCut";
TString DataFolder = "merged_419_LHC18h1_420_LHC17h3_421_LHC17d1_LF_OnlyRec_DalitzGammaGamma_widerSecSecPrefilter_0-0.35MassCut";
// TString DataFolder = "merged_429_LHC18h1_430_LHC17d1_LF_431_LHC17h3_DalitzGammaGamma_withoutPrefilter_withMassCut0-0.35";

// TString DataFolder2 = "merged_403_LHC18h1_child3_404_LHC17h3_405_LHC17d1_LF_407_LHC18h1_child1+2_Dalitz_withPreFilter_withMasscut0-0.35_lowerSplitLevel";
// TString DataFolder2 = "merged_412_LHC18h1_413_LHC17h3_414_LHC17d1_LF_OnlyRec_DalitzGamma_withPrefilter_0.1-0.2MassCut";
// TString DataFolder2 = "merged_426_LHC17h3_427_LHC17d1_LF_428_LHC18h1_OnlyRec_DalitzGammaGamma_withPrefilters_withMassCut0.1-0.2";
TString DataFolder2 = "merged_429_LHC18h1_430_LHC17d1_LF_431_LHC17h3_DalitzGammaGamma_withoutPrefilter_withMassCut0-0.35";

TString pathToFirstRootFile  = Form("Plots/%s/PrimSecPairing/FourPair/FourPairJPID_sum1_pt75_sec_kV0List/TemplateFit/Nrec_FourAnyPartPair1_FinalState_UndefinedMother/Signals/",DataFolder.Data());
TString pathToSecondRootFile = Form("Plots/%s/SecSecPairing/FourPair/FourPairpairkV0ListSecSec/TemplateFit/Nrec_FourAnyPartPair1_Secondary_UndefinedMother/Signals/",DataFolder.Data());
// used for merged_352
// TString pathToSecondRootFile = Form("Plots/%s/SecSecPairing/FourPair/FourPairtrack_V0_standardListSecSec/TemplateFit/Nrec_FourAnyPartPair1_Secondary_UndefinedMother/Signals/",DataFolder.Data());

TString path2ToFirstRootFile  = Form("Plots/%s/PrimSecPairing/FourPair/FourPairJPID_sum1_pt75_sec_kV0List/TemplateFit/Nrec_FourAnyPartPair1_FinalState_UndefinedMother/Signals/",DataFolder2.Data());
// TString path2ToSecondRootFile = Form("Plots/%s/SecSecPairing/FourPair/FourPairpairkV0ListSecSec/TemplateFit/Nrec_FourAnyPartPair1_Secondary_UndefinedMother/Signals/",DataFolder2.Data());

void PlotTwoHistograms() {
  if(doSoB)          fFileHistNames = "SignalOverBackgroundEta_pt0:0.3_rebin20;SignalOverBackgroundEta_pt0.3:0.5_rebin20;SignalOverBackgroundEta_pt0.5:0.8_rebin20;SignalOverBackgroundEta_pt0.8:1_rebin20;SignalOverBackgroundEta_pt1:2_rebin20;SignalOverBackgroundEta_pt2:4_rebin20";
  // if(doSoB)          fFileHistNames = "SignalOverBackgroundEta_pt0:0.3_rebin50;SignalOverBackgroundEta_pt0.3:0.5_rebin50;SignalOverBackgroundEta_pt0.5:0.8_rebin50;SignalOverBackgroundEta_pt0.8:1_rebin50;SignalOverBackgroundEta_pt1:2_rebin50;SignalOverBackgroundEta_pt2:4_rebin50";
  if(doSignificance) fFileHistNames = "SignificanceEta_pt0:0.3_rebin20;SignificanceEta_pt0.3:0.5_rebin20;SignificanceEta_pt0.5:0.8_rebin20;SignificanceEta_pt0.8:1_rebin20;SignificanceEta_pt1:2_rebin20;SignificanceEta_pt2:4_rebin20";
  // if(doSignificance) fFileHistNames = "SignificanceEta_pt0:0.3_rebin50;SignificanceEta_pt0.3:0.5_rebin50;SignificanceEta_pt0.5:0.8_rebin50;SignificanceEta_pt0.8:1_rebin50;SignificanceEta_pt1:2_rebin50;SignificanceEta_pt2:4_rebin50";


  gSystem->Exec(Form("mkdir -p Plots/%s/Comparison_Dalitz_GammaGamma/",DataFolder.Data()));
  TObjArray *arrFileHistNames =  fFileHistNames.Tokenize(";");
  const unsigned int nHists=arrFileHistNames->GetEntriesFast();

                                                                                if(DoDebug) Printf("%d", __LINE__);

  TCanvas *cSignalCanvas  = new TCanvas("cSignalCanvas","",800,800);
  cSignalCanvas->DivideSquare(nHists);
  TGaxis::SetMaxDigits(3);

  for (size_t iHist = 0; iHist < nHists; iHist++) {
    TFile *fFile1 = new TFile(Form("%s%s.root",pathToFirstRootFile.Data() ,arrFileHistNames->At(iHist)->GetName()));
    TFile *fFile2 = new TFile(Form("%s%s.root",pathToSecondRootFile.Data(),arrFileHistNames->At(iHist)->GetName()));
    TFile * fFile3;
    // TFile * fFile4;

    if(CompDiffTrains){
      fFile3 = new TFile(Form("%s%s.root",path2ToFirstRootFile.Data(),arrFileHistNames->At(iHist)->GetName()));
      // fFile4 = new TFile(Form("%s%s.root",path2ToSecondRootFile.Data(),arrFileHistNames->At(iHist)->GetName()));
    }

    TH1D* hist1;
    TH1D* hist2;

    TH1D* hist3;
    TH1D* hist4;

    if(doSoB){
      hist1 = (TH1D*) fFile1->Get("SignalOverBackgroundEta");
      hist2 = (TH1D*) fFile2->Get("SignalOverBackgroundEta");
      if(CompDiffTrains){
        hist3 = (TH1D*) fFile3->Get("SignalOverBackgroundEta");
        // hist4 = (TH1D*) fFile4->Get("SignalOverBackgroundEta");
      }
    }
    if (doSignificance) {
      hist1 = (TH1D*) fFile1->Get("SignificanceEta");
      hist2 = (TH1D*) fFile2->Get("SignificanceEta");
      if(CompDiffTrains){
        hist3 = (TH1D*) fFile3->Get("SignificanceEta");
        // hist4 = (TH1D*) fFile4->Get("SignificanceEta");
      }
    }

    hist1->SetLineColor(kBlue);   hist1->SetLineWidth(2);
    hist2->SetLineColor(kRed);    hist2->SetLineWidth(2);

    if(doSoB) hist1->GetXaxis()->SetTitleOffset(1.4);
    if(doSoB) hist2->GetXaxis()->SetTitleOffset(1.4);
    if(doSignificance) hist1->GetXaxis()->SetTitleOffset(1.2);
    if(doSignificance) hist2->GetXaxis()->SetTitleOffset(1.2);
    hist1->GetYaxis()->SetTitleOffset(1.2);
    hist2->GetYaxis()->SetTitleOffset(1.2);

    if(CompDiffTrains){
      hist3->SetLineColor(kGreen+2);   hist3->SetLineWidth(2); hist3->SetLineStyle(2);
      // hist4->SetLineColor(kViolet);    hist4->SetLineWidth(2); hist4->SetLineStyle(2);
      hist3->SetTitle(" ");  hist3->GetXaxis()->SetTitleOffset(1.2);
      // hist4->SetTitle(" ");  hist4->GetXaxis()->SetTitleOffset(1.2);
    }
                                                                                if(DoDebug) Printf("%d", __LINE__);

    // Double_t yAxis1Hist = hist1->GetYaxis()->GetXmax();
    Double_t yAxis1Hist = hist1->GetMaximum();
    // Double_t yAxis2Hist = hist2->GetYaxis()->GetXmax();
    Double_t yAxis2Hist;
    if (!CompDiffTrains) yAxis2Hist = hist2->GetMaximum();
    if (CompDiffTrains)  yAxis2Hist = hist3->GetMaximum();
                                                                                if(DoDebug) Printf("%d", __LINE__);
    Double_t x_low; Double_t y_low; Double_t x_max; Double_t y_max;
    if (!CompDiffTrains) {x_low = 0.15; y_low = 0.75; x_max = 0.39; y_max = 0.89;}
    if (CompDiffTrains)  {x_low = 0.12; y_low = 0.65; x_max = 0.39; y_max = 0.89;}
    auto legend = new TLegend(x_low,y_low,x_max,y_max);
    if (!CompDiffTrains) {
      TLegendEntry *entry1=legend->AddEntry(hist1,"Dalitz Rec","l");
      TLegendEntry *entry2=legend->AddEntry(hist2,"#gamma #gamma Rec","l");
    }
    if(CompDiffTrains && !compWandWOPreFilter){
      TLegendEntry *entry1=legend->AddEntry(hist1,"#splitline{Dalitz Rec}{0.0 < m^{prim}_{ee} < 0.35 #frac{GeV}{c^{2}}}","l");
      // TLegendEntry *entry1=legend->AddEntry(hist1,"Dalitz Rec with m^{prim}_{ee} 0.0-0.35","l");
      // TLegendEntry *entry3=legend->AddEntry(hist3,"Dalitz Rec with m^{prim}_{ee} 0.1-0.2","l");
      TLegendEntry *entry3=legend->AddEntry(hist3,"#splitline{Dalitz Rec}{0.1 < m^{prim}_{ee} < 0.2 #frac{GeV}{c^{2}}}","l");
      // TLegendEntry *entry4=legend->AddEntry(hist4,"#gamma #gamma Rec","l");
    }
    else if(CompDiffTrains && compWandWOPreFilter){
      TLegendEntry *entry1=legend->AddEntry(hist1,"Dalitz Rec with prefilter","l");
      TLegendEntry *entry3=legend->AddEntry(hist3,"Dalitz Rec without prefilter","l");

    }
    legend->SetBorderSize(0);
    legend->SetFillColorAlpha(0, 0.0);
    if (!CompDiffTrains) legend->SetTextSize(0.05);
    if (CompDiffTrains) legend->SetTextSize(0.035);
                                                                                if(DoDebug) Printf("%d", __LINE__);

    TString HistName = arrFileHistNames->At(iHist)->GetName();
    if (doSoB)           HistName = HistName.ReplaceAll("SignalOverBackgroundEta_pt","").Data();
    if (doSignificance)  HistName = HistName.ReplaceAll("SignificanceEta_pt","").Data();
    HistName = HistName.ReplaceAll(Form("_%s",rebinCase.Data()),"").Data();
    TObjArray* arrayPtIntervall = HistName.Tokenize(":");

                                                                                if(DoDebug) Printf("%d", __LINE__);

    auto legendInfo = new TLegend(0.5,0.65,0.85,0.89);
    TLegendEntry *entry5=legendInfo->AddEntry("collisionSystem"   ,"pp, #sqrt{s} = 13 TeV, |#eta_{e}| < 0.8","");
    TLegendEntry *entry6=legendInfo->AddEntry("SinglePtPrim","#font[12]{p}_{T,e}^{prim} > 0.075 GeV/#font[12]{c}","");
    TLegendEntry *entry7=legendInfo->AddEntry("SinglePtSec","#font[12]{p}_{T,e}^{sec} > 0.02 GeV/#font[12]{c}","");
    // TLegendEntry *entry7=legendInfo->AddEntry("PairPt"  ,"#font[12]{p}_{T,ee}^{prim} > 0.075 GeV/#font[12]{c} , #font[12]{p}_{T,ee}^{sec} > 0.02 GeV/#font[12]{c}","");
    TLegendEntry *entry8=legendInfo->AddEntry("FourPairPt"   ,Form("%s < #font[12]{p}_{T,eeee} < %s  GeV/#font[12]{c}",arrayPtIntervall->At(0)->GetName() ,arrayPtIntervall->At(1)->GetName()),"");
    legendInfo->SetBorderSize(0);
    legendInfo->SetFillColorAlpha(0, 0.0);
    legendInfo->SetTextSize(0.04);

    cSignalCanvas->cd(iHist+1);
                                                                                if(DoDebug) Printf("%d", __LINE__);
    if(!CompDiffTrains) {
      if (yAxis1Hist > yAxis2Hist)   hist1 -> Draw();    hist2 -> Draw("same");     legend -> Draw("same");   legendInfo->Draw("same");
      if (yAxis2Hist > yAxis1Hist)   hist2 -> Draw();    hist1 -> Draw("same");     legend -> Draw("same");   legendInfo->Draw("same");
    }
    if(CompDiffTrains){
      if (yAxis1Hist > yAxis2Hist)   hist1 -> Draw();    if(CompDiffTrains) {hist3->Draw("same");}    legend -> Draw("same");   legendInfo->Draw("same");
      if (yAxis2Hist > yAxis1Hist)   hist3 -> Draw();    hist1 -> Draw("same");  if(CompDiffTrains) {hist3->Draw("same"); /*hist4->Draw("same");*/}    legend -> Draw("same");   legendInfo->Draw("same");
    }
  }
  // cSignalCanvas -> SaveAs(Form("Plots/%s/Comparison_PrimSec_GammaGamma/%s.pdf",DataFolder.Data(),arrFileHistNames->At(iHist)->GetName()));

  TString CompareCase;  if(!CompDiffTrains) CompareCase = "DDYY"; if(CompDiffTrains && !compWandWOPreFilter) CompareCase = "DD0.0-0.35&0.1-0.2"; if(CompDiffTrains && compWandWOPreFilter) CompareCase = "DDw&woPrefilter";
  if (output_PDF) {
    if(doSoB)          cSignalCanvas -> SaveAs(Form("Plots/%s/Comparison_Dalitz_GammaGamma/AllSignalOverBackground_Eta_%s%s.pdf",DataFolder.Data(),rebinCase.Data(),CompareCase.Data()));
    if(doSignificance) cSignalCanvas -> SaveAs(Form("Plots/%s/Comparison_Dalitz_GammaGamma/AllSignificance_Eta_%s%s.pdf",DataFolder.Data(),rebinCase.Data(),CompareCase.Data()));
  }
  if (output_PNG) {
    if(doSoB)          cSignalCanvas -> SaveAs(Form("Plots/%s/Comparison_Dalitz_GammaGamma/AllSignalOverBackground_Eta_%s%s.png",DataFolder.Data(),rebinCase.Data(),CompareCase.Data()));
    if(doSignificance) cSignalCanvas -> SaveAs(Form("Plots/%s/Comparison_Dalitz_GammaGamma/AllSignificance_Eta_%s%s.png",DataFolder.Data(),rebinCase.Data(),CompareCase.Data()));
  }

}
