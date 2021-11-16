//////////////////////////////////////////////////////////////////
// AnalyzeOscSpectra.C
// Mike Wallbank, June 2020
// 
// First attempt at writing a DUNE CAFAna analysis!
// Analyze basic Far Detector spectra to test the effects of
// oscillation in the framework.
//////////////////////////////////////////////////////////////////

// framework
#include "CAFAna/Analysis/Calcs.h"
#include "CAFAna/Core/LoadFromFile.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Prediction/IPrediction.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Prediction/NDPredictionSterile.h"
#include "CAFAna/Prediction/FDPredictionSterile.h"
#include "OscLib/OscCalcSterile.h"

// stl

// ROOT
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMath.h"

using namespace ana;

// ---------------------------------------------------------------
void DrawPrediction(IPrediction* pred,
		    osc::OscCalcSterile* calc3f, osc::OscCalcSterile* calc4f,
		    double pot, TFile* outFile, std::string label);

// ---------------------------------------------------------------
void AnalyzeOscSpectra(std::string inFileName) {

  // Open input file and extract predictions etc
  std::string fileName = inFileName.substr(0, inFileName.find(".root"));
  std::string saveName = fileName + "Ana.root";
  TFile* inFile = new TFile(inFileName.c_str(), "READ");
  IPrediction* fd_pred = PredictionNoExtrap::LoadFrom(inFile->GetDirectory("fd_pred"), "fd_pred").release();
  IPrediction* nd_pred = NDPredictionSterile::LoadFrom(inFile->GetDirectory("nd_pred"), "nd_pred").release();
  IPrediction* fd_pred_bsm = FDPredictionSterile::LoadFrom(inFile->GetDirectory("fd_pred_bsm"), "fd_pred_bsm").release();

  // Make an oscillation calculator
  osc::OscCalcSterile* calc3f = DefaultSterileCalc(4);
  calc3f->SetNFlavors(4);
  osc::OscCalcSterile* calc4f = DefaultSterileCalc(4);
  calc4f->SetNFlavors(4);

  // Set angles
  calc3f->SetDm(2, 7.53e-5);
  calc3f->SetDm(3, 7.53e-5 + 2.44e-3);
  calc3f->SetDm(4, 0);
  calc3f->SetAngle(1, 2, 0.588);
  calc3f->SetAngle(1, 3, 0);//0.154);
  calc3f->SetAngle(2, 3, 0.722);
  calc3f->SetAngle(1, 4, 0);
  calc3f->SetAngle(2, 4, 0);
  calc3f->SetAngle(3, 4, 0);
  calc3f->SetDelta(1, 3, 0);//1.37*M_PI);

  calc4f->SetDm(2, 7.53e-5);
  calc4f->SetDm(3, 7.53e-5 + 2.44e-3);
  calc4f->SetDm(4, 1.0);
  calc4f->SetAngle(1, 2, 0.588);
  calc4f->SetAngle(1, 3, 0);//0.154);
  calc4f->SetAngle(2, 3, 0.722);
  calc4f->SetAngle(1, 4, 0);//0.16);
  calc4f->SetAngle(2, 4, 0.2);//10*TMath::DegToRad());
  calc4f->SetAngle(3, 4, 0.6);//35*TMath::DegToRad());
  calc4f->SetDelta(1, 3, 0);//1.37*M_PI);
  calc4f->SetDelta(2, 4, 0);//1.5*M_PI);
  calc4f->SetDelta(1, 4, 0);

  // Make up some exposures
  double pot = 1e20;

  // Out file
  TFile* outFile = new TFile(saveName.c_str(), "RECREATE");

  DrawPrediction(fd_pred, calc3f, calc4f, pot, outFile, "FD");
  DrawPrediction(nd_pred, calc3f, calc4f, pot, outFile, "ND");
  DrawPrediction(fd_pred_bsm, calc3f, calc4f, pot, outFile, "FDBSM");

  outFile->Close();

  return;

}

// ---------------------------------------------------------------
void DrawPrediction(IPrediction* pred,
		    osc::OscCalcSterile* calc3f, osc::OscCalcSterile* calc4f,
		    double pot, TFile* outFile, std::string label) {

  // Make some predictions
  TH1* all = pred->Predict(calc3f).ToTH1(pot);
  TH1* numu = pred->PredictComponent(calc3f, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* nue = pred->PredictComponent(calc3f, Flavors::kAllNuE, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* nutau = pred->PredictComponent(calc3f, Flavors::kAllNuTau, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* nc = pred->PredictComponent(calc3f, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(pot);

  TH1* all_st = pred->Predict(calc4f).ToTH1(pot);
  TH1* numu_st = pred->PredictComponent(calc4f, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* nue_st = pred->PredictComponent(calc4f, Flavors::kAllNuE, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* nutau_st = pred->PredictComponent(calc4f, Flavors::kAllNuTau, Current::kCC, Sign::kBoth).ToTH1(pot);
  TH1* nc_st = pred->PredictComponent(calc4f, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(pot);

  // Draw!
  TCanvas* canv = new TCanvas("canv", "", 800, 600);
  TLegend* leg_comp = new TLegend(0.65, 0.5, 0.85, 0.8);
  TLegend* leg_model = new TLegend(0.65, 0.3, 0.85, 0.45);

  all->SetLineColor(kBlack);
  all->Draw("hist");
  leg_comp->AddEntry(all, "All", "l");
  numu->SetLineColor(kRed);
  numu->Draw("hist same");
  leg_comp->AddEntry(numu, "#nu_{#mu}", "l");
  nue->SetLineColor(kGreen+2);
  nue->Draw("hist same");
  leg_comp->AddEntry(nue, "#nu_{e}", "l");
  nutau->SetLineColor(kMagenta+2);
  nutau->Draw("hist same");
  leg_comp->AddEntry(nutau, "#nu_{#tau}", "l");
  nc->SetLineColor(kAzure+2);
  nc->Draw("hist same");
  leg_comp->AddEntry(nc, "NC", "l");
  leg_comp->Draw();
  outFile->cd();
  canv->Write(Form("%sOscSpectra3f", label.c_str()));

  canv->cd();
  canv->Clear();
  leg_comp->Clear();
  all_st->SetLineColor(kBlack);
  all_st->Draw("hist");
  leg_comp->AddEntry(all, "All", "l");
  numu_st->SetLineColor(kRed);
  numu_st->Draw("hist same");
  leg_comp->AddEntry(numu, "#nu_{#mu}", "l");
  nue_st->SetLineColor(kGreen+2);
  nue_st->Draw("hist same");
  leg_comp->AddEntry(nue, "#nu_{e}", "l");
  nutau_st->SetLineColor(kMagenta+2);
  nutau_st->Draw("hist same");
  leg_comp->AddEntry(nutau, "#nu_{#tau}", "l");
  nc_st->SetLineColor(kAzure+2);
  nc_st->Draw("hist same");
  leg_comp->AddEntry(nc, "NC", "l");
  leg_comp->Draw();
  outFile->cd();
  canv->Write(Form("%sOscSpectra4f", label.c_str()));

  canv->cd();
  canv->Clear();
  leg_comp->Clear();
  all->Draw("hist");
  all_st->SetLineStyle(kDashed);
  all_st->Draw("hist same");
  leg_comp->AddEntry(all, "All", "l");
  leg_model->AddEntry(all, "3-flavor", "l");
  leg_model->AddEntry(all_st, "4-flavor", "l");
  numu->Draw("hist same");
  numu_st->SetLineStyle(kDashed);
  numu_st->Draw("hist same");
  leg_comp->AddEntry(numu, "#nu_{#mu}", "l");
  nue->Draw("hist same");
  nue_st->SetLineStyle(kDashed);
  nue_st->Draw("hist same");
  leg_comp->AddEntry(nue, "#nu_{e}", "l");
  nutau->Draw("hist same");
  nutau_st->SetLineStyle(kDashed);
  nutau_st->Draw("hist same");
  leg_comp->AddEntry(nutau, "#nu_{#tau}", "l");
  nc->Draw("hist same");
  nc_st->SetLineStyle(kDashed);
  nc_st->Draw("hist same");
  leg_comp->AddEntry(nc, "NC", "l");
  leg_comp->Draw();
  leg_model->Draw();
  outFile->cd();
  canv->Write(Form("%sOscSpectra", label.c_str()));

  delete canv;

  return;

}
