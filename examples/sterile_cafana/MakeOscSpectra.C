//////////////////////////////////////////////////////////////////
// MakeOscSpectra.C
// Mike Wallbank, June 2020
// 
// First attempt at writing a DUNE CAFAna analysis!
// Make some basic Far Detector spectra to test the effects of
// oscillation in the framework.
//////////////////////////////////////////////////////////////////

// framework
#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Prediction/IPrediction.h"
#include "CAFAna/Prediction/PredictionGenerator.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Prediction/NDPredictionSterile.h"
#include "CAFAna/Prediction/FDPredictionSterile.h"
#include "StandardRecord/StandardRecord.h"

// stl
#include <vector>
#include <string>

// ROOT
#include "TFile.h"

using namespace ana;

// ---------------------------------------------------------------
void MakeOscSpectra() {

  SpectrumLoader* fd_nonswap_loader
    = new SpectrumLoader("/pnfs/dune/persistent/users/LBL_TDR/CAFs/v4/FD_FHC_nonswap.root");
  SpectrumLoader* nd_nonswap_loader
    = new SpectrumLoader("/pnfs/dune/persistent/users/LBL_TDR/CAFs/v4/ND_FHC_FV_00.root");

  Loaders* loaders = new Loaders();
  loaders->SetLoaderPath("/pnfs/dune/persistent/users/LBL_TDR/CAFs/v4/ND_FHC_FV_00.root",
			 caf::kNEARDET, Loaders::DataMC::kMC, Loaders::SwappingConfig::kNonSwap);
  loaders->SetLoaderPath("/pnfs/dune/persistent/users/LBL_TDR/CAFs/v4/FD_FHC_nonswap.root",
			 caf::kFARDET, Loaders::DataMC::kMC, Loaders::SwappingConfig::kNonSwap);
  loaders->SetLoaderPath("/pnfs/dune/persistent/users/LBL_TDR/CAFs/v4/FD_FHC_nueswap.root",
			 caf::kFARDET, Loaders::DataMC::kMC, Loaders::SwappingConfig::kNueSwap);
  loaders->SetLoaderPath("/pnfs/dune/persistent/users/LBL_TDR/CAFs/v4/FD_FHC_tauswap.root",
			 caf::kFARDET, Loaders::DataMC::kMC, Loaders::SwappingConfig::kNuTauSwap);

  const Var kTrueE([](const caf::StandardRecord* sr)
		   {
		     return sr->Ev;
		   });

  Spectrum* nd_nonswap = new Spectrum("True Energy (GeV)", Binning::Simple(100,0,10), *nd_nonswap_loader, kTrueE, kNoCut);
  Spectrum* fd_nonswap = new Spectrum("True Energy (GeV)", Binning::Simple(100,0,10), *fd_nonswap_loader, kTrueE, kNoCut);

  // FD
  PredictionNoExtrap* fd_pred = new PredictionNoExtrap(*loaders, "True Energy (GeV)", Binning::Simple(100,0,10), kTrueE, kNoCut);
  FDPredictionGenerator* fd_predgen = new FDPredictionGenerator(HistAxis("True Energy (GeV)", Binning::Simple(100,0,10), kTrueE), kNoCut);
  IPrediction* fd_pred_bsm = fd_predgen->Generate(*loaders).release();

  // ND
  loaders->SetND(true);
  NDPredictionGenerator* nd_predgen = new NDPredictionGenerator(HistAxis("True Energy (GeV)", Binning::Simple(100,0,10), kTrueE), kNoCut);
  IPrediction* nd_pred = nd_predgen->Generate(*loaders).release();

  fd_nonswap_loader->Go();
  loaders->Go();

  TFile* outFile = new TFile("OscSpectra.root", "RECREATE");
  nd_nonswap->SaveTo(outFile->mkdir("nd_nonswap"), "nd_nonswap");
  fd_nonswap->SaveTo(outFile->mkdir("fd_nonswap"), "fd_nonswap");
  fd_pred->SaveTo(outFile->mkdir("fd_pred"), "fd_pred");
  fd_pred_bsm->SaveTo(outFile->mkdir("fd_pred_bsm"), "fd_pred_bsm");
  nd_pred->SaveTo(outFile->mkdir("nd_pred"), "nd_pred");
  outFile->Close();
  delete outFile;

  return;

}
