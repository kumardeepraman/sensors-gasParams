// EICrecon/src/detectors/DRICH/DRICH.cc

// Copyright (C) 2022, 2023, Christopher Dilks, Luigi Dello Stritto
// Subject to the terms in the LICENSE file found in the top-level directory.

// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024, Dmitry Kalinkin

#include <edm4eic/EDM4eicVersion.h>
#include <DD4hep/DetElement.h>
#include <DD4hep/Detector.h>
#include <Evaluator/DD4hepUnits.h>
#include <JANA/JApplication.h>
#include <JANA/JApplicationFwd.h>
#include <JANA/Utils/JTypeInfo.h>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <variant>
#include <vector>

// algorithm configurations
#include "algorithms/digi/PhotoMultiplierHitDigiConfig.h"
#include "algorithms/pid/IrtCherenkovParticleIDConfig.h"
#include "algorithms/pid/MergeParticleIDConfig.h"
#include "algorithms/tracking/TrackPropagationConfig.h"
#include "extensions/jana/JOmniFactoryGeneratorT.h"
// factories
#include "factories/digi/PhotoMultiplierHitDigi_factory.h"
#include "factories/pid/IrtCherenkovParticleID_factory.h"
#include "factories/pid/MergeCherenkovParticleID_factory.h"
#include "factories/pid/MergeTrack_factory.h"
#include "factories/pid/RichTrack_factory.h"
#include "services/geometry/richgeo/ActsGeo.h"
#include "services/geometry/richgeo/RichGeo.h"
#include "services/geometry/richgeo/RichGeo_service.h"

extern "C" {
void InitPlugin(JApplication* app) {
  InitJANAPlugin(app);

  using namespace eicrecon;

  // configuration parameters ///////////////////////////////////////////////

  // digitization
  PhotoMultiplierHitDigiConfig digi_cfg;
  digi_cfg.detectorName    = "DRICH";
  digi_cfg.readoutClass    = "DRICHHits";
  digi_cfg.hitTimeWindow   = 20.0;     // [ns]
  digi_cfg.timeResolution  = 1 / 16.0; // [ns]
  digi_cfg.speMean         = 80.0;
  digi_cfg.speError        = 16.0;
  digi_cfg.pedMean         = 200.0;
  digi_cfg.pedError        = 3.0;
  digi_cfg.enablePixelGaps = true;
  digi_cfg.safetyFactor    = 0.7;
  digi_cfg.enableNoise     = false;
  digi_cfg.noiseRate       = 20000;             // [Hz]
  digi_cfg.noiseTimeWindow = 20.0 * dd4hep::ns; // [ns]
  digi_cfg.quantumEfficiency.clear();

  // sipm 75 um
  digi_cfg.quantumEfficiency = {// wavelength units are [nm]
    {285, 0.10}, {295, 0.23}, {306, 0.28}, {316, 0.31}, {326, 0.33},
    {337, 0.35}, {347, 0.36}, {358, 0.36}, {368, 0.37}, {378, 0.39},
    {389, 0.43}, {399, 0.45}, {410, 0.47}, {420, 0.48}, {430, 0.49},
    {441, 0.50}, {451, 0.50}, {462, 0.50}, {472, 0.50}, {483, 0.50},
    {493, 0.49}, {503, 0.48}, {514, 0.47}, {524, 0.46}, {535, 0.44},
    {545, 0.42}, {555, 0.41}, {566, 0.39}, {576, 0.37}, {587, 0.36},
    {597, 0.34}, {607, 0.33}, {618, 0.31}, {628, 0.29}, {639, 0.27},
    {649, 0.26}, {660, 0.24}, {670, 0.23}, {680, 0.21}, {691, 0.20},
    {701, 0.19}, {712, 0.18}, {722, 0.17}, {732, 0.16}, {743, 0.15},
    {753, 0.14}, {764, 0.13}, {774, 0.12}, {784, 0.11}, {795, 0.10},
    {805, 0.10}, {816, 0.09}, {826, 0.08}, {836, 0.08}, {847, 0.07},
    {857, 0.06}, {868, 0.06}, {878, 0.05}, {1000, 0.00}
  };
 
  // Track propagation
  TrackPropagationConfig aerogel_track_cfg;
  TrackPropagationConfig gas_track_cfg;

  // get RICH geo service
  auto richGeoSvc       = app->GetService<RichGeo_service>();
  const auto* dd4hepGeo = richGeoSvc->GetDD4hepGeo();
  if (dd4hepGeo->world().children().contains("DRICH")) {
    const auto* actsGeo          = richGeoSvc->GetActsGeo("DRICH");
    auto aerogel_tracking_planes = actsGeo->TrackingPlanes(richgeo::kAerogel, 5);
    auto aerogel_track_point_cut = actsGeo->TrackPointCut(richgeo::kAerogel);
    auto gas_tracking_planes     = actsGeo->TrackingPlanes(richgeo::kGas, 10);
    auto gas_track_point_cut     = actsGeo->TrackPointCut(richgeo::kGas);
    auto filter_surface          = gas_tracking_planes.back();
    // track propagation to each radiator
    aerogel_track_cfg.filter_surfaces.push_back(filter_surface);
    aerogel_track_cfg.target_surfaces = aerogel_tracking_planes;
    aerogel_track_cfg.track_point_cut = aerogel_track_point_cut;
    gas_track_cfg.filter_surfaces.push_back(filter_surface);
    gas_track_cfg.target_surfaces = gas_tracking_planes;
    gas_track_cfg.track_point_cut = gas_track_point_cut;
  }

  // IRT PID
  IrtCherenkovParticleIDConfig irt_cfg;
  // - refractive index interpolation
  irt_cfg.numRIndexBins = 100;
  // - aerogel
  irt_cfg.radiators.insert({"Aerogel", RadiatorConfig{}});
  irt_cfg.radiators.at("Aerogel").referenceRIndex = 1.0260;
  irt_cfg.radiators.at("Aerogel").attenuation     = 48; // [mm]
  irt_cfg.radiators.at("Aerogel").smearingMode    = "gaussian";
  irt_cfg.radiators.at("Aerogel").smearing        = 2e-3; // [radians]
  // - gas
  irt_cfg.radiators.insert({"Gas", RadiatorConfig{}});
  irt_cfg.radiators.at("Gas").referenceRIndex = 1.00076;
  irt_cfg.radiators.at("Gas").attenuation     = 0; // [mm]
  irt_cfg.radiators.at("Gas").smearingMode    = "gaussian";
  irt_cfg.radiators.at("Gas").smearing        = 5e-3; // [radians]
  // - PDG list
  irt_cfg.pdgList.insert(irt_cfg.pdgList.end(), {11, 211, 321, 2212});
  // - cheat modes
  irt_cfg.cheatPhotonVertex = false;
  irt_cfg.cheatTrueRadiator = false;

  // Merge PID from radiators
  MergeParticleIDConfig merge_cfg;
  merge_cfg.mergeMode = MergeParticleIDConfig::kAddWeights;

  // wiring between factories and data ///////////////////////////////////////

  // digitization
  app->Add(new JOmniFactoryGeneratorT<PhotoMultiplierHitDigi_factory>("DRICHRawHits",
                                                                      {"EventHeader", "DRICHHits"},
                                                                      {"DRICHRawHits",
#if EDM4EIC_BUILD_VERSION >= EDM4EIC_VERSION(8, 7, 0)
                                                                       "DRICHRawHitsLinks",
#endif
                                                                       "DRICHRawHitsAssociations"},
                                                                      digi_cfg, app));

  // charged particle tracks
  app->Add(new JOmniFactoryGeneratorT<RichTrack_factory>(
      "DRICHAerogelTracks",
      {"CentralCKFTracks", "CentralCKFActsTrackStates", "CentralCKFActsTracks"},
      {"DRICHAerogelTracks"}, aerogel_track_cfg, app));
  app->Add(new JOmniFactoryGeneratorT<RichTrack_factory>(
      "DRICHGasTracks", {"CentralCKFTracks", "CentralCKFActsTrackStates", "CentralCKFActsTracks"},
      {"DRICHGasTracks"}, gas_track_cfg, app));

  app->Add(new JOmniFactoryGeneratorT<MergeTrack_factory>("DRICHMergedTracks",
                                                          {"DRICHAerogelTracks", "DRICHGasTracks"},
                                                          {"DRICHMergedTracks"}, {}, app));

  // PID algorithm
  app->Add(new JOmniFactoryGeneratorT<IrtCherenkovParticleID_factory>(
      "DRICHIrtCherenkovParticleID",
      {"DRICHAerogelTracks", "DRICHGasTracks", "DRICHMergedTracks", "DRICHRawHits",
       "DRICHRawHitsAssociations"},
      {"DRICHAerogelIrtCherenkovParticleID", "DRICHGasIrtCherenkovParticleID"}, irt_cfg, app));

  // merge aerogel and gas PID results
  app->Add(new JOmniFactoryGeneratorT<MergeCherenkovParticleID_factory>(
      "DRICHMergedIrtCherenkovParticleID",
      {"DRICHAerogelIrtCherenkovParticleID", "DRICHGasIrtCherenkovParticleID"},
      {"DRICHMergedIrtCherenkovParticleID"}, merge_cfg, app));

  // clang-format on
}
}
