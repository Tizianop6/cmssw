#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "FWCore/Utilities/interface/isFinite.h"

// TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "TFile.h"

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "DataFormats/Common/interface/ValidHandle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/GeantUnits.h"
#include "DataFormats/Math/interface/angle_units.h"
#include "DataFormats/ForwardDetId/interface/ETLDetId.h"
#include "DataFormats/ForwardDetId/interface/BTLDetId.h"

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"
#include "DataFormats/Common/interface/RefProd.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "Geometry/Records/interface/MTDDigiGeometryRecord.h"
#include "Geometry/Records/interface/MTDTopologyRcd.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDTopology.h"
#include "Geometry/MTDCommonData/interface/MTDTopologyMode.h"
#include "Geometry/MTDGeometryBuilder/interface/MTDGeometry.h"
#include "Geometry/MTDGeometryBuilder/interface/ProxyMTDTopology.h"
#include "Geometry/MTDGeometryBuilder/interface/RectangularMTDTopology.h"

#include "RecoMTD/DetLayers/interface/MTDDetLayerGeometry.h"
#include "RecoMTD/Records/interface/MTDRecoGeometryRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoMTD/DetLayers/interface/MTDDetLayerGeometry.h"
#include "RecoMTD/DetLayers/interface/MTDTrayBarrelLayer.h"
#include "RecoMTD/DetLayers/interface/MTDDetTray.h"
#include "RecoMTD/DetLayers/interface/MTDSectorForwardDoubleLayer.h"
#include "RecoMTD/DetLayers/interface/MTDDetSector.h"
#include "RecoMTD/Records/interface/MTDRecoGeometryRecord.h"

#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/KalmanUpdators/interface/Chi2MeasurementEstimator.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementError.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
#include "DataFormats/FTLRecHit/interface/FTLRecHitCollections.h"
#include "DataFormats/FTLRecHit/interface/FTLClusterCollections.h"
#include "DataFormats/TrackerRecHit2D/interface/MTDTrackingRecHit.h"

#include "DataFormats/Common/interface/OneToMany.h"
#include "DataFormats/Common/interface/AssociationMap.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/Associations/interface/MtdSimLayerClusterToTPAssociatorBaseImpl.h"
#include "SimDataFormats/CaloAnalysis/interface/MtdSimLayerCluster.h"
#include "SimDataFormats/Associations/interface/MtdRecoClusterToSimLayerClusterAssociationMap.h"
#include "SimDataFormats/Associations/interface/MtdSimLayerClusterToRecoClusterAssociationMap.h"

#include "CLHEP/Units/PhysicalConstants.h"
#include "MTDHit.h"

class MtdTracksValidation : public DQMEDAnalyzer {
public:
  explicit MtdTracksValidation(const edm::ParameterSet&);
  ~MtdTracksValidation() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void ClearVectors();
  void MakeBranches();

private:
  void bookHistograms(DQMStore::IBooker&, edm::Run const&, edm::EventSetup const&) override;

  void analyze(const edm::Event&, const edm::EventSetup&) override;

  const std::pair<bool, bool> checkAcceptance(
      const reco::Track&, const edm::Event&, const edm::EventSetup&, size_t&, float&, float&, float&, float&);

  const bool trkTPSelLV(const TrackingParticle&);
  const bool trkTPSelAll(const TrackingParticle&);
  const bool trkRecSel(const reco::TrackBase&);
  const bool trkRecSelLowPt(const reco::TrackBase&);
  const edm::Ref<std::vector<TrackingParticle>>* getMatchedTP(const reco::TrackBaseRef&);
  void isParticle(const reco::TrackRef&,
    const edm::ValueMap<float>&,
    const edm::ValueMap<float>&,
    const edm::ValueMap<float>&,
    const edm::ValueMap<float>&,
    const edm::ValueMap<float>&,
    unsigned int&,
    bool&,
    bool&,
    bool&,
    bool&);


  const unsigned long int uniqueId(const uint32_t x, const EncodedEventId& y) {
    const uint64_t a = static_cast<uint64_t>(x);
    const uint64_t b = static_cast<uint64_t>(y.rawId());

    if (x < y.rawId())
      return (b << 32) | a;
    else
      return (a << 32) | b;
  }

  bool isETL(const double eta) const { return (std::abs(eta) > trackMinEtlEta_) && (std::abs(eta) < trackMaxEtlEta_); }

  void fillTrackClusterMatchingHistograms(MonitorElement* me1,
                                          MonitorElement* me2,
                                          MonitorElement* me3,
                                          MonitorElement* me4,
                                          MonitorElement* me5,
                                          float var1,
                                          float var2,
                                          float var3,
                                          float var4,
                                          float var5,
                                          bool flag);

  // ------------ member data ------------

  const std::string folder_;
  const bool optionalPlots_;
  const float trackMaxPt_;
  const float trackMaxBtlEta_;
  const float trackMinEtlEta_;
  const float trackMaxEtlEta_;

  static constexpr double simUnit_ = 1e9;                // sim time in s while reco time in ns
  static constexpr double etacutGEN_ = 4.;               // |eta| < 4;
  static constexpr double etacutREC_ = 3.;               // |eta| < 3;
  static constexpr double pTcutBTL_ = 0.7;               // PT > 0.7 GeV
  static constexpr double pTcutETL_ = 0.2;               // PT > 0.2 GeV
  static constexpr double depositBTLthreshold_ = 1;      // threshold for energy deposit in BTL cell [MeV]
  static constexpr double depositETLthreshold_ = 0.001;  // threshold for energy deposit in ETL cell [MeV]
  static constexpr double rBTL_ = 110.0;
  static constexpr double zETL_ = 290.0;
  static constexpr double etaMatchCut_ = 0.05;
  static constexpr double cluDRradius_ = 0.05;  // to cluster rechits around extrapolated track

  const reco::RecoToSimCollection* r2s_;

  edm::EDGetTokenT<reco::TrackCollection> GenRecTrackToken_;
  edm::EDGetTokenT<reco::TrackCollection> RecTrackToken_;

  edm::EDGetTokenT<TrackingParticleCollection> trackingParticleCollectionToken_;
  edm::EDGetTokenT<reco::SimToRecoCollection> simToRecoAssociationToken_;
  edm::EDGetTokenT<reco::RecoToSimCollection> recoToSimAssociationToken_;
  edm::EDGetTokenT<reco::TPToSimCollectionMtd> tp2SimAssociationMapToken_;
  edm::EDGetTokenT<MtdRecoClusterToSimLayerClusterAssociationMap> r2sAssociationMapToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> probPiToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> probKToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> probPToken_;

  edm::EDGetTokenT<FTLRecHitCollection> btlRecHitsToken_;
  edm::EDGetTokenT<FTLRecHitCollection> etlRecHitsToken_;
  edm::EDGetTokenT<FTLClusterCollection> btlRecCluToken_;
  edm::EDGetTokenT<FTLClusterCollection> etlRecCluToken_;

  edm::EDGetTokenT<edm::ValueMap<int>> trackAssocToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> pathLengthToken_;

  edm::EDGetTokenT<edm::ValueMap<float>> tmtdToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> SigmatmtdToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0SrcToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> Sigmat0SrcToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0PidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> Sigmat0PidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> t0SafePidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> Sigmat0SafePidToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> SigmaTofPiToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> SigmaTofKToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> SigmaTofPToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> TofPiToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> TofKToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> TofPToken_;
  
  edm::EDGetTokenT<edm::ValueMap<float>> trackMVAQualToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> outermostHitPositionToken_;

  edm::ESGetToken<MTDGeometry, MTDDigiGeometryRecord> mtdgeoToken_;
  edm::ESGetToken<MTDTopology, MTDTopologyRcd> mtdtopoToken_;
  edm::ESGetToken<MTDDetLayerGeometry, MTDRecoGeometryRecord> mtdlayerToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magfieldToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> builderToken_;

  MonitorElement* meBTLTrackRPTime_;
  MonitorElement* meBTLTrackEtaTot_;
  MonitorElement* meBTLTrackPhiTot_;
  MonitorElement* meBTLTrackPtTot_;
  MonitorElement* meBTLTrackEtaMtd_;
  MonitorElement* meBTLTrackPhiMtd_;
  MonitorElement* meBTLTrackPtMtd_;
  MonitorElement* meBTLTrackPtRes_;

  MonitorElement* meETLTrackRPTime_;
  MonitorElement* meETLTrackEtaTot_;
  MonitorElement* meETLTrackPhiTot_;
  MonitorElement* meETLTrackPtTot_;
  MonitorElement* meETLTrackEtaMtd_;
  MonitorElement* meETLTrackPhiMtd_;
  MonitorElement* meETLTrackPtMtd_;
  MonitorElement* meETLTrackEta2Mtd_;
  MonitorElement* meETLTrackPhi2Mtd_;
  MonitorElement* meETLTrackPt2Mtd_;
  MonitorElement* meETLTrackPtRes_;

  MonitorElement* meETLTrackEtaTotLowPt_[2];
  MonitorElement* meETLTrackEtaMtdLowPt_[2];
  MonitorElement* meETLTrackEta2MtdLowPt_[2];

  MonitorElement* meBTLTrackMatchedTPEtaTot_;
  MonitorElement* meBTLTrackMatchedTPPtTot_;
  MonitorElement* meBTLTrackMatchedTPEtaMtd_;
  MonitorElement* meBTLTrackMatchedTPPtMtd_;
  MonitorElement* meETLTrackMatchedTPEtaTot_;
  MonitorElement* meETLTrackMatchedTPPtTot_;
  MonitorElement* meETLTrackMatchedTPEtaMtd_;
  MonitorElement* meETLTrackMatchedTPPtMtd_;
  MonitorElement* meETLTrackMatchedTPEta2Mtd_;
  MonitorElement* meETLTrackMatchedTPPt2Mtd_;
  MonitorElement* meETLTrackMatchedTPEtaMtdCorrect_;
  MonitorElement* meETLTrackMatchedTPPtMtdCorrect_;

  MonitorElement* meTracktmtd_;
  MonitorElement* meTrackt0Src_;
  MonitorElement* meTrackSigmat0Src_;
  MonitorElement* meTrackt0Pid_;
  MonitorElement* meTrackSigmat0Pid_;
  MonitorElement* meTrackt0SafePid_;
  MonitorElement* meTrackSigmat0SafePid_;
  MonitorElement* meTrackNumHits_;
  MonitorElement* meTrackNumHitsNT_;
  MonitorElement* meTrackMVAQual_;
  MonitorElement* meTrackPathLenghtvsEta_;
  MonitorElement* meTrackOutermostHitR_;
  MonitorElement* meTrackOutermostHitZ_;

  MonitorElement* meTrackSigmaTof_[3];
  MonitorElement* meTrackSigmaTofvsP_[3];

  MonitorElement* meBTLTrackMatchedTPPtResMtd_;
  MonitorElement* meETLTrackMatchedTPPtResMtd_;
  MonitorElement* meETLTrackMatchedTP2PtResMtd_;
  MonitorElement* meBTLTrackMatchedTPPtRatioGen_;
  MonitorElement* meETLTrackMatchedTPPtRatioGen_;
  MonitorElement* meETLTrackMatchedTP2PtRatioGen_;
  MonitorElement* meBTLTrackMatchedTPPtRatioMtd_;
  MonitorElement* meETLTrackMatchedTPPtRatioMtd_;
  MonitorElement* meETLTrackMatchedTP2PtRatioMtd_;
  MonitorElement* meBTLTrackMatchedTPPtResvsPtMtd_;
  MonitorElement* meETLTrackMatchedTPPtResvsPtMtd_;
  MonitorElement* meETLTrackMatchedTP2PtResvsPtMtd_;
  MonitorElement* meBTLTrackMatchedTPDPtvsPtGen_;
  MonitorElement* meETLTrackMatchedTPDPtvsPtGen_;
  MonitorElement* meETLTrackMatchedTP2DPtvsPtGen_;
  MonitorElement* meBTLTrackMatchedTPDPtvsPtMtd_;
  MonitorElement* meETLTrackMatchedTPDPtvsPtMtd_;
  MonitorElement* meETLTrackMatchedTP2DPtvsPtMtd_;

  MonitorElement* meTrackResTot_;
  MonitorElement* meTrackPullTot_;
  MonitorElement* meTrackResTotvsMVAQual_;
  MonitorElement* meTrackPullTotvsMVAQual_;

  MonitorElement* meTrackMatchedTPPtTotLV_;
  MonitorElement* meTrackMatchedTPEtaTotLV_;
  MonitorElement* meExtraPtMtd_;
  MonitorElement* meExtraPtEtl2Mtd_;
  MonitorElement* meExtraEtaMtd_;
  MonitorElement* meExtraEtaEtl2Mtd_;
  MonitorElement* meExtraPhiAtBTL_;
  MonitorElement* meExtraPhiAtBTLmatched_;
  MonitorElement* meExtraBTLeneInCone_;
  MonitorElement* meExtraMTDfailExtenderEta_;
  MonitorElement* meExtraMTDfailExtenderPt_;

  // ====== Trak-cluster matching based on MC truth
  // - BTL TPmtd Direct, TPmtd Other, TPnomtd
  MonitorElement* meBTLTrackMatchedTPmtdDirectEta_;
  MonitorElement* meBTLTrackMatchedTPmtdDirectPt_;

  MonitorElement* meBTLTrackMatchedTPmtdOtherEta_;
  MonitorElement* meBTLTrackMatchedTPmtdOtherPt_;

  MonitorElement* meBTLTrackMatchedTPnomtdEta_;
  MonitorElement* meBTLTrackMatchedTPnomtdPt_;

  // - BTL TPmtd Direct hits: correct, wrong, missing association in MTD
  MonitorElement* meBTLTrackMatchedTPmtdDirectCorrectAssocEta_;
  MonitorElement* meBTLTrackMatchedTPmtdDirectCorrectAssocPt_;
  MonitorElement* meBTLTrackMatchedTPmtdDirectCorrectAssocMVAQual_;
  MonitorElement* meBTLTrackMatchedTPmtdDirectCorrectAssocTimeRes_;
  MonitorElement* meBTLTrackMatchedTPmtdDirectCorrectAssocTimePull_;

  MonitorElement* meBTLTrackMatchedTPmtdDirectWrongAssocEta_;
  MonitorElement* meBTLTrackMatchedTPmtdDirectWrongAssocPt_;
  MonitorElement* meBTLTrackMatchedTPmtdDirectWrongAssocMVAQual_;
  MonitorElement* meBTLTrackMatchedTPmtdDirectWrongAssocTimeRes_;
  MonitorElement* meBTLTrackMatchedTPmtdDirectWrongAssocTimePull_;

  MonitorElement* meBTLTrackMatchedTPmtdDirectNoAssocEta_;
  MonitorElement* meBTLTrackMatchedTPmtdDirectNoAssocPt_;

  // - BTL TPmtd "other" hits: correct, wrong, missing association in MTD
  MonitorElement* meBTLTrackMatchedTPmtdOtherCorrectAssocEta_;
  MonitorElement* meBTLTrackMatchedTPmtdOtherCorrectAssocPt_;
  MonitorElement* meBTLTrackMatchedTPmtdOtherCorrectAssocMVAQual_;
  MonitorElement* meBTLTrackMatchedTPmtdOtherCorrectAssocTimeRes_;
  MonitorElement* meBTLTrackMatchedTPmtdOtherCorrectAssocTimePull_;

  MonitorElement* meBTLTrackMatchedTPmtdOtherWrongAssocEta_;
  MonitorElement* meBTLTrackMatchedTPmtdOtherWrongAssocPt_;
  MonitorElement* meBTLTrackMatchedTPmtdOtherWrongAssocMVAQual_;
  MonitorElement* meBTLTrackMatchedTPmtdOtherWrongAssocTimeRes_;
  MonitorElement* meBTLTrackMatchedTPmtdOtherWrongAssocTimePull_;

  MonitorElement* meBTLTrackMatchedTPmtdOtherNoAssocEta_;
  MonitorElement* meBTLTrackMatchedTPmtdOtherNoAssocPt_;

  // - BTL TPnomtd but a reco cluster is associated
  MonitorElement* meBTLTrackMatchedTPnomtdAssocEta_;
  MonitorElement* meBTLTrackMatchedTPnomtdAssocPt_;
  MonitorElement* meBTLTrackMatchedTPnomtdAssocMVAQual_;
  MonitorElement* meBTLTrackMatchedTPnomtdAssocTimeRes_;
  MonitorElement* meBTLTrackMatchedTPnomtdAssocTimePull_;

  // - ETL: one, two o no sim hits
  MonitorElement* meETLTrackMatchedTPmtd1Eta_;  // -- sim hit in >=1 etl disk
  MonitorElement* meETLTrackMatchedTPmtd1Pt_;
  MonitorElement* meETLTrackMatchedTPmtd2Eta_;  // -- sim hits in 2 etl disks
  MonitorElement* meETLTrackMatchedTPmtd2Pt_;
  MonitorElement* meETLTrackMatchedTPnomtdEta_;  // -- no sim hits in etl
  MonitorElement* meETLTrackMatchedTPnomtdPt_;

  // - ETL >=1 sim hit: each correct, at least one wrong, each sim hit missing reco association
  MonitorElement* meETLTrackMatchedTPmtd1CorrectAssocEta_;
  MonitorElement* meETLTrackMatchedTPmtd1CorrectAssocPt_;
  MonitorElement* meETLTrackMatchedTPmtd1CorrectAssocMVAQual_;
  MonitorElement* meETLTrackMatchedTPmtd1CorrectAssocTimeRes_;
  MonitorElement* meETLTrackMatchedTPmtd1CorrectAssocTimePull_;

  MonitorElement* meETLTrackMatchedTPmtd1WrongAssocEta_;
  MonitorElement* meETLTrackMatchedTPmtd1WrongAssocPt_;
  MonitorElement* meETLTrackMatchedTPmtd1WrongAssocMVAQual_;
  MonitorElement* meETLTrackMatchedTPmtd1WrongAssocTimeRes_;
  MonitorElement* meETLTrackMatchedTPmtd1WrongAssocTimePull_;

  MonitorElement* meETLTrackMatchedTPmtd1NoAssocEta_;
  MonitorElement* meETLTrackMatchedTPmtd1NoAssocPt_;
  MonitorElement* meETLTrackMatchedTPmtd1NoAssocMVAQual_;
  MonitorElement* meETLTrackMatchedTPmtd1NoAssocTimeRes_;
  MonitorElement* meETLTrackMatchedTPmtd1NoAssocTimePull_;

  // - ETL - 2 sim hits: both correct, at least one wrong or one missing, both missing reco association
  MonitorElement* meETLTrackMatchedTPmtd2CorrectAssocEta_;
  MonitorElement* meETLTrackMatchedTPmtd2CorrectAssocPt_;
  MonitorElement* meETLTrackMatchedTPmtd2CorrectAssocMVAQual_;
  MonitorElement* meETLTrackMatchedTPmtd2CorrectAssocTimeRes_;
  MonitorElement* meETLTrackMatchedTPmtd2CorrectAssocTimePull_;

  MonitorElement* meETLTrackMatchedTPmtd2WrongAssocEta_;
  MonitorElement* meETLTrackMatchedTPmtd2WrongAssocPt_;
  MonitorElement* meETLTrackMatchedTPmtd2WrongAssocMVAQual_;
  MonitorElement* meETLTrackMatchedTPmtd2WrongAssocTimeRes_;
  MonitorElement* meETLTrackMatchedTPmtd2WrongAssocTimePull_;

  MonitorElement* meETLTrackMatchedTPmtd2NoAssocEta_;
  MonitorElement* meETLTrackMatchedTPmtd2NoAssocPt_;
  MonitorElement* meETLTrackMatchedTPmtd2NoAssocMVAQual_;
  MonitorElement* meETLTrackMatchedTPmtd2NoAssocTimeRes_;
  MonitorElement* meETLTrackMatchedTPmtd2NoAssocTimePull_;

  // - ETL - no sim hits, but reco hit associated to the track
  MonitorElement* meETLTrackMatchedTPnomtdAssocEta_;
  MonitorElement* meETLTrackMatchedTPnomtdAssocPt_;
  MonitorElement* meETLTrackMatchedTPnomtdAssocMVAQual_;
  MonitorElement* meETLTrackMatchedTPnomtdAssocTimeRes_;
  MonitorElement* meETLTrackMatchedTPnomtdAssocTimePull_;
  TTree* dump_tree; 
  //vectors for dump in TTree:

  std::vector<float> xsim_vec;
  std::vector<float> ysim_vec;
  std::vector<float> zsim_vec;
  std::vector<float> xPCA_vec;
  std::vector<float> yPCA_vec;
  std::vector<float> zPCA_vec;
  std::vector<float> tsim_vec;
  std::vector<bool> isBTL_vec;
  std::vector<bool> isETL_vec;
  std::vector<float> outerZ_vec;
  std::vector<float> outerR_vec;
  std::vector<bool> isTPmtdDirectCorrectBTL_vec, isTPmtdOtherCorrectBTL_vec,isTPmtdDirectBTL_vec, isTPmtdOtherBTL_vec;
  std::vector<bool> isTPmtdCorrectETLD1_vec, isTPmtdCorrectETLD2_vec;
  std::vector<bool> ETLdisc1_vec, ETLdisc2_vec;
  std::vector<bool> isTPmtdETLD1_vec, isTPmtdETLD2_vec;
  std::vector<bool> trackMatchedtoTPnosimyesReco_vec;
  std::vector<bool> trkTPSelLV_vec;
  std::vector<float> TrackRes_vec, TrackPull_vec,sigmat0Safe_vec,t0safe_vec;
  std::vector<float> simPt_vec, simEta_vec, simPhi_vec,simPdgId_vec;
  std::vector<float> recoP_vec, recoPt_vec,recoPhi_vec, recoEta_vec;
  std::vector<float> outermostHitPosition_vec,mtdQualMVA_vec;
  std::vector<float> pathLength_vec;
  std::vector<int> ndof_vec,noPIDtype_vec;
  std::vector<bool> isPi_vec, isP_vec, isK_vec, noPID_vec;
  std::vector<bool> withMTD_vec;

  std::vector<float> par_curvature_vec, par_curvatureErr_vec;
  std::vector<float> par_d0_vec, par_d0Err_vec;
  std::vector<float> par_dz_vec, par_dzErr_vec;
  std::vector<float> par_phi_vec, par_phiErr_vec;
  std::vector<float> par_theta_vec, par_thetaErr_vec;
  std::vector<float> par_curvature_theta_cov_vec;
  std::vector<float> par_curvature_phi_cov_vec;
  std::vector<float> recoEta_err_vec, recoPhi_err_vec,recoP_err_vec;
  std::vector<float> sigmatmtd_vec;
  std::vector<float> tmtd_vec;
  std::vector<float> sigmat0_vec;
  std::vector<float> sigmat0safePID_vec;
  std::vector<float> tosafe_PID_vec, t0_PID_vec, sigmaTOFPi_vec, sigmaTOFK_vec, sigmaTOFP_vec;
  std::vector<float> t0Src_vec, tofK_vec, tofPi_vec, tofP_vec;
  std::vector<float> simCluster_time_vec;
  std::vector<float> simCluster_time_ETLD1_vec;
  std::vector<float> simCluster_time_ETLD2_vec;
  std::vector<float> MCTOF_vec;
  std::vector<float> MCTOF_ETLD1_vec;
  std::vector<float> MCTOF_ETLD2_vec;
  std::vector<float> simP_vec;
  std::vector<float> par_dsz_vec, par_dszErr_vec, par_dsz_lambda_cov_vec, par_dsz_curvature_cov_vec;
  std::vector<float> zHit_vec;
  std::vector<float> vx_vec,vy_vec,vz_vec;
  std::vector<float> TPstatus_vec, TPnoOfSimTrks_vec, isSimTrkPrimary_vec;
  std::vector<float> pl_SC_vec, pl_SC_earliesthit_vec, time_SC_vec, time_SC_earliesthit_vec;
  std::vector<float> pl_SC_ETLD1_vec, pl_SC_ETLD1_earliesthit_vec, time_SC_ETLD1_vec, time_SC_ETLD1_earliesthit_vec;
  std::vector<float> p_SC_vec, p_SC_earliesthit_vec;
  std::vector<float> p_SC_ETLD1_vec, p_SC_ETLD1_earliesthit_vec;
  std::vector<float> n_simHits_vec;
  

  
};

void MtdTracksValidation::ClearVectors(){
  xsim_vec.clear();
  ysim_vec.clear();
  zsim_vec.clear();
  xPCA_vec.clear();
  yPCA_vec.clear();
  zPCA_vec.clear();
  isBTL_vec.clear();
  isETL_vec.clear();
  isTPmtdDirectCorrectBTL_vec.clear(); 
  isTPmtdOtherCorrectBTL_vec.clear();
  isTPmtdDirectBTL_vec.clear();
  isTPmtdOtherBTL_vec.clear();
  isTPmtdCorrectETLD1_vec.clear(); 
  isTPmtdCorrectETLD2_vec.clear();
  ETLdisc1_vec.clear();
  ETLdisc2_vec.clear();
  isTPmtdETLD1_vec.clear();
  isTPmtdETLD2_vec.clear();
  trackMatchedtoTPnosimyesReco_vec.clear();
  trkTPSelLV_vec.clear();
  TrackRes_vec.clear();
  TrackPull_vec.clear();
  sigmat0Safe_vec.clear();
  t0safe_vec.clear();
  simPt_vec.clear();
  simEta_vec.clear();
  simPhi_vec.clear();
  simPdgId_vec.clear();
  recoP_vec.clear(); 
  recoPt_vec.clear();
  recoPhi_vec.clear(); 
  recoEta_vec.clear();
  outermostHitPosition_vec.clear();
  mtdQualMVA_vec.clear();
  pathLength_vec.clear();
  ndof_vec.clear();
  isPi_vec.clear();
  isP_vec.clear();
  isK_vec.clear();
  noPID_vec.clear();
  outerZ_vec.clear();
  outerR_vec.clear();
  withMTD_vec.clear();
  par_curvature_vec.clear(); 
  par_curvatureErr_vec.clear();
  par_d0_vec.clear(); 
  par_d0Err_vec.clear();
  par_dz_vec.clear(); 
  par_dzErr_vec.clear();
  par_phi_vec.clear();
  par_phiErr_vec.clear();
  par_theta_vec.clear();
  par_thetaErr_vec.clear();
  par_curvature_theta_cov_vec.clear();
  par_curvature_phi_cov_vec.clear();
  recoEta_err_vec.clear();
  recoPhi_err_vec.clear();
  recoP_err_vec.clear();
  sigmatmtd_vec.clear();
  tmtd_vec.clear();
  sigmat0_vec.clear();
  sigmat0safePID_vec.clear();
  tosafe_PID_vec.clear();
  t0_PID_vec.clear();
  sigmaTOFPi_vec.clear();
  sigmaTOFK_vec.clear();
  sigmaTOFP_vec.clear();
  t0Src_vec.clear();
  tofK_vec.clear();
  tofPi_vec.clear();
  tofP_vec.clear();
  simCluster_time_vec.clear();
  simCluster_time_ETLD1_vec.clear();
  simCluster_time_ETLD2_vec.clear();
  MCTOF_vec.clear();
  MCTOF_ETLD1_vec.clear();
  MCTOF_ETLD2_vec.clear();
  tsim_vec.clear();
  simP_vec.clear();
  par_dsz_vec.clear();
  par_dszErr_vec.clear();
  par_dsz_lambda_cov_vec.clear();
  zHit_vec.clear();
  vx_vec.clear();
  vy_vec.clear();
  vz_vec.clear();
  par_dsz_curvature_cov_vec.clear();

  TPstatus_vec.clear();
  TPnoOfSimTrks_vec.clear();
  isSimTrkPrimary_vec.clear();
  pl_SC_vec.clear();
  pl_SC_earliesthit_vec.clear();
  time_SC_vec.clear();
  time_SC_earliesthit_vec.clear();
  pl_SC_ETLD1_vec.clear();
  pl_SC_ETLD1_earliesthit_vec.clear();
  p_SC_vec.clear();
  p_SC_earliesthit_vec.clear();
  p_SC_ETLD1_vec.clear();
  p_SC_ETLD1_earliesthit_vec.clear();
  
  time_SC_ETLD1_vec.clear();
  time_SC_ETLD1_earliesthit_vec.clear();
  n_simHits_vec.clear();
  


}

// ------------ constructor and destructor --------------
MtdTracksValidation::MtdTracksValidation(const edm::ParameterSet& iConfig)
    : folder_(iConfig.getParameter<std::string>("folder")),
      optionalPlots_(iConfig.getParameter<bool>("optionalPlots")),
      trackMaxPt_(iConfig.getParameter<double>("trackMaximumPt")),
      trackMaxBtlEta_(iConfig.getParameter<double>("trackMaximumBtlEta")),
      trackMinEtlEta_(iConfig.getParameter<double>("trackMinimumEtlEta")),
      trackMaxEtlEta_(iConfig.getParameter<double>("trackMaximumEtlEta")) {
  GenRecTrackToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("inputTagG"));
  RecTrackToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("inputTagT"));
  trackingParticleCollectionToken_ =
      consumes<TrackingParticleCollection>(iConfig.getParameter<edm::InputTag>("SimTag"));
  simToRecoAssociationToken_ =
      consumes<reco::SimToRecoCollection>(iConfig.getParameter<edm::InputTag>("TPtoRecoTrackAssoc"));
  recoToSimAssociationToken_ =
      consumes<reco::RecoToSimCollection>(iConfig.getParameter<edm::InputTag>("TPtoRecoTrackAssoc"));
  tp2SimAssociationMapToken_ =
      consumes<reco::TPToSimCollectionMtd>(iConfig.getParameter<edm::InputTag>("tp2SimAssociationMapTag"));
  r2sAssociationMapToken_ = consumes<MtdRecoClusterToSimLayerClusterAssociationMap>(
      iConfig.getParameter<edm::InputTag>("r2sAssociationMapTag"));
  btlRecHitsToken_ = consumes<FTLRecHitCollection>(iConfig.getParameter<edm::InputTag>("btlRecHits"));
  etlRecHitsToken_ = consumes<FTLRecHitCollection>(iConfig.getParameter<edm::InputTag>("etlRecHits"));
  btlRecCluToken_ = consumes<FTLClusterCollection>(iConfig.getParameter<edm::InputTag>("recCluTagBTL"));
  etlRecCluToken_ = consumes<FTLClusterCollection>(iConfig.getParameter<edm::InputTag>("recCluTagETL"));
  trackAssocToken_ = consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("trackAssocSrc"));
  pathLengthToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("pathLengthSrc"));
  probPiToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("probPi"));
  probKToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("probK"));
  probPToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("probP"));
  tmtdToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("tmtd"));
  SigmatmtdToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmatmtd"));
  t0SrcToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0Src"));
  Sigmat0SrcToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0Src"));
  t0PidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0PID"));
  Sigmat0PidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0PID"));
  t0SafePidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("t0SafePID"));
  Sigmat0SafePidToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmat0SafePID"));
  SigmaTofPiToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmaTofPi"));
  SigmaTofKToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmaTofK"));
  SigmaTofPToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("sigmaTofP"));
  TofPiToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("TofPi"));
  TofKToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("TofK"));
  TofPToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("TofP"));
  
  trackMVAQualToken_ = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trackMVAQual"));
  outermostHitPositionToken_ =
      consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("outermostHitPositionSrc"));
  mtdgeoToken_ = esConsumes<MTDGeometry, MTDDigiGeometryRecord>();
  mtdtopoToken_ = esConsumes<MTDTopology, MTDTopologyRcd>();
  mtdlayerToken_ = esConsumes<MTDDetLayerGeometry, MTDRecoGeometryRecord>();
  magfieldToken_ = esConsumes<MagneticField, IdealMagneticFieldRecord>();
  builderToken_ = esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"));
  MakeBranches();
}

void MtdTracksValidation::MakeBranches(){
  edm::Service<TFileService> fs;
  dump_tree = fs->make<TTree>( "treeTrackValidation", "treeTrackValidation1" );
  dump_tree->Branch("trackPulls",&TrackPull_vec);
  dump_tree->Branch("trackRes",&TrackRes_vec);
  dump_tree->Branch("recoP",&recoP_vec);
  dump_tree->Branch("recoPt",&recoPt_vec);
  dump_tree->Branch("trackEta",&recoEta_vec);
  dump_tree->Branch("trackPhi",&recoPhi_vec);
  dump_tree->Branch("simPt",&simPt_vec);
  dump_tree->Branch("simEta",&simEta_vec);
  dump_tree->Branch("simPhi",&simPhi_vec);
  dump_tree->Branch("simPdgId",&simPdgId_vec);
  dump_tree->Branch("mtdQualMVA",&mtdQualMVA_vec);
  dump_tree->Branch("t0safe",&t0safe_vec);
  dump_tree->Branch("sigmat0Safe",&sigmat0Safe_vec);
  dump_tree->Branch("simPdgId",&simPdgId_vec);
  dump_tree->Branch("isPi",&isPi_vec);
  dump_tree->Branch("isP",&isP_vec);
  dump_tree->Branch("isK",&isK_vec);
  dump_tree->Branch("noPID",&noPID_vec);
  dump_tree->Branch("noPIDtype",&noPIDtype_vec);
  
  //dump_tree->Branch("matchCategory",&matchCategory_vec);
  dump_tree->Branch("xsim",&xsim_vec);
  dump_tree->Branch("ysim",&ysim_vec);
  dump_tree->Branch("zsim",&zsim_vec); 
  dump_tree->Branch("xPCA",&xPCA_vec); 
  dump_tree->Branch("yPCA",&yPCA_vec);  
  dump_tree->Branch("zPCA",&zPCA_vec);
  dump_tree->Branch("tsim",&tsim_vec);
  dump_tree->Branch("isBTL",&isBTL_vec);
  dump_tree->Branch("isETL",&isETL_vec);
  dump_tree->Branch("isTPmtdDirectCorrectBTL",&isTPmtdDirectCorrectBTL_vec);
  dump_tree->Branch("isTPmtdOtherCorrectBTL",&isTPmtdOtherCorrectBTL_vec);
  dump_tree->Branch("isTPmtdDirectBTL",&isTPmtdDirectBTL_vec);
  dump_tree->Branch("isTPmtdOtherBTL",&isTPmtdOtherBTL_vec);
  dump_tree->Branch("isTPmtdCorrectETLD1",&isTPmtdCorrectETLD1_vec);
  dump_tree->Branch("isTPmtdCorrectETLD2",&isTPmtdCorrectETLD2_vec);
  dump_tree->Branch("ETLdisc1",&ETLdisc1_vec);
  dump_tree->Branch("ETLdisc2",&ETLdisc2_vec);
  dump_tree->Branch("isTPmtdETLD1",&isTPmtdETLD1_vec);
  dump_tree->Branch("isTPmtdETLD2",&isTPmtdETLD2_vec);
  dump_tree->Branch("trackMatchedtoTPnosimyesReco",&trackMatchedtoTPnosimyesReco_vec);
  dump_tree->Branch("trkTPSelLV",&trkTPSelLV_vec);
  dump_tree->Branch("pathLength",&pathLength_vec);
  dump_tree->Branch("ndof",&ndof_vec);
  dump_tree->Branch("outermostHitPosition",&outermostHitPosition_vec);
  dump_tree->Branch("outerZ",&outerZ_vec);
  dump_tree->Branch("outerR",&outerR_vec);
  dump_tree->Branch("withMTD",&withMTD_vec);
    
  dump_tree->Branch("par_curvature",&par_curvature_vec);
  dump_tree->Branch("par_curvatureErr",&par_curvatureErr_vec);
  dump_tree->Branch("par_d0",&par_d0_vec);
  dump_tree->Branch("par_d0Err",&par_d0Err_vec);
  dump_tree->Branch("par_dz",&par_dz_vec);
  dump_tree->Branch("par_dzErr",&par_dzErr_vec);
  dump_tree->Branch("par_phi",&par_phi_vec);
  dump_tree->Branch("par_phiErr",&par_phiErr_vec);
  dump_tree->Branch("par_theta",&par_theta_vec);
  dump_tree->Branch("par_thetaErr",&par_thetaErr_vec);
  dump_tree->Branch("par_curvature_theta_cov",&par_curvature_theta_cov_vec);
  dump_tree->Branch("par_curvature_phi_cov",&par_curvature_phi_cov_vec);
  dump_tree->Branch("par_dsz_curvature_cov",&par_dsz_curvature_cov_vec);
  dump_tree->Branch("recoEta_err",&recoEta_err_vec);
  dump_tree->Branch("recoPhi_err",&recoPhi_err_vec);
  dump_tree->Branch("recoP_err",&recoP_err_vec);
  dump_tree->Branch("sigmatmtd_vec",&sigmatmtd_vec);
  dump_tree->Branch("tmtd",&tmtd_vec);
  dump_tree->Branch("sigmat0",&sigmat0_vec);
  dump_tree->Branch("sigmat0safePID",&sigmat0safePID_vec);
  dump_tree->Branch("tosafe_PID",&tosafe_PID_vec);
  dump_tree->Branch("t0_PID",&t0_PID_vec);
  dump_tree->Branch("sigmaTOFPi",&sigmaTOFPi_vec);
  dump_tree->Branch("sigmaTOFK",&sigmaTOFK_vec);
  dump_tree->Branch("sigmaTOFP",&sigmaTOFP_vec);
  dump_tree->Branch("t0Src",&t0Src_vec);
  dump_tree->Branch("tofK",&tofK_vec);
  dump_tree->Branch("tofPi",&tofPi_vec);
  dump_tree->Branch("tofP",&tofP_vec);
  dump_tree->Branch("simCluster_time",&simCluster_time_vec);
  dump_tree->Branch("simCluster_time_ETLD1",&simCluster_time_ETLD1_vec);
  dump_tree->Branch("simCluster_time_ETLD2",&simCluster_time_ETLD2_vec);
  dump_tree->Branch("MCTOF",&MCTOF_vec);
  dump_tree->Branch("MCTOF_ETLD1",&MCTOF_ETLD1_vec);
  dump_tree->Branch("MCTOF_ETLD2",&MCTOF_ETLD2_vec);
  dump_tree->Branch("simP",&simP_vec);
  dump_tree->Branch("par_dsz",&par_dsz_vec);
  dump_tree->Branch("par_dszErr",&par_dszErr_vec);
  dump_tree->Branch("par_dsz_lambda_cov",&par_dsz_lambda_cov_vec);
  dump_tree->Branch("zHit",&zHit_vec);
  dump_tree->Branch("vx",&vx_vec);
  dump_tree->Branch("vy",&vy_vec);
  dump_tree->Branch("vz",&vz_vec);
  
  dump_tree->Branch("TPnoOfSimTrks_vec",&TPnoOfSimTrks_vec);
  dump_tree->Branch("TPstatus",&TPstatus_vec);
  dump_tree->Branch("isSimTrkPrimary_vec",&isSimTrkPrimary_vec);
  
  dump_tree->Branch("pl_SC",&pl_SC_vec);
  dump_tree->Branch("pl_SC_earliesthit",&pl_SC_earliesthit_vec);
  dump_tree->Branch("time_SC",&time_SC_vec);
  dump_tree->Branch("time_SC_earliesthit",&time_SC_earliesthit_vec);
  dump_tree->Branch("pl_SC_ETLD1",&pl_SC_ETLD1_vec);
  dump_tree->Branch("pl_SC_ETLD1_earliesthit",&pl_SC_ETLD1_earliesthit_vec);

  dump_tree->Branch("time_SC_ETLD1",&time_SC_ETLD1_vec);
  dump_tree->Branch("time_SC_ETLD1_earliesthit",&time_SC_ETLD1_earliesthit_vec);
  
  dump_tree->Branch("p_SC",&p_SC_vec);
  dump_tree->Branch("p_SC_earliesthit",&p_SC_earliesthit_vec);
  dump_tree->Branch("p_SC_ETLD1",&p_SC_ETLD1_vec);
  dump_tree->Branch("p_SC_ETLD1_earliesthit",&p_SC_ETLD1_earliesthit_vec);

  dump_tree->Branch("n_simHits",&n_simHits_vec);
  
  
}

MtdTracksValidation::~MtdTracksValidation() {}

// ------------ method called for each event  ------------
void MtdTracksValidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace geant_units::operators;
  using namespace std;

  auto GenRecTrackHandle = makeValid(iEvent.getHandle(GenRecTrackToken_));

  auto btlRecCluHandle = makeValid(iEvent.getHandle(btlRecCluToken_));
  auto etlRecCluHandle = makeValid(iEvent.getHandle(etlRecCluToken_));

  std::unordered_map<uint32_t, MTDHit> m_btlHits;
  std::unordered_map<uint32_t, MTDHit> m_etlHits;
  std::unordered_map<uint32_t, std::set<unsigned long int>> m_btlTrkPerCell;
  std::unordered_map<uint32_t, std::set<unsigned long int>> m_etlTrkPerCell;
  const auto& tp2SimAssociationMap = iEvent.get(tp2SimAssociationMapToken_);
  const auto& r2sAssociationMap = iEvent.get(r2sAssociationMapToken_);

  const auto& tMtd = iEvent.get(tmtdToken_);
  const auto& SigmatMtd = iEvent.get(SigmatmtdToken_);
  const auto& t0Src = iEvent.get(t0SrcToken_);
  const auto& Sigmat0Src = iEvent.get(Sigmat0SrcToken_);
  const auto& t0Pid = iEvent.get(t0PidToken_);
  const auto& Sigmat0Pid = iEvent.get(Sigmat0PidToken_);
  const auto& t0Safe = iEvent.get(t0SafePidToken_);
  const auto& Sigmat0Safe = iEvent.get(Sigmat0SafePidToken_);
  const auto& SigmaTofPi = iEvent.get(SigmaTofPiToken_);
  const auto& SigmaTofK = iEvent.get(SigmaTofKToken_);
  const auto& SigmaTofP = iEvent.get(SigmaTofPToken_);
  const auto& TofPi = iEvent.get(TofPiToken_);
  const auto& TofK = iEvent.get(TofKToken_);
  const auto& TofP = iEvent.get(TofPToken_);
  
  const auto& mtdQualMVA = iEvent.get(trackMVAQualToken_);
  const auto& trackAssoc = iEvent.get(trackAssocToken_);
  const auto& pathLength = iEvent.get(pathLengthToken_);
  const auto& outermostHitPosition = iEvent.get(outermostHitPositionToken_);
  const auto& probPi = iEvent.get(probPiToken_);
  const auto& probK = iEvent.get(probKToken_);
  const auto& probP = iEvent.get(probPToken_);

  auto recoToSimH = makeValid(iEvent.getHandle(recoToSimAssociationToken_));
  r2s_ = recoToSimH.product();

  unsigned int index = 0;

  // --- Loop over all RECO tracks ---
  for (const auto& trackGen : *GenRecTrackHandle) {
    const reco::TrackRef trackref(iEvent.getHandle(GenRecTrackToken_), index);
    index++;

    if (trackAssoc[trackref] == -1) {
      LogWarning("mtdTracks") << "Extended track not associated";
      continue;
    }

    const reco::TrackRef mtdTrackref = reco::TrackRef(iEvent.getHandle(RecTrackToken_), trackAssoc[trackref]);
    const reco::Track& track = *mtdTrackref;

    bool isBTL = false;
    bool isETL = false;
    bool ETLdisc1 = false;
    bool ETLdisc2 = false;
    bool twoETLdiscs = false;
    bool noCrack = std::abs(trackGen.eta()) < trackMaxBtlEta_ || std::abs(trackGen.eta()) > trackMinEtlEta_;
    float MCTOF = std::numeric_limits<float>::max();
    float MCTOF_ETLD1 = std::numeric_limits<float>::max();
    float MCTOF_ETLD2 = std::numeric_limits<float>::max();
    float mtdhitposz = -1.;
        
    if (trkRecSel(trackGen)) {
      meTracktmtd_->Fill(tMtd[trackref]);
      if (std::round(SigmatMtd[trackref] - Sigmat0Pid[trackref]) != 0) {
        LogWarning("mtdTracks")
            << "TimeError associated to refitted track is different from TimeError stored in tofPID "
               "sigmat0 ValueMap: this should not happen";
      }

      meTrackt0Src_->Fill(t0Src[trackref]);
      meTrackSigmat0Src_->Fill(Sigmat0Src[trackref]);

      meTrackt0Pid_->Fill(t0Pid[trackref]);
      meTrackSigmat0Pid_->Fill(Sigmat0Pid[trackref]);
      meTrackt0SafePid_->Fill(t0Safe[trackref]);
      meTrackSigmat0SafePid_->Fill(std::log10(std::max(Sigmat0Safe[trackref], 0.001f)));
      meTrackMVAQual_->Fill(mtdQualMVA[trackref]);

      meTrackSigmaTof_[0]->Fill(SigmaTofPi[trackref] * 1e3);  //save as ps
      meTrackSigmaTof_[1]->Fill(SigmaTofK[trackref] * 1e3);
      meTrackSigmaTof_[2]->Fill(SigmaTofP[trackref] * 1e3);
      meTrackSigmaTofvsP_[0]->Fill(trackGen.p(), SigmaTofPi[trackref] * 1e3);
      meTrackSigmaTofvsP_[1]->Fill(trackGen.p(), SigmaTofK[trackref] * 1e3);
      meTrackSigmaTofvsP_[2]->Fill(trackGen.p(), SigmaTofP[trackref] * 1e3);

      meTrackPathLenghtvsEta_->Fill(std::abs(trackGen.eta()), pathLength[trackref]);
      bool MTDEtlZnegD1 = false;
      bool MTDEtlZnegD2 = false;
      bool MTDEtlZposD1 = false;
      bool MTDEtlZposD2 = false;
      std::vector<edm::Ref<edmNew::DetSetVector<FTLCluster>, FTLCluster>> recoClustersRefs;

      if (std::abs(trackGen.eta()) < trackMaxBtlEta_) {
        // --- all BTL tracks (with and without hit in MTD) ---
        meBTLTrackEtaTot_->Fill(std::abs(trackGen.eta()));
        meBTLTrackPhiTot_->Fill(trackGen.phi());
        meBTLTrackPtTot_->Fill(trackGen.pt());

        bool MTDBtl = false;
        int numMTDBtlvalidhits = 0;
        for (const auto hit : track.recHits()) {
          if (hit->isValid() == false)
            continue;
          MTDDetId Hit = hit->geographicalId();
          if ((Hit.det() == 6) && (Hit.subdetId() == 1) && (Hit.mtdSubDetector() == 1)) {
            MTDBtl = true;
            numMTDBtlvalidhits++;
            const auto* mtdhit = static_cast<const MTDTrackingRecHit*>(hit);
            const auto& hitCluster = mtdhit->mtdCluster();
            if (hitCluster.size() != 0) {
              auto recoClusterRef = edmNew::makeRefTo(btlRecCluHandle, &hitCluster);
              recoClustersRefs.push_back(recoClusterRef);
            }
          }
        }
        meTrackNumHits_->Fill(numMTDBtlvalidhits);

        // --- keeping only tracks with last hit in MTD ---
        if (MTDBtl == true) {
          isBTL = true;
          meBTLTrackEtaMtd_->Fill(std::abs(trackGen.eta()));
          meBTLTrackPhiMtd_->Fill(trackGen.phi());
          meBTLTrackPtMtd_->Fill(trackGen.pt());
          meBTLTrackRPTime_->Fill(track.t0());
          meBTLTrackPtRes_->Fill((trackGen.pt() - track.pt()) / trackGen.pt());
        }
        if (isBTL && Sigmat0Safe[trackref] < 0.) {
          meTrackNumHitsNT_->Fill(numMTDBtlvalidhits);
        }
      }  //loop over (geometrical) BTL tracks

      else {
        // --- all ETL tracks (with and without hit in MTD) ---
        meETLTrackEtaTot_->Fill(std::abs(trackGen.eta()));
        meETLTrackPhiTot_->Fill(trackGen.phi());
        meETLTrackPtTot_->Fill(trackGen.pt());

        int numMTDEtlvalidhits = 0;
        for (const auto hit : track.recHits()) {
          if (hit->isValid() == false)
            continue;
          MTDDetId Hit = hit->geographicalId();
          if ((Hit.det() == 6) && (Hit.subdetId() == 1) && (Hit.mtdSubDetector() == 2)) {
            isETL = true;
            ETLDetId ETLHit = hit->geographicalId();

            const auto* mtdhit = static_cast<const MTDTrackingRecHit*>(hit);
            const auto& hitCluster = mtdhit->mtdCluster();
            if (hitCluster.size() != 0) {
              auto recoClusterRef = edmNew::makeRefTo(etlRecCluHandle, &hitCluster);
              recoClustersRefs.push_back(recoClusterRef);
            }

            if ((ETLHit.zside() == -1) && (ETLHit.nDisc() == 1)) {
              MTDEtlZnegD1 = true;
              meETLTrackRPTime_->Fill(track.t0());
              meETLTrackPtRes_->Fill((trackGen.pt() - track.pt()) / trackGen.pt());
              mtdhitposz = mtdhit->globalPosition().z();
              numMTDEtlvalidhits++;
            }
            if ((ETLHit.zside() == -1) && (ETLHit.nDisc() == 2)) {
              MTDEtlZnegD2 = true;
              meETLTrackRPTime_->Fill(track.t0());
              meETLTrackPtRes_->Fill((trackGen.pt() - track.pt()) / trackGen.pt());
              numMTDEtlvalidhits++;
            }
            if ((ETLHit.zside() == 1) && (ETLHit.nDisc() == 1)) {
              MTDEtlZposD1 = true;
              meETLTrackRPTime_->Fill(track.t0());
              meETLTrackPtRes_->Fill((trackGen.pt() - track.pt()) / trackGen.pt());
              mtdhitposz = mtdhit->globalPosition().z();
            
              numMTDEtlvalidhits++;
            }
            if ((ETLHit.zside() == 1) && (ETLHit.nDisc() == 2)) {
              MTDEtlZposD2 = true;
              meETLTrackRPTime_->Fill(track.t0());
              meETLTrackPtRes_->Fill((trackGen.pt() - track.pt()) / trackGen.pt());
              numMTDEtlvalidhits++;
            }
          }
        }
        meTrackNumHits_->Fill(-numMTDEtlvalidhits);
        if (isETL && Sigmat0Safe[trackref] < 0.) {
          meTrackNumHitsNT_->Fill(-numMTDEtlvalidhits);
        }

        // --- keeping only tracks with last hit in MTD ---
        ETLdisc1 = (MTDEtlZnegD1 || MTDEtlZposD1);
        ETLdisc2 = (MTDEtlZnegD2 || MTDEtlZposD2);
        twoETLdiscs =
            ((MTDEtlZnegD1 == true) && (MTDEtlZnegD2 == true)) || ((MTDEtlZposD1 == true) && (MTDEtlZposD2 == true));
        if (ETLdisc1 || ETLdisc2) {
          meETLTrackEtaMtd_->Fill(std::abs(trackGen.eta()));
          meETLTrackPhiMtd_->Fill(trackGen.phi());
          meETLTrackPtMtd_->Fill(trackGen.pt());
          if (twoETLdiscs) {
            meETLTrackEta2Mtd_->Fill(std::abs(trackGen.eta()));
            meETLTrackPhi2Mtd_->Fill(trackGen.phi());
            meETLTrackPt2Mtd_->Fill(trackGen.pt());
          }
        }
      }

      if (isBTL)
        meTrackOutermostHitR_->Fill(outermostHitPosition[trackref]);
      if (isETL)
        meTrackOutermostHitZ_->Fill(std::abs(outermostHitPosition[trackref]));

      LogDebug("MtdTracksValidation") << "Track p/pt = " << trackGen.p() << " " << trackGen.pt() << " eta "
                                      << trackGen.eta() << " BTL " << isBTL << " ETL " << isETL << " 2disks "
                                      << twoETLdiscs;

      // == TrackingParticle based matching
      const reco::TrackBaseRef trkrefb(trackref);
      auto tp_info = getMatchedTP(trkrefb);
      
      if (tp_info != nullptr && trkTPSelAll(**tp_info)) {
        // -- pT resolution plots
        if (optionalPlots_) {
          if (trackGen.pt() < trackMaxPt_) {
            if (isBTL) {
              meBTLTrackMatchedTPPtResMtd_->Fill(std::abs(track.pt() - (*tp_info)->pt()) /
                                                 std::abs(trackGen.pt() - (*tp_info)->pt()));
              meBTLTrackMatchedTPPtRatioGen_->Fill(trackGen.pt() / (*tp_info)->pt());
              meBTLTrackMatchedTPPtRatioMtd_->Fill(track.pt() / (*tp_info)->pt());
              meBTLTrackMatchedTPPtResvsPtMtd_->Fill(
                  (*tp_info)->pt(),
                  std::abs(track.pt() - (*tp_info)->pt()) / std::abs(trackGen.pt() - (*tp_info)->pt()));
              meBTLTrackMatchedTPDPtvsPtGen_->Fill((*tp_info)->pt(),
                                                   (trackGen.pt() - (*tp_info)->pt()) / (*tp_info)->pt());
              meBTLTrackMatchedTPDPtvsPtMtd_->Fill((*tp_info)->pt(),
                                                   (track.pt() - (*tp_info)->pt()) / (*tp_info)->pt());
            }
            if (isETL && !twoETLdiscs && (std::abs(trackGen.eta()) > trackMinEtlEta_) &&
                (std::abs(trackGen.eta()) < trackMaxEtlEta_)) {
              meETLTrackMatchedTPPtResMtd_->Fill(std::abs(track.pt() - (*tp_info)->pt()) /
                                                 std::abs(trackGen.pt() - (*tp_info)->pt()));
              meETLTrackMatchedTPPtRatioGen_->Fill(trackGen.pt() / (*tp_info)->pt());
              meETLTrackMatchedTPPtRatioMtd_->Fill(track.pt() / (*tp_info)->pt());
              meETLTrackMatchedTPPtResvsPtMtd_->Fill(
                  (*tp_info)->pt(),
                  std::abs(track.pt() - (*tp_info)->pt()) / std::abs(trackGen.pt() - (*tp_info)->pt()));
              meETLTrackMatchedTPDPtvsPtGen_->Fill((*tp_info)->pt(),
                                                   (trackGen.pt() - (*tp_info)->pt()) / ((*tp_info)->pt()));
              meETLTrackMatchedTPDPtvsPtMtd_->Fill((*tp_info)->pt(),
                                                   (track.pt() - (*tp_info)->pt()) / ((*tp_info)->pt()));
            }
            if (isETL && twoETLdiscs) {
              meETLTrackMatchedTP2PtResMtd_->Fill(std::abs(track.pt() - (*tp_info)->pt()) /
                                                  std::abs(trackGen.pt() - (*tp_info)->pt()));
              meETLTrackMatchedTP2PtRatioGen_->Fill(trackGen.pt() / (*tp_info)->pt());
              meETLTrackMatchedTP2PtRatioMtd_->Fill(track.pt() / (*tp_info)->pt());
              meETLTrackMatchedTP2PtResvsPtMtd_->Fill(
                  (*tp_info)->pt(),
                  std::abs(track.pt() - (*tp_info)->pt()) / std::abs(trackGen.pt() - (*tp_info)->pt()));
              meETLTrackMatchedTP2DPtvsPtGen_->Fill((*tp_info)->pt(),
                                                    (trackGen.pt() - (*tp_info)->pt()) / ((*tp_info)->pt()));
              meETLTrackMatchedTP2DPtvsPtMtd_->Fill((*tp_info)->pt(),
                                                    (track.pt() - (*tp_info)->pt()) / ((*tp_info)->pt()));
            }
          }
        }

        // -- Track matched to TP: all and with last hit in MTD
        if (std::abs(trackGen.eta()) < trackMaxBtlEta_) {
          meBTLTrackMatchedTPEtaTot_->Fill(std::abs(trackGen.eta()));
          meBTLTrackMatchedTPPtTot_->Fill(trackGen.pt());
          if (isBTL) {
            meBTLTrackMatchedTPEtaMtd_->Fill(std::abs(trackGen.eta()));
            meBTLTrackMatchedTPPtMtd_->Fill(trackGen.pt());
          }
        } else {
          meETLTrackMatchedTPEtaTot_->Fill(std::abs(trackGen.eta()));
          meETLTrackMatchedTPPtTot_->Fill(trackGen.pt());
          if (isETL) {
            meETLTrackMatchedTPEtaMtd_->Fill(std::abs(trackGen.eta()));
            meETLTrackMatchedTPPtMtd_->Fill(trackGen.pt());
            if (twoETLdiscs) {
              meETLTrackMatchedTPEta2Mtd_->Fill(std::abs(trackGen.eta()));
              meETLTrackMatchedTPPt2Mtd_->Fill(trackGen.pt());
            }
          }
        }

        if (noCrack) {
          if (trkTPSelLV(**tp_info)) {
            meTrackMatchedTPEtaTotLV_->Fill(std::abs(trackGen.eta()));
            meTrackMatchedTPPtTotLV_->Fill(trackGen.pt());
          }
        }

        bool hasTime = false;
        double tsim = (*tp_info)->parentVertex()->position().t() * simUnit_;
        double dT(-9999.);
        double pullT(-9999.);
        if (Sigmat0Safe[trackref] != -1.) {
          dT = t0Safe[trackref] - tsim;
          pullT = dT / Sigmat0Safe[trackref];
          hasTime = true;
        }
        TrackRes_vec.push_back(dT);
        TrackPull_vec.push_back(pullT);
        sigmat0Safe_vec.push_back(Sigmat0Safe[trackref]);
        t0safe_vec.push_back(t0Safe[trackref]);
        trkTPSelLV_vec.push_back(trkTPSelLV(**tp_info));
        simPt_vec.push_back((*tp_info)->pt());
        simEta_vec.push_back((*tp_info)->eta());
        simPhi_vec.push_back((*tp_info)->phi());
        simPdgId_vec.push_back((*tp_info)->pdgId());
        recoP_vec.push_back(trackGen.p());
        recoPt_vec.push_back(trackGen.pt());
        recoPhi_vec.push_back(trackGen.phi());
        recoEta_vec.push_back(trackGen.eta());
        xsim_vec.push_back((*tp_info)->parentVertex()->position().x() * simUnit_);
        ysim_vec.push_back((*tp_info)->parentVertex()->position().y() * simUnit_);
        zsim_vec.push_back((*tp_info)->parentVertex()->position().z() * simUnit_);
        xPCA_vec.push_back((*trackref).vx());
        yPCA_vec.push_back((*trackref).vy());
        zPCA_vec.push_back((*trackref).vz());
        tsim_vec.push_back(tsim);
        sigmatmtd_vec.push_back(SigmatMtd[trackref]);
        tmtd_vec.push_back(tMtd[trackref]);
        sigmat0_vec.push_back(Sigmat0Src[trackref]);
        sigmat0safePID_vec.push_back(Sigmat0Safe[trackref]);
        outerZ_vec.push_back((*trackref).outerZ());
        //std::cout << "outerZ: " << (*trackref).outerZ() << std::endl;
        outerR_vec.push_back((*trackref).outerRadius());
        outermostHitPosition_vec.push_back(outermostHitPosition[trackref]);
        mtdQualMVA_vec.push_back(mtdQualMVA[trackref]);
        isBTL_vec.push_back(isBTL);
        isETL_vec.push_back(isETL);
        pathLength_vec.push_back(pathLength[trackref]);
        ndof_vec.push_back(trackGen.ndof());
                // ==  PID
        unsigned int no_PIDtype = 0;
        bool no_PID, is_Pi, is_K, is_P;
        isParticle(trackref, Sigmat0Src, Sigmat0Safe, probPi, probK, probP, no_PIDtype, no_PID, is_Pi, is_K, is_P);
        isPi_vec.push_back(is_Pi);
        isP_vec.push_back(is_P);
        isK_vec.push_back(is_K);
        noPIDtype_vec.push_back(no_PIDtype);
        noPID_vec.push_back(no_PID);
        par_curvature_vec.push_back(trackGen.parameter(0));
        par_curvatureErr_vec.push_back(pow(trackGen.covariance(0,0),0.5));
        par_d0_vec.push_back(trackGen.parameter(3));
        par_d0Err_vec.push_back(pow(trackGen.covariance(3,3),0.5));
        par_dz_vec.push_back(trackGen.parameter(4));
        par_dzErr_vec.push_back(pow(trackGen.covariance(4,4),0.5));
        par_phi_vec.push_back(trackGen.parameter(2));
        par_phiErr_vec.push_back(pow(trackGen.covariance(2,2),0.5));
        par_theta_vec.push_back(trackGen.parameter(1));
        par_thetaErr_vec.push_back(pow(trackGen.covariance(1,1),0.5));
        par_dsz_vec.push_back(trackGen.parameter(4));
        par_dszErr_vec.push_back(pow(trackGen.covariance(4,4),0.5));
        par_dsz_lambda_cov_vec.push_back(trackGen.covariance(4, 1));
        par_dsz_curvature_cov_vec.push_back(trackGen.covariance(4,0));
        recoEta_err_vec.push_back(trackGen.etaError());
        recoPhi_err_vec.push_back(trackGen.phiError());
        //recoP_err_vec.push_back(trackGen.pError());
        par_curvature_theta_cov_vec.push_back(trackGen.covariance(0, 1));
        par_curvature_phi_cov_vec.push_back(trackGen.covariance(0, 2));
        tosafe_PID_vec.push_back(t0Safe[trackref]);
        t0_PID_vec.push_back(t0Pid[trackref]);
        sigmaTOFPi_vec.push_back(SigmaTofPi[trackref]);
        sigmaTOFK_vec.push_back(SigmaTofK[trackref]);
        sigmaTOFP_vec.push_back(SigmaTofP[trackref]);
        t0Src_vec.push_back(t0Src[trackref]);
        tofK_vec.push_back(TofK[trackref]);
        tofPi_vec.push_back(TofPi[trackref]);
        tofP_vec.push_back(TofP[trackref]);
        simP_vec.push_back((*tp_info)->p());
        zHit_vec.push_back(mtdhitposz);
        
        TPstatus_vec.push_back((*tp_info)->status());
        TPnoOfSimTrks_vec.push_back((*tp_info)->g4Tracks().size());
        if((*tp_info)->g4Tracks().size()==1){
          isSimTrkPrimary_vec.push_back((*tp_info)->g4Tracks()[0].isPrimary());
        }else{
          isSimTrkPrimary_vec.push_back(-10); 
        }

        //vx_vec.push_back()
        
    
        

        // ==  MC truth matching
        bool isTPmtdDirectBTL = false, isTPmtdOtherBTL = false, isTPmtdDirectCorrectBTL = false,
             isTPmtdOtherCorrectBTL = false, isTPmtdETLD1 = false, isTPmtdETLD2 = false, isTPmtdCorrectETLD1 = false,
             isTPmtdCorrectETLD2 = false;

        auto simClustersRefsIt = tp2SimAssociationMap.find(*tp_info);
        const bool withMTD = (simClustersRefsIt != tp2SimAssociationMap.end());
        withMTD_vec.push_back(withMTD);
        // If there is a mtdSimLayerCluster from the tracking particle
        if (withMTD) {
          // -- Get the refs to MtdSimLayerClusters associated to the TP
          std::vector<edm::Ref<MtdSimLayerClusterCollection>> simClustersRefs;
          for (const auto& ref : simClustersRefsIt->val) {
            simClustersRefs.push_back(ref);
            MTDDetId mtddetid = ref->detIds_and_rows().front().first;
            if (mtddetid.mtdSubDetector() == 2) {
              ETLDetId detid(mtddetid.rawId());
              if (detid.nDisc() == 1)
                isTPmtdETLD1 = true;
              if (detid.nDisc() == 2)
                isTPmtdETLD2 = true;
            }
          }

          // === BTL
          // -- Sort BTL sim clusters by time
          std::vector<edm::Ref<MtdSimLayerClusterCollection>>::iterator directSimClusIt;
          if (std::abs(trackGen.eta()) < trackMaxBtlEta_ && !simClustersRefs.empty()) {
            std::sort(simClustersRefs.begin(), simClustersRefs.end(), [](const auto& a, const auto& b) {
              return a->simLCTime() < b->simLCTime();
            });
            // Find the first direct hit in time
            directSimClusIt = std::find_if(simClustersRefs.begin(), simClustersRefs.end(), [](const auto& simCluster) {
              MTDDetId mtddetid = simCluster->detIds_and_rows().front().first;
              return (mtddetid.mtdSubDetector() == 1 && simCluster->trackIdOffset() == 0);
            });
            // Check if TP has direct or other sim cluster for BTL
            for (const auto& simClusterRef : simClustersRefs) {
              if (directSimClusIt != simClustersRefs.end() && simClusterRef == *directSimClusIt) {
                isTPmtdDirectBTL = true;
              } else if (simClusterRef->trackIdOffset() != 0) {
                isTPmtdOtherBTL = true;
              }
            }
          }
          // ==  Check if the track-cluster association is correct: Track->RecoClus->SimClus == Track->TP->SimClus
          float maxClEnETL = -10;
          float maxClEnBTL = -10;
          bool earliesthit = false;
          float pl_SC, pl_SC_earliesthit, time_SC, time_SC_earliesthit, pl_SC_ETLD1, pl_SC_ETLD1_earliesthit, time_SC_ETLD1, time_SC_ETLD1_earliesthit;
          float p_SC, p_SC_earliesthit, p_SC_ETLD1, p_SC_ETLD1_earliesthit;
          float n_simHits=-1.;
          
          for (const auto& recClusterRef : recoClustersRefs) {
            if (recClusterRef.isNonnull()) {
              auto itp = r2sAssociationMap.equal_range(recClusterRef);
              if (itp.first != itp.second) {
                auto& simClustersRefs_RecoMatch = (*itp.first).second;

                for (const auto& simClusterRef_RecoMatch : simClustersRefs_RecoMatch) {
                  // Check if simClusterRef_RecoMatch  exists in SimClusters
                  auto simClusterIt =
                      std::find(simClustersRefs.begin(), simClustersRefs.end(), simClusterRef_RecoMatch);

                  // SimCluster found in SimClusters
                  if (simClusterIt != simClustersRefs.end()) {
                    if (isBTL) {
                      if (directSimClusIt != simClustersRefs.end() && simClusterRef_RecoMatch == *directSimClusIt) {
                        if (isTPmtdDirectCorrectBTL){
                          std::cout << "Problem: two assocs of direct in BTL";
                          /*if ((*simClusterIt)->simLCEnergy() > maxClEnBTL){
                            MCTOF = (*simClusterIt)->simLCTime();
                            maxClEnBTL = (*simClusterIt)->simLCEnergy();
                          }
                        }else{
                          MCTOF = simClusterRef_RecoMatch->simLCTime();
                          maxClEnBTL =  (*simClusterIt)->simLCEnergy();*/
                        }
                        if (earliesthit){
                          auto hits_and_times = (*simClusterIt)->hits_and_times();
                          for (uint32_t icounter = 0; icounter < hits_and_times.size(); icounter++) {
                            MCTOF = std::min(MCTOF, hits_and_times[icounter].second);
                          }
                        }else{
                          if (maxClEnBTL > 0 ){
                            MCTOF = -10;
                          }else{
                            n_simHits = (*simClusterIt)->hits_and_times().size();
                            if ((*simClusterIt)->simLCEnergy() > maxClEnBTL){
                              MCTOF = (*simClusterIt)->simLCTime();
                              maxClEnBTL = (*simClusterIt)->simLCEnergy();
                              pl_SC = (*simClusterIt)->simLCPL();
                              pl_SC_earliesthit = (*simClusterIt)->simLCearliestPL();
                              p_SC = (*simClusterIt)->simLCP();
                              p_SC_earliesthit = (*simClusterIt)->simLCearliestP();
                              
                              time_SC = (*simClusterIt)->simLCTime();
                              time_SC_earliesthit = (*simClusterIt)->simLCearliestTime();
                        }
                      }
                      }
                    
                        isTPmtdDirectCorrectBTL = true;
                        
                      } else if (simClusterRef_RecoMatch->trackIdOffset() != 0) {
                        isTPmtdOtherCorrectBTL = true;                        
                      }
                    }
                    if (isETL) {
                      MTDDetId mtddetid = (*simClusterIt)->detIds_and_rows().front().first;
                      ETLDetId detid(mtddetid.rawId());
                      if (detid.nDisc() == 1){
                        if (isTPmtdCorrectETLD1){
                          std::cout << "Problem: two assocs of correct D1 in ETL";
                          
                        }
                        isTPmtdCorrectETLD1 = true;
                        if(earliesthit){
                          auto hits_and_times = (*simClusterIt)->hits_and_times();
                          for (uint32_t icounter = 0; icounter < hits_and_times.size(); icounter++) {
                            MCTOF_ETLD1 = std::min(MCTOF_ETLD1, hits_and_times[icounter].second);
                          }
                        }else{
                          if (maxClEnETL > 0 ){
                            MCTOF_ETLD1 = -10;
                          }else{
                          if ((*simClusterIt)->simLCEnergy() > maxClEnETL){
                            MCTOF_ETLD1 = (*simClusterIt)->simLCTime();
                            maxClEnETL = (*simClusterIt)->simLCEnergy();
                            pl_SC_ETLD1 = (*simClusterIt)->simLCPL();
                            pl_SC_ETLD1_earliesthit = (*simClusterIt)->simLCearliestPL();
                            p_SC_ETLD1 = (*simClusterIt)->simLCP();
                            p_SC_ETLD1_earliesthit = (*simClusterIt)->simLCearliestP();
                            time_SC_ETLD1 = (*simClusterIt)->simLCTime();
                            time_SC_ETLD1_earliesthit = (*simClusterIt)->simLCearliestTime();
                        
                          }
                        }
                        }
                        
                        
                      }
                      if (detid.nDisc() == 2){
                        isTPmtdCorrectETLD2 = true;
                        MCTOF_ETLD2 = (*simClusterIt)->simLCTime();
                      }
                      
                    }

                  }
                  }
              }
            }
          }  /// end loop over reco clusters associated to this track.
          
          isTPmtdDirectCorrectBTL_vec.push_back(isTPmtdDirectCorrectBTL);
          isTPmtdOtherCorrectBTL_vec.push_back(isTPmtdOtherCorrectBTL);
          isTPmtdDirectBTL_vec.push_back(isTPmtdDirectBTL);
          isTPmtdOtherBTL_vec.push_back(isTPmtdOtherBTL);
          isTPmtdCorrectETLD1_vec.push_back(isTPmtdCorrectETLD1);
          isTPmtdCorrectETLD2_vec.push_back(isTPmtdCorrectETLD2);
          ETLdisc1_vec.push_back(ETLdisc1);
          ETLdisc2_vec.push_back(ETLdisc2);
          isTPmtdETLD1_vec.push_back(isTPmtdETLD1);
          isTPmtdETLD2_vec.push_back(isTPmtdETLD2);
          if (isTPmtdDirectCorrectBTL) {
            simCluster_time_vec.push_back(MCTOF);
            pl_SC_vec.push_back(pl_SC);
            pl_SC_earliesthit_vec.push_back(pl_SC_earliesthit);
            p_SC_vec.push_back(p_SC);
            p_SC_earliesthit_vec.push_back(p_SC_earliesthit);
            time_SC_vec.push_back(time_SC);
            time_SC_earliesthit_vec.push_back(time_SC_earliesthit);
            n_simHits_vec.push_back(n_simHits);
          }else{
            simCluster_time_vec.push_back(-99.);  // if no sim cluster associated, fill with -9999
            pl_SC_vec.push_back(-99.);
            pl_SC_earliesthit_vec.push_back(-99.);
            p_SC_vec.push_back(-99.);
            p_SC_earliesthit_vec.push_back(-99.);
            time_SC_vec.push_back(-99.);
            time_SC_earliesthit_vec.push_back(-99.);
            n_simHits_vec.push_back(-1.);
          }

          if  (isTPmtdCorrectETLD1){
            simCluster_time_ETLD1_vec.push_back(MCTOF_ETLD1);
            pl_SC_ETLD1_vec.push_back(pl_SC_ETLD1);
            pl_SC_ETLD1_earliesthit_vec.push_back(pl_SC_ETLD1_earliesthit);
            p_SC_ETLD1_vec.push_back(p_SC_ETLD1);
            p_SC_ETLD1_earliesthit_vec.push_back(p_SC_ETLD1_earliesthit);
            time_SC_ETLD1_vec.push_back(time_SC_ETLD1);
            time_SC_ETLD1_earliesthit_vec.push_back(time_SC_ETLD1_earliesthit);
            }else{
            simCluster_time_ETLD1_vec.push_back(-99.);  // if no sim cluster associated, fill with -9999
            pl_SC_ETLD1_vec.push_back(-99.);
            pl_SC_ETLD1_earliesthit_vec.push_back(-99.);
            p_SC_ETLD1_vec.push_back(-99.);
            p_SC_ETLD1_earliesthit_vec.push_back(-99.);
            time_SC_ETLD1_vec.push_back(-99.);
            time_SC_ETLD1_earliesthit_vec.push_back(-99.);
          }

          if (isTPmtdCorrectETLD2) {
            simCluster_time_ETLD2_vec.push_back(MCTOF_ETLD2);              
          } else{
            simCluster_time_ETLD2_vec.push_back(-99.);
          }

          // == BTL
          if (std::abs(trackGen.eta()) < trackMaxBtlEta_) {
            // -- Track matched to TP with sim hit in MTD
            if (isTPmtdDirectBTL) {
              meBTLTrackMatchedTPmtdDirectEta_->Fill(std::abs(trackGen.eta()));
              meBTLTrackMatchedTPmtdDirectPt_->Fill(trackGen.pt());
            } else if (isTPmtdOtherBTL) {
              meBTLTrackMatchedTPmtdOtherEta_->Fill(std::abs(trackGen.eta()));
              meBTLTrackMatchedTPmtdOtherPt_->Fill(trackGen.pt());
            }
            //-- Track matched to TP with sim hit in MTD, with associated reco cluster
            if (isBTL) {
              if (isTPmtdDirectBTL) {
                // -- Track matched to TP with sim hit (direct), correctly associated reco cluster
                if (isTPmtdDirectCorrectBTL) {
                  fillTrackClusterMatchingHistograms(meBTLTrackMatchedTPmtdDirectCorrectAssocEta_,
                                                     meBTLTrackMatchedTPmtdDirectCorrectAssocPt_,
                                                     meBTLTrackMatchedTPmtdDirectCorrectAssocMVAQual_,
                                                     meBTLTrackMatchedTPmtdDirectCorrectAssocTimeRes_,
                                                     meBTLTrackMatchedTPmtdDirectCorrectAssocTimePull_,
                                                     std::abs(trackGen.eta()),
                                                     trackGen.pt(),
                                                     mtdQualMVA[trackref],
                                                     dT,
                                                     pullT,
                                                     hasTime);
                }
                // -- Track matched to TP with sim hit (direct), incorrectly associated reco cluster
                else {
                  fillTrackClusterMatchingHistograms(meBTLTrackMatchedTPmtdDirectWrongAssocEta_,
                                                     meBTLTrackMatchedTPmtdDirectWrongAssocPt_,
                                                     meBTLTrackMatchedTPmtdDirectWrongAssocMVAQual_,
                                                     meBTLTrackMatchedTPmtdDirectWrongAssocTimeRes_,
                                                     meBTLTrackMatchedTPmtdDirectWrongAssocTimePull_,
                                                     std::abs(trackGen.eta()),
                                                     trackGen.pt(),
                                                     mtdQualMVA[trackref],
                                                     dT,
                                                     pullT,
                                                     hasTime);
                }
              }
              // -- Track matched to TP with sim hit (other), correctly associated reco cluster
              else if (isTPmtdOtherBTL) {
                if (isTPmtdOtherCorrectBTL) {
                  fillTrackClusterMatchingHistograms(meBTLTrackMatchedTPmtdOtherCorrectAssocEta_,
                                                     meBTLTrackMatchedTPmtdOtherCorrectAssocPt_,
                                                     meBTLTrackMatchedTPmtdOtherCorrectAssocMVAQual_,
                                                     meBTLTrackMatchedTPmtdOtherCorrectAssocTimeRes_,
                                                     meBTLTrackMatchedTPmtdOtherCorrectAssocTimePull_,
                                                     std::abs(trackGen.eta()),
                                                     trackGen.pt(),
                                                     mtdQualMVA[trackref],
                                                     dT,
                                                     pullT,
                                                     hasTime);
                }
                // -- Track matched to TP with sim hit (other), incorrectly associated reco cluster
                else {
                  fillTrackClusterMatchingHistograms(meBTLTrackMatchedTPmtdOtherWrongAssocEta_,
                                                     meBTLTrackMatchedTPmtdOtherWrongAssocPt_,
                                                     meBTLTrackMatchedTPmtdOtherWrongAssocMVAQual_,
                                                     meBTLTrackMatchedTPmtdOtherWrongAssocTimeRes_,
                                                     meBTLTrackMatchedTPmtdOtherWrongAssocTimePull_,
                                                     std::abs(trackGen.eta()),
                                                     trackGen.pt(),
                                                     mtdQualMVA[trackref],
                                                     dT,
                                                     pullT,
                                                     hasTime);
                }
              }
            }
            // -- Track matched to TP with sim hit in MTD, missing associated reco cluster
            else {
              if (isTPmtdDirectBTL) {
                meBTLTrackMatchedTPmtdDirectNoAssocEta_->Fill(std::abs(trackGen.eta()));
                meBTLTrackMatchedTPmtdDirectNoAssocPt_->Fill(trackGen.pt());
              } else if (isTPmtdOtherBTL) {
                meBTLTrackMatchedTPmtdOtherNoAssocEta_->Fill(std::abs(trackGen.eta()));
                meBTLTrackMatchedTPmtdOtherNoAssocPt_->Fill(trackGen.pt());
              }
            }
          }  // == end BTL
          // == ETL
          else {
            // -- Track matched to TP with reco hits (one or two) correctly matched
            if ((ETLdisc1 && isTPmtdCorrectETLD1) || (ETLdisc2 && isTPmtdCorrectETLD2)) {
              meETLTrackMatchedTPEtaMtdCorrect_->Fill(std::abs(trackGen.eta()));
              meETLTrackMatchedTPPtMtdCorrect_->Fill(trackGen.pt());
            }
            // -- Track matched to TP with sim hit in one etl layer
            if (isTPmtdETLD1 || isTPmtdETLD2) {  // at least one hit (D1 or D2 or both)
              meETLTrackMatchedTPmtd1Eta_->Fill(std::abs(trackGen.eta()));
              meETLTrackMatchedTPmtd1Pt_->Fill(trackGen.pt());
            }
            // -- Track matched to TP with sim hits in both etl layers (D1 and D2)
            if (isTPmtdETLD1 && isTPmtdETLD2) {
              meETLTrackMatchedTPmtd2Eta_->Fill(std::abs(trackGen.eta()));
              meETLTrackMatchedTPmtd2Pt_->Fill(trackGen.pt());
            }
            if (isETL) {
              // -- Track matched to TP with sim hit in >=1 etl layer
              if (isTPmtdETLD1 || isTPmtdETLD2) {
                // - each hit is correctly associated to the track
                if ((isTPmtdETLD1 && !isTPmtdETLD2 && ETLdisc1 && isTPmtdCorrectETLD1) ||
                    (isTPmtdETLD2 && !isTPmtdETLD1 && ETLdisc2 && isTPmtdCorrectETLD2) ||
                    (isTPmtdETLD1 && isTPmtdETLD2 && ETLdisc1 && ETLdisc2 && isTPmtdCorrectETLD1 &&
                     isTPmtdCorrectETLD2)) {
                  fillTrackClusterMatchingHistograms(meETLTrackMatchedTPmtd1CorrectAssocEta_,
                                                     meETLTrackMatchedTPmtd1CorrectAssocPt_,
                                                     meETLTrackMatchedTPmtd1CorrectAssocMVAQual_,
                                                     meETLTrackMatchedTPmtd1CorrectAssocTimeRes_,
                                                     meETLTrackMatchedTPmtd1CorrectAssocTimePull_,
                                                     std::abs(trackGen.eta()),
                                                     trackGen.pt(),
                                                     mtdQualMVA[trackref],
                                                     dT,
                                                     pullT,
                                                     hasTime);
                }
                // - at least one reco hit is incorrectly associated or, if two sim hits, one reco hit is missing
                else if ((isTPmtdETLD1 && !isTPmtdCorrectETLD1) || (isTPmtdETLD2 && !isTPmtdCorrectETLD2)) {
                  fillTrackClusterMatchingHistograms(meETLTrackMatchedTPmtd1WrongAssocEta_,
                                                     meETLTrackMatchedTPmtd1WrongAssocPt_,
                                                     meETLTrackMatchedTPmtd1WrongAssocMVAQual_,
                                                     meETLTrackMatchedTPmtd1WrongAssocTimeRes_,
                                                     meETLTrackMatchedTPmtd1WrongAssocTimePull_,
                                                     std::abs(trackGen.eta()),
                                                     trackGen.pt(),
                                                     mtdQualMVA[trackref],
                                                     dT,
                                                     pullT,
                                                     hasTime);
                }
              }
              // -- Track matched to TP with sim hits in both etl layers (D1 and D2)
              if (isTPmtdETLD1 && isTPmtdETLD2) {
                // - each hit correctly associated to the track
                if (ETLdisc1 && ETLdisc2 && isTPmtdCorrectETLD1 && isTPmtdCorrectETLD2) {
                  fillTrackClusterMatchingHistograms(meETLTrackMatchedTPmtd2CorrectAssocEta_,
                                                     meETLTrackMatchedTPmtd2CorrectAssocPt_,
                                                     meETLTrackMatchedTPmtd2CorrectAssocMVAQual_,
                                                     meETLTrackMatchedTPmtd2CorrectAssocTimeRes_,
                                                     meETLTrackMatchedTPmtd2CorrectAssocTimePull_,
                                                     std::abs(trackGen.eta()),
                                                     trackGen.pt(),
                                                     mtdQualMVA[trackref],
                                                     dT,
                                                     pullT,
                                                     hasTime);
                }
                // - at least one reco hit incorrectly associated or one hit missing
                else if ((ETLdisc1 || ETLdisc2) && (!isTPmtdCorrectETLD1 || !isTPmtdCorrectETLD2)) {
                  fillTrackClusterMatchingHistograms(meETLTrackMatchedTPmtd2WrongAssocEta_,
                                                     meETLTrackMatchedTPmtd2WrongAssocPt_,
                                                     meETLTrackMatchedTPmtd2WrongAssocMVAQual_,
                                                     meETLTrackMatchedTPmtd2WrongAssocTimeRes_,
                                                     meETLTrackMatchedTPmtd2WrongAssocTimePull_,
                                                     std::abs(trackGen.eta()),
                                                     trackGen.pt(),
                                                     mtdQualMVA[trackref],
                                                     dT,
                                                     pullT,
                                                     hasTime);
                }
              }
            }
            // -- Missing association with reco hits in MTD
            else {
              // -- Track matched to TP with sim hit in >=1 etl layers, no reco hits associated to the track
              if (isTPmtdETLD1 || isTPmtdETLD2) {
                meETLTrackMatchedTPmtd1NoAssocEta_->Fill(std::abs(trackGen.eta()));
                meETLTrackMatchedTPmtd1NoAssocPt_->Fill(trackGen.pt());
              }
              // -- Track matched to TP with sim hit in 2 etl layers, no reco hits associated to the track
              if (isTPmtdETLD1 && isTPmtdETLD2) {
                meETLTrackMatchedTPmtd2NoAssocEta_->Fill(std::abs(trackGen.eta()));
                meETLTrackMatchedTPmtd2NoAssocPt_->Fill(trackGen.pt());
              }
            }
          }  // == end ETL
        }  // --- end "withMTD"
        else{
          isTPmtdDirectCorrectBTL_vec.push_back(false);
          isTPmtdOtherCorrectBTL_vec.push_back(false);
          isTPmtdDirectBTL_vec.push_back(false);
          isTPmtdOtherBTL_vec.push_back(false);
          isTPmtdCorrectETLD1_vec.push_back(false);
          isTPmtdCorrectETLD2_vec.push_back(false);
          ETLdisc1_vec.push_back(false);
          ETLdisc2_vec.push_back(false);
          isTPmtdETLD1_vec.push_back(false);
          isTPmtdETLD2_vec.push_back(false);
          simCluster_time_vec.push_back(-99.);  // if no sim cluster associated, fill with -9999
          simCluster_time_ETLD1_vec.push_back(-99.);  // if no sim cluster associated, fill with -9999
          simCluster_time_ETLD2_vec.push_back(-99.);
          pl_SC_ETLD1_vec.push_back(-99.);
          pl_SC_ETLD1_earliesthit_vec.push_back(-99.);
          p_SC_ETLD1_vec.push_back(-99.);
          p_SC_ETLD1_earliesthit_vec.push_back(-99.);
          
          time_SC_ETLD1_vec.push_back(-99.);
          time_SC_ETLD1_earliesthit_vec.push_back(-99.);
          pl_SC_vec.push_back(-99.);
          pl_SC_earliesthit_vec.push_back(-99.);
          p_SC_vec.push_back(-99.);
          p_SC_earliesthit_vec.push_back(-99.);
          n_simHits_vec.push_back(-1.);
          
          time_SC_vec.push_back(-99.);
          time_SC_earliesthit_vec.push_back(-99.);
          


        }
        MCTOF_vec.push_back(MCTOF-tsim);
        MCTOF_ETLD1_vec.push_back(MCTOF_ETLD1-tsim);
        MCTOF_ETLD2_vec.push_back(MCTOF_ETLD2-tsim);


        // - Track matched to TP without sim hit in MTD, but with reco cluster associated
        // - BTL
        if (std::abs(trackGen.eta()) < trackMaxBtlEta_) {
          if (!isTPmtdDirectBTL && !isTPmtdOtherBTL) {
            meBTLTrackMatchedTPnomtdEta_->Fill(std::abs(trackGen.eta()));
            meBTLTrackMatchedTPnomtdPt_->Fill(trackGen.pt());
            if (isBTL) {
              trackMatchedtoTPnosimyesReco_vec.push_back(true);
              fillTrackClusterMatchingHistograms(meBTLTrackMatchedTPnomtdAssocEta_,
                                                 meBTLTrackMatchedTPnomtdAssocPt_,
                                                 meBTLTrackMatchedTPnomtdAssocMVAQual_,
                                                 meBTLTrackMatchedTPnomtdAssocTimeRes_,
                                                 meBTLTrackMatchedTPnomtdAssocTimePull_,
                                                 std::abs(trackGen.eta()),
                                                 trackGen.pt(),
                                                 mtdQualMVA[trackref],
                                                 dT,
                                                 pullT,
                                                 hasTime);
            }else{
              trackMatchedtoTPnosimyesReco_vec.push_back(false);
            }
          }else{
            trackMatchedtoTPnosimyesReco_vec.push_back(false);
          }
        }
        // - ETL
        else if (!isTPmtdETLD1 && !isTPmtdETLD2) {
          meETLTrackMatchedTPnomtdEta_->Fill(std::abs(trackGen.eta()));
          meETLTrackMatchedTPnomtdPt_->Fill(trackGen.pt());
          if (isETL) {
            trackMatchedtoTPnosimyesReco_vec.push_back(true);
            fillTrackClusterMatchingHistograms(meETLTrackMatchedTPnomtdAssocEta_,
                                               meETLTrackMatchedTPnomtdAssocPt_,
                                               meETLTrackMatchedTPnomtdAssocMVAQual_,
                                               meETLTrackMatchedTPnomtdAssocTimeRes_,
                                               meETLTrackMatchedTPnomtdAssocTimePull_,
                                               std::abs(trackGen.eta()),
                                               trackGen.pt(),
                                               mtdQualMVA[trackref],
                                               dT,
                                               pullT,
                                               hasTime);
          }else{
            trackMatchedtoTPnosimyesReco_vec.push_back(false);
          }
        }else{
          trackMatchedtoTPnosimyesReco_vec.push_back(false);
        }

        // == Time pull and detailed extrapolation check only on tracks associated to TP from signal event
        if (!trkTPSelLV(**tp_info)) {
          continue;
        }
        size_t nlayers(0);
        float extrho(0.);
        float exteta(0.);
        float extphi(0.);
        float selvar(0.);
        auto accept = checkAcceptance(trackGen, iEvent, iSetup, nlayers, extrho, exteta, extphi, selvar);
        if (accept.first && std::abs(exteta) < trackMaxBtlEta_) {
          meExtraPhiAtBTL_->Fill(angle_units::operators::convertRadToDeg(extphi));
          meExtraBTLeneInCone_->Fill(selvar);
        }
        if (accept.second) {
          if (std::abs(exteta) < trackMaxBtlEta_) {
            meExtraPhiAtBTLmatched_->Fill(angle_units::operators::convertRadToDeg(extphi));
          }
          if (noCrack) {
            meExtraPtMtd_->Fill(trackGen.pt());
            if (nlayers == 2) {
              meExtraPtEtl2Mtd_->Fill(trackGen.pt());
            }
          }
          meExtraEtaMtd_->Fill(std::abs(trackGen.eta()));
          if (nlayers == 2) {
            meExtraEtaEtl2Mtd_->Fill(std::abs(trackGen.eta()));
          }
          if (accept.first && accept.second && !(isBTL || isETL)) {
            edm::LogInfo("MtdTracksValidation")
                << "MtdTracksValidation: extender fail in " << iEvent.id().run() << " " << iEvent.id().event()
                << " pt= " << trackGen.pt() << " eta= " << trackGen.eta();
            meExtraMTDfailExtenderEta_->Fill(std::abs(trackGen.eta()));
            if (noCrack) {
              meExtraMTDfailExtenderPt_->Fill(trackGen.pt());
            }
          }
        }  // detailed extrapolation check

        // time res and time pull
        if (Sigmat0Safe[trackref] != -1.) {
          if (isBTL || isETL) {
            meTrackResTot_->Fill(dT);
            meTrackPullTot_->Fill(pullT);
            meTrackResTotvsMVAQual_->Fill(mtdQualMVA[trackref], dT);
            meTrackPullTotvsMVAQual_->Fill(mtdQualMVA[trackref], pullT);
          }
        }  // time res and time pull
      }  // TP matching
    }  // trkRecSel

    // ETL tracks with low pt (0.2 < Pt [GeV] < 0.7)
    if (trkRecSelLowPt(trackGen)) {
      if ((std::abs(trackGen.eta()) > trackMinEtlEta_) && (std::abs(trackGen.eta()) < trackMaxEtlEta_)) {
        if (trackGen.pt() < 0.45) {
          meETLTrackEtaTotLowPt_[0]->Fill(std::abs(trackGen.eta()));
        } else {
          meETLTrackEtaTotLowPt_[1]->Fill(std::abs(trackGen.eta()));
        }
      }
      bool MTDEtlZnegD1 = false;
      bool MTDEtlZnegD2 = false;
      bool MTDEtlZposD1 = false;
      bool MTDEtlZposD2 = false;
      for (const auto hit : track.recHits()) {
        if (hit->isValid() == false)
          continue;
        MTDDetId Hit = hit->geographicalId();
        if ((Hit.det() == 6) && (Hit.subdetId() == 1) && (Hit.mtdSubDetector() == 2)) {
          isETL = true;
          ETLDetId ETLHit = hit->geographicalId();
          if ((ETLHit.zside() == -1) && (ETLHit.nDisc() == 1)) {
            MTDEtlZnegD1 = true;
          }
          if ((ETLHit.zside() == -1) && (ETLHit.nDisc() == 2)) {
            MTDEtlZnegD2 = true;
          }
          if ((ETLHit.zside() == 1) && (ETLHit.nDisc() == 1)) {
            MTDEtlZposD1 = true;
          }
          if ((ETLHit.zside() == 1) && (ETLHit.nDisc() == 2)) {
            MTDEtlZposD2 = true;
          }
        }
      }
      if ((trackGen.eta() < -trackMinEtlEta_) && (trackGen.eta() > -trackMaxEtlEta_)) {
        twoETLdiscs = (MTDEtlZnegD1 == true) && (MTDEtlZnegD2 == true);
      }
      if ((trackGen.eta() > trackMinEtlEta_) && (trackGen.eta() < trackMaxEtlEta_)) {
        twoETLdiscs = (MTDEtlZposD1 == true) && (MTDEtlZposD2 == true);
      }
      if (isETL && (std::abs(trackGen.eta()) > trackMinEtlEta_) && (std::abs(trackGen.eta()) < trackMaxEtlEta_)) {
        if (trackGen.pt() < 0.45) {
          meETLTrackEtaMtdLowPt_[0]->Fill(std::abs(trackGen.eta()));
        } else {
          meETLTrackEtaMtdLowPt_[1]->Fill(std::abs(trackGen.eta()));
        }
      }
      if (isETL && twoETLdiscs) {
        if (trackGen.pt() < 0.45) {
          meETLTrackEta2MtdLowPt_[0]->Fill(std::abs(trackGen.eta()));
        } else {
          meETLTrackEta2MtdLowPt_[1]->Fill(std::abs(trackGen.eta()));
        }
      }
    }  // trkRecSelLowPt

  }  // RECO tracks loop
  dump_tree->Fill();
  ClearVectors();
}

const std::pair<bool, bool> MtdTracksValidation::checkAcceptance(const reco::Track& track,
                                                                 const edm::Event& iEvent,
                                                                 edm::EventSetup const& iSetup,
                                                                 size_t& nlayers,
                                                                 float& extrho,
                                                                 float& exteta,
                                                                 float& extphi,
                                                                 float& selvar) {
  bool isMatched(false);
  nlayers = 0;
  extrho = 0.;
  exteta = -999.;
  extphi = -999.;
  selvar = 0.;

  auto geometryHandle = iSetup.getTransientHandle(mtdgeoToken_);
  const MTDGeometry* geom = geometryHandle.product();
  auto topologyHandle = iSetup.getTransientHandle(mtdtopoToken_);
  const MTDTopology* topology = topologyHandle.product();

  auto layerHandle = iSetup.getTransientHandle(mtdlayerToken_);
  const MTDDetLayerGeometry* layerGeo = layerHandle.product();

  auto magfieldHandle = iSetup.getTransientHandle(magfieldToken_);
  const MagneticField* mfield = magfieldHandle.product();

  auto ttrackBuilder = iSetup.getTransientHandle(builderToken_);

  auto tTrack = ttrackBuilder->build(track);
  TrajectoryStateOnSurface tsos = tTrack.outermostMeasurementState();
  float theMaxChi2 = 500.;
  float theNSigma = 10.;
  std::unique_ptr<MeasurementEstimator> theEstimator =
      std::make_unique<Chi2MeasurementEstimator>(theMaxChi2, theNSigma);
  SteppingHelixPropagator prop(mfield, anyDirection);

  auto btlRecHitsHandle = makeValid(iEvent.getHandle(btlRecHitsToken_));
  auto etlRecHitsHandle = makeValid(iEvent.getHandle(etlRecHitsToken_));

  edm::LogVerbatim("MtdTracksValidation")
      << "MtdTracksValidation: extrapolating track, pt= " << track.pt() << " eta= " << track.eta();

  //try BTL
  bool inBTL = false;
  float eneSum(0.);
  const std::vector<const DetLayer*>& layersBTL = layerGeo->allBTLLayers();
  for (const DetLayer* ilay : layersBTL) {
    std::pair<bool, TrajectoryStateOnSurface> comp = ilay->compatible(tsos, prop, *theEstimator);
    if (!comp.first)
      continue;
    if (!inBTL) {
      inBTL = true;
      extrho = comp.second.globalPosition().perp();
      exteta = comp.second.globalPosition().eta();
      extphi = comp.second.globalPosition().phi();
      edm::LogVerbatim("MtdTracksValidation") << "MtdTracksValidation: extrapolation at BTL surface, rho= " << extrho
                                              << " eta= " << exteta << " phi= " << extphi;
    }
    std::vector<DetLayer::DetWithState> compDets = ilay->compatibleDets(tsos, prop, *theEstimator);
    for (const auto& detWithState : compDets) {
      const auto& det = detWithState.first;

      // loop on compatible rechits and check energy in a fixed size cone around the extrapolation point

      edm::LogVerbatim("MtdTracksValidation")
          << "MtdTracksValidation: DetId= " << det->geographicalId().rawId()
          << " gp= " << detWithState.second.globalPosition().x() << " " << detWithState.second.globalPosition().y()
          << " " << detWithState.second.globalPosition().z() << " rho= " << detWithState.second.globalPosition().perp()
          << " eta= " << detWithState.second.globalPosition().eta()
          << " phi= " << detWithState.second.globalPosition().phi();

      for (const auto& recHit : *btlRecHitsHandle) {
        BTLDetId detId = recHit.id();
        DetId geoId = detId.geographicalId(MTDTopologyMode::crysLayoutFromTopoMode(topology->getMTDTopologyMode()));
        const MTDGeomDet* thedet = geom->idToDet(geoId);
        if (thedet == nullptr)
          throw cms::Exception("MtdTracksValidation") << "GeographicalID: " << std::hex << geoId.rawId() << " ("
                                                      << detId.rawId() << ") is invalid!" << std::dec << std::endl;
        if (geoId == det->geographicalId()) {
          const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(thedet->topology());
          const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());

          Local3DPoint local_point(0., 0., 0.);
          local_point = topo.pixelToModuleLocalPoint(local_point, detId.row(topo.nrows()), detId.column(topo.nrows()));
          const auto& global_point = thedet->toGlobal(local_point);
          edm::LogVerbatim("MtdTracksValidation")
              << "MtdTracksValidation: Hit id= " << detId.rawId() << " ene= " << recHit.energy()
              << " dr= " << reco::deltaR(global_point, detWithState.second.globalPosition());
          if (reco::deltaR(global_point, detWithState.second.globalPosition()) < cluDRradius_) {
            eneSum += recHit.energy();
            //extrho = detWithState.second.globalPosition().perp();
            //exteta = detWithState.second.globalPosition().eta();
            //extphi = detWithState.second.globalPosition().phi();
          }
        }
      }
    }
    if (eneSum > depositBTLthreshold_) {
      nlayers++;
      selvar = eneSum;
      isMatched = true;
      edm::LogVerbatim("MtdTracksValidation")
          << "MtdTracksValidation: BTL matched, energy= " << eneSum << " #layers= " << nlayers;
    }
  }
  if (inBTL) {
    return std::make_pair(inBTL, isMatched);
  }

  //try ETL
  bool inETL = false;
  const std::vector<const DetLayer*>& layersETL = layerGeo->allETLLayers();
  for (const DetLayer* ilay : layersETL) {
    size_t hcount(0);
    const BoundDisk& disk = static_cast<const MTDSectorForwardDoubleLayer*>(ilay)->specificSurface();
    const double diskZ = disk.position().z();
    if (tsos.globalPosition().z() * diskZ < 0)
      continue;  // only propagate to the disk that's on the same side
    std::pair<bool, TrajectoryStateOnSurface> comp = ilay->compatible(tsos, prop, *theEstimator);
    if (!comp.first)
      continue;
    if (!inETL) {
      inETL = true;
      extrho = comp.second.globalPosition().perp();
      exteta = comp.second.globalPosition().eta();
      extphi = comp.second.globalPosition().phi();
    }
    edm::LogVerbatim("MtdTracksValidation") << "MtdTracksValidation: extrapolation at ETL surface, rho= " << extrho
                                            << " eta= " << exteta << " phi= " << extphi;
    std::vector<DetLayer::DetWithState> compDets = ilay->compatibleDets(tsos, prop, *theEstimator);
    for (const auto& detWithState : compDets) {
      const auto& det = detWithState.first;

      // loop on compatible rechits and check hits in a fixed size cone around the extrapolation point

      edm::LogVerbatim("MtdTracksValidation")
          << "MtdTracksValidation: DetId= " << det->geographicalId().rawId()
          << " gp= " << detWithState.second.globalPosition().x() << " " << detWithState.second.globalPosition().y()
          << " " << detWithState.second.globalPosition().z() << " rho= " << detWithState.second.globalPosition().perp()
          << " eta= " << detWithState.second.globalPosition().eta()
          << " phi= " << detWithState.second.globalPosition().phi();

      for (const auto& recHit : *etlRecHitsHandle) {
        ETLDetId detId = recHit.id();
        DetId geoId = detId.geographicalId();
        const MTDGeomDet* thedet = geom->idToDet(geoId);
        if (thedet == nullptr)
          throw cms::Exception("MtdTracksValidation") << "GeographicalID: " << std::hex << geoId.rawId() << " ("
                                                      << detId.rawId() << ") is invalid!" << std::dec << std::endl;
        if (geoId == det->geographicalId()) {
          const ProxyMTDTopology& topoproxy = static_cast<const ProxyMTDTopology&>(thedet->topology());
          const RectangularMTDTopology& topo = static_cast<const RectangularMTDTopology&>(topoproxy.specificTopology());

          Local3DPoint local_point(topo.localX(recHit.row()), topo.localY(recHit.column()), 0.);
          const auto& global_point = thedet->toGlobal(local_point);
          edm::LogVerbatim("MtdTracksValidation")
              << "MtdTracksValidation: Hit id= " << detId.rawId() << " time= " << recHit.time()
              << " dr= " << reco::deltaR(global_point, detWithState.second.globalPosition());
          if (reco::deltaR(global_point, detWithState.second.globalPosition()) < cluDRradius_) {
            hcount++;
            if (hcount == 1) {
              //extrho = detWithState.second.globalPosition().perp();
              //exteta = detWithState.second.globalPosition().eta();
              //extphi = detWithState.second.globalPosition().phi();
            }
          }
        }
      }
    }
    if (hcount > 0) {
      nlayers++;
      selvar = (float)hcount;
      isMatched = true;
      edm::LogVerbatim("MtdTracksValidation")
          << "MtdTracksValidation: ETL matched, counts= " << hcount << " #layers= " << nlayers;
    }
  }

  if (!inBTL && !inETL) {
    edm::LogVerbatim("MtdTracksValidation")
        << "MtdTracksValidation: track not extrapolating to MTD: pt= " << track.pt() << " eta= " << track.eta()
        << " phi= " << track.phi() << " vz= " << track.vz()
        << " vxy= " << std::sqrt(track.vx() * track.vx() + track.vy() * track.vy());
  }
  return std::make_pair(inETL, isMatched);
}

// ------------ method for histogram booking ------------
void MtdTracksValidation::bookHistograms(DQMStore::IBooker& ibook, edm::Run const& run, edm::EventSetup const& iSetup) {
  ibook.setCurrentFolder(folder_);

  // histogram booking
  meBTLTrackRPTime_ = ibook.book1D("TrackBTLRPTime", "Track t0 with respect to R.P.;t0 [ns]", 100, -1, 3);
  meBTLTrackEtaTot_ = ibook.book1D("TrackBTLEtaTot", "Eta of tracks (Tot);#eta_{RECO}", 30, 0., 1.5);
  meBTLTrackPhiTot_ = ibook.book1D("TrackBTLPhiTot", "Phi of tracks (Tot);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meBTLTrackPtTot_ = ibook.book1D("TrackBTLPtTot", "Pt of tracks (Tot);pt_{RECO} [GeV]", 50, 0, 10);
  meBTLTrackEtaMtd_ = ibook.book1D("TrackBTLEtaMtd", "Eta of tracks (Mtd);#eta_{RECO}", 30, 0., 1.5);
  meBTLTrackPhiMtd_ = ibook.book1D("TrackBTLPhiMtd", "Phi of tracks (Mtd);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meBTLTrackPtMtd_ = ibook.book1D("TrackBTLPtMtd", "Pt of tracks (Mtd);pt_{RECO} [GeV]", 50, 0, 10);
  meBTLTrackPtRes_ =
      ibook.book1D("TrackBTLPtRes", "Track pT resolution  ;pT_{Gentrack}-pT_{MTDtrack}/pT_{Gentrack} ", 100, -0.1, 0.1);
  meETLTrackRPTime_ = ibook.book1D("TrackETLRPTime", "Track t0 with respect to R.P.;t0 [ns]", 100, -1, 3);
  meETLTrackEtaTot_ = ibook.book1D("TrackETLEtaTot", "Eta of tracks (Tot);#eta_{RECO}", 30, 1.5, 3.0);
  meETLTrackPhiTot_ = ibook.book1D("TrackETLPhiTot", "Phi of tracks (Tot);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meETLTrackPhiTot_ = ibook.book1D("TrackETLPhiTot", "Phi of tracks (Tot);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meETLTrackPtTot_ = ibook.book1D("TrackETLPtTot", "Pt of tracks (Tot);pt_{RECO} [GeV]", 50, 0, 10);

  meETLTrackEtaTotLowPt_[0] =
      ibook.book1D("TrackETLEtaTotLowPt0", "Eta of tracks, 0.2 < pt < 0.45 (Tot);#eta_{RECO}", 30, 1.5, 3.0);
  meETLTrackEtaTotLowPt_[1] =
      ibook.book1D("TrackETLEtaTotLowPt1", "Eta of tracks, 0.45 < pt < 0.7 (Tot);#eta_{RECO}", 30, 1.5, 3.0);

  meETLTrackEtaMtd_ = ibook.book1D("TrackETLEtaMtd", "Eta of tracks (Mtd);#eta_{RECO}", 30, 1.5, 3.0);
  meETLTrackEtaMtdLowPt_[0] =
      ibook.book1D("TrackETLEtaMtdLowPt0", "Eta of tracks, 0.2 < pt < 0.45 (Mtd);#eta_{RECO}", 30, 1.5, 3.0);
  meETLTrackEtaMtdLowPt_[1] =
      ibook.book1D("TrackETLEtaMtdLowPt1", "Eta of tracks, 0.45 < pt < 0.7 (Mtd);#eta_{RECO}", 30, 1.5, 3.0);
  meETLTrackEta2MtdLowPt_[0] =
      ibook.book1D("TrackETLEta2MtdLowPt0", "Eta of tracks, 0.2 < pt < 0.45 (Mtd 2 hit);#eta_{RECO}", 30, 1.5, 3.0);
  meETLTrackEta2MtdLowPt_[1] =
      ibook.book1D("TrackETLEta2MtdLowPt1", "Eta of tracks, 0.45 < pt < 0.7 (Mtd 2 hit);#eta_{RECO}", 30, 1.5, 3.0);

  meETLTrackPhiMtd_ = ibook.book1D("TrackETLPhiMtd", "Phi of tracks (Mtd);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meETLTrackPtMtd_ = ibook.book1D("TrackETLPtMtd", "Pt of tracks (Mtd);pt_{RECO} [GeV]", 50, 0, 10);
  meETLTrackEta2Mtd_ = ibook.book1D("TrackETLEta2Mtd", "Eta of tracks (Mtd 2 hit);#eta_{RECO}", 30, 1.5, 3.0);
  meETLTrackPhi2Mtd_ = ibook.book1D("TrackETLPhi2Mtd", "Phi of tracks (Mtd 2 hit);#phi_{RECO} [rad]", 100, -3.2, 3.2);
  meETLTrackPt2Mtd_ = ibook.book1D("TrackETLPt2Mtd", "Pt of tracks (Mtd 2 hit);pt_{RECO} [GeV]", 50, 0, 10);

  meETLTrackPtRes_ =
      ibook.book1D("TrackETLPtRes", "Track pT resolution;pT_{Gentrack}-pT_{MTDtrack}/pT_{Gentrack} ", 100, -0.1, 0.1);

  meTracktmtd_ = ibook.book1D("Tracktmtd", "Track time from TrackExtenderWithMTD;tmtd [ns]", 150, 1, 16);
  meTrackt0Src_ = ibook.book1D("Trackt0Src", "Track time from TrackExtenderWithMTD;t0Src [ns]", 100, -1.5, 1.5);
  meTrackSigmat0Src_ =
      ibook.book1D("TrackSigmat0Src", "Time Error from TrackExtenderWithMTD; #sigma_{t0Src} [ns]", 100, 0, 0.1);

  meTrackt0Pid_ = ibook.book1D("Trackt0Pid", "Track t0 as stored in TofPid;t0 [ns]", 100, -1, 1);
  meTrackSigmat0Pid_ = ibook.book1D("TrackSigmat0Pid", "Sigmat0 as stored in TofPid; #sigma_{t0} [ns]", 100, 0, 0.1);
  meTrackt0SafePid_ = ibook.book1D("Trackt0SafePID", "Track t0 Safe as stored in TofPid;t0 [ns]", 100, -1, 1);
  meTrackSigmat0SafePid_ = ibook.book1D(
      "TrackSigmat0SafePID", "Log10(Sigmat0 Safe) as stored in TofPid; Log10(#sigma_{t0} [ns])", 80, -3, 1);
  meTrackNumHits_ = ibook.book1D("TrackNumHits", "Number of valid MTD hits per track ; Number of hits", 10, -5, 5);
  meTrackNumHitsNT_ = ibook.book1D(
      "TrackNumHitsNT", "Number of valid MTD hits per track no time associated; Number of hits", 10, -5, 5);
  meTrackMVAQual_ = ibook.book1D("TrackMVAQual", "Track MVA Quality as stored in Value Map ; MVAQual", 100, -1, 1);
  meTrackPathLenghtvsEta_ = ibook.bookProfile(
      "TrackPathLenghtvsEta", "MTD Track pathlength vs MTD track Eta;|#eta|;Pathlength", 100, 0, 3.2, 100.0, 400.0, "S");

  meTrackOutermostHitR_ = ibook.book1D("TrackOutermostHitR", "Track outermost hit position R; R[cm]", 40, 0, 120.);
  meTrackOutermostHitZ_ = ibook.book1D("TrackOutermostHitZ", "Track outermost hit position Z; z[cm]", 100, 0, 300.);

  meTrackSigmaTof_[0] =
      ibook.book1D("TrackSigmaTof_Pion", "Sigma(TOF) for pion hypothesis; #sigma_{t0} [ps]", 10, 0, 5);
  meTrackSigmaTof_[1] =
      ibook.book1D("TrackSigmaTof_Kaon", "Sigma(TOF) for kaon hypothesis; #sigma_{t0} [ps]", 25, 0, 25);
  meTrackSigmaTof_[2] =
      ibook.book1D("TrackSigmaTof_Proton", "Sigma(TOF) for proton hypothesis; #sigma_{t0} [ps]", 50, 0, 50);

  meTrackSigmaTofvsP_[0] = ibook.bookProfile("TrackSigmaTofvsP_Pion",
                                             "Sigma(TOF) for pion hypothesis vs p; p [GeV]; #sigma_{t0} [ps]",
                                             20,
                                             0,
                                             10.,
                                             0,
                                             50.,
                                             "S");
  meTrackSigmaTofvsP_[1] = ibook.bookProfile("TrackSigmaTofvsP_Kaon",
                                             "Sigma(TOF) for kaon hypothesis vs p; p [GeV]; #sigma_{t0} [ps]",
                                             20,
                                             0,
                                             10.,
                                             0,
                                             50.,
                                             "S");
  meTrackSigmaTofvsP_[2] = ibook.bookProfile("TrackSigmaTofvsP_Proton",
                                             "Sigma(TOF) for proton hypothesis vs p; p [GeV]; #sigma_{t0} [ps]",
                                             20,
                                             0,
                                             10.,
                                             0,
                                             50.,
                                             "S");

  meExtraPtMtd_ =
      ibook.book1D("ExtraPtMtd", "Pt of tracks associated to LV extrapolated to hits; track pt [GeV] ", 110, 0., 11.);
  meExtraPtEtl2Mtd_ = ibook.book1D("ExtraPtEtl2Mtd",
                                   "Pt of tracks associated to LV extrapolated to hits, 2 ETL layers; track pt [GeV] ",
                                   110,
                                   0.,
                                   11.);
  meExtraEtaMtd_ =
      ibook.book1D("ExtraEtaMtd", "Eta of tracks associated to LV extrapolated to hits; track eta ", 66, 0., 3.3);
  meExtraEtaEtl2Mtd_ = ibook.book1D(
      "ExtraEtaEtl2Mtd", "Eta of tracks associated to LV extrapolated to hits, 2 ETL layers; track eta ", 66, 0., 3.3);
  meTrackMatchedTPEtaTotLV_ =
      ibook.book1D("MatchedTPEtaTotLV", "Eta of tracks associated to LV matched to TP; track eta ", 66, 0., 3.3);
  meTrackMatchedTPPtTotLV_ =
      ibook.book1D("MatchedTPPtTotLV", "Pt of tracks associated to LV matched to TP; track pt [GeV] ", 110, 0., 11.);

  meBTLTrackMatchedTPEtaTot_ =
      ibook.book1D("BTLTrackMatchedTPEtaTot", "Eta of tracks matched to TP; track eta ", 30, 0., 1.5);
  meBTLTrackMatchedTPEtaMtd_ =
      ibook.book1D("BTLTrackMatchedTPEtaMtd", "Eta of tracks matched to TP with time; track eta ", 30, 0., 1.5);
  meBTLTrackMatchedTPPtTot_ =
      ibook.book1D("BTLTrackMatchedTPPtTot", "Pt of tracks matched to TP; track pt [GeV] ", 50, 0., 10.);
  meBTLTrackMatchedTPPtMtd_ =
      ibook.book1D("BTLTrackMatchedTPPtMtd", "Pt of tracks matched to TP with time; track pt [GeV] ", 50, 0., 10.);
  meETLTrackMatchedTPEtaTot_ =
      ibook.book1D("ETLTrackMatchedTPEtaTot", "Eta of tracks matched to TP; track eta ", 30, 1.5, 3.0);
  meETLTrackMatchedTPEtaMtd_ = ibook.book1D(
      "ETLTrackMatchedTPEtaMtd", "Eta of tracks matched to TP with time (>=1 ETL hit); track eta ", 30, 1.5, 3.0);
  meETLTrackMatchedTPEtaMtdCorrect_ =
      ibook.book1D("ETLTrackMatchedTPEtaMtdCorrect",
                   "Eta of tracks matched to TP with time (>=1 ETL hit), correct reco match; track eta ",
                   30,
                   1.5,
                   3.0);
  meETLTrackMatchedTPEta2Mtd_ = ibook.book1D(
      "ETLTrackMatchedTPEta2Mtd", "Eta of tracks matched to TP with time (2 ETL hits); track eta ", 30, 1.5, 3.0);
  meETLTrackMatchedTPPtTot_ =
      ibook.book1D("ETLTrackMatchedTPPtTot", "Pt of tracks matched to TP; track pt [GeV] ", 50, 0., 10.);
  meETLTrackMatchedTPPtMtd_ = ibook.book1D(
      "ETLTrackMatchedTPPtMtd", "Pt of tracks matched to TP with time (>=1 ETL hit); track pt [GeV] ", 50, 0., 10.);
  meETLTrackMatchedTPPtMtdCorrect_ =
      ibook.book1D("ETLTrackMatchedTPPtMtdCorrect",
                   "Pt of tracks matched to TP with time (>=1 ETL hit), correct reco match; track pt [GeV] ",
                   50,
                   0.,
                   10.);
  meETLTrackMatchedTPPt2Mtd_ = ibook.book1D(
      "ETLTrackMatchedTPPt2Mtd", "Pt of tracks matched to TP with time (2 ETL hits); track pt [GeV] ", 50, 0., 10.);

  if (optionalPlots_) {
    meBTLTrackMatchedTPPtResMtd_ = ibook.book1D(
        "TrackMatchedTPBTLPtResMtd",
        "Pt resolution of tracks matched to TP-BTL hit  ;|pT_{MTDtrack}-pT_{truth}|/|pT_{Gentrack}-pT_{truth}| ",
        100,
        0.,
        4.);
    meETLTrackMatchedTPPtResMtd_ = ibook.book1D(
        "TrackMatchedTPETLPtResMtd",
        "Pt resolution of tracks matched to TP-ETL hit  ;|pT_{MTDtrack}-pT_{truth}|/|pT_{Gentrack}-pT_{truth}| ",
        100,
        0.,
        4.);
    meETLTrackMatchedTP2PtResMtd_ = ibook.book1D(
        "TrackMatchedTPETL2PtResMtd",
        "Pt resolution of tracks matched to TP-ETL 2hits  ;|pT_{MTDtrack}-pT_{truth}|/|pT_{Gentrack}-pT_{truth}| ",
        100,
        0.,
        4.);
    meBTLTrackMatchedTPPtRatioGen_ = ibook.book1D(
        "TrackMatchedTPBTLPtRatioGen", "Pt ratio of Gentracks (BTL)  ;pT_{Gentrack}/pT_{truth} ", 100, 0.9, 1.1);
    meETLTrackMatchedTPPtRatioGen_ = ibook.book1D(
        "TrackMatchedTPETLPtRatioGen", "Pt ratio of Gentracks (ETL 1hit)  ;pT_{Gentrack}/pT_{truth} ", 100, 0.9, 1.1);
    meETLTrackMatchedTP2PtRatioGen_ = ibook.book1D(
        "TrackMatchedTPETL2PtRatioGen", "Pt ratio of Gentracks (ETL 2hits)  ;pT_{Gentrack}/pT_{truth} ", 100, 0.9, 1.1);
    meBTLTrackMatchedTPPtRatioMtd_ =
        ibook.book1D("TrackMatchedTPBTLPtRatioMtd",
                     "Pt ratio of tracks matched to TP-BTL hit  ;pT_{MTDtrack}/pT_{truth} ",
                     100,
                     0.9,
                     1.1);
    meETLTrackMatchedTPPtRatioMtd_ =
        ibook.book1D("TrackMatchedTPETLPtRatioMtd",
                     "Pt ratio of tracks matched to TP-ETL hit  ;pT_{MTDtrack}/pT_{truth} ",
                     100,
                     0.9,
                     1.1);
    meETLTrackMatchedTP2PtRatioMtd_ =
        ibook.book1D("TrackMatchedTPETL2PtRatioMtd",
                     "Pt ratio of tracks matched to TP-ETL 2hits  ;pT_{MTDtrack}/pT_{truth} ",
                     100,
                     0.9,
                     1.1);
    meBTLTrackMatchedTPPtResvsPtMtd_ =
        ibook.bookProfile("TrackMatchedTPBTLPtResvsPtMtd",
                          "Pt resolution of tracks matched to TP-BTL hit vs Pt;pT_{truth} "
                          "[GeV];|pT_{MTDtrack}-pT_{truth}|/|pT_{Gentrack}-pT_{truth}| ",
                          20,
                          0.7,
                          10.,
                          0.,
                          4.,
                          "s");
    meETLTrackMatchedTPPtResvsPtMtd_ =
        ibook.bookProfile("TrackMatchedTPETLPtResvsPtMtd",
                          "Pt resolution of tracks matched to TP-ETL hit vs Pt;pT_{truth} "
                          "[GeV];|pT_{MTDtrack}-pT_{truth}|/|pT_{Gentrack}-pT_{truth}| ",
                          20,
                          0.7,
                          10.,
                          0.,
                          4.,
                          "s");
    meETLTrackMatchedTP2PtResvsPtMtd_ =
        ibook.bookProfile("TrackMatchedTPETL2PtResvsPtMtd",
                          "Pt resolution of tracks matched to TP-ETL 2hits Pt pT;pT_{truth} "
                          "[GeV];|pT_{MTDtrack}-pT_{truth}|/|pT_{Gentrack}-pT_{truth}| ",
                          20,
                          0.7,
                          10.,
                          0.,
                          4.,
                          "s");
    meBTLTrackMatchedTPDPtvsPtGen_ = ibook.bookProfile(
        "TrackMatchedTPBTLDPtvsPtGen",
        "Pt relative difference of Gentracks (BTL) vs Pt;pT_{truth} [GeV];pT_{Gentrack}-pT_{truth}/pT_{truth} ",
        20,
        0.7,
        10.,
        -0.1,
        0.1,
        "s");
    meETLTrackMatchedTPDPtvsPtGen_ = ibook.bookProfile(
        "TrackMatchedTPETLDPtvsPtGen",
        "Pt relative difference of Gentracks (ETL 1hit) vs Pt;pT_{truth} [GeV];pT_{Gentrack}-pT_{truth}/pT_{truth} ",
        20,
        0.7,
        10.,
        -0.1,
        0.1,
        "s");
    meETLTrackMatchedTP2DPtvsPtGen_ = ibook.bookProfile(
        "TrackMatchedTPETL2DPtvsPtGen",
        "Pt relative difference  of Gentracks (ETL 2hits) vs Pt;pT_{truth} [GeV];pT_{Gentrack}-pT_{truth}/pT_{truth} ",
        20,
        0.7,
        10.,
        -0.1,
        0.1,
        "s");
    meBTLTrackMatchedTPDPtvsPtMtd_ = ibook.bookProfile("TrackMatchedTPBTLDPtvsPtMtd",
                                                       "Pt relative difference of tracks matched to TP-BTL hits vs "
                                                       "Pt;pT_{truth} [GeV];pT_{MTDtrack}-pT_{truth}/pT_{truth} ",
                                                       20,
                                                       0.7,
                                                       10.,
                                                       -0.1,
                                                       0.1,
                                                       "s");
    meETLTrackMatchedTPDPtvsPtMtd_ = ibook.bookProfile("TrackMatchedTPETLDPtvsPtMtd",
                                                       "Pt relative difference of tracks matched to TP-ETL hits vs "
                                                       "Pt;pT_{truth} [GeV];pT_{MTDtrack}-pT_{truth}/pT_{truth} ",
                                                       20,
                                                       0.7,
                                                       10.,
                                                       -0.1,
                                                       0.1,
                                                       "s");
    meETLTrackMatchedTP2DPtvsPtMtd_ = ibook.bookProfile("TrackMatchedTPETL2DPtvsPtMtd",
                                                        "Pt relative difference of tracks matched to TP-ETL 2hits vs "
                                                        "Pt;pT_{truth} [GeV];pT_{MTDtrack}-pT_{truth}/pT_{truth} ",
                                                        20,
                                                        0.7,
                                                        10.,
                                                        -0.1,
                                                        0.1,
                                                        "s");
  }  // end optional plots

  meTrackResTot_ = ibook.book1D(
      "TrackRes", "t_{rec} - t_{sim} for LV associated tracks matched to TP; t_{rec} - t_{sim} [ns] ", 120, -0.15, 0.15);
  meTrackPullTot_ = ibook.book1D(
      "TrackPull", "Pull for LV associated tracks matched to TP; (t_{rec}-t_{sim})/#sigma_{t}", 50, -5., 5.);
  meTrackResTotvsMVAQual_ = ibook.bookProfile(
      "TrackResvsMVA",
      "t_{rec} - t_{sim} for LV associated tracks matched to TP vs MVA Quality; MVAQual; t_{rec} - t_{sim} [ns] ",
      100,
      -1.,
      1.,
      -0.15,
      0.15,
      "s");
  meTrackPullTotvsMVAQual_ = ibook.bookProfile(
      "TrackPullvsMVA",
      "Pull for LV associated tracks matched to TP vs MVA Quality; MVAQual; (t_{rec}-t_{sim})/#sigma_{t}",
      100,
      -1.,
      1.,
      -5.,
      5.,
      "s");

  meExtraPhiAtBTL_ = ibook.book1D(
      "ExtraPhiAtBTL", "Phi at BTL surface of extrapolated tracks associated to LV; phi [deg]", 720, -180., 180.);
  meExtraPhiAtBTLmatched_ =
      ibook.book1D("ExtraPhiAtBTLmatched",
                   "Phi at BTL surface of extrapolated tracks associated to LV matched with BTL hits; phi [deg]",
                   720,
                   -180.,
                   180.);
  meExtraBTLeneInCone_ =
      ibook.book1D("ExtraBTLeneInCone",
                   "BTL reconstructed energy in cone arounnd extrapolated track associated to LV; E [MeV]",
                   100,
                   0.,
                   50.);
  meExtraMTDfailExtenderEta_ =
      ibook.book1D("ExtraMTDfailExtenderEta",
                   "Eta of tracks associated to LV extrapolated to MTD with no track extender match to hits; track eta",
                   66,
                   0.,
                   3.3);
  ;
  meExtraMTDfailExtenderPt_ = ibook.book1D(
      "ExtraMTDfailExtenderPt",
      "Pt of tracks associated to LV extrapolated to MTD with no track extender match to hits; track pt [GeV] ",
      110,
      0.,
      11.);

  // Book the histograms for track-hit matching based on MC truth
  meBTLTrackMatchedTPmtdDirectEta_ =
      ibook.book1D("BTLTrackMatchedTPmtdDirectEta",
                   "Eta of tracks matched to TP with sim hit in MTD (direct);#eta_{RECO}",
                   30,
                   0.,
                   1.5);
  meBTLTrackMatchedTPmtdDirectPt_ =
      ibook.book1D("BTLTrackMatchedTPmtdDirectPt",
                   "Pt of tracks matched to TP with sim hit in MTD (direct); track pt [GeV]",
                   50,
                   0.,
                   10.);

  meBTLTrackMatchedTPmtdOtherEta_ = ibook.book1D("BTLTrackMatchedTPmtdOtherEta",
                                                 "Eta of tracks matched to TP with sim hit in MTD (other);#eta_{RECO}",
                                                 30,
                                                 0.,
                                                 1.5);
  meBTLTrackMatchedTPmtdOtherPt_ =
      ibook.book1D("BTLTrackMatchedTPmtdOtherPt",
                   "Pt of tracks matched to TP with sim hit in MTD (other); track pt [GeV]",
                   50,
                   0.,
                   10.);

  meBTLTrackMatchedTPnomtdEta_ = ibook.book1D(
      "BTLTrackMatchedTPnomtdEta", "Eta of tracks matched to TP w/o sim hit in MTD;#eta_{RECO}", 30, 0., 1.5);
  meBTLTrackMatchedTPnomtdPt_ = ibook.book1D(
      "BTLTrackMatchedTPnomtdPt", "Pt of tracks matched to TP w/o sim hit in MTD; track pt [GeV]", 50, 0., 10.);

  meBTLTrackMatchedTPmtdDirectCorrectAssocEta_ = ibook.book1D(
      "BTLTrackMatchedTPmtdDirectCorrectAssocEta",
      "Eta of tracks matched to TP with sim hit in MTD (direct), correct track-MTD association;#eta_{RECO}",
      30,
      0.,
      1.5);
  meBTLTrackMatchedTPmtdDirectCorrectAssocPt_ = ibook.book1D(
      "BTLTrackMatchedTPmtdDirectCorrectAssocPt",
      "Pt of tracks matched to TP with sim hit in MTD (direct) - correct track-MTD association;track pt [GeV]",
      50,
      0.,
      10.);
  meBTLTrackMatchedTPmtdDirectCorrectAssocMVAQual_ = ibook.book1D(
      "BTLTrackMatchedTPmtdDirectCorrectAssocMVAQual",
      "MVA of tracks matched to TP with sim hit in MTD (direct) - correct track-MTD association; MVA score",
      100,
      -1.,
      1.);
  meBTLTrackMatchedTPmtdDirectCorrectAssocTimeRes_ =
      ibook.book1D("BTLTrackMatchedTPmtdDirectCorrectAssocTimeRes",
                   "Time resolution of tracks matched to TP with sim hit in MTD (direct) - correct track-MTD "
                   "association; t_{rec} - t_{sim} [ns] ",
                   120,
                   -0.15,
                   0.15);
  meBTLTrackMatchedTPmtdDirectCorrectAssocTimePull_ =
      ibook.book1D("BTLTrackMatchedTPmtdDirectCorrectAssocTimePull",
                   "Time pull of tracks matched to TP with sim hit in MTD (direct) - correct track-MTD association; "
                   "(t_{rec}-t_{sim})/#sigma_{t}",
                   50,
                   -5.,
                   5.);

  meBTLTrackMatchedTPmtdDirectWrongAssocEta_ =
      ibook.book1D("BTLTrackMatchedTPmtdDirectWrongAssocEta",
                   "Eta of tracks matched to TP with sim hit in MTD (direct) - wrong track-MTD association;#eta_{RECO}",
                   30,
                   0.,
                   1.5);
  meBTLTrackMatchedTPmtdDirectWrongAssocPt_ = ibook.book1D(
      "BTLTrackMatchedTPmtdDirectWrongAssocPt",
      "Pt of tracks matched to TP with sim hit in MTD (direct) - wrong track-MTD association;track pt [GeV]",
      50,
      0.,
      10.);
  meBTLTrackMatchedTPmtdDirectWrongAssocMVAQual_ =
      ibook.book1D("BTLTrackMatchedTPmtdDirectWrongAssocMVAQual",
                   "MVA of tracks matched to TP with sim hit in MTD (direct) - wrong track-MTD association; MVA score",
                   100,
                   -1.,
                   1.);
  meBTLTrackMatchedTPmtdDirectWrongAssocTimeRes_ =
      ibook.book1D("BTLTrackMatchedTPmtdDirectWrongAssocTimeRes",
                   "Time resolution of tracks matched to TP with sim hit in MTD (direct) - wrong track-MTD "
                   "association; t_{rec} - t_{sim} [ns] ",
                   120,
                   -0.15,
                   0.15);
  meBTLTrackMatchedTPmtdDirectWrongAssocTimePull_ =
      ibook.book1D("BTLTrackMatchedTPmtdDirectWrongAssocTimePull",
                   "Time pull of tracks matched to TP with sim hit in MTD (direct) - wrong track-MTD association; "
                   "(t_{rec}-t_{sim})/#sigma_{t}",
                   50,
                   -5.,
                   5.);

  meBTLTrackMatchedTPmtdDirectNoAssocEta_ =
      ibook.book1D("BTLTrackMatchedTPmtdDirectNoAssocEta",
                   "Eta of tracks matched to TP with sim hit in MTD (direct) - no track-MTD association;#eta_{RECO}",
                   30,
                   0.,
                   1.5);
  meBTLTrackMatchedTPmtdDirectNoAssocPt_ =
      ibook.book1D("BTLTrackMatchedTPmtdDirectNoAssocPt",
                   "Pt of tracks matched to TP with sim hit in MTD (direct) - no track-MTD association;track pt [GeV]",
                   50,
                   0.,
                   10.);

  meBTLTrackMatchedTPmtdOtherCorrectAssocEta_ = ibook.book1D(
      "BTLTrackMatchedTPmtdOtherCorrectAssocEta",
      "Eta of tracks matched to TP with sim hit in MTD (direct), correct track-MTD association;#eta_{RECO}",
      30,
      0.,
      1.5);
  meBTLTrackMatchedTPmtdOtherCorrectAssocPt_ = ibook.book1D(
      "BTLTrackMatchedTPmtdOtherCorrectAssocPt",
      "Pt of tracks matched to TP with sim hit in MTD (direct) - correct track-MTD association;track pt [GeV]",
      50,
      0.,
      10.);
  meBTLTrackMatchedTPmtdOtherCorrectAssocMVAQual_ = ibook.book1D(
      "BTLTrackMatchedTPmtdOtherCorrectAssocMVAQual",
      "MVA of tracks matched to TP with sim hit in MTD (direct) - correct track-MTD association; MVA score",
      100,
      -1.,
      1.);
  meBTLTrackMatchedTPmtdOtherCorrectAssocTimeRes_ =
      ibook.book1D("BTLTrackMatchedTPmtdOtherCorrectAssocTimeRes",
                   "Time resolution of tracks matched to TP with sim hit in MTD (direct) - correct track-MTD "
                   "association; t_{rec} - t_{sim} [ns] ",
                   120,
                   -0.15,
                   0.15);
  meBTLTrackMatchedTPmtdOtherCorrectAssocTimePull_ =
      ibook.book1D("BTLTrackMatchedTPmtdOtherCorrectAssocTimePull",
                   "Time pull of tracks matched to TP with sim hit in MTD (direct) - correct track-MTD association; "
                   "(t_{rec}-t_{sim})/#sigma_{t}",
                   50,
                   -5.,
                   5.);

  meBTLTrackMatchedTPmtdOtherWrongAssocEta_ =
      ibook.book1D("BTLTrackMatchedTPmtdOtherWrongAssocEta",
                   "Eta of tracks matched to TP with sim hit in MTD (direct) - wrong track-MTD association;#eta_{RECO}",
                   30,
                   0.,
                   1.5);
  meBTLTrackMatchedTPmtdOtherWrongAssocPt_ = ibook.book1D(
      "BTLTrackMatchedTPmtdOtherWrongAssocPt",
      "Pt of tracks matched to TP with sim hit in MTD (direct) - wrong track-MTD association;track pt [GeV]",
      50,
      0.,
      10.);
  meBTLTrackMatchedTPmtdOtherWrongAssocMVAQual_ =
      ibook.book1D("BTLTrackMatchedTPmtdOtherWrongAssocMVAQual",
                   "MVA of tracks matched to TP with sim hit in MTD (direct) - wrong track-MTD association; MVA score",
                   100,
                   -1.,
                   1.);
  meBTLTrackMatchedTPmtdOtherWrongAssocTimeRes_ =
      ibook.book1D("BTLTrackMatchedTPmtdOtherWrongAssocTimeRes",
                   "Time resolution of tracks matched to TP with sim hit in MTD (direct) - wrong track-MTD "
                   "association; t_{rec} - t_{sim} [ns] ",
                   120,
                   -0.15,
                   0.15);
  meBTLTrackMatchedTPmtdOtherWrongAssocTimePull_ =
      ibook.book1D("BTLTrackMatchedTPmtdOtherWrongAssocTimePull",
                   "Time pull of tracks matched to TP with sim hit in MTD (direct) - wrong track-MTD association; "
                   "(t_{rec}-t_{sim})/#sigma_{t}",
                   50,
                   -5.,
                   5.);

  meBTLTrackMatchedTPmtdOtherNoAssocEta_ =
      ibook.book1D("BTLTrackMatchedTPmtdOtherNoAssocEta",
                   "Eta of tracks matched to TP with sim hit in MTD (direct) - no track-MTD association;#eta_{RECO}",
                   30,
                   0.,
                   1.5);
  meBTLTrackMatchedTPmtdOtherNoAssocPt_ =
      ibook.book1D("BTLTrackMatchedTPmtdOtherNoAssocPt",
                   "Pt of tracks matched to TP with sim hit in MTD (direct) - no track-MTD association;track pt [GeV]",
                   50,
                   0.,
                   10.);

  meBTLTrackMatchedTPnomtdAssocEta_ =
      ibook.book1D("BTLTrackMatchedTPnomtdAssocEta",
                   "Eta of tracks matched to TP w/o sim hit in MTD, with associated reco cluster;#eta_{RECO}",
                   30,
                   0.,
                   1.5);
  meBTLTrackMatchedTPnomtdAssocPt_ =
      ibook.book1D("BTLTrackMatchedTPnomtdAssocPt",
                   "Pt of tracks matched to TP w/o sim hit in MTD, with associated reco cluster;track pt [GeV]",
                   50,
                   0.,
                   10.);
  meBTLTrackMatchedTPnomtdAssocMVAQual_ =
      ibook.book1D("BTLTrackMatchedTPnomtdAssocMVAQual",
                   "MVA of tracks matched to TP w/o sim hit in MTD, with associated reco cluster; MVA score",
                   100,
                   -1.,
                   1.);
  meBTLTrackMatchedTPnomtdAssocTimeRes_ = ibook.book1D("BTLTrackMatchedTPnomtdAssocTimeRes",
                                                       "Time resolution of tracks matched to TP w/o sim hit in MTD, "
                                                       "with associated reco cluster; t_{rec} - t_{sim} [ns] ",
                                                       120,
                                                       -0.15,
                                                       0.15);
  meBTLTrackMatchedTPnomtdAssocTimePull_ = ibook.book1D("BTLTrackMatchedTPnomtdAssocTimePull",
                                                        "Time pull of tracks matched to TP w/o sim hit in MTD, with "
                                                        "associated reco cluster; (t_{rec}-t_{sim})/#sigma_{t}",
                                                        50,
                                                        -5.,
                                                        5.);

  meETLTrackMatchedTPmtd1Eta_ = ibook.book1D("ETLTrackMatchedTPmtd1Eta",
                                             "Eta of tracks matched to TP with sim hit in MTD (>= 1 hit);#eta_{RECO}",
                                             30,
                                             1.5,
                                             3.0);
  meETLTrackMatchedTPmtd1Pt_ = ibook.book1D("ETLTrackMatchedTPmtd1Pt",
                                            "Pt of tracks matched to TP with sim hit in MTD (>= 1 hit); track pt [GeV]",
                                            50,
                                            0.,
                                            10.);

  meETLTrackMatchedTPmtd2Eta_ = ibook.book1D(
      "ETLTrackMatchedTPmtd2Eta", "Eta of tracks matched to TP with sim hit in MTD (2 hits);#eta_{RECO}", 30, 1.5, 3.0);
  meETLTrackMatchedTPmtd2Pt_ = ibook.book1D(
      "ETLTrackMatchedTPmtd2Pt", "Pt of tracks matched to TP with sim hit in MTD (2 hits); track pt [GeV]", 50, 0., 10.);

  meETLTrackMatchedTPnomtdEta_ = ibook.book1D(
      "ETLTrackMatchedTPnomtdEta", "Eta of tracks matched to TP w/o sim hit in MTD;#eta_{RECO}", 30, 1.5, 3.0);
  meETLTrackMatchedTPnomtdPt_ = ibook.book1D(
      "ETLTrackMatchedTPnomtdPt", "Pt of tracks matched to TP w/o sim hit in MTD; track pt [GeV]", 50, 0., 10.);

  meETLTrackMatchedTPmtd1CorrectAssocEta_ = ibook.book1D(
      "ETLTrackMatchedTPmtd1CorrectAssocEta",
      "Eta of tracks matched to TP with sim hit in MTD (>= 1 hit), correct track-MTD association;#eta_{RECO}",
      30,
      1.5,
      3.0);
  meETLTrackMatchedTPmtd1CorrectAssocPt_ = ibook.book1D(
      "ETLTrackMatchedTPmtd1CorrectAssocPt",
      "Pt of tracks matched to TP with sim hit in MTD (>= 1 hit) - correct track-MTD association;track pt [GeV]",
      50,
      0.,
      10.);
  meETLTrackMatchedTPmtd1CorrectAssocMVAQual_ = ibook.book1D(
      "ETLTrackMatchedTPmtd1CorrectAssocMVAQual",
      "MVA of tracks matched to TP with sim hit in MTD (>= 1 hit) - correct track-MTD association; MVA score",
      100,
      -1.,
      1.);
  meETLTrackMatchedTPmtd1CorrectAssocTimeRes_ =
      ibook.book1D("ETLTrackMatchedTPmtd1CorrectAssocTimeRes",
                   "Time resolution of tracks matched to TP with sim hit in MTD (>= 1 hit) - correct track-MTD "
                   "association; t_{rec} - t_{sim} [ns] ",
                   120,
                   -0.15,
                   0.15);
  meETLTrackMatchedTPmtd1CorrectAssocTimePull_ =
      ibook.book1D("ETLTrackMatchedTPmtd1CorrectAssocTimePull",
                   "Time pull of tracks matched to TP with sim hit in MTD (>= 1 hit) - correct track-MTD association; "
                   "(t_{rec}-t_{sim})/#sigma_{t}",
                   50,
                   -5.,
                   5.);

  meETLTrackMatchedTPmtd2CorrectAssocEta_ = ibook.book1D(
      "ETLTrackMatchedTPmtd2CorrectAssocEta",
      "Eta of tracks matched to TP with sim hit in MTD (2 hits), correct track-MTD association;#eta_{RECO}",
      30,
      1.5,
      3.0);
  meETLTrackMatchedTPmtd2CorrectAssocPt_ = ibook.book1D(
      "ETLTrackMatchedTPmtd2CorrectAssocPt",
      "Pt of tracks matched to TP with sim hit in MTD (2 hits) - correct track-MTD association;track pt [GeV]",
      50,
      0.,
      10.);
  meETLTrackMatchedTPmtd2CorrectAssocMVAQual_ = ibook.book1D(
      "ETLTrackMatchedTPmtd2CorrectAssocMVAQual",
      "MVA of tracks matched to TP with sim hit in MTD (2 hits) - correct track-MTD association; MVA score",
      100,
      -1.,
      1.);
  meETLTrackMatchedTPmtd2CorrectAssocTimeRes_ =
      ibook.book1D("ETLTrackMatchedTPmtd2CorrectAssocTimeRes",
                   "Time resolution of tracks matched to TP with sim hit in MTD (2 hits) - correct track-MTD "
                   "association; t_{rec} - t_{sim} [ns] ",
                   120,
                   -0.15,
                   0.15);
  meETLTrackMatchedTPmtd2CorrectAssocTimePull_ =
      ibook.book1D("ETLTrackMatchedTPmtd2CorrectAssocTimePull",
                   "Time pull of tracks matched to TP with sim hit in MTD (2 hits) - correct track-MTD association; "
                   "(t_{rec}-t_{sim})/#sigma_{t}",
                   50,
                   -5.,
                   5.);

  meETLTrackMatchedTPmtd1WrongAssocEta_ = ibook.book1D(
      "ETLTrackMatchedTPmtd1WrongAssocEta",
      "Eta of tracks matched to TP with sim hit in MTD (>= 1 hit), wrong track-MTD association;#eta_{RECO}",
      30,
      1.5,
      3.0);
  meETLTrackMatchedTPmtd1WrongAssocPt_ = ibook.book1D(
      "ETLTrackMatchedTPmtd1WrongAssocPt",
      "Pt of tracks matched to TP with sim hit in MTD (>= 1 hit) - wrong track-MTD association;track pt [GeV]",
      50,
      0.,
      10.);
  meETLTrackMatchedTPmtd1WrongAssocMVAQual_ = ibook.book1D(
      "ETLTrackMatchedTPmtd1WrongAssocMVAQual",
      "MVA of tracks matched to TP with sim hit in MTD (>= 1 hit) - wrong track-MTD association; MVA score",
      100,
      -1.,
      1.);
  meETLTrackMatchedTPmtd1WrongAssocTimeRes_ =
      ibook.book1D("ETLTrackMatchedTPmtd1WrongAssocTimeRes",
                   "Time resolution of tracks matched to TP with sim hit in MTD (>= 1 hit) - wrong track-MTD "
                   "association; t_{rec} - t_{sim} [ns] ",
                   120,
                   -0.15,
                   0.15);
  meETLTrackMatchedTPmtd1WrongAssocTimePull_ =
      ibook.book1D("ETLTrackMatchedTPmtd1WrongAssocTimePull",
                   "Time pull of tracks matched to TP with sim hit in MTD (>= 1 hit) - wrong track-MTD association; "
                   "(t_{rec}-t_{sim})/#sigma_{t}",
                   50,
                   -5.,
                   5.);

  meETLTrackMatchedTPmtd2WrongAssocEta_ =
      ibook.book1D("ETLTrackMatchedTPmtd2WrongAssocEta",
                   "Eta of tracks matched to TP with sim hit in MTD (2 hits), wrong track-MTD association;#eta_{RECO}",
                   30,
                   1.5,
                   3.0);
  meETLTrackMatchedTPmtd2WrongAssocPt_ = ibook.book1D(
      "ETLTrackMatchedTPmtd2WrongAssocPt",
      "Pt of tracks matched to TP with sim hit in MTD (2 hits) - wrong track-MTD association;track pt [GeV]",
      50,
      0.,
      10.);
  meETLTrackMatchedTPmtd2WrongAssocMVAQual_ =
      ibook.book1D("ETLTrackMatchedTPmtd2WrongAssocMVAQual",
                   "MVA of tracks matched to TP with sim hit in MTD (2 hits) - wrong track-MTD association; MVA score",
                   100,
                   -1.,
                   1.);
  meETLTrackMatchedTPmtd2WrongAssocTimeRes_ =
      ibook.book1D("ETLTrackMatchedTPmtd2WrongAssocTimeRes",
                   "Time resolution of tracks matched to TP with sim hit in MTD (2 hits) - wrong track-MTD "
                   "association; t_{rec} - t_{sim} [ns] ",
                   120,
                   -0.15,
                   0.15);
  meETLTrackMatchedTPmtd2WrongAssocTimePull_ =
      ibook.book1D("ETLTrackMatchedTPmtd2WrongAssocTimePull",
                   "Time pull of tracks matched to TP with sim hit in MTD (2 hits) - wrong track-MTD association; "
                   "(t_{rec}-t_{sim})/#sigma_{t}",
                   50,
                   -5.,
                   5.);

  meETLTrackMatchedTPmtd1NoAssocEta_ = ibook.book1D(
      "ETLTrackMatchedTPmtd1NoAssocEta",
      "Eta of tracks matched to TP with sim hit in MTD (>= 1 hit), missing track-MTD association;#eta_{RECO}",
      30,
      1.5,
      3.0);
  meETLTrackMatchedTPmtd1NoAssocPt_ = ibook.book1D(
      "ETLTrackMatchedTPmtd1NoAssocPt",
      "Pt of tracks matched to TP with sim hit in MTD (>= 1 hit) - missing track-MTD association;track pt [GeV]",
      50,
      0.,
      10.);

  meETLTrackMatchedTPmtd2NoAssocEta_ = ibook.book1D(
      "ETLTrackMatchedTPmtd2NoAssocEta",
      "Eta of tracks matched to TP with sim hit in MTD (2 hits), missing track-MTD association;#eta_{RECO}",
      30,
      1.5,
      3.0);
  meETLTrackMatchedTPmtd2NoAssocPt_ = ibook.book1D(
      "ETLTrackMatchedTPmtd2NoAssocPt",
      "Pt of tracks matched to TP with sim hit in MTD (2 hits) - missing track-MTD association;track pt [GeV]",
      50,
      0.,
      10.);
  meETLTrackMatchedTPnomtdAssocEta_ =
      ibook.book1D("ETLTrackMatchedTPnomtdAssocEta",
                   "Eta of tracks matched to TP w/o sim hit in MTD, with associated reco cluster;#eta_{RECO}",
                   30,
                   1.5,
                   3.0);
  meETLTrackMatchedTPnomtdAssocPt_ =
      ibook.book1D("ETLTrackMatchedTPnomtdAssocPt",
                   "Pt of tracks matched to TP w/o sim hit in MTD, with associated reco cluster;track pt [GeV]",
                   50,
                   0.,
                   10.);
  meETLTrackMatchedTPnomtdAssocMVAQual_ =
      ibook.book1D("ETLTrackMatchedTPnomtdAssocMVAQual",
                   "MVA of tracks matched to TP w/o sim hit in MTD, with associated reco cluster; MVA score",
                   100,
                   -1.,
                   1.);
  meETLTrackMatchedTPnomtdAssocTimeRes_ = ibook.book1D("ETLTrackMatchedTPnomtdAssocTimeRes",
                                                       "Time resolution of tracks matched to TP w/o sim hit in MTD, "
                                                       "with associated reco cluster; t_{rec} - t_{sim} [ns] ",
                                                       120,
                                                       -0.15,
                                                       0.15);
  meETLTrackMatchedTPnomtdAssocTimePull_ = ibook.book1D("ETLTrackMatchedTPnomtdAssocTimePull",
                                                        "Time pull of tracks matched to TP w/o sim hit in MTD, with "
                                                        "associated reco cluster; (t_{rec}-t_{sim})/#sigma_{t}",
                                                        50,
                                                        -5.,
                                                        5.);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

void MtdTracksValidation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;

  desc.add<std::string>("folder", "MTD/Tracks");
  desc.add<bool>("optionalPlots", false);
  desc.add<edm::InputTag>("inputTagG", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("inputTagT", edm::InputTag("trackExtenderWithMTD"));
  desc.add<edm::InputTag>("inputTagV", edm::InputTag("offlinePrimaryVertices4D"));
  desc.add<edm::InputTag>("inputTagH", edm::InputTag("generatorSmeared"));
  desc.add<edm::InputTag>("SimTag", edm::InputTag("mix", "MergedTrackTruth"));
  desc.add<edm::InputTag>("TPtoRecoTrackAssoc", edm::InputTag("trackingParticleRecoTrackAsssociation"));
  desc.add<edm::InputTag>("tp2SimAssociationMapTag", edm::InputTag("mtdSimLayerClusterToTPAssociation"));
  desc.add<edm::InputTag>("r2sAssociationMapTag", edm::InputTag("mtdRecoClusterToSimLayerClusterAssociation"));
  desc.add<edm::InputTag>("btlRecHits", edm::InputTag("mtdRecHits", "FTLBarrel"));
  desc.add<edm::InputTag>("etlRecHits", edm::InputTag("mtdRecHits", "FTLEndcap"));
  desc.add<edm::InputTag>("recCluTagBTL", edm::InputTag("mtdClusters", "FTLBarrel"));
  desc.add<edm::InputTag>("recCluTagETL", edm::InputTag("mtdClusters", "FTLEndcap"));
  desc.add<edm::InputTag>("tmtd", edm::InputTag("trackExtenderWithMTD:generalTracktmtd"));
  desc.add<edm::InputTag>("sigmatmtd", edm::InputTag("trackExtenderWithMTD:generalTracksigmatmtd"));
  desc.add<edm::InputTag>("t0Src", edm::InputTag("trackExtenderWithMTD:generalTrackt0"));
  desc.add<edm::InputTag>("sigmat0Src", edm::InputTag("trackExtenderWithMTD:generalTracksigmat0"));
  desc.add<edm::InputTag>("trackAssocSrc", edm::InputTag("trackExtenderWithMTD:generalTrackassoc"))
      ->setComment("Association between General and MTD Extended tracks");
  desc.add<edm::InputTag>("pathLengthSrc", edm::InputTag("trackExtenderWithMTD:generalTrackPathLength"));
  desc.add<edm::InputTag>("t0SafePID", edm::InputTag("tofPID:t0safe"));
  desc.add<edm::InputTag>("sigmat0SafePID", edm::InputTag("tofPID:sigmat0safe"));
  desc.add<edm::InputTag>("sigmat0PID", edm::InputTag("tofPID:sigmat0"));
  desc.add<edm::InputTag>("t0PID", edm::InputTag("tofPID:t0"));
  desc.add<edm::InputTag>("TofPi", edm::InputTag("trackExtenderWithMTD:generalTrackTofPi"));
  desc.add<edm::InputTag>("TofK", edm::InputTag("trackExtenderWithMTD:generalTrackTofK"));
  desc.add<edm::InputTag>("TofP", edm::InputTag("trackExtenderWithMTD:generalTrackTofP"));
  desc.add<edm::InputTag>("sigmaTofPi", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofPi"));
  desc.add<edm::InputTag>("sigmaTofK", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofK"));
  desc.add<edm::InputTag>("sigmaTofP", edm::InputTag("trackExtenderWithMTD:generalTrackSigmaTofP"));
  desc.add<edm::InputTag>("trackMVAQual", edm::InputTag("mtdTrackQualityMVA:mtdQualMVA"));
  desc.add<edm::InputTag>("outermostHitPositionSrc",
                          edm::InputTag("trackExtenderWithMTD:generalTrackOutermostHitPosition"));
  desc.add<double>("trackMaximumPt", 12.);  // [GeV]
  desc.add<double>("trackMaximumBtlEta", 1.5);
  desc.add<double>("trackMinimumEtlEta", 1.6);
  desc.add<double>("trackMaximumEtlEta", 3.);
  desc.add<edm::InputTag>("probPi", edm::InputTag("tofPID:probPi"));
  desc.add<edm::InputTag>("probK", edm::InputTag("tofPID:probK"));
  desc.add<edm::InputTag>("probP", edm::InputTag("tofPID:probP"));

  descriptions.add("mtdTracksValid", desc);

}

const bool MtdTracksValidation::trkTPSelLV(const TrackingParticle& tp) {
  bool match = (tp.status() != 1) ? false : true;
  return match;
}

const bool MtdTracksValidation::trkTPSelAll(const TrackingParticle& tp) {
  bool match = false;

  auto x_pv = tp.parentVertex()->position().x();
  auto y_pv = tp.parentVertex()->position().y();
  auto z_pv = tp.parentVertex()->position().z();

  auto r_pv = std::sqrt(x_pv * x_pv + y_pv * y_pv);

  match = tp.charge() != 0 && std::abs(tp.eta()) < etacutGEN_ && tp.pt() > pTcutBTL_ && r_pv < rBTL_ &&
          std::abs(z_pv) < zETL_;
  return match;
}

const bool MtdTracksValidation::trkRecSel(const reco::TrackBase& trk) {
  bool match = false;
  match = std::abs(trk.eta()) <= etacutREC_ && trk.pt() > pTcutBTL_;
  return match;
}

const bool MtdTracksValidation::trkRecSelLowPt(const reco::TrackBase& trk) {
  bool match = false;
  match = std::abs(trk.eta()) <= etacutREC_ && trk.pt() > pTcutETL_ && trk.pt() < pTcutBTL_;
  return match;
}

const edm::Ref<std::vector<TrackingParticle>>* MtdTracksValidation::getMatchedTP(const reco::TrackBaseRef& recoTrack) {
  auto found = r2s_->find(recoTrack);

  // reco track not matched to any TP
  if (found == r2s_->end())
    return nullptr;

  //matched TP equal to any TP associated to in time events
  for (const auto& tp : found->val) {
    if (tp.first->eventId().bunchCrossing() == 0)
      return &tp.first;
  }

  // reco track not matched to any TP from vertex
  return nullptr;
}

void MtdTracksValidation::fillTrackClusterMatchingHistograms(MonitorElement* me1,
                                                             MonitorElement* me2,
                                                             MonitorElement* me3,
                                                             MonitorElement* me4,
                                                             MonitorElement* me5,
                                                             float var1,
                                                             float var2,
                                                             float var3,
                                                             float var4,
                                                             float var5,
                                                             bool flag) {
  me1->Fill(var1);
  me2->Fill(var2);
  if (flag) {
    me3->Fill(var3);
    me4->Fill(var4);
    me5->Fill(var5);
  }
}


void MtdTracksValidation::isParticle(const reco::TrackRef& recoTrack,
  const edm::ValueMap<float>& sigmat0,
  const edm::ValueMap<float>& sigmat0Safe,
  const edm::ValueMap<float>& probPi,
  const edm::ValueMap<float>& probK,
  const edm::ValueMap<float>& probP,
  unsigned int& no_PIDtype,
  bool& no_PID,
  bool& is_Pi,
  bool& is_K,
  bool& is_P) {
no_PIDtype = 0;
no_PID = false;
is_Pi = false;
is_K = false;
is_P = false;
if (probPi[recoTrack] == -1) {
no_PIDtype = 1;
} else if (edm::isNotFinite(probPi[recoTrack])) {
no_PIDtype = 2;
} else if (probPi[recoTrack] == 1 && probK[recoTrack] == 0 && probP[recoTrack] == 0 &&
sigmat0[recoTrack] < sigmat0Safe[recoTrack]) {
no_PIDtype = 3;
}
no_PID = no_PIDtype > 0;
is_Pi = !no_PID && 1. - probPi[recoTrack] < 0.75;
is_K = !no_PID && !is_Pi && probK[recoTrack] > probP[recoTrack];
is_P = !no_PID && !is_Pi && !is_K;
}

DEFINE_FWK_MODULE(MtdTracksValidation);
