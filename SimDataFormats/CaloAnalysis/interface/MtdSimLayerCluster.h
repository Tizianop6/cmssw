// Author: Aurora Perego, Fabio Cossutti - aurora.perego@cern.ch, fabio.cossutti@ts.infn.it
// Date: 05/2023

#ifndef SimDataFormats_CaloAnalysis_MtdSimLayerCluster_h
#define SimDataFormats_CaloAnalysis_MtdSimLayerCluster_h

#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "SimDataFormats/CaloAnalysis/interface/MtdSimCluster.h"
#include <vector>

class MtdSimLayerCluster : public MtdSimCluster {
  friend std::ostream &operator<<(std::ostream &s, MtdSimLayerCluster const &tp);

public:
  MtdSimLayerCluster();
  MtdSimLayerCluster(const SimTrack &simtrk);
  MtdSimLayerCluster(EncodedEventId eventID, uint32_t particleID);  // for PU

  // destructor
  ~MtdSimLayerCluster();

  /** @brief computes the time of the cluster */
  float computeClusterTime() {
    simLC_time_ = 0.;
    float tot_en = 0.;
    for (uint32_t i = 0; i < times_.size(); i++) {
      if (pathlengths_[i]>100){
        simLC_time_ += times_[i] * energies_[i];
        tot_en += energies_[i];
      }
    }
    if (tot_en != 0.)
      simLC_time_ = simLC_time_ / tot_en;
    return simLC_time_;
  }

 

  /** @brief computes the time of the cluster */
  float computeClusterPL() {
    simLC_PL_ = 0.;
    float tot_en = 0.;
    
    for (uint32_t i = 0; i < pathlengths_.size(); i++) {
      if(pathlengths_[i]>100){
        simLC_PL_ += pathlengths_[i] * energies_[i];
        tot_en += energies_[i];      
      }
      else{
        std::cout << "Discarded: ";
      }
      std::cout << "MTDSimLayerCluster: pathlengths_[" << i << "] = " << pathlengths_[i] << " energies_[" << i << "] = " << energies_[i] << " time = "<< times_[i]<< " p= "<< ps_[i] <<std::endl;
      }
    if (tot_en != 0.)
      simLC_PL_ = simLC_PL_ / tot_en;
    std::cout << "--> MTDSimLayerCluster weighted average PL: " << simLC_PL_ << std::endl;
    
    return simLC_PL_;
  }
  
   /** @brief computes the time of the cluster */
   float computeClusterP() {
    simLC_P_ = 0.;
    float tot_en = 0.;
    
    for (uint32_t i = 0; i < pathlengths_.size(); i++) {
      if(pathlengths_[i]>100){
        simLC_P_ += ps_[i] * energies_[i];
        tot_en += energies_[i];      
      }
      }
    if (tot_en != 0.)
    simLC_P_ = simLC_P_ / tot_en;
    std::cout << "--> MTDSimLayerCluster weighted average P: " << simLC_P_ << std::endl;
    
    return simLC_P_;
  }



  
  /** @brief computes the time of the cluster */
  float computeClusterearliestPL() {
    simLC_earliestPL_ = 0;
    float simLC_time_temp = std::numeric_limits<float>::max();
    for (uint32_t i = 0; i < times_.size(); i++) {
      if (pathlengths_[i]>100){
        if (simLC_time_temp > times_[i]){
          simLC_time_temp = times_[i];
          simLC_earliestPL_ = pathlengths_[i];
        }
    }
    }
    std::cout << "--> MTDSimLayerCluster Earliest PL: " << simLC_earliestPL_ << std::endl;
    return simLC_earliestPL_;
  }


  /** @brief computes the time of the cluster */
  float computeClusterearliestP() {
    simLC_earliestP_ = 0;
    float simLC_time_temp = std::numeric_limits<float>::max();
    for (uint32_t i = 0; i < times_.size(); i++) {
      if (pathlengths_[i]>100){
        if (simLC_time_temp > times_[i]){
          simLC_time_temp = times_[i];
          simLC_earliestP_ = ps_[i];
        }
    }
    }
    std::cout << "--> MTDSimLayerCluster Earliest P: " << simLC_earliestP_ << std::endl;
    return simLC_earliestP_;
  }
  

  
  float computeClusterEarliestTime() {
    
    simLC_earliesttime_ = std::numeric_limits<float>::max();
    for (uint32_t i = 0; i < times_.size(); i++) {
      if (pathlengths_[i]>100){
        if (simLC_earliesttime_ > times_[i]) {
          simLC_earliesttime_ = times_[i];
        }
      }
    }

    return simLC_earliesttime_;
  }

  /** @brief computes the energy of the cluster */
  void addCluEnergy(float energy) { simLC_energy_ = energy; }

  /** @brief computes the position of the cluster */
  void addCluLocalPos(LocalPoint pos) { simLC_pos_ = pos; }

  /** @brief add the index of the simcluster */
  void addCluIndex(const uint32_t index) { seedId_ = index; }

  /** @brief returns the time of the cluster */
  float simLCTime() const { return simLC_time_; }
  
  /** @brief returns the time of the cluster */
  float simLCearliestTime() const { return simLC_earliesttime_; }
  

  /** @brief returns the time of the cluster */
  float simLCPL() const { return simLC_PL_; }

  /** @brief returns the time of the cluster */
  float simLCearliestPL() const { return simLC_earliestPL_; }

  /** @brief returns the time of the cluster */
  float simLCP() const { return simLC_P_; }

  /** @brief returns the time of the cluster */
  float simLCearliestP() const { return simLC_earliestP_; }
    

  /** @brief returns the local position of the cluster */
  LocalPoint simLCPos() const { return simLC_pos_; }

  /** @brief returns the accumulated sim energy in the cluster */
  float simLCEnergy() const { return simLC_energy_; }

  uint32_t seedId() const { return seedId_; }

private:
  // id of the simCluster it comes from
  uint32_t seedId_;

  float simLC_time_{0.f};
  float simLC_earliesttime_{0.f};
  float simLC_PL_{0.f};
  float simLC_earliestPL_{0.f};
  float simLC_P_{0.f};
  float simLC_earliestP_{0.f};
  
  float simLC_energy_{0.f};
  LocalPoint simLC_pos_;
};

#endif
