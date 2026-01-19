import FWCore.ParameterSet.Config as cms

from RecoMTD.TrackExtender.PropagatorWithMaterialForMTD_cfi import *
from RecoMTD.TrackExtender.trackExtenderWithMTDBase_cfi import *

trackExtenderWithMTD = trackExtenderWithMTDBase.clone()

trackExtenderWithMTDrelaxedTimeChi2 = trackExtenderWithMTDBase.clone(etlTimeChi2Cut=-1., btlTimeChi2Cut=-1.)
