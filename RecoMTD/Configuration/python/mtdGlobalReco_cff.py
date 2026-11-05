import FWCore.ParameterSet.Config as cms

from RecoMTD.TrackExtender.trackExtenderWithMTD_cfi import *
from RecoMTD.TimingIDTools.mtdTrackQualityMVA_cfi import *

fastTimingGlobalRecoTask = cms.Task(trackExtenderWithMTD,trackExtenderWithMTDrelaxedTimeChi2,mtdTrackQualityMVA)
#fastTimingGlobalRecoTask2 = cms.Task(,mtdTrackQualityMVA)

fastTimingGlobalReco = cms.Sequence(fastTimingGlobalRecoTask)
