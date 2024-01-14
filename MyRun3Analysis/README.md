Dir:
alice/O2Physics/PWGCF/Flow/Tasks/
In CMakeList:
o2physics_add_dpl_workflow(flow-pbpb-task
                    SOURCES FlowPbPbTask.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::GFWCore
                    COMPONENT_NAME Analysis)

executable name:
o2-analysis-cf-flow-pbpb-task


data file:
AO2D_LHC23zzh_apass1.root: /alice/data/2023/LHC23zzh/544122/apass1/0420/o2_ctf_run00544122_orbit0140388224_tf0000136787_epn102/001/AO2D.root
