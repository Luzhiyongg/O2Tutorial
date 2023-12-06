Dir:
alice/O2Physics/PWGCF/Flow/Tasks/
In CMakeList:
o2physics_add_dpl_workflow(flow-pbpb-task
                    SOURCES FlowPbPbTask.cxx
                    PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore O2Physics::GFWCore
                    COMPONENT_NAME Analysis)

executable name:
o2-analysis-cf-flow-pbpb-task