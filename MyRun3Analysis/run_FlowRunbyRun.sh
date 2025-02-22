opt="-b --configuration json://configuration_flow_runbyrun.json"
    aodfile="AO2D_LHC23zzh_apass4_544122.root"   
    o2-analysis-track-propagation ${opt} |
    o2-analysis-tracks-extra-v002-converter ${opt} | 
    o2-analysis-timestamp ${opt} | 
    o2-analysis-multiplicity-table ${opt} | 
    o2-analysis-trackselection ${opt} | 
    o2-analysis-centrality-table ${opt} |
    o2-analysis-event-selection ${opt} | 
    o2-analysis-cf-flow-runby-run ${opt}  --aod-file ${aodfile} --aod-memory-rate-limit 429496320
