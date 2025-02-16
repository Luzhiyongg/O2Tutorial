opt="-b --configuration json://configuration_flowPtEfficiency.json"
    aodfile="AO2D_LHC24k2.root"   
    o2-analysis-track-propagation ${opt} |
    o2-analysis-tracks-extra-v002-converter ${opt} | 
    o2-analysis-timestamp ${opt} | 
    o2-analysis-trackselection ${opt} | 
    o2-analysis-event-selection ${opt} | 
    o2-analysis-cf-flow-pt-efficiency ${opt}  --aod-file ${aodfile} --aod-memory-rate-limit 429496320
