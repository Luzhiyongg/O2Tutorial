opt="-b --configuration json://configuration_run3.json"
    o2-analysis-track-propagation ${opt} |
    o2-analysis-timestamp ${opt} | 
    o2-analysis-multiplicity-table ${opt} | 
    o2-analysis-trackselection ${opt} | 
    o2-analysis-centrality-table ${opt} |
    o2-analysis-event-selection ${opt} | 
    o2-analysis-cf-flow-pbpb-task ${opt}  --aod-file @testaod.txt --aod-memory-rate-limit 429496320
