LOGFILE="log.txt"
OPTION="-b --configuration json://Config.json"
## Tree creator
# o2-analysis-hf-tree-creator-ds-to-k-k-pi -b $OPTION | o2-analysis-hf-candidate-selector-ds-to-k-k-pi -b $OPTION | o2-analysis-hf-candidate-creator-3prong -b $OPTION | o2-analysis-bc-converter -b $OPTION | o2-analysis-pid-tpc-full -b $OPTION | o2-analysis-pid-tof-full -b $OPTION | o2-analysis-pid-tof-base -b $OPTION | o2-analysis-multiplicity-table -b $OPTION | o2-analysis-tracks-extra-converter -b $OPTION | o2-analysis-hf-track-index-skim-creator -b $OPTION | o2-analysis-pid-tpc-base -b $OPTION | o2-analysis-timestamp -b $OPTION | o2-analysis-centrality-table -b $OPTION | o2-analysis-ft0-corrected-table -b $OPTION | o2-analysis-event-selection -b $OPTION | o2-analysis-track-propagation -b $OPTION | o2-analysis-collision-converter -b $OPTION | o2-analysis-trackselection -b $OPTION | o2-analysis-zdc-converter -b $OPTION | o2-analysis-track-to-collision-associator -b $OPTION --aod-file @input_data.txt --aod-writer-json OutputDirector.json --shm-segment-size 3000000000 > "$LOGFILE" 2>&1

## Task without centrality
#o2-analysis-hf-task-ds -b $OPTION | o2-analysis-hf-candidate-selector-ds-to-k-k-pi -b $OPTION | o2-analysis-hf-candidate-creator-3prong -b $OPTION | o2-analysis-bc-converter -b $OPTION | o2-analysis-pid-tpc-full -b $OPTION | o2-analysis-pid-tof-full -b $OPTION | o2-analysis-pid-tof-base -b $OPTION | o2-analysis-multiplicity-table -b $OPTION | o2-analysis-tracks-extra-converter -b $OPTION | o2-analysis-hf-track-index-skim-creator -b $OPTION | o2-analysis-pid-tpc-base -b $OPTION | o2-analysis-timestamp -b $OPTION | o2-analysis-ft0-corrected-table -b $OPTION | o2-analysis-event-selection -b $OPTION | o2-analysis-track-propagation -b $OPTION | o2-analysis-collision-converter -b $OPTION | o2-analysis-trackselection -b $OPTION | o2-analysis-zdc-converter -b $OPTION | o2-analysis-track-to-collision-associator -b $OPTION --aod-file @input_data.txt --aod-writer-json OutputDirector.json --shm-segment-size 3000000000 > "$LOGFILE" 2>&1

## Task with centrality
#o2-analysis-hf-task-ds -b $OPTION | o2-analysis-hf-candidate-selector-ds-to-k-k-pi -b $OPTION | o2-analysis-hf-candidate-creator-3prong -b $OPTION | o2-analysis-bc-converter -b $OPTION | o2-analysis-pid-tpc-full -b $OPTION | o2-analysis-pid-tof-full -b $OPTION | o2-analysis-pid-tof-base -b $OPTION | o2-analysis-multiplicity-table -b $OPTION | o2-analysis-tracks-extra-converter -b $OPTION | o2-analysis-hf-track-index-skim-creator -b $OPTION | o2-analysis-pid-tpc-base -b $OPTION | o2-analysis-timestamp -b $OPTION | o2-analysis-centrality-table -b $OPTION | o2-analysis-ft0-corrected-table -b $OPTION | o2-analysis-event-selection -b $OPTION | o2-analysis-track-propagation -b $OPTION | o2-analysis-trackselection -b $OPTION | o2-analysis-track-to-collision-associator -b $OPTION --aod-file @input_data_LHC23j4.txt --aod-writer-json OutputDirector.json --shm-segment-size 3000000000 > "$LOGFILE" 2>&1

## Task PbPb MC
#o2-analysis-hf-task-ds -b $OPTION | o2-analysis-hf-candidate-selector-ds-to-k-k-pi -b $OPTION | o2-analysis-hf-candidate-creator-3prong -b $OPTION | o2-analysis-pid-tpc-full -b $OPTION | o2-analysis-pid-tof-full -b $OPTION | o2-analysis-pid-tof-base -b $OPTION | o2-analysis-multiplicity-table -b $OPTION | o2-analysis-hf-track-index-skim-creator -b $OPTION | o2-analysis-pid-tpc-base -b $OPTION | o2-analysis-timestamp -b $OPTION | o2-analysis-centrality-table -b $OPTION | o2-analysis-ft0-corrected-table -b $OPTION | o2-analysis-event-selection -b $OPTION | o2-analysis-track-propagation -b $OPTION | o2-analysis-trackselection -b $OPTION | o2-analysis-track-to-collision-associator -b $OPTION --aod-file @input_data_LHC23k6c.txt --aod-writer-json OutputDirector.json --shm-segment-size 3000000000 > "$LOGFILE" 2>&1

## Task PbPb
o2-analysis-hf-task-ds -b $OPTION | o2-analysis-hf-candidate-selector-ds-to-k-k-pi -b $OPTION | o2-analysis-hf-candidate-creator-3prong -b $OPTION | o2-analysis-pid-tpc-full -b $OPTION | o2-analysis-pid-tof-full -b $OPTION | o2-analysis-pid-tof-base -b $OPTION | o2-analysis-multiplicity-table -b $OPTION | o2-analysis-hf-track-index-skim-creator -b $OPTION | o2-analysis-pid-tpc-base -b $OPTION | o2-analysis-timestamp -b $OPTION | o2-analysis-centrality-table -b $OPTION | o2-analysis-ft0-corrected-table -b $OPTION | o2-analysis-event-selection -b $OPTION | o2-analysis-track-propagation -b $OPTION | o2-analysis-trackselection -b $OPTION | o2-analysis-track-to-collision-associator -b $OPTION --aod-file @input_data_PbPb.txt --aod-writer-json OutputDirector.json --shm-segment-size 3000000000 > "$LOGFILE" 2>&1

# report status
rc=$?
if [ $rc -eq 0 ]; then
  echo "No problems!"
else
  echo "Error: Exit code $rc"
  echo "Check the log file $LOGFILE"
  exit $rc
fi