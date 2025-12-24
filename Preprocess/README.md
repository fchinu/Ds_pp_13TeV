# Pre-processing Script
This utility prepares large `AnalysisResults.root` files by performing pT-bin splitting, background skimming, and MC reweighting.

## Usage
Run the script by providing the path to your YAML configuration:

```bash
python3 pre_process.py config.yml --workers 4
```

## Note:
For MC, we only use it to apply weights to be used for efficiency evaluation.