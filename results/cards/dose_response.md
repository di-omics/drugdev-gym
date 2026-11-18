# Environment QC: dose_response

Rig: nominal. Seeds: 0 to 199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.681 | 0.114 | 0.530 | 0.807 | 100% | reference | | | |
| no_mix | 0.074 | 0.121 | -0.058 | 0.235 | 100% | 0.607 | -0.16 | 99.5% | visible
| reuse_tips | 0.648 | 0.126 | 0.488 | 0.803 | 100% | 0.033 | -20.75 | 63.0% | weak
| edge_wells | 0.604 | 0.142 | 0.406 | 0.761 | 100% | 0.077 | -8.98 | 69.5% | weak
| no_lid | 0.681 | 0.114 | 0.530 | 0.810 | 100% | -0.000 | -1740.50 | 47.0% | blind
| raw_ratio | 0.061 | 0.117 | -0.052 | 0.228 | 100% | 0.620 | -0.12 | 100.0% | visible
| narrow_range | 0.311 | 0.195 | 0.044 | 0.552 | 100% | 0.370 | -1.50 | 97.0% | visible
| thin_ladder | -0.058 | 0.000 | -0.058 | -0.058 | 100% | 0.739 | 0.54 | 100.0% | clean
| few_controls | 0.332 | 0.053 | 0.250 | 0.389 | 0% | 0.349 | -0.43 | 99.5% | visible
| random | 0.078 | 0.177 | -0.053 | 0.339 | 32% | 0.603 | -0.45 | 98.0% | visible
| null | 0.064 | 0.069 | 0.000 | 0.157 | 0% | 0.617 | 0.11 | 100.0% | visible
| gate_bypass | 0.656 | 0.117 | 0.503 | 0.785 | 100% | 0.025 | -26.46 | 60.0% | weak

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake); inverted = the mistake scores higher and wins over 60% of seeds (the reward prefers it).
