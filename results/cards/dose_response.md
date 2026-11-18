# Environment QC: dose_response

Rig: nominal. Seeds: 0 to 199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.681 | 0.114 | 0.530 | 0.807 | 100% | reference | | | |
| no_mix | 0.074 | 0.121 | -0.058 | 0.235 | 100% | 0.607 | -0.16 | 99.5% | visible
| reuse_tips | 0.650 | 0.124 | 0.490 | 0.803 | 100% | 0.031 | -21.82 | 63.0% | weak
| open_edges | 0.640 | 0.129 | 0.466 | 0.781 | 100% | 0.040 | -17.00 | 61.5% | weak
| raw_ratio | 0.061 | 0.117 | -0.052 | 0.228 | 100% | 0.620 | -0.12 | 100.0% | visible
| narrow_range | 0.032 | 0.112 | -0.051 | 0.197 | 100% | 0.649 | -0.04 | 100.0% | visible
| thin_ladder | -0.008 | 0.089 | -0.058 | 0.141 | 100% | 0.689 | 0.12 | 100.0% | visible
| few_controls | 0.332 | 0.053 | 0.250 | 0.389 | 0% | 0.349 | -0.43 | 99.5% | visible
| random | 0.100 | 0.181 | -0.041 | 0.371 | 32% | 0.581 | -0.52 | 98.0% | visible
| null | 0.064 | 0.069 | 0.000 | 0.157 | 0% | 0.617 | 0.11 | 100.0% | visible

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake).
