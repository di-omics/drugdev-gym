# Environment QC: stability

Rig: nominal. Seeds: 0 to 199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.916 | 0.019 | 0.889 | 0.934 | 100% | reference | | | |
| two_point | 0.626 | 0.223 | 0.310 | 0.937 | 100% | 0.290 | -1.51 | 86% | weak
| early_only | 0.759 | 0.158 | 0.573 | 0.921 | 100% | 0.157 | -2.40 | 90% | visible
| no_internal_standard | 0.496 | 0.192 | 0.352 | 0.872 | 22% | 0.420 | -0.51 | 99% | visible
| no_mix | 0.450 | 0.207 | 0.224 | 0.780 | 42% | 0.466 | -0.45 | 100% | visible
| open_edges | 0.726 | 0.142 | 0.539 | 0.904 | 100% | 0.189 | -1.55 | 98% | visible
| single_replicate | 0.452 | 0.012 | 0.434 | 0.464 | 0% | 0.464 | 0.80 | 100% | clean
| no_lloq | 0.697 | 0.184 | 0.497 | 0.919 | 100% | 0.218 | -1.79 | 88% | weak
| random | 0.320 | 0.157 | 0.133 | 0.446 | 7% | 0.595 | 0.11 | 99% | visible
| null | 0.085 | 0.084 | 0.000 | 0.214 | 0% | 0.830 | 0.63 | 100% | clean

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake).
