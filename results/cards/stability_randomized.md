# Environment QC: stability

Rig: randomized. Seeds: 0 to 199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.905 | 0.070 | 0.884 | 0.934 | 98% | reference | | | |
| two_point | 0.620 | 0.227 | 0.309 | 0.937 | 98% | 0.286 | -2.13 | 85.0% | weak
| early_only | 0.749 | 0.168 | 0.552 | 0.923 | 98% | 0.157 | -3.57 | 90.0% | visible
| no_internal_standard | 0.498 | 0.194 | 0.352 | 0.874 | 23% | 0.407 | -0.94 | 98.5% | visible
| no_mix | 0.458 | 0.231 | 0.204 | 0.818 | 44% | 0.447 | -1.02 | 99.0% | visible
| open_edges | 0.747 | 0.149 | 0.559 | 0.909 | 98% | 0.159 | -3.15 | 97.0% | visible
| single_replicate | 0.452 | 0.012 | 0.435 | 0.464 | 0% | 0.453 | 0.45 | 98.0% | visible
| no_lloq | 0.689 | 0.189 | 0.434 | 0.919 | 98% | 0.216 | -2.61 | 87.5% | weak
| random | 0.320 | 0.162 | 0.129 | 0.447 | 7% | 0.585 | -0.19 | 98.5% | visible
| null | 0.085 | 0.084 | 0.000 | 0.214 | 0% | 0.820 | 0.44 | 100.0% | visible

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake).
