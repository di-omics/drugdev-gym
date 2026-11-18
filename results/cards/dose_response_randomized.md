# Environment QC: dose_response

Rig: randomized. Seeds: 0 to 199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.659 | 0.127 | 0.498 | 0.799 | 100% | reference | | | |
| no_mix | 0.092 | 0.149 | -0.058 | 0.309 | 100% | 0.566 | -0.46 | 99.5% | visible
| reuse_tips | 0.629 | 0.136 | 0.454 | 0.795 | 100% | 0.029 | -25.83 | 59.5% | blind
| open_edges | 0.606 | 0.148 | 0.392 | 0.768 | 100% | 0.052 | -14.77 | 61.0% | weak
| raw_ratio | 0.070 | 0.125 | -0.052 | 0.237 | 100% | 0.588 | -0.28 | 99.5% | visible
| narrow_range | 0.089 | 0.145 | -0.051 | 0.285 | 100% | 0.570 | -0.43 | 100.0% | visible
| thin_ladder | 0.062 | 0.172 | -0.058 | 0.270 | 100% | 0.597 | -0.50 | 99.5% | visible
| few_controls | 0.323 | 0.059 | 0.237 | 0.389 | 0% | 0.335 | -0.67 | 99.0% | visible
| random | 0.098 | 0.176 | -0.044 | 0.346 | 32% | 0.560 | -0.63 | 98.0% | visible
| null | 0.064 | 0.069 | 0.000 | 0.157 | 0% | 0.595 | 0.01 | 99.5% | visible

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake).
