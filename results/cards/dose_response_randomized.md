# Environment QC: dose_response

Rig: randomized. Seeds: 0 to 199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.658 | 0.129 | 0.498 | 0.799 | 100% | reference | | | |
| no_mix | 0.092 | 0.149 | -0.058 | 0.309 | 100% | 0.566 | -0.48 | 99.5% | visible
| reuse_tips | 0.626 | 0.141 | 0.444 | 0.795 | 100% | 0.032 | -23.93 | 61.5% | weak
| edge_wells | 0.569 | 0.160 | 0.373 | 0.769 | 100% | 0.089 | -8.76 | 69.5% | weak
| no_lid | 0.658 | 0.130 | 0.508 | 0.799 | 100% | -0.000 | -3179.86 | 44.8% | blind
| raw_ratio | 0.070 | 0.125 | -0.052 | 0.237 | 100% | 0.588 | -0.30 | 99.5% | visible
| narrow_range | 0.297 | 0.187 | 0.054 | 0.531 | 100% | 0.361 | -1.63 | 96.0% | visible
| thin_ladder | 0.022 | 0.170 | -0.058 | 0.264 | 100% | 0.636 | -0.41 | 99.5% | visible
| few_controls | 0.323 | 0.060 | 0.237 | 0.389 | 0% | 0.335 | -0.70 | 99.0% | visible
| random | 0.080 | 0.175 | -0.054 | 0.330 | 32% | 0.578 | -0.58 | 98.0% | visible
| null | 0.064 | 0.069 | 0.000 | 0.157 | 0% | 0.595 | 0.00 | 99.5% | visible
| gate_bypass | 0.638 | 0.130 | 0.449 | 0.776 | 100% | 0.020 | -37.67 | 58.5% | blind

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake); inverted = the mistake scores higher and wins over 60% of seeds (the reward prefers it).
