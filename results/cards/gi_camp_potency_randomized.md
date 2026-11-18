# Environment QC: gi_camp_potency

Rig: randomized. Seeds: 0 to 199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.600 | 0.176 | 0.359 | 0.794 | 94% | reference | | | |
| forskolin_too_high | 0.606 | 0.169 | 0.361 | 0.794 | 96% | -0.006 | -164.68 | 42.2% | blind
| skip_titration | 0.567 | 0.206 | 0.258 | 0.815 | 86% | 0.033 | -33.45 | 55.5% | blind
| no_forskolin | -0.058 | 0.000 | -0.058 | -0.058 | 0% | 0.658 | 0.20 | 99.8% | visible
| no_ibmx | 0.428 | 0.269 | 0.064 | 0.739 | 62% | 0.172 | -6.74 | 75.8% | weak
| raw_ratio | 0.250 | 0.154 | 0.060 | 0.459 | 94% | 0.350 | -1.83 | 96.0% | visible
| free_plateau | 0.541 | 0.172 | 0.293 | 0.745 | 94% | 0.059 | -16.55 | 64.0% | weak
| too_many_cells | 0.469 | 0.271 | 0.041 | 0.787 | 85% | 0.131 | -9.26 | 67.2% | weak
| no_mix | 0.061 | 0.148 | -0.058 | 0.279 | 92% | 0.540 | -0.80 | 99.0% | visible
| random | 0.004 | 0.113 | -0.060 | 0.095 | 21% | 0.596 | -0.45 | 99.0% | visible
| null | 0.065 | 0.065 | 0.000 | 0.151 | 0% | 0.535 | -0.35 | 98.5% | visible

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake); inverted = the mistake scores higher and wins over 60% of seeds (the reward prefers it).
