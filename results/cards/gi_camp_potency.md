# Environment QC: gi_camp_potency

Rig: nominal. Seeds: 0 to 199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.627 | 0.153 | 0.432 | 0.799 | 96% | reference | | | |
| forskolin_too_high | 0.630 | 0.149 | 0.442 | 0.803 | 97% | -0.003 | -261.85 | 43.8% | blind
| skip_titration | 0.600 | 0.187 | 0.319 | 0.803 | 88% | 0.027 | -36.43 | 53.5% | blind
| no_forskolin | -0.058 | 0.000 | -0.058 | -0.058 | 0% | 0.685 | 0.33 | 99.8% | visible
| no_ibmx | 0.461 | 0.258 | 0.095 | 0.749 | 66% | 0.166 | -6.44 | 78.8% | weak
| raw_ratio | 0.261 | 0.150 | 0.080 | 0.449 | 96% | 0.366 | -1.48 | 97.5% | visible
| free_plateau | 0.558 | 0.153 | 0.346 | 0.738 | 96% | 0.068 | -12.39 | 69.5% | weak
| too_many_cells | 0.494 | 0.266 | 0.043 | 0.786 | 86% | 0.132 | -8.47 | 70.0% | weak
| no_mix | 0.040 | 0.102 | -0.058 | 0.183 | 97% | 0.587 | -0.30 | 99.5% | visible
| random | 0.012 | 0.133 | -0.061 | 0.161 | 23% | 0.615 | -0.39 | 98.5% | visible
| null | 0.065 | 0.065 | 0.000 | 0.151 | 0% | 0.562 | -0.16 | 99.5% | visible

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake); inverted = the mistake scores higher and wins over 60% of seeds (the reward prefers it).
