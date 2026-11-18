# Environment QC: cytotox_window

Rig: randomized. Seeds: 0 to 199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.812 | 0.094 | 0.713 | 0.910 | 100% | reference | | | |
| unmatched_dmso | 0.395 | 0.042 | 0.345 | 0.441 | 0% | 0.417 | 0.02 | 99.0% | visible
| short_exposure | 0.667 | 0.188 | 0.429 | 0.892 | 100% | 0.145 | -4.85 | 82.0% | weak
| edge_wells | 0.719 | 0.165 | 0.484 | 0.886 | 91% | 0.093 | -7.34 | 77.0% | weak
| no_kill_control | 0.392 | 0.046 | 0.334 | 0.440 | 0% | 0.420 | -0.00 | 99.5% | visible
| no_interference_control | 0.719 | 0.163 | 0.467 | 0.916 | 100% | 0.093 | -7.29 | 62.0% | weak
| nominal_above_solubility | 0.546 | 0.212 | 0.277 | 0.814 | 100% | 0.266 | -2.45 | 93.8% | visible
| exact_not_bound | 0.033 | 0.257 | -0.282 | 0.373 | 100% | 0.779 | -0.35 | 100.0% | visible
| random | 0.127 | 0.193 | -0.051 | 0.345 | 24% | 0.685 | -0.26 | 100.0% | visible
| null | -0.500 | 0.000 | -0.500 | -0.500 | 0% | 1.312 | 0.78 | 100.0% | clean

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake); inverted = the mistake scores higher and wins over 60% of seeds (the reward prefers it).
