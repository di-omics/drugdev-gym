# Environment QC: cytotox_window

Rig: nominal. Seeds: 0 to 199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.825 | 0.075 | 0.726 | 0.910 | 100% | reference | | | |
| unmatched_dmso | 0.399 | 0.038 | 0.352 | 0.438 | 0% | 0.426 | 0.21 | 100.0% | visible
| short_exposure | 0.678 | 0.184 | 0.448 | 0.893 | 100% | 0.147 | -4.27 | 80.2% | weak
| edge_wells | 0.765 | 0.093 | 0.636 | 0.881 | 100% | 0.060 | -7.43 | 78.5% | weak
| no_kill_control | 0.396 | 0.041 | 0.345 | 0.441 | 0% | 0.429 | 0.19 | 100.0% | visible
| no_interference_control | 0.724 | 0.160 | 0.480 | 0.912 | 100% | 0.101 | -6.01 | 67.0% | weak
| nominal_above_solubility | 0.552 | 0.211 | 0.291 | 0.831 | 100% | 0.272 | -2.15 | 93.2% | visible
| exact_not_bound | 0.039 | 0.262 | -0.284 | 0.372 | 100% | 0.786 | -0.28 | 100.0% | visible
| random | 0.130 | 0.198 | -0.051 | 0.352 | 24% | 0.694 | -0.18 | 100.0% | visible
| null | -0.500 | 0.000 | -0.500 | -0.500 | 0% | 1.325 | 0.83 | 100.0% | clean

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake); inverted = the mistake scores higher and wins over 60% of seeds (the reward prefers it).
