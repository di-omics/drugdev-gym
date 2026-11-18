# Environment QC: assay_window

Rig: randomized. Seeds: 0 to 199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.587 | 0.096 | 0.472 | 0.697 | 95% | reference | | | |
| edge_wells | 0.523 | 0.108 | 0.395 | 0.653 | 80% | 0.064 | -8.53 | 77.0% | weak
| no_lid | 0.587 | 0.096 | 0.471 | 0.697 | 94% | 0.000 | -3166.88 | 59.5% | blind
| raw_ratio_zprime | 0.585 | 0.096 | 0.462 | 0.699 | 95% | 0.002 | -255.74 | 51.5% | blind
| max_cells | 0.352 | 0.307 | -0.014 | 0.701 | 48% | 0.234 | -4.15 | 63.5% | weak
| min_cells | 0.339 | 0.297 | -0.010 | 0.691 | 44% | 0.247 | -3.77 | 72.5% | weak
| kit_default | 0.465 | 0.231 | -0.011 | 0.700 | 57% | 0.122 | -7.02 | 66.0% | weak
| random | 0.350 | 0.270 | -0.011 | 0.668 | 41% | 0.237 | -3.63 | 78.0% | weak

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake); inverted = the mistake scores higher and wins over 60% of seeds (the reward prefers it).
