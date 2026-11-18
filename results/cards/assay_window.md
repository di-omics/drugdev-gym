# Environment QC: assay_window

Rig: nominal. Seeds: 0 to 199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.599 | 0.079 | 0.509 | 0.695 | 99% | reference | | | |
| edge_wells | 0.545 | 0.085 | 0.442 | 0.639 | 90% | 0.055 | -7.96 | 80.5% | weak
| no_lid | 0.599 | 0.079 | 0.507 | 0.697 | 99% | 0.000 | -3142.27 | 59.5% | blind
| raw_ratio_zprime | 0.600 | 0.079 | 0.509 | 0.695 | 98% | -0.000 | -4046.01 | 49.0% | blind
| max_cells | 0.364 | 0.307 | -0.014 | 0.699 | 51% | 0.235 | -3.91 | 64.0% | weak
| min_cells | 0.352 | 0.296 | -0.010 | 0.683 | 48% | 0.247 | -3.54 | 73.5% | weak
| kit_default | 0.482 | 0.221 | 0.078 | 0.717 | 60% | 0.118 | -6.61 | 67.5% | weak
| random | 0.368 | 0.270 | -0.011 | 0.672 | 45% | 0.232 | -3.51 | 75.5% | weak

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake); inverted = the mistake scores higher and wins over 60% of seeds (the reward prefers it).
