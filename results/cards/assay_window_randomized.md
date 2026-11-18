# Environment QC: assay_window

Rig: randomized. Seeds: 0 to 199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.584 | 0.106 | 0.470 | 0.697 | 94% | reference | | | |
| open_edges | 0.546 | 0.111 | 0.433 | 0.666 | 87% | 0.037 | -16.41 | 69.7% | weak
| raw_ratio_zprime | 0.560 | 0.122 | 0.425 | 0.700 | 88% | 0.024 | -28.09 | 55.8% | blind
| max_cells | 0.393 | 0.301 | -0.014 | 0.712 | 53% | 0.190 | -5.42 | 59.0% | weak
| min_cells | 0.339 | 0.297 | -0.010 | 0.691 | 44% | 0.244 | -3.96 | 72.5% | weak
| kit_default | 0.467 | 0.229 | 0.011 | 0.700 | 57% | 0.117 | -7.60 | 65.5% | weak
| random | 0.358 | 0.268 | -0.010 | 0.677 | 42% | 0.225 | -3.97 | 77.0% | weak

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake).
