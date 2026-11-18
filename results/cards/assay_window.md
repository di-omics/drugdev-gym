# Environment QC: assay_window

Rig: nominal. Seeds: 0 to 199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.600 | 0.078 | 0.510 | 0.695 | 99% | reference | | | |
| open_edges | 0.575 | 0.097 | 0.474 | 0.670 | 94% | 0.025 | -20.31 | 64.5% | weak
| raw_ratio_zprime | 0.573 | 0.109 | 0.458 | 0.686 | 92% | 0.027 | -19.51 | 59.2% | blind
| max_cells | 0.409 | 0.297 | -0.014 | 0.719 | 56% | 0.191 | -4.90 | 60.5% | weak
| min_cells | 0.352 | 0.296 | -0.010 | 0.683 | 48% | 0.248 | -3.53 | 73.5% | weak
| kit_default | 0.487 | 0.220 | 0.108 | 0.717 | 61% | 0.113 | -6.90 | 67.0% | weak
| random | 0.379 | 0.268 | -0.010 | 0.678 | 46% | 0.221 | -3.69 | 74.5% | weak

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake).
