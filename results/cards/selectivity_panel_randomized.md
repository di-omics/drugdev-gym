# Environment QC: selectivity_panel

Rig: randomized. Seeds: 0 to 199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.605 | 0.131 | 0.402 | 0.757 | 100% | reference | | | |
| on_target_range | 0.087 | 0.168 | -0.097 | 0.318 | 100% | 0.518 | -0.73 | 99.5% | visible
| single_point | 0.244 | 0.198 | -0.054 | 0.502 | 100% | 0.361 | -1.73 | 95.0% | visible
| single_point_as_curve | 0.110 | 0.305 | -0.277 | 0.505 | 100% | 0.496 | -1.64 | 91.5% | visible
| no_controls | -0.077 | 0.027 | -0.085 | -0.085 | 0% | 0.682 | 0.31 | 100.0% | visible
| no_standards | 0.109 | 0.175 | -0.089 | 0.356 | 100% | 0.496 | -0.84 | 100.0% | visible
| on_target_only | -0.097 | 0.001 | -0.097 | -0.097 | 100% | 0.702 | 0.44 | 100.0% | visible
| raw_ratio | 0.344 | 0.225 | 0.054 | 0.639 | 100% | 0.261 | -3.09 | 84.5% | weak
| claim_past_range | 0.092 | 0.308 | -0.325 | 0.466 | 100% | 0.514 | -1.56 | 94.0% | visible
| inactive_as_number | 0.407 | 0.185 | 0.144 | 0.653 | 100% | 0.199 | -3.76 | 84.5% | weak
| random | -0.545 | 0.038 | -0.591 | -0.511 | 47% | 1.150 | 0.56 | 100.0% | clean
| null | -0.500 | 0.000 | -0.500 | -0.500 | 0% | 1.105 | 0.65 | 100.0% | clean

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake); inverted = the mistake scores higher and wins over 60% of seeds (the reward prefers it).
