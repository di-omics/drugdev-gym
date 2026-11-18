# Environment QC: selectivity_panel

Rig: nominal. Seeds: 0 to 199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.621 | 0.119 | 0.439 | 0.755 | 100% | reference | | | |
| on_target_range | 0.093 | 0.172 | -0.097 | 0.330 | 100% | 0.528 | -0.65 | 99.5% | visible
| single_point | 0.250 | 0.198 | -0.054 | 0.519 | 100% | 0.371 | -1.56 | 95.5% | visible
| single_point_as_curve | 0.115 | 0.308 | -0.288 | 0.528 | 100% | 0.506 | -1.53 | 92.0% | visible
| no_controls | -0.077 | 0.027 | -0.085 | -0.085 | 0% | 0.698 | 0.37 | 100.0% | visible
| no_standards | 0.112 | 0.175 | -0.089 | 0.358 | 100% | 0.509 | -0.73 | 100.0% | visible
| on_target_only | -0.097 | 0.000 | -0.097 | -0.097 | 100% | 0.718 | 0.50 | 100.0% | clean
| raw_ratio | 0.338 | 0.227 | 0.058 | 0.618 | 100% | 0.283 | -2.66 | 86.5% | weak
| claim_past_range | 0.102 | 0.315 | -0.324 | 0.473 | 100% | 0.520 | -1.50 | 94.5% | visible
| inactive_as_number | 0.416 | 0.180 | 0.190 | 0.632 | 100% | 0.205 | -3.36 | 85.8% | weak
| random | -0.545 | 0.038 | -0.591 | -0.511 | 47% | 1.166 | 0.60 | 100.0% | clean
| null | -0.500 | 0.000 | -0.500 | -0.500 | 0% | 1.121 | 0.68 | 100.0% | clean

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake); inverted = the mistake scores higher and wins over 60% of seeds (the reward prefers it).
