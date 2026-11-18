# Environment QC: lead_opt

Rig: nominal. Seeds: 0 to 199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.816 | 0.080 | 0.703 | 0.870 | n/a | reference | | | |
| upstream_no_mix | 0.785 | 0.110 | 0.647 | 0.870 | n/a | 0.031 | -17.52 | 63.3% | weak
| upstream_raw_ratio | 0.812 | 0.083 | 0.718 | 0.872 | n/a | 0.004 | -123.77 | 34.0% | blind
| upstream_two_point | 0.517 | 0.418 | -0.085 | 0.915 | n/a | 0.299 | -3.99 | 61.0% | weak
| upstream_no_internal_standard | 0.805 | 0.085 | 0.696 | 0.870 | n/a | 0.011 | -45.86 | 54.0% | blind
| no_model | 0.631 | 0.137 | 0.449 | 0.804 | n/a | 0.185 | -2.51 | 89.5% | weak
| one_round | 0.551 | 0.220 | 0.254 | 0.807 | n/a | 0.265 | -2.39 | 91.5% | visible
| null | 0.215 | 0.275 | 0.000 | 0.642 | n/a | 0.601 | -0.77 | 95.5% | visible

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake); inverted = the mistake scores higher and wins over 60% of seeds (the reward prefers it).
