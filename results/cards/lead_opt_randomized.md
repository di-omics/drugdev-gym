# Environment QC: lead_opt

Rig: randomized. Seeds: 0 to 199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.810 | 0.078 | 0.697 | 0.870 | n/a | reference | | | |
| upstream_no_mix | 0.776 | 0.120 | 0.636 | 0.870 | n/a | 0.034 | -16.23 | 60.7% | weak
| upstream_raw_ratio | 0.813 | 0.074 | 0.723 | 0.872 | n/a | -0.003 | -151.36 | 31.0% | inverted
| upstream_two_point | 0.523 | 0.413 | -0.085 | 0.915 | n/a | 0.287 | -4.13 | 56.5% | weak
| upstream_no_internal_standard | 0.810 | 0.080 | 0.709 | 0.870 | n/a | -0.000 | -4070.67 | 50.8% | blind
| no_model | 0.628 | 0.143 | 0.448 | 0.809 | n/a | 0.182 | -2.64 | 88.2% | weak
| one_round | 0.554 | 0.215 | 0.281 | 0.793 | n/a | 0.256 | -2.43 | 92.0% | visible
| null | 0.215 | 0.275 | 0.000 | 0.642 | n/a | 0.595 | -0.78 | 95.5% | visible

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake); inverted = the mistake scores higher and wins over 60% of seeds (the reward prefers it).
