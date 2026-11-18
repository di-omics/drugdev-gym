# Environment QC: lead_opt

Rig: nominal. Seeds: 0 to 199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.816 | 0.078 | 0.703 | 0.870 | n/a | reference | | | |
| sloppy_measurement | 0.542 | 0.417 | -0.083 | 0.917 | n/a | 0.274 | -4.41 | 57% | weak
| no_model | 0.632 | 0.137 | 0.449 | 0.809 | n/a | 0.184 | -2.50 | 89% | weak
| one_round | 0.551 | 0.220 | 0.254 | 0.807 | n/a | 0.266 | -2.36 | 92% | visible
| null | 0.215 | 0.275 | 0.000 | 0.642 | n/a | 0.602 | -0.76 | 96% | visible

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake).
