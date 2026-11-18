# Environment QC: lead_opt

Rig: randomized. Seeds: 0 to 199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.813 | 0.074 | 0.698 | 0.870 | n/a | reference | | | |
| sloppy_measurement | 0.535 | 0.408 | -0.083 | 0.917 | n/a | 0.278 | -4.20 | 57% | weak
| no_model | 0.631 | 0.138 | 0.449 | 0.809 | n/a | 0.182 | -2.49 | 88% | weak
| one_round | 0.554 | 0.215 | 0.281 | 0.793 | n/a | 0.259 | -2.34 | 92% | visible
| null | 0.215 | 0.275 | 0.000 | 0.642 | n/a | 0.598 | -0.75 | 96% | visible

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake).
