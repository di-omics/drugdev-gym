# Environment QC: timecourse_review

Rig: nominal. Seeds: 10000 to 10199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.882 | 0.247 | 0.400 | 1.000 | n/a | reference | | | |
| always_accept | 0.750 | 0.291 | 0.400 | 1.000 | n/a | 0.132 | -11.24 | 61.0% | weak
| always_repeat | 0.594 | 0.324 | 0.400 | 1.000 | n/a | 0.288 | -4.95 | 74.0% | weak
| random | 0.360 | 0.361 | 0.000 | 1.000 | n/a | 0.522 | -2.50 | 86.5% | weak
| null | 0.520 | 0.430 | 0.000 | 1.000 | n/a | 0.362 | -4.61 | 67.0% | weak
| constant_repeat | 0.364 | 0.240 | 0.000 | 0.600 | n/a | 0.518 | -1.82 | 96.7% | visible

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake); inverted = the mistake scores higher and wins over 60% of seeds (the reward prefers it).
