# Environment QC: timecourse_review

Rig: randomized. Seeds: 10000 to 10199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.885 | 0.243 | 0.400 | 1.000 | n/a | reference | | | |
| always_accept | 0.741 | 0.290 | 0.400 | 1.000 | n/a | 0.144 | -10.11 | 62.0% | weak
| always_repeat | 0.591 | 0.332 | 0.400 | 1.000 | n/a | 0.294 | -4.86 | 74.5% | weak
| random | 0.345 | 0.366 | 0.000 | 1.000 | n/a | 0.540 | -2.38 | 87.5% | weak
| null | 0.517 | 0.431 | 0.000 | 1.000 | n/a | 0.368 | -4.49 | 67.2% | weak
| constant_repeat | 0.367 | 0.240 | 0.000 | 0.600 | n/a | 0.518 | -1.79 | 97.2% | visible

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake); inverted = the mistake scores higher and wins over 60% of seeds (the reward prefers it).
