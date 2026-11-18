# Environment QC: timecourse_review

Rig: randomized. Seeds: 10000 to 10199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.876 | 0.249 | 0.400 | 1.000 | n/a | reference | | | |
| always_accept | 0.500 | 0.451 | 0.000 | 1.000 | n/a | 0.376 | -4.58 | 61.5% | weak
| always_repeat | 0.314 | 0.372 | 0.000 | 1.000 | n/a | 0.562 | -2.31 | 91.8% | visible
| random | 0.360 | 0.333 | 0.000 | 0.600 | n/a | 0.516 | -2.39 | 88.8% | weak

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake); inverted = the mistake scores higher and wins over 60% of seeds (the reward prefers it).
