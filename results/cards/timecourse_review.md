# Environment QC: timecourse_review

Rig: nominal. Seeds: 10000 to 10199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.875 | 0.251 | 0.400 | 1.000 | n/a | reference | | | |
| always_accept | 0.500 | 0.451 | 0.000 | 1.000 | n/a | 0.375 | -4.61 | 61.3% | weak
| always_repeat | 0.314 | 0.368 | 0.000 | 1.000 | n/a | 0.561 | -2.31 | 91.5% | visible
| random | 0.354 | 0.330 | 0.000 | 0.600 | n/a | 0.521 | -2.35 | 89.2% | weak

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake); inverted = the mistake scores higher and wins over 60% of seeds (the reward prefers it).
