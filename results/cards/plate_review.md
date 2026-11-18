# Environment QC: plate_review

Rig: nominal. Seeds: 10000 to 10199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.957 | 0.168 | 1.000 | 1.000 | n/a | reference | | | |
| always_accept | 0.389 | 0.475 | 0.000 | 1.000 | n/a | 0.568 | -2.39 | 77.0% | weak
| always_repeat | 0.403 | 0.346 | 0.000 | 0.640 | n/a | 0.554 | -1.78 | 92.5% | visible
| random | 0.371 | 0.319 | 0.000 | 0.600 | n/a | 0.586 | -1.49 | 94.0% | visible

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake); inverted = the mistake scores higher and wins over 60% of seeds (the reward prefers it).
