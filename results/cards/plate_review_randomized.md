# Environment QC: plate_review

Rig: randomized. Seeds: 10000 to 10199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.913 | 0.215 | 0.600 | 1.000 | n/a | reference | | | |
| always_accept | 0.356 | 0.454 | 0.000 | 1.000 | n/a | 0.557 | -2.60 | 78.2% | weak
| always_repeat | 0.436 | 0.331 | 0.000 | 0.600 | n/a | 0.477 | -2.43 | 86.0% | weak
| random | 0.368 | 0.308 | 0.000 | 0.600 | n/a | 0.545 | -1.88 | 90.5% | visible

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake); inverted = the mistake scores higher and wins over 60% of seeds (the reward prefers it).
