# Environment QC: plate_review

Rig: randomized. Seeds: 10000 to 10199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.930 | 0.185 | 0.600 | 1.000 | n/a | reference | | | |
| always_accept | 0.681 | 0.311 | 0.400 | 1.000 | n/a | 0.249 | -4.97 | 70.8% | weak
| always_repeat | 0.639 | 0.335 | 0.400 | 1.000 | n/a | 0.291 | -4.36 | 74.2% | weak
| random | 0.391 | 0.343 | 0.000 | 0.600 | n/a | 0.539 | -1.94 | 89.5% | weak
| null | 0.505 | 0.473 | 0.000 | 1.000 | n/a | 0.425 | -3.64 | 71.8% | weak
| constant_repeat | 0.463 | 0.191 | 0.400 | 0.600 | n/a | 0.467 | -1.41 | 95.2% | visible

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake); inverted = the mistake scores higher and wins over 60% of seeds (the reward prefers it).
