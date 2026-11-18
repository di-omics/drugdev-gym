# Environment QC: plate_review

Rig: nominal. Seeds: 10000 to 10199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.943 | 0.173 | 0.600 | 1.000 | n/a | reference | | | |
| always_accept | 0.697 | 0.303 | 0.400 | 1.000 | n/a | 0.246 | -4.80 | 70.5% | weak
| always_repeat | 0.643 | 0.332 | 0.400 | 1.000 | n/a | 0.300 | -4.05 | 75.0% | weak
| random | 0.397 | 0.338 | 0.000 | 0.600 | n/a | 0.546 | -1.80 | 90.5% | visible
| null | 0.511 | 0.478 | 0.000 | 1.000 | n/a | 0.432 | -3.52 | 71.5% | weak
| constant_repeat | 0.457 | 0.176 | 0.400 | 0.600 | n/a | 0.486 | -1.16 | 97.2% | visible

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake); inverted = the mistake scores higher and wins over 60% of seeds (the reward prefers it).
