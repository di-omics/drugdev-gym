# Environment QC: gi_fluid_stability

Rig: randomized. Seeds: 0 to 199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.816 | 0.155 | 0.399 | 0.891 | 88% | reference | | | |
| acid_quench_in_sgf | 0.793 | 0.167 | 0.391 | 0.885 | 86% | 0.023 | -40.05 | 80.8% | weak
| acid_quench_in_sgf_warm_overnight | 0.242 | 0.133 | 0.143 | 0.340 | 4% | 0.575 | -0.50 | 100.0% | visible
| stale_enzyme | 0.406 | 0.068 | 0.383 | 0.409 | 2% | 0.410 | -0.63 | 88.5% | weak
| no_activity_control | 0.390 | 0.016 | 0.370 | 0.408 | 0% | 0.426 | -0.21 | 96.0% | visible
| late_first_point | 0.731 | 0.175 | 0.384 | 0.879 | 87% | 0.085 | -10.60 | 88.0% | weak
| no_true_time_zero | 0.799 | 0.155 | 0.393 | 0.880 | 88% | 0.017 | -53.34 | 75.5% | weak
| no_internal_standard | 0.353 | 0.058 | 0.307 | 0.383 | 1% | 0.463 | -0.38 | 99.5% | visible
| ignore_enzyme_decay | 0.777 | 0.149 | 0.381 | 0.867 | 88% | 0.039 | -22.33 | 94.5% | visible
| random | 0.155 | 0.079 | 0.058 | 0.251 | 2% | 0.661 | -0.06 | 100.0% | visible
| null | 0.081 | 0.060 | 0.006 | 0.165 | 0% | 0.735 | 0.12 | 100.0% | visible

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake); inverted = the mistake scores higher and wins over 60% of seeds (the reward prefers it).
