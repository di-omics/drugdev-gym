# Environment QC: gi_fluid_stability

Rig: nominal. Seeds: 0 to 199 (200 episodes per policy, same hidden truth across policies).

| policy | mean reward | sd | p10 | p90 | QC gate pass | gap to expert | Z' vs expert | paired win rate | verdict |
|---|---|---|---|---|---|---|---|---|---|
| expert | 0.834 | 0.133 | 0.822 | 0.893 | 92% | reference | | | |
| acid_quench_in_sgf | 0.800 | 0.159 | 0.391 | 0.885 | 88% | 0.033 | -25.20 | 81.2% | weak
| acid_quench_in_sgf_warm_overnight | 0.240 | 0.133 | 0.143 | 0.341 | 4% | 0.593 | -0.35 | 100.0% | visible
| stale_enzyme | 0.404 | 0.059 | 0.383 | 0.408 | 2% | 0.430 | -0.34 | 92.0% | visible
| no_activity_control | 0.390 | 0.016 | 0.371 | 0.408 | 0% | 0.443 | -0.01 | 99.5% | visible
| late_first_point | 0.746 | 0.158 | 0.515 | 0.880 | 90% | 0.088 | -8.94 | 87.0% | weak
| no_true_time_zero | 0.811 | 0.140 | 0.789 | 0.880 | 91% | 0.022 | -35.59 | 78.5% | weak
| no_internal_standard | 0.353 | 0.058 | 0.308 | 0.383 | 1% | 0.481 | -0.19 | 99.5% | visible
| ignore_enzyme_decay | 0.794 | 0.129 | 0.772 | 0.868 | 92% | 0.040 | -18.65 | 93.5% | visible
| random | 0.155 | 0.079 | 0.060 | 0.253 | 2% | 0.678 | 0.06 | 100.0% | visible
| null | 0.081 | 0.060 | 0.006 | 0.165 | 0% | 0.753 | 0.23 | 100.0% | visible

Verdicts: clean = Z' >= 0.5 between the two reward distributions; visible = the expert wins on >= 90% of seeds; weak = the expert wins on 60 to 90% of seeds, or the mean gap is >= 0.1 (punished on average, not reliably per episode); blind = neither (the reward cannot see this mistake); inverted = the mistake scores higher and wins over 60% of seeds (the reward prefers it).
