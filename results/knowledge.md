# Expert rules: what the rewards can teach

24 rules, 33 rule-trap pairs. The reward teaches 13 of them (the expert wins on at least 90% of seeds), 15 only on average, is blind to 5, and rewards the mistake in 0.
24 of the 24 rules are drafts seeded from general bench practice and have not been reviewed by the expert yet.

| env | rule | status | trap | reward teaches it | expert wins | reward gap |
|---|---|---|---|---|---|---|
| dose_response | Mix after every transfer in a serial dilution. | draft | no_mix | yes | 99.5% | 0.607 |
| dose_response | Change tips between ladder steps. | draft | reuse_tips | on average | 63.0% | 0.033 |
| dose_response | Keep assay wells off the outer ring. | draft | edge_wells | on average | 69.5% | 0.077 |
| dose_response | Lid the plate during the 37 C incubation. | draft | no_lid | no | 47.0% | -0.000 |
| dose_response | Convert the ratio to analyte through the standard curve before fitting potency. | draft | raw_ratio | yes | 100.0% | 0.620 |
| dose_response | Dose wide enough to see both plateaus when you do not know the potency. | draft | narrow_range | yes | 97.0% | 0.370 |
| dose_response | A ladder well has to hold everything you will take out of it, plus the dead volume. | draft | thin_ladder | yes | 100.0% | 0.739 |
| dose_response | Put enough negative and positive controls on every plate to compute a Z-prime. | draft | few_controls | yes | 99.5% | 0.349 |
| stability | Quench into internal standard and read the analyte over it. | draft | no_internal_standard | yes | 99.0% | 0.420 |
| stability | Mix the matrix after spiking, before the time-zero sample. | draft | no_mix | yes | 100.0% | 0.466 |
| stability | Lid the incubation plate and stay off the outer ring for hours-long time courses. | draft | no_lid | on average | 69.0% | 0.011 |
| stability | Lid the incubation plate and stay off the outer ring for hours-long time courses. | draft | edge_wells | on average | 66.0% | 0.007 |
| stability | Lid the incubation plate and stay off the outer ring for hours-long time courses. | draft | no_lid_on_edge_wells | yes | 98.0% | 0.189 |
| stability | Spread time points so both a five-minute and a five-hour half-life leave several usable points. | draft | two_point | on average | 85.5% | 0.290 |
| stability | Spread time points so both a five-minute and a five-hour half-life leave several usable points. | draft | early_only | yes | 90.0% | 0.157 |
| stability | Run at least two incubation wells per compound. | draft | single_replicate | yes | 100.0% | 0.464 |
| stability | Drop points below the limit of quantification before fitting. | draft | no_lloq | on average | 87.5% | 0.218 |
| assay_window | Titrate cell density; do not assume more cells or the kit default is right. | draft | max_cells | on average | 64.0% | 0.235 |
| assay_window | Titrate cell density; do not assume more cells or the kit default is right. | draft | min_cells | on average | 73.5% | 0.247 |
| assay_window | Titrate cell density; do not assume more cells or the kit default is right. | draft | kit_default | on average | 67.5% | 0.118 |
| assay_window | Judge the assay window after back-calculation, not on the raw ratio. | draft | raw_ratio_zprime | no | 49.0% | -0.000 |
| assay_window | Validate with the lid on and off the outer ring, the way the screen will run. | draft | no_lid | no | 59.5% | 0.000 |
| assay_window | Validate with the lid on and off the outer ring, the way the screen will run. | draft | edge_wells | on average | 80.5% | 0.055 |
| lead_opt | Do not rank analogs on a sloppy assay. | draft | upstream_no_mix | on average | 63.3% | 0.031 |
| lead_opt | Do not rank analogs on a sloppy assay. | draft | upstream_raw_ratio | no | 34.0% | 0.004 |
| lead_opt | Do not rank analogs on a sloppy assay. | draft | upstream_two_point | on average | 61.0% | 0.299 |
| lead_opt | Do not rank analogs on a sloppy assay. | draft | upstream_no_internal_standard | no | 54.0% | 0.011 |
| lead_opt | Use what each round taught you to choose the next analogs. | draft | no_model | on average | 89.5% | 0.185 |
| lead_opt | Run more than one design-make-test round before nominating. | draft | one_round | yes | 91.5% | 0.265 |
| plate_review | Check the controls, the replicates and the reference compound before you trust a potency. | draft | always_accept | on average | 77.0% | 0.568 |
| plate_review | Do not repeat a plate that passes its checks. | draft | always_repeat | yes | 92.5% | 0.554 |
| timecourse_review | Check time zero, the replicates and the control compound before you trust a half-life. | draft | always_accept | on average | 61.3% | 0.375 |
| timecourse_review | Do not repeat a time course that passes its checks. | draft | always_repeat | yes | 91.5% | 0.561 |

## What a learner rediscovered from reward alone

Cross-entropy learner on dose_response, one run per training rig (see learning.md). A rule counts as rediscovered if the learner's final design follows it.

| rule | reward teaches it | learned on nominal | learned on randomized |
|---|---|---|---|
| Mix after every transfer in a serial dilution. | yes | yes | yes |
| Change tips between ladder steps. | on average | yes | yes |
| Keep assay wells off the outer ring. | on average | no | no |
| Lid the plate during the 37 C incubation. | no | no | no |
| Convert the ratio to analyte through the standard curve before fitting potency. | yes | yes | yes |
| Dose wide enough to see both plateaus when you do not know the potency. | yes | yes | yes |
| A ladder well has to hold everything you will take out of it, plus the dead volume. | yes | yes | yes |
| Put enough negative and positive controls on every plate to compute a Z-prime. | yes | yes | yes |

Rediscovered: learned on nominal 6 of 8, learned on randomized 6 of 8.
