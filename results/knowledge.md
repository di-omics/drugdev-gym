# Expert rules: what the rewards can teach

46 rules, 58 rule-trap pairs. The reward teaches 26 of them (the expert wins on at least 90% of seeds), 25 only on average, is blind to 7, and rewards the mistake in 0.
46 of the 46 rules are drafts seeded from general bench practice and have not been reviewed by the expert yet.

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
| plate_review | Check the controls, the replicates and the reference compound before you trust a potency. | draft | always_accept | on average | 70.5% | 0.246 |
| plate_review | Do not repeat a plate that passes its checks. | draft | always_repeat | on average | 75.0% | 0.300 |
| timecourse_review | Check time zero, the replicates and the control compound before you trust a half-life. | draft | always_accept | on average | 61.0% | 0.132 |
| timecourse_review | Do not repeat a time course that passes its checks. | draft | always_repeat | on average | 74.0% | 0.288 |
| gi_fluid_stability | Stop gastric samples by neutralising them, not with the acid stop used for intestinal fluid or plasma. | draft | acid_quench_in_sgf | on average | 81.2% | 0.033 |
| gi_fluid_stability | Stop gastric samples by neutralising them, not with the acid stop used for intestinal fluid or plasma. | draft | acid_quench_in_sgf_warm_overnight | yes | 100.0% | 0.593 |
| gi_fluid_stability | Dissolve the pepsin and the pancreatin fresh for every run. | draft | stale_enzyme | yes | 92.0% | 0.430 |
| gi_fluid_stability | Run a control substrate with a known half-life in every fluid, every run. | draft | no_activity_control | yes | 99.5% | 0.443 |
| gi_fluid_stability | Take the first samples within the first few minutes. | draft | late_first_point | on average | 87.0% | 0.088 |
| gi_fluid_stability | Make time zero by quenching the fluid before the compound meets the enzyme. | draft | no_true_time_zero | on average | 78.5% | 0.022 |
| gi_fluid_stability | Quench into internal standard and read the peptide over it. | draft | no_internal_standard | yes | 99.5% | 0.481 |
| gi_fluid_stability | Fit against what the enzyme could still do, not against the clock. | draft | ignore_enzyme_decay | yes | 93.5% | 0.040 |
| gi_camp_potency | Run the agonist plate with forskolin near its EC80, not at a saturating dose. | draft | forskolin_too_high | no | 43.8% | -0.003 |
| gi_camp_potency | Titrate forskolin on your own cells before the agonist plate; do not run it on handbook conditions. | draft | skip_titration | no | 53.5% | 0.027 |
| gi_camp_potency | A Gi agonist plate needs forskolin in every cell well. | draft | no_forskolin | yes | 99.8% | 0.685 |
| gi_camp_potency | Put IBMX in the stimulation buffer. | draft | no_ibmx | on average | 78.8% | 0.166 |
| gi_camp_potency | Convert the ratio to cAMP through the standard curve before fitting potency. | draft | raw_ratio | yes | 97.5% | 0.366 |
| gi_camp_potency | Pin the top of every inhibition curve on the forskolin-only controls. | draft | free_plateau | on average | 69.5% | 0.068 |
| gi_camp_potency | Titrate the cell number so the forskolin control lands inside the standard curve; more cells is not more signal. | draft | too_many_cells | on average | 70.0% | 0.132 |
| gi_camp_potency | Mix after every transfer in the agonist dilution. | draft | no_mix | yes | 99.5% | 0.587 |
| selectivity_panel | Dose the off-target plates up to the highest concentration the assay allows, not over the on-target range. | draft | on_target_range | yes | 99.5% | 0.528 |
| selectivity_panel | When nothing is seen, say 'at least' up to the top dose you tested and no further. | draft | claim_past_range | yes | 94.5% | 0.520 |
| selectivity_panel | A single concentration tells you whether something happens; run a curve wherever you need a potency. | draft | single_point | yes | 95.5% | 0.371 |
| selectivity_panel | A single concentration tells you whether something happens; run a curve wherever you need a potency. | draft | single_point_as_curve | yes | 92.0% | 0.506 |
| selectivity_panel | Put negative and positive controls on every off-target plate. | draft | no_controls | yes | 100.0% | 0.698 |
| selectivity_panel | Spend the plates on the panel: every off-target needs data before any window can be stated. | draft | on_target_only | yes | 100.0% | 0.718 |
| selectivity_panel | Read every target's plate through its own standard curve before you compare potencies. | draft | raw_ratio | on average | 86.5% | 0.283 |
| selectivity_panel | Read every target's plate through its own standard curve before you compare potencies. | draft | no_standards | yes | 100.0% | 0.509 |
| selectivity_panel | Report an inactive off-target as greater than the top dose, never as the top dose. | draft | inactive_as_number | on average | 85.8% | 0.205 |

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
