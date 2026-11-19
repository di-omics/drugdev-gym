# Campaign: what would price each mistake on a real rig

62 rule-trap pairs count. Sorted by the kind of choice the mistake changes:

- 15 can be priced from complete records of correct runs. 6 read the same raw data differently and 9 leave data off the run. Of the 15, 6 hold on a plain complete record and 9 need the record to keep something specific, listed below.
- 14 are a point on a grid. 10 are a setting of the assay, which moves between days and cell batches, so the grid has to be repeated. 4 are a choice of concentrations: 2 are the top well of the expert's own ladder and hold on a plain record today, and 2 need a second ladder until every well is dosed on its own.
- 14 change something physical on the rig: 10 how the liquid or the sample is handled, 4 where the wells sit. Direct measurement comes first, paired runs second.
- 6 leave out run QC. A correct run has nothing for a control to catch, so these show their value only on a run that carries the fault the control exists for.
- 13 are not plate questions: 3 reporting conventions, 4 assay mistakes made inside a campaign, 2 campaign strategy, 4 judgment on finished data.

The kinds are drafts the expert has not read. The simulator checks what each trap changes, by 4 rules. A reading or a convention asks the rig for exactly what the expert asks for; the check cannot tell those 2 apart, so that split is the author's call. A judgment is made where nothing is run. An upstream trap changes the same assay field, to the same value, as the assay trap it names. Any other trap changes a design field of its own kind, or drops a run of that kind, and every other field it changes is explained at the end. That is a check on the simulator's actions. Whether the named data would measure the cost on a bench is a judgment. The note on each pair comes from an automated review, with one reviewer prompted to read as a bench scientist. Nothing here is measured on hardware.

## By kind

| group | kind | what prices it | pairs | reward teaches it today | only on average | blind |
|---|---|---|---|---|---|---|
| from complete records of correct runs | reading | a complete raw record of a correct run, read again | 6 | 3 | 3 | 0 |
| from complete records of correct runs | data | a record of a correct run designed for it, with the same data left out | 9 | 6 | 3 | 0 |
| from grids | dose | a dense grid of concentrations where every well is dosed on its own; on a serial ladder, a second ladder | 4 | 4 | 0 | 0 |
| from grids | condition | a grid of assay conditions, repeated over days and cell batches | 10 | 1 | 6 | 3 |
| on the rig | handling | the physical quantity it changes, measured directly on the rig; then a few same-day paired runs | 10 | 8 | 2 | 0 |
| on the rig | position | uniformity plates and a split layout on one plate | 4 | 0 | 4 | 0 |
| only when something is wrong | qc | the same item left out of runs that carry a fault or a liability, and how often faults happen | 6 | 5 | 1 | 0 |
| not a plate question | convention | nothing on a plate: a rule about what may be claimed | 3 | 2 | 1 | 0 |
| not a plate question | upstream | the size of the error from its own assay environment; what it costs a campaign, only from a campaign | 4 | 0 | 2 | 2 |
| not a plate question | strategy | a whole campaign of experiments, not a plate | 2 | 1 | 1 | 0 |
| not a plate question | judgment | review cases with a planted answer, not a new run | 4 | 0 | 4 | 0 |

## By environment

| environment | from complete records of correct runs | from grids | on the rig | only when something is wrong | not a plate question |
|---|---|---|---|---|---|
| dose_response | 1 | 1 | 4 | 1 | 0 |
| stability | 4 | 0 | 3 | 1 | 0 |
| assay_window | 0 | 4 | 1 | 0 | 0 |
| lead_opt | 0 | 0 | 0 | 0 | 6 |
| plate_review | 0 | 0 | 0 | 0 | 2 |
| timecourse_review | 0 | 0 | 0 | 0 | 2 |
| gi_fluid_stability | 4 | 0 | 3 | 1 | 0 |
| gi_camp_potency | 2 | 5 | 1 | 0 | 0 |
| selectivity_panel | 3 | 3 | 0 | 1 | 2 |
| cytotox_window | 1 | 1 | 2 | 2 | 1 |

## From complete records of correct runs

A record prices one of these only if it kept, as raw values, what the mistake leaves out or reads differently.

| environment | trap | kind | what it changes | the record also needs |
|---|---|---|---|---|
| dose_response | raw_ratio | data | standards | nothing more |
| stability | no_internal_standard | data | internal_standard | a fixed injection order: the price holds for that order only |
| stability | two_point | data | timepoints_min | one well per time point or a large incubation, and the real quench time of every sample |
| stability | early_only | data | timepoints_min | one well per time point or a large incubation, and the real quench time of every sample |
| stability | no_lloq | reading | nothing the rig is asked for | blanks and low calibrators, so the limit of quantification is measured |
| gi_fluid_stability | late_first_point | data | timepoints_min | one well per time point or a large incubation, and the real quench time of every sample |
| gi_fluid_stability | no_true_time_zero | data | prequenched_t0 | both zeros on the run: the pre-quenched one and the earliest sample after the spike, with its real quench time |
| gi_fluid_stability | no_internal_standard | data | internal_standard | a fixed injection order: the price holds for that order only |
| gi_fluid_stability | ignore_enzyme_decay | reading | nothing the rig is asked for | the activity of each fluid lot followed over the incubation: a series of its own, since the control substrate is gone too soon |
| gi_camp_potency | raw_ratio | reading | nothing the rig is asked for | nothing more |
| gi_camp_potency | free_plateau | reading | nothing the rig is asked for | nothing more |
| selectivity_panel | on_target_only | data | -run_plate | nothing more |
| selectivity_panel | raw_ratio | reading | nothing the rig is asked for | nothing more |
| selectivity_panel | no_standards | data | standards | nothing more |
| cytotox_window | nominal_above_solubility | reading | nothing the rig is asked for | a turbidity read or a solubility measurement in the same medium |

A stability record that holds every time-point design here samples at 0, 5, 10, 15, 20, 30, 60, 120, 240, 360 min. That is more samples than one incubation well gives.

A gi_fluid_stability record that holds every time-point design here samples at 0, 2, 5, 10, 15, 20, 30, 45, 60, 90, 120, 240 min. That is more samples than one incubation well gives.

## From grids

| environment | grid over | prices |
|---|---|---|
| dose_response | dilution_factor, n_points, top_conc_nM | narrow_range |
| assay_window | cells_k, stim_min | kit_default, max_cells, min_cells, raw_ratio_zprime |
| gi_camp_potency | cells_k, forskolin_uM, ibmx | forskolin_too_high, no_forskolin, no_ibmx, skip_titration, too_many_cells |
| selectivity_panel | top_conc_nM, n_points | on_target_range, single_point, single_point_as_curve |
| cytotox_window | exposure_h | short_exposure |

A setting chosen for one day's cells is the answer for that day. The price of not titrating is how far the best setting moves between days and cell batches, so a condition grid is repeated, and a setting chosen from a pilot can move a second setting with it. A dose grid is a lookup only where every well is dosed on its own. Every environment here doses by serial ladder.

- dose_response/narrow_range: on a serial ladder this is a second ladder, run beside the first.
- assay_window/raw_ratio_zprime: what the trap changes is how the pilots are read; the condition it commits follows from that, and on some seeds it is the expert's.
- selectivity_panel/on_target_range: on a serial ladder this is a second ladder, run beside the first.
- selectivity_panel/single_point: it is the top well of the expert's own ladder, so a plain record of the expert's plate holds it today.
- selectivity_panel/single_point_as_curve: it is the top well of the expert's own ladder, so a plain record of the expert's plate holds it today; it also misreports the hit, which no rig field shows.
- cytotox_window/short_exposure: staggered dosing or a readout that does not end the plate, and a slow-acting reference compound.

## On the rig

Handling: measure the physical quantity the mistake changes, then confirm with paired runs, the run done right and the run done wrong on the same day. Position is an effect inside a plate.

| environment | trap | kind | reward teaches it today | direct measurement first | plate wells per run | LC-MS samples per run |
|---|---|---|---|---|---|---|
| dose_response | no_mix | handling | yes | a dye ladder with and without mixing | 60 | 0 |
| dose_response | reuse_tips | handling | on average | the ladder wells of a sticky reference compound read by LC-MS, with fresh and reused tips; a dye does not stick to a tip | 60 | 0 |
| dose_response | edge_wells | position | on average | uniformity plates and a split inner and outer layout | 60 | 0 |
| dose_response | thin_ladder | handling | yes | a dye ladder at the small and the large transfer volume | 60 | 0 |
| stability | no_mix | handling | yes | time-zero samples alone, mixed and unmixed | 0 | 48 |
| stability | edge_wells | position | on average | a stable marker compound in inner and outer wells of the incubation plate; there is no plate read here, so the uniformity study does not cover it | 0 | 48 |
| stability | no_lid_on_edge_wells | handling | yes | plate weights or a stable marker in the matrix; and the expert's answer first, since the lid alone is ruled moot | 0 | 48 |
| assay_window | edge_wells | position | on average | uniformity plates and a split inner and outer layout | 216 | 0 |
| gi_fluid_stability | acid_quench_in_sgf | handling | on average | each gastric sample split into the acid and the base quench | 0 | 128 |
| gi_fluid_stability | acid_quench_in_sgf_warm_overnight | handling | yes | the same split, injected again after the hold, with the temperature of the held plate | 0 | 128 |
| gi_fluid_stability | stale_enzyme | handling | yes | a control substrate in fresh and aged fluid, the age and storage stated | 0 | 128 |
| gi_camp_potency | no_mix | handling | yes | a dye ladder with and without mixing | 86 | 0 |
| cytotox_window | unmatched_dmso | handling | yes | a DMSO tolerance curve on the cells; here the mistake is punished by a plate-map rule, not by a simulated effect | 96 | 0 |
| cytotox_window | edge_wells | position | on average | uniformity plates and a split inner and outer layout | 84 | 0 |

If each of the 10 handling mistakes were also confirmed by paired runs, with the expert run beside it in each of the 5 environments, the assumed numbers give about 8 same-day pairs per arm: 120 runs, 4832 plate wells and 5248 LC-MS samples. The 8 is for a shift of 0.25 log units in the mean when one run scatters by 0.15 within a day (two-sided 0.05, power 0.8, a paired t-test by a normal approximation that can fall 1 pair short). Both numbers are assumptions, so this is a starting size and not a bound. A smaller shift, a mistake that mostly widens the scatter or fails a run now and then, or an arm that often returns no fitted value needs more, and 10 comparisons read at 0.05 each will show a false difference somewhere more often than that. Run sizes are one pass of each simulated design in the 96-well format (seeds 0, 1, 2).

## Only when something is wrong

| environment | trap | kind | what it changes | needs |
|---|---|---|---|---|
| dose_response | few_controls | qc | controls_per_side | a run with a planted cell, reader or reagent fault; a ladder mistake does not touch the control wells, so the handling arms cannot show it |
| stability | single_replicate | qc | replicates | runs with an unmixed spike, to see whether one replicate still flags them |
| gi_fluid_stability | no_activity_control | qc | activity_control | an aged or off-specification enzyme lot among the runs |
| selectivity_panel | no_controls | qc | controls_per_side | a run with a planted cell, reader or reagent fault; a ladder mistake does not touch the control wells, so the handling arms cannot show it |
| cytotox_window | no_kill_control | qc | kill_controls | a run with a planted cell, reader or reagent fault; a ladder mistake does not touch the control wells, so the handling arms cannot show it |
| cytotox_window | no_interference_control | qc | interference_control | a reference compound that interferes with the readout |

## Not a plate question

| environment | trap | kind | what it changes | priced by |
|---|---|---|---|---|
| lead_opt | upstream_no_mix | upstream | analogs, dose_response.mix | the same mistake as dose_response/no_mix: the size of the error from its own assay environment; what it costs a campaign, only from a campaign |
| lead_opt | upstream_raw_ratio | upstream | analogs, dose_response.standards | the same mistake as dose_response/raw_ratio: the size of the error from its own assay environment; what it costs a campaign, only from a campaign |
| lead_opt | upstream_two_point | upstream | analogs, stability.timepoints_min | the same mistake as stability/two_point: the size of the error from its own assay environment; what it costs a campaign, only from a campaign |
| lead_opt | upstream_no_internal_standard | upstream | analogs, stability.internal_standard | the same mistake as stability/no_internal_standard: the size of the error from its own assay environment; what it costs a campaign, only from a campaign |
| lead_opt | no_model | strategy | analogs | a whole campaign of experiments, not a plate |
| lead_opt | one_round | strategy | -make_test | a whole campaign of experiments, not a plate |
| plate_review | always_accept | judgment | nothing the rig is asked for | review cases with a planted answer, not a new run |
| plate_review | always_repeat | judgment | nothing the rig is asked for | review cases with a planted answer, not a new run |
| timecourse_review | always_accept | judgment | nothing the rig is asked for | review cases with a planted answer, not a new run |
| timecourse_review | always_repeat | judgment | nothing the rig is asked for | review cases with a planted answer, not a new run |
| selectivity_panel | claim_past_range | convention | nothing the rig is asked for | nothing on a plate: a rule about what may be claimed |
| selectivity_panel | inactive_as_number | convention | nothing the rig is asked for | nothing on a plate: a rule about what may be claimed |
| cytotox_window | exact_not_bound | convention | nothing the rig is asked for | nothing on a plate: a rule about what may be claimed |

## Every pair

| environment | trap | kind | what it changes | differs from the expert on | reward teaches it today |
|---|---|---|---|---|---|
| dose_response | no_mix | handling | mix | 3 of 3 seeds | yes |
| dose_response | reuse_tips | handling | new_tips | 3 of 3 seeds | on average |
| dose_response | edge_wells | position | avoid_edges | 3 of 3 seeds | on average |
| dose_response | raw_ratio | data | standards | 3 of 3 seeds | yes |
| dose_response | narrow_range | dose | dilution_factor, n_points, top_conc_nM, transfer_ul | 3 of 3 seeds | yes |
| dose_response | thin_ladder | handling | transfer_ul | 3 of 3 seeds | yes |
| dose_response | few_controls | qc | controls_per_side | 3 of 3 seeds | yes |
| stability | no_internal_standard | data | internal_standard | 3 of 3 seeds | yes |
| stability | no_mix | handling | mix_after_spike | 3 of 3 seeds | yes |
| stability | edge_wells | position | avoid_edges | 3 of 3 seeds | on average |
| stability | no_lid_on_edge_wells | handling | avoid_edges, lid | 3 of 3 seeds | yes |
| stability | two_point | data | timepoints_min | 3 of 3 seeds | on average |
| stability | early_only | data | timepoints_min | 3 of 3 seeds | yes |
| stability | single_replicate | qc | replicates | 3 of 3 seeds | yes |
| stability | no_lloq | reading | nothing the rig is asked for | no seed | on average |
| assay_window | max_cells | condition | -run_pilot, cells_k, stim_min | 3 of 3 seeds | on average |
| assay_window | min_cells | condition | -run_pilot, cells_k, stim_min | 3 of 3 seeds | on average |
| assay_window | kit_default | condition | -run_pilot, cells_k, stim_min | 3 of 3 seeds | on average |
| assay_window | raw_ratio_zprime | condition | cells_k, stim_min | 3 of 3 seeds | no |
| assay_window | edge_wells | position | avoid_edges | 3 of 3 seeds | on average |
| lead_opt | upstream_no_mix | upstream | analogs, dose_response.mix | 3 of 3 seeds | on average |
| lead_opt | upstream_raw_ratio | upstream | analogs, dose_response.standards | 3 of 3 seeds | no |
| lead_opt | upstream_two_point | upstream | analogs, stability.timepoints_min | 3 of 3 seeds | on average |
| lead_opt | upstream_no_internal_standard | upstream | analogs, stability.internal_standard | 3 of 3 seeds | no |
| lead_opt | no_model | strategy | analogs | 3 of 3 seeds | on average |
| lead_opt | one_round | strategy | -make_test | 3 of 3 seeds | yes |
| plate_review | always_accept | judgment | nothing the rig is asked for | no seed | on average |
| plate_review | always_repeat | judgment | nothing the rig is asked for | no seed | on average |
| timecourse_review | always_accept | judgment | nothing the rig is asked for | no seed | on average |
| timecourse_review | always_repeat | judgment | nothing the rig is asked for | no seed | on average |
| gi_fluid_stability | acid_quench_in_sgf | handling | sgf_quench | 3 of 3 seeds | on average |
| gi_fluid_stability | acid_quench_in_sgf_warm_overnight | handling | hold_min, hold_temp_c, sgf_quench | 3 of 3 seeds | yes |
| gi_fluid_stability | stale_enzyme | handling | fresh_enzyme | 3 of 3 seeds | yes |
| gi_fluid_stability | no_activity_control | qc | activity_control | 3 of 3 seeds | yes |
| gi_fluid_stability | late_first_point | data | timepoints_min | 3 of 3 seeds | on average |
| gi_fluid_stability | no_true_time_zero | data | prequenched_t0 | 3 of 3 seeds | on average |
| gi_fluid_stability | no_internal_standard | data | internal_standard | 3 of 3 seeds | yes |
| gi_fluid_stability | ignore_enzyme_decay | reading | nothing the rig is asked for | no seed | yes |
| gi_camp_potency | forskolin_too_high | condition | cells_k, forskolin_uM | 3 of 3 seeds | no |
| gi_camp_potency | skip_titration | condition | -run_titration, cells_k, forskolin_uM | 3 of 3 seeds | no |
| gi_camp_potency | no_forskolin | condition | cells_k, forskolin_uM | 3 of 3 seeds | yes |
| gi_camp_potency | no_ibmx | condition | cells_k, forskolin_uM, ibmx | 3 of 3 seeds | on average |
| gi_camp_potency | raw_ratio | reading | nothing the rig is asked for | no seed | yes |
| gi_camp_potency | free_plateau | reading | nothing the rig is asked for | no seed | on average |
| gi_camp_potency | too_many_cells | condition | cells_k | 2 of 3 seeds | on average |
| gi_camp_potency | no_mix | handling | mix | 3 of 3 seeds | yes |
| selectivity_panel | on_target_range | dose | top_conc_nM | 3 of 3 seeds | yes |
| selectivity_panel | claim_past_range | convention | nothing the rig is asked for | no seed | yes |
| selectivity_panel | single_point | dose | n_points | 3 of 3 seeds | yes |
| selectivity_panel | single_point_as_curve | dose | n_points | 3 of 3 seeds | yes |
| selectivity_panel | no_controls | qc | controls_per_side | 3 of 3 seeds | yes |
| selectivity_panel | on_target_only | data | -run_plate | 3 of 3 seeds | yes |
| selectivity_panel | raw_ratio | reading | nothing the rig is asked for | no seed | on average |
| selectivity_panel | no_standards | data | standards | 3 of 3 seeds | yes |
| selectivity_panel | inactive_as_number | convention | nothing the rig is asked for | no seed | on average |
| cytotox_window | unmatched_dmso | handling | match_dmso | 3 of 3 seeds | yes |
| cytotox_window | short_exposure | condition | exposure_h | 3 of 3 seeds | on average |
| cytotox_window | edge_wells | position | avoid_edges | 3 of 3 seeds | on average |
| cytotox_window | no_kill_control | qc | kill_controls | 3 of 3 seeds | yes |
| cytotox_window | no_interference_control | qc | interference_control | 3 of 3 seeds | on average |
| cytotox_window | nominal_above_solubility | reading | nothing the rig is asked for | no seed | yes |
| cytotox_window | exact_not_bound | convention | nothing the rig is asked for | no seed | yes |

## Fields a trap changes outside its own kind

Where a trap changes fields of 2 kinds, this is the author's choice of which one is the mistake.

- dose_response/narrow_range (dose) also changes transfer_ul: keeps the wells of the 2-fold ladder above the dead volume.
- stability/no_lid_on_edge_wells (handling) also changes avoid_edges: the mistake is the 2 together.

## Kinds nothing checks yet

No trap changes these design fields, so their kind is stated and unchecked. lead_opt is left out: it nests the 2 assay designs.

- assay_window: conditions
- dose_response: replicates
- gi_fluid_stability: fluids, replicates, sif_quench
- gi_camp_potency: controls_per_side, dilution_factor, n_points, replicates, standards, stim_min, top_conc_nM, transfer_ul
- selectivity_panel: avoid_edges, compounds, dilution_factor, lid, mix, new_tips, replicates, target, transfer_ul
- cytotox_window: cells_k, dilution_factor, lid, n_points, replicates, top_conc_uM, vehicle_controls
- whole runs no trap drops or adds: commit, run

## Ruled moot by the expert

Left out of every count above.

| environment | trap | kind |
|---|---|---|
| dose_response | no_lid | handling |
| stability | no_lid | handling |
| assay_window | no_lid | handling |
