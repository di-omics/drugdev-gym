# Assumptions

Nothing in this repository has touched hardware, cells, or compounds. Every number below is
an assumed order-of-magnitude default chosen to be plausible, not a measurement. The point
of listing them is that each one is a named dial that real QC data can replace.

## Physical layer (src/drugdevgym/rig.py, RigProfile)

| parameter | default | what it stands for | how to measure it for real |
|---|---|---|---|
| cv_floor, cv_k | 1%, 4% x uL | pipetting CV = floor + k / volume | gravimetric or dye QC across volumes -> `RigProfile.from_qc` |
| bias_k | -3% x uL | small volumes under-deliver | same QC table |
| dead_volume_ul | 15 uL | volume a tip cannot reach in the source well | aspirate-to-air test per labware |
| carryover | 1% x stickiness | aspirate fraction left in a tip | blank after a high standard, per compound class |
| unmixed_recovery | 0.7 +/- 0.08 | solute fraction recovered from an unmixed well | dye ladder with and without mix steps |
| evap_edge / evap_inner | 4 / 0.5 uL per h | open plate at 37 C | weigh or dye-read a plate over time |
| lid_factor | 0.1 | fraction of evaporation a lid lets through | same, lidded |
| edge_cell_loss | 15% | weaker cell response on the outer ring | uniformity plate of controls |
| reader_cv, reader_floor | 2%, 30 ratio units | proportional and additive reader noise | repeated reads of one plate |
| row_gradient | 4% A to H | signal drift across rows | uniformity plate |

## Assay, compound and LC-MS constants

These are assumed too, and they set how large each mistake looks.

| constant | value | what it stands for | how to measure it for real |
|---|---|---|---|
| kit half-displacement k | 20 to 60 nM per seed | analyte at half signal in the competitive immunoassay | the kit standard curve |
| kit ratio window | 600 to 6000 | signal at saturating and zero analyte | the kit standard curve |
| standard series | 1 to 3000 nM, 8 points | on-plate standard curve | kit insert |
| basal and stimulated analyte | 0.1 to 0.3 x k and 3 to 8 x k | where the cell response sits on the curve; this ratio sets the size of the raw-ratio error | controls back-calculated on a real plate |
| cell-number CV | 4 to 8% (4 to 10% in assay_window) | well-to-well cell variation | uniformity plate |
| compound pEC50, Emax, Hill | 5.5 to 9.5, 0.6 to 1.0, 0.8 to 1.3 | hidden potency truth | not measurable in advance; set from public potency distributions |
| compound stickiness | log-normal, median 2 | multiplier on tip carryover | carryover test per compound class |
| cell-line productivity (assay_window) | 0.02 to 5 x k per 1000 cells | how much analyte a cell line makes | cell titration |
| half-life | 5 to 2000 min, log-uniform | hidden stability truth | set from public stability distributions |
| MS response variability | 18% log-scale | lumped sample-to-sample variability without an internal standard (matrix effect, recovery, drift), not injection repeatability alone | repeat injections and repeat extractions |
| internal-standard noise | 3% | residual after normalisation | same |
| detector noise floor | 0.4% of the time-zero area | sets the limit of quantification; the analysis default drops points under 2% | blank injections |
| cAMP kit half-displacement (gi_camp_potency) | 2.5 to 7.5 nM per seed | cAMP in the read well at half signal; the ratio window is the 600 to 6000 above, identical on every seed, so it can be memorized (see Known simplifications) | the kit cAMP standard curve |
| cAMP standard series | 0.04 to 712 nM, 8 points, 4-fold | on-plate cAMP standard curve, as final concentration in the read well | kit insert |
| Gi well volumes | 30 uL stimulation buffer over the plated cells, 10 uL agonist, 10 uL forskolin, 50 uL detection mix; 300 uL ladder wells | one 96-well stimulation and its lysis and detection step; the cells are plated the day before at the chosen number per well, and one 30 uL fill stands for the plating and the buffer exchange | the protocol in use |
| forskolin EC50 | 0.5 to 10 uM, log-uniform | hidden potency of forskolin on the cell line | forskolin titration read through the standard curve |
| forskolin slope | 1 | Hill slope of the forskolin response | the same titration |
| forskolin titration series | 100 uM top, half-log, 8 points, 40 uL ladder transfers | the fixed series of a titration plate, and the most forskolin a design can ask for | solubility and vehicle tolerance of the cells |
| cAMP output | 0.5 to 8 fmol per 1000 cells per min, log-uniform | hidden synthesis rate at saturating forskolin (15 to 240 pmol per million cells in 30 min) | cell titration at saturating forskolin with IBMX |
| unstimulated cyclase | 1 to 5% of the forskolin maximum | hidden cAMP synthesis with no forskolin; not inhibited by the agonist | wells without forskolin |
| PDE activity | 0.1 to 0.5 per min | hidden first-order cAMP breakdown with no inhibitor (half-life 1.4 to 7 min); sets how far the signal collapses without IBMX | cAMP time course with and without IBMX |
| IBMX | 500 uM in the stimulation, Ki 20 uM | one lumped competitive PDE inhibitor; leaves 1 / 26 of the PDE activity | IBMX titration on the cell line |
| full-agonist inhibition | 0.7 to 0.95 | hidden fraction of the forskolin response the reference agonist removes; sets the control window | forskolin controls with and without reference agonist, back-calculated |
| agonist efficacy | 0.6 to 1.0 of the reference | hidden; Emax of a test agonist is this times the full-agonist inhibition | not measurable in advance |
| Gi agonist pEC50, Hill, stickiness, cell CV | 5.5 to 9.5, 0.8 to 1.3, log-normal median 2, 4 to 8% | same ranges as dose_response | as in dose_response |
| Gi QC gate | Z' >= 0.5 on the raw ratio, 3 or more controls a side, every potency plate | when the reward is halved; a reported pEC50 must lie in 0 to 14 | the acceptance rule of the real screen |
| Gi expert operating point | forskolin control at 3 x the kit midpoint; titration at 10k cells; 10 uM forskolin if the titration cannot be fitted | where the scripted expert puts the assay: a policy choice made before the first card, not physics. The expert's curve fit is a change made after the first card: the first expert fitted free plateaus and scored 0.57 on 40 seeds; pinning the top plateau on the forskolin-only controls was added after seeing that number and moves the expert from about 0.56 to 0.63. The earlier analysis stays on the card as free_plateau | compare operating points on real plates |
| Gi skip_titration conditions | 10 uM forskolin, 30k cells | the handbook conditions of the trap that runs no titration plate. Added after the first card because the review asked for the trap; the two values are the ones the review named and were not adjusted. The verdict depends on the cell number: at 30k the trap loses 0.04 and is blind, at 10k (the expert's own fallback) it loses 0.19 and is weak | a cell titration on the real line |
| reward tolerance | 0.5 log units | error at which credit reaches zero | minimum significant ratio of the real assay |
| GI fluid pH (gi_fluid_stability) | SGF 1.2, SIF 6.8 | USP simulated gastric and intestinal fluid; the SGF value is also the reference pH for acid hydrolysis | pH meter on the prepared fluid |
| GI enzymatic half-life | SGF 1 to 3000 min, SIF 0.5 to 1000 min, log-uniform | hidden truth: digestion half-life at specification pepsin or pancreatin activity | set from public peptide SGF and SIF stability data |
| acid-labile compounds | 30% of compounds; non-enzymatic half-life 30 to 600 min at pH 1.2 and 37 C, log-uniform | hidden truth: loss in gastric acid with no enzyme, first order in [H+] | incubate in SGF without pepsin |
| enzyme activity half-life at 37 C | pepsin in SGF 24 h, pancreatin in SIF 4 h | first-order loss of proteolytic activity in the fluid (self-digestion), during the run and, slowed by the cold, in storage | control substrate spiked into fluid aged at 37 C for 0 to 4 h |
| Q10 | 2 | every rate (digestion, enzyme loss, acid hydrolysis) halves per 10 C of cooling; sets what a cold quench plate and a fridge are worth | the same rates at 4, 25 and 37 C |
| stale enzyme preparation | 24 h at 4 C | the batch that was not made for this run; with the two rows above it keeps 93% (pepsin) and 66% (pancreatin) of its activity | activity control on the aged batch |
| fresh-lot activity spread | 10% log-scale | lot potency, weighing and dissolution of a fresh preparation around specification | activity assay across preparations |
| residual activity after quench | pepsin: acid 1, base 0, organic 0.02, none 1; pancreatin: acid 0.01, base 1, organic 0.05, none 1 | fraction of the enzyme still working in the quenched sample, before the five-fold dilution and the cold are counted | control substrate spiked into quenched blank fluid and followed for a day |
| quenched-sample pH | SGF aliquot: 9 after base, otherwise 1.9; SIF aliquot: 2.5 after acid, 9 after base, otherwise 6.8 | how fast an acid-labile compound keeps hydrolysing while it waits | micro pH electrode on quenched blanks |
| activity control substrate | half-life 15 min at specification activity, both fluids | reference substrate that reads the enzyme activity of a run | validation runs with fresh, assayed enzyme |
| activity acceptance window | 0.8 to 1.25 | control-read activity outside this fails the run's QC | historical spread of the control substrate |
| handling time | 1 min | spike, mix and return for the first aliquot: the earliest later time point, and the real age of a time zero that was not pre-quenched | time the deck method |
| injection cycle | 3 min per sample, plate order, SGF plate first | how long each quenched sample waits in the queue after the hold | the LC-MS method and sequence table |
| GI quantifiable range and horizon | 1 to 960 min; incubation up to 240 min | reported and true half-lives are clipped to the range before scoring | assay validation with slow and fast reference peptides |
| fresh enzyme preparation | 20 cost units per fluid | what dissolving enzyme for one run costs, next to 2 per LC-MS sample | enzyme price and bench time |
| selectivity_panel: plate budget | 6 plates | one plate per target (one on-target, four off-targets) plus one spare for a repeat. With four compounds the budget never binds: all four fit on one 10-point curve plate per target, so the budget only forbids an uncovered panel (two wasted plates leave an off-target untested) and a full screen-then-curve plan (9 plates), and one curve plate per target is the only sensible plan. On the nominal rig the expert used exactly 5 plates on all 400 seeds measured and its gate never failed, so the spare was never used. Eight compounds would make the budget bind; that is a different environment and was not built | the campaign's real plate and cell budget |
| selectivity_panel: cost reference | 600 | 100 per budgeted plate, the same proportion as dose_response; the expert's five plates cost about 0.10 of reward | a real price list |
| selectivity_panel: on-target pEC50 | 7.0 to 9.0, uniform | hidden potency of a lead at its own target | the project's potency table |
| selectivity_panel: off-target margin | 3.0 to 8.0 log units below the on-target, uniform, independent per off-target | most off-targets are weak or inactive: 74% of compound and off-target pairs sit below the labware limit | panel data on a real series |
| selectivity_panel: liability rate and margin | 0.3; 0.0 to 2.0 log units below the on-target | a minority of compounds carry one close off-target | the same panel data |
| selectivity_panel: off-target efficacy, slope and cell lines | same ranges as the on-target (Emax 0.6 to 1.0, Hill 0.8 to 1.3 per target; each cell line draws its own unstimulated and stimulated level and cell CV; one kit lot for all five) | every off-target hit is an agonist read like the on-target | full curves and controls on each off-target line |
| selectivity_panel: labware top concentration | 100000 nM, a pEC50 floor of 4.0 | the highest final dose, so the widest window an episode can show is on-target pEC50 minus 4 | solubility and cell health at the top dose |
| selectivity_panel: tolerance | 0.5 log units of window error | the repository's reward tolerance, applied to the window and to how far a bound sits from the labware limit | minimum significant ratio of the two assays combined |
| selectivity_panel: support margin | 0.5 log units | how far a claimed window may exceed the tested-range limit before it counts as unsupported | a reporting policy, to set with the project |
| selectivity_panel: claim penalty | -0.5 per compound, and no floor on the mean | an unsupported or false claim of selectivity is worse than no answer, and a report made of such claims scores below no report. A failed QC gate halves credit, never a penalty. The floor at 0 was removed after the first card | a reporting policy |
| selectivity_panel: false-claim margin | 1.0 log unit | how far a claimed window may exceed the true window before it scores -0.5 instead of 0, whatever range was tested. Added after the first card at the size the review proposed; see What the rewards do and do not claim | how much overstated selectivity a project can absorb |
| selectivity_panel: plate QC | 3 controls a side, Z' 0.5 | the dose_response gate, per plate. Each compound and target pair rests on one plate: its QC-passing plate with the highest top dose, or its highest-dose plate of any kind when none passes. Support and gate are read from that plate, so a follow-up plate can add support and never takes it away from a pair that already has a QC-passing plate. Changed after the first card, which used the latest plate of each pair | n/a |
| selectivity_panel: accepted report range | -12 to 12 log units | input validation only, wider than any possible window | n/a |
| ACTIVE_FRACTION (analysis.estimate_potency) | 0.5 of the control window | a response under this at every dose is called "greater than the top dose" instead of fitted | bench convention; check it against partial agonists |
| POTENCY_FIT_STEP (analysis.estimate_potency) | 0.05 log units | grid resolution of the midpoint search; 0.02 moved the expert's mean reward by 0.003 on 150 seeds | numerical, not physical |
| selectivity_panel: reference protocol and the sizes of its mistakes | the expert's plates are the dose_response expert protocol with the top dose set to 3000 nM on-target and 100000 nM off-target. Traps: off-target top dose 3000 nM; a bound stated 1.0 log unit past the tested range; a single point is one well per compound at the top dose; no standard curve on the off-target plates. One labelled interaction: the single point plus a report that reads an active well as a potency at the dose tested. The null claims "at least 2.0" untested; the random policy draws its plates from SP_CHOICES in policies.py and claims a window drawn uniformly from 0 to 5, bound or number at random, for every compound | how large each reference mistake is | review with the expert |
| cytotox_window labware: STOCK_MM, CELLS_UL, DOSE_UL, FINAL_UL, INTER_DMSO_UL, INTER_MEDIUM_UL, REAGENT_UL | 10 mM; 90, 10, 100, 10, 90, 100 uL | stock in neat DMSO, and the volumes of the two 1:10 steps (DMSO solution into medium, dosing solution onto cells) and of the detection reagent | your own plate map and worklist |
| INTO_DOSING, INTO_ASSAY, MATCHED_DMSO_PCT, TOP_LIMIT_UM | 10, 10, 1%, 100 uM | derived from those volumes: the two dilution steps, the DMSO in every well of a matched run (and in the top well of any run), and the highest dose a 10 mM stock reaches at 1% DMSO. That dose is the dose limit and also the ceiling the reward scores a lower bound against. Between the second and the third pass the limit was 300 uM at up to 3% DMSO; see the last section | arithmetic; confirm the DMSO in a well by weight |
| LADDER_KEEP_UL, MEDIUM_LADDER_KEEP_UL | 45 and 70 uL | volume left in a ladder well after the transfer on, sized to cover what is taken from it later plus the dead volume, up to a dead volume of about 35 uL for the matched ladder and 45 uL for the unmatched one. About 1 randomized rig in 130 draws more; there the matched ladder under-delivers silently, the reference design passes its gate and its quality falls, and the unmatched route does not (found in re-review) | aspirate-to-air test per labware |
| CELL_FREE_REF_WELLS | 2 per plate | cell-free ATP wells with vehicle, the reference of the interference control | layout choice |
| EXPOSURES_H, REFERENCE_EXPOSURE_H | 24, 48 or 72 h; 72 h | allowed exposure times, and the exposure the true CC50 is defined at | protocol choice |
| DOUBLING_H | 24 h | cell doubling time | growth curve of the cell line |
| CONFLUENT_CELLS | 50,000 per well | cells a well holds at confluence; growth stops there, which hides toxicity in over-seeded wells | cell titration read after 72 h |
| CELL_CV | 4 to 8% per seed | well-to-well cell-number variation | uniformity plate |
| LUM_PER_CELL, LUM_BACKGROUND, ATP_STANDARD_LUM | 1 count per live cell, 100 counts, 20,000 counts | luminescence scale: signal per cell, reagent background, the cell-free ATP standard | cell titration, blank wells and an ATP standard on your reader |
| DMSO_IC50_PCT, DMSO_HILL | 2% v/v, slope 3 | DMSO alone: final cell number halves at 2%; about 11% is lost at 1% and 1.5% at 0.5% | DMSO titration on the cell line, read at 72 h |
| OSMOTIC_HALF, OSMOTIC_HILL | 0.4, slope 3 | evaporation: a 40% rise in medium concentration halves the final cell number | wells with a known volume removed, or medium made hypertonic |
| SLOW_EXPONENT, P_SLOW | 2; half of the compounds | slow compounds: CC50(t) = CC50(72 h) x (72 / t)^2, so 2.25x at 48 h and 9x at 24 h (a pure growth-arrest compound with a 24 h doubling time would give 2x at 48 h and no CC50 at 24 h) | reference compounds read at 24, 48 and 72 h |
| P_LIABILITY, LOG_CC50_LIABLE_UM, LOG_CC50_CLEAN_UM | 0.35; 1 to 100 uM; 100 to 3160 uM, log-uniform | hidden CC50 truth: a minority inside the dose range, the rest above it | not measurable in advance; set from public cytotoxicity distributions of your compound class |
| VIABILITY_HILL | 1.0 to 2.5 | slope of the viability curve | fits of reference cytotoxic compounds |
| LOG_SOLUBILITY_UM | 10 to 1000 uM, log-uniform | hidden solubility limit in the assay well | kinetic solubility (nephelometry) in the assay medium |
| TURBIDITY_BLANK, TURBIDITY_PER_UM | 500 units, 20 units per uM | turbidity pre-read: a clear well, and the rise per uM of compound out of solution | nephelometer blank and a compound of known solubility |
| P_INTERFERENCE, LOG_LUC_IC50_UM | 0.15; 1 to 50 uM, log-uniform | share of compounds that inhibit the luciferase readout, and their IC50 against it | cell-free ATP spike test; public luciferase counter-screens |
| LOG_EC50_UM | 1 nM to 1 uM, log-uniform | the functional potency given to the agent | the potency assay |
| STICKY_MEDIAN, STICKY_SD | 2, log-scale sd 0.8 | multiplier on tip carryover, the same draw as dose_response | carryover test per compound class |
| TOLERANCE, PENALTY_FLOOR | 0.5 log units, -0.5 | cytotox_window reward: error at which credit reaches zero, and the lowest score for one compound (a false alarm, a missed liability, or a bound with no run behind it) | minimum significant ratio of the real assay; the floor is a judgment about how bad a false alarm or a missed liability is |
| HEDGE_CAP | 0.5 | cytotox_window reward: the most a lower bound can score for a compound whose true CC50 the run dosed past. Added after the first card, following review; the value is the reviewer's example | a judgment about what a true statement that leaves a liability unreported is worth |
| ZPRIME_GATE, MIN_CONTROLS | 0.5, 3 wells | cytotox_window QC gate, the measured term: Z' between vehicle and full-kill luminescence on every plate, and the wells of each of the two controls a plate needs | convention |
| cytotox_window QC gate: solvent-matched vehicle | no number | the second gate term, a convention and not a measurement: the vehicle wells of a plate must carry the nominal DMSO percentage of the top dose on that plate. It is read from the plate map, never from the luminescence, so a run with match_dmso false fails it by rule. Added after the first card, in the third pass; see the last section | the solvent-control requirement of the cytotoxicity guideline you work to |
| CC50_REPORT_RANGE_UM | 1e-6 to 1e6 uM | range a reported CC50 is accepted in | not physical |
| PLATE_COST, WELL_COST, TIP_COST, cost_reference | 10, 1, 0.1; 600 | cytotox_window costs, the same round numbers as dose_response | your own consumables prices |
| TURBIDITY_SD_MULT, INTERFERENCE_FLAG, MIN_READOUT_FRACTION (analysis.py) | 4 blank SDs, 0.8, 0.05 | thresholds of the scripted CC50 analysis: when a dose counts as precipitated, when a compound counts as hitting the readout, and when a dose is unreadable. Set before the first measurement and not changed after it | blank scatter and control precision of your reader |

## Planted deviations in the review environments

Each is one step done wrong or one condition gone wrong; the sizes are assumed.

| deviation | what is changed | size |
|---|---|---|
| unmixed ladder, unmixed spike | the mix step is skipped | rig's own unmixed recovery |
| starved ladder | ladder transfer volume | 7 uL instead of 25 uL |
| tip carryover | one tip reused down each ladder | rig's own carryover |
| edge effect | cell response on the outer ring | 50% lost instead of 15% |
| reader drift | signal gradient from row A to row H | 35% instead of 4% |
| weak cells | stimulated analyte level | 1.3 to 2.0 x basal |
| evaporation | lid off, incubation wells on the outer ring | rig's own evaporation |
| failed internal standard | internal-standard CV | 30% instead of 3% |
| short aliquots | matrix volume per incubation well | 105 uL instead of 200 uL |

Review labels use a tolerance of 0.5 log units for pEC50 and 0.3 log units (2-fold) for half-life.
In time-course cases each compound gets its own MS response factor (log-scale sd 0.5) and each run
its own internal-standard level (sd 0.3), so absolute peak area carries no information.

## Known simplifications

- A well is either mixed or not; unmixed sampling applies the same recovery to every solute.
  Standing 20 minutes mixes a well by diffusion.
- Tip carryover does not wash out when the tip mixes the destination well.
- A fresh tip still keeps its carryover fraction, so sticky compounds lose material at every
  ladder step even with perfect technique. That is deliberate (adsorption), but the size is assumed.
- Bulk dispenses (cells, buffer, detection mix) have a flat 2% CV and always mix.
- Short aspirations are silent. A rig with pressure monitoring or liquid-level detection would flag them.
- Cell response is a Hill curve with no time dependence beyond the assay-window accumulation term;
  no cytotoxicity, no receptor desensitisation, no compound solubility limit.
- Detection is one competitive TR-FRET immunoassay with slope 1.
- Stability is single-phase first-order decay; no protein binding, no non-specific loss to plastic,
  no matrix lot effect. LC-MS is one lumped log-normal response factor per sample plus a noise floor.
- The lead-optimization landscape is additive with one pairwise term, 625 analogs. Real SAR is rougher.
- GI fluid stability treats each fluid as one protease activity, first order in enzyme and in substrate (far below Km),
  with no product inhibition and one Q10 for every rate. A quench acts at once and evenly; a base quench never harms the
  compound; peptides do not stick or precipitate when the pH changes. A pre-quenched time zero is an aliquot drawn at
  reaction time zero. Samples are injected in plate order at a fixed cycle time, so the wait of each sample is known.
  A compound's acid loss counts toward its SGF half-life. The MS noise model is the stability one, the same for every quench.
  Its hidden draws are salted with the environment's name, so it shares none with another environment on the same seed.
- In GI fluid stability, which way a quench that does not stop the enzyme biases a half-life comes from the injection
  order and from the expert's activity scaling as simulated, not from bench data. Cold and injected the same day, an
  acid quench in gastric fluid moves only peptides with half-lives under 3 min, which read about 0.06 log units too
  long. Left at 25 C overnight, the control substrate is digested away before injection, the run fails its gate, and the
  activity scaling then pulls half-lives over an hour about 0.27 log units too short (without that scaling they read
  0.02 too long). No rule states a direction, and none should until an expert has confirmed one.

- In gi_camp_potency the Gi pathway is noncompetitive with forskolin: the agonist removes the same fraction of
  the forskolin response at any forskolin concentration. In real cells saturating forskolin overrides part of
  the inhibition and shifts the apparent potency. That surmountable term is missing, so saturating forskolin is
  not punished here. The reward slightly prefers it: forskolin_too_high scores 0.003 above the expert on 4000
  seeds (7 standard errors) and the expert wins 44% of pairs. Most of that comes from the one seed in six where
  the expert's cells are already at the 80k limit and more forskolin is the only signal left (0.009 there, 0.002
  elsewhere). The card reads the trap as blind on most 200-seed windows and as inverted on about one in twenty
  (1 of 20 here, 1 of 27 in review). A surmountable term tried in memory (inhibition scaled by
  1 - 0.3 x drive^4) made the trap visible but took the expert to 0.51, so it needs its own constants and its
  own card and is not in the simulator.
- Skipping the forskolin titration is nearly free for the same reason. With the forskolin level harmless, the
  titration only earns the cell number, and 30k cells is a fair guess for most cell lines in these ranges, so
  skip_titration loses 0.04 (the expert wins 55% of pairs) and the card shows it as blind (weak on 2 of 20
  windows of 200 seeds).
- The kit ratio window is the same on every seed, as in dose_response, so a policy can convert ratio to cAMP
  with a memorized window instead of on-plate standards. One plate at 100 uM forskolin and 30k cells with no
  titration and no standards scores 0.60 against the expert's 0.63, at lower cost. Drawing the window per seed
  within 25% of nominal was tried in memory and cost that policy only 0.01, so the constant window was kept.
  What the card measures under the standard-curve rule is fitting the raw ratio, not skipping the standards.
- gi_camp_potency cells are plated the day before, so the rig's edge_cell_loss applies to its outer ring. That
  term is why a design that spills onto the outer ring scores below one that stays on the inner 60: 9 points
  in duplicate on 92 wells reach accuracy 0.60 against 0.71 for the expert's singles, and 0.74 with the term
  switched off.
- cAMP there is constant synthesis against first-order breakdown: no feedback activation of PDEs, no
  desensitisation, no export. Without IBMX the signal is a low steady state, not a peak that then falls, so it
  is small but not time dependent. no_ibmx is only weak on the card because the same policy answers the low
  signal with more cells: with the expert's cell number held it scores 0.24 against 0.64 and the expert wins
  96% of pairs.
- cAMP scales linearly with cell number, forskolin and agonist act at once, the reference agonist is always
  saturating, and IBMX-insensitive PDEs are ignored.
- gi_camp_potency plates always use fresh tips and a lid and fill the inner 60 wells first; those choices are
  taught in dose_response.
- Costs are round numbers that only set the relative price of wells, tips, samples, and synthesis.
- The four original environments share their first hidden draws on a given seed (the review
  environments do not). Do not reuse a seed across environments where one's revealed truth could
  inform another.
- gi_camp_potency is salted like the review environments. Its first card was measured before the salt was
  added, on other hidden draws from the same ranges, so every number quoted for it here was measured again
  afterwards.
- selectivity_panel: a plate carries one cell line. There is no solubility limit, no cytotoxicity and no
  assay interference at 100 uM, which is what limits the top dose of a real panel. Off-target margins are
  independent across the four off-targets and unrelated to structure. Nothing responds more than the
  reference agonist and nothing less than 0.6 of it, so there are no antagonists, no silent binders and no
  weak partial agonists; analysis.estimate_potency leans on both (see its docstring).
- selectivity_panel: the reward judges support on the nominal top concentration of the plate a pair rests on.
  A plate the compound never reached (a starved ladder) still counts as tested. A token plate with controls at
  the same top dose also satisfies the gate for a pair whose real plate failed QC: one off-target curve plate with
  2 controls a side scores 0.27, and adding a single-point token plate with 3 controls lifts it to 0.63 (seeds
  60000 to 60199, found in re-review). One spare plate covers one target, and a plate with 3 or more controls
  never failed QC for the expert, so a rule-following agent gains nothing; it is the same kind of hole as
  gate_bypass in dose_response and it is not registered as a policy. What catches the false
  "inactive" it produces is the false-claim rule, and only where the truth is more than a log unit below the
  claim. A 5 uL ladder transfer on the off-target plates (measured, not registered as a trap) scores 0.09
  against 0.63 for the expert on seeds 40000 to 40199 and the expert wins 98% of pairs, because its bounds on
  close liabilities score -0.5; its false bound on a weak off-target within a log unit of the labware limit
  still scores 0.
- selectivity_panel: a true bound earns credit only within 0.5 log units of the labware limit. That is a rule,
  not physics: nothing in the simulator limits the top dose, so dosing to 100 uM is free and always right.
  Off-target plates topped at 50 uM score 0.39, at 30 uM 0.24 and at 10 uM 0.15, against 0.63 for the expert
  (seeds 40000 to 40199). A true bound from 30 uM or below earns the same as no answer.
- selectivity_panel: a pair with a QC-failing plate at a high top dose and a passing plate at a lower one rests
  on the passing plate, whichever was run first. Adding that passing low-dose plate moves the pair from "supported
  to the high dose, reward halved" to "gate passes, supported to the low dose", and a claim left at the high
  dose then scores -0.5. The spec states the rule; the remedy is the repeat at the high dose.
- cytotox_window: viability is one Hill curve times independent factors for DMSO and for
  evaporation. No synergy between DMSO stress and compound, no time dependence of the DMSO effect.
- cytotox_window: slow or fast is a two-class power law in exposure time. No recovery, no
  cytostatic plateau. The true CC50 is the 72 h one, whatever exposure the agent runs.
- cytotox_window: growth is exponential up to a hard cap, and luminescence is proportional to live
  cells and does not depend on the well volume.
- cytotox_window: solubility acts in the assay well only. The 10x dosing solutions transfer in full
  even above the limit, DMSO does not raise the limit, and precipitate neither redissolves when the
  reagent is added nor changes the luminescence.
- cytotox_window: the cell-free control sees exactly the luciferase inhibition the cell wells see.
  On a real plate the two agree only roughly.
- cytotox_window: the rig loses compound to the tip at every transfer, in neat DMSO the same as in
  medium. Scripted CC50s therefore read high, by 0.06 log units on average on the nominal rig.
  Real adsorption is mostly from aqueous solution, so the DMSO ladder is treated too harshly.
- cytotox_window: a 10 mM stock reaches 100 uM at 1% DMSO and the dose limit is set there, so no
  design can put more than 1% DMSO in a well. For a lower top dose the DMSO working solution is the
  stock diluted in DMSO, a free step that is not pipetted on the rig, so every run puts 1% DMSO in
  its top well, a 30 uM run included. A bench would often run low doses at less.
- cytotox_window: in that range an unmatched vehicle lowers the top well by about 11% and cannot
  produce a false CC50. The reward does not see the mistake through the data. Skipping the DMSO
  match is punished by a gate convention (the vehicle wells must carry the DMSO of the top dose),
  not by a simulated effect. The last section gives the numbers with the convention switched off.
- cytotox_window: DMSO and compound act as independent factors, so a matched vehicle leaves the
  normalised curve unbiased by construction. Real cells under solvent stress can respond to a
  compound differently. That is not simulated.
- cytotox_window: the edge effect is the rig's evaporation (its open-plate rate times the lid
  factor, not a humidified-incubator measurement) plus its edge cell loss.

## What the rewards do and do not claim

- Truth is known by construction, so every reward is exactly verifiable. That is the whole value of
  simulation here, and also its limit: a policy that scores well has learned this simulator.
- "Expert" means a scripted policy that follows standard bench practice. It is a reference point,
  not an upper bound. The experts of the four original environments are not tuned. The review
  judges are: their thresholds and their list of serious causes were set on seeds 0 to 399, so
  they are reported on other seeds.
- The gi_camp_potency expert's protocol was fixed before its first card, but its curve fit was not: the
  pinned top plateau was chosen after the first measurement (see the Gi expert operating point row). That
  choice was made on low seeds; on seeds 40000 to 40199 the expert scores 0.65 against 0.63 on seeds 0 to 199,
  and every trap gets the same verdict.
- The transfer experiment in results/learning.md is sim-to-sim: train on one simulated rig, test on
  a different simulated rig. It says nothing yet about a physical robot.
- The gi_fluid_stability expert is a reference, not this reward's optimum. Its eight time points are over-sampled
  for this reward: [0, 3, 30, 240] with the same analysis is as accurate and cheaper, scores 0.846 against 0.819 and
  wins 89% of 300 fresh seeds. The expert is also halved on about 9% of seeds by things one run cannot avoid (a fresh
  lot whose activity lands outside the window, or two time-zero replicates that disagree), which is why none of its
  traps separates cleanly and stale_enzyme sits on the 90% line.
- What the gi_fluid_stability reward sees of each rule. Stale enzyme loses through the QC gate, not through the
  answer: the control substrate corrects the half-lives (accuracy 0.960 against 0.967) and reads the intestinal fluid at
  0.65, so 98% of runs fail the activity window. Stale pepsin on its own is not punished: it keeps 93% of its
  activity, and in gastric-only runs yesterday's batch scores about 0.01 below a fresh one and wins 70% of seeds
  because it is cheaper. A late first point, a time zero taken after the spike and the acid quench in gastric fluid
  are punished on average only, each through peptides with half-lives under about 3 min. The acid quench is ruinous
  only together with a warm overnight wait; that pair is on the card as a labelled interaction
  (acid_quench_in_sgf_warm_overnight), not as a single mistake. ignore_enzyme_decay turns one analysis flag that acts
  on the compound fit and on the control read alike; the two partly cancel, and the loss is the compound fit reading
  intestinal half-lives over an hour 0.10 log units too long.
- A gi_fluid_stability pair that was not run counts at half, whatever the rest of the run did. Before that rule a
  12-sample token run (one compound in gastric fluid, stale pepsin) lifted constant guesses for the five other pairs
  from 0.08 to 0.26; it now scores 0.20, and what it gains over 0.08 is the one pair it measured, less its cost.
- gi_fluid_stability changed after its first card. Before review, as its build report states: the fit applies the
  limit of quantification per time point instead of per injection (expert 0.827 to 0.836 on seeds 0 to 199), and a
  single-replicate trap was replaced by ignore_enzyme_decay. After an adversarial review: hidden draws were salted,
  observations became copies of the data, the half credit for a pair that was not run was added, the acid quench left
  warm overnight was registered as an interaction, worklist steps name the hold temperature and the age of the
  enzyme, and the rules on fresh enzyme and enzyme decay were reworded to the mechanism that was measured. The review
  fixes changed no physical constant and nothing in the expert's protocol or analysis. The card and the measurements
  quoted above are from after these changes; the two build-report numbers are from before the salt.
- The selectivity_panel expert's analysis was changed after its first measurement. As first written,
  analysis.estimate_potency fitted a free 4-parameter curve to every compound and only then asked whether
  the top dose was active; on seeds 0 to 199 that expert scored 0.531. Three choices were then added, on
  seeds 1000 to 1199: (1) no potency is reported unless the response reaches half the control window, because
  the foot of an incomplete curve was being fitted with a midpoint inside the tested range; (2) "reaches" means
  at any dose, not only the top one, because noise at the top dose was turning complete partial-agonist curves
  into bounds; (3) the fitted plateau is held at or below the positive control (fit_4pl_capped), because the
  free fit of an incomplete curve put the plateau several control windows high and the midpoint past the top
  dose. All three are ordinary bench conventions, but they were chosen after looking, and the third relies on
  the simulator having no response above the reference agonist. No simulator constant, reward constant or trap
  size was changed in that round; the reward rule changed later, after review (the next items).
- selectivity_panel changed after its first card. That card, on seeds 0 to 199, read expert 0.621, null 0.000,
  five traps visible or clean and two weak. An adversarial review then found four holes in the reward, and the
  fixes are the ones it proposed:
  (1) Quality is no longer floored at 0, and a failed gate halves credit but not a penalty. Before, a report
  made only of unsupported claims scored the same as no report (0.000 with nothing run); it now scores -0.500.
  (2) Support and gate rest, per compound and target pair, on the QC-passing plate with the highest top dose
  instead of the latest plate. Before, a sixth plate that re-tested one off-target at 3 uM took the expert from
  0.625 to 0.009 with its report unchanged; it now costs only its price (0.605).
  (3) A claim more than 1.0 log unit above the true window scores -0.5. Before, calling a compound selective
  against the agent's own off-target curve cost nothing. The labware-limit bound claimed for every compound
  without reading the off-target data scored 0.30 on token screens; it now scores 0.12 there and 0.06 on full
  curves.
  (4) A compound never tested on the on-target has no supported window. Before, four token off-target screens
  and "exactly 3.5" for every compound scored 0.18 with the gate passing; they now score -0.52.
  The "before" numbers are the review's and the "after" numbers were measured again, all on seeds 40000 to 40199.
- The 1.0 in (3) is a constant chosen after a measurement, and not by the builder. The review sized it on
  seeds 40000 to 40199 so that no honest claim there is touched: the expert's largest overstatement on those
  seeds is +0.85 and on_target_range's is +0.93. It is not clear of the expert everywhere. On seeds 50000 to
  51999 it penalises 6 of 7994 expert claims and on seeds 1000 to 1199 1 of 800. Every one is an off-target
  partial agonist (Emax 0.60 to 0.65, pEC50 5.0 to 5.4): six were called inactive by the expert's
  half-of-the-control-window rule, the limit already listed for analysis.estimate_potency, and one was fitted
  0.9 log units too weak. Those claims scored 0 before and -0.5 now. The margin was not adjusted after these
  counts.
- In the same round the card gained two rows and the random policy a report. no_standards leaves the standard
  curve off the off-target plates, the one-field form of the raw-ratio mistake. single_point_as_curve is the
  single point plus a report that reads an active well as a potency at the dose tested; it is two changes and
  is registered as a labelled interaction, like the one in gi_fluid_stability. The random policy used to report
  nothing (0 claims in 800 on the review's seeds), so its row measured plate cost only; it now claims a window
  for every compound. Two rules were reworded to the mechanism that was measured (off-target dose range and
  analyte units). Nothing in the expert's protocol or analysis, no simulator constant and no trap size changed
  in this round. After it, on seeds 0 to 199: expert 0.621, null -0.500, seven trap rows visible or clean (six
  single mistakes and the labelled interaction), two weak, none blind or inverted.
- The selectivity_panel null row is the untested claim and scores -0.5. An empty report scores 0.0 without
  running a plate, and that, not the null row, is the floor a do-nothing agent can reach.
- What the selectivity_panel reward sees of each rule. Two rows are weak because of how the mistake acts, not
  because the reward misses it. raw_ratio shifts the on-target and the off-target fits by almost the same
  amount (0.55 and 0.57 log units on seeds 40000 to 40199), so measured windows barely move; most of the loss
  is in compounds with nothing measurable (mean score 0.10 against 0.77 for the expert), where 123 of 268 get a
  number instead of a bound and 101 claims (94 bounds and 7 numbers) are pushed past the tested range by the
  on-target shift alone. That cancellation holds
  only because all five cell lines sit on the same part of one kit's standard curve. inactive_as_number ties
  the expert on 30 of those 200 seeds (on 21 of them no compound gets a bound) and beats it on 13, each time
  through a compound whose true window lies within 0.35 log units of the labware limit, where the limit
  written as a number is nearly right. Four of the visible or clean rows are punished by a rule and not by
  the simulated data: no_controls and on_target_only score about -0.08 and -0.10 because the expert's own
  analysis refuses to state a window without a control window or a full panel, so they pay the plate cost and
  report next to nothing; claim_past_range is punished by the support rule; and on_target_range by the rule
  that a bound earns credit only near the labware limit (570 of its 597 bounds are true and all 597 score 0).
  single_point and no_standards lose through the data: the first cannot give a potency for an active well,
  the second reads every off-target curve 0.57 log units too potent. A policy that skips the judgment but not all of the work still earns
  something: single-point screens, one follow-up curve and guessed windows for the remaining hits score 0.45
  against 0.63 for the expert on seeds 40000 to 40199, and the expert wins 74% of pairs. Only overstated
  selectivity is penalised, so a guessed low window on a single-point hit earns guess credit: "exactly 0.75" for
  every active well scores 0.32 against 0.24 for no answer and is never worse (seeds 60000 to 60199, found in
  re-review). The single_point row is the honest form of that plan, not its best-scoring report.

## cytotox_window: what was fixed first, and what changed after the first card

Before the first card:

- The magnitudes above, the first expert design (2 wells per dose, 8 points, 100 uM top dose, matched DMSO, 72 h,
  4,000 cells, 5 vehicle and 5 full-kill wells, cell-free control, inner wells, lid) and the three analysis
  thresholds were chosen before the first measurement. The thresholds and the scripted analysis have not changed.

After the first card, before review (from the build report):

- One simulator step. The unmatched protocol took its top dosing solution as a ready-made reagent, which gave it
  one transfer fewer per dose than the matched protocol. It is now pipetted from the DMSO stock on the rig, so every
  dose goes through the same number of transfers either way. The unmatched_dmso row moved from a gap of -0.008 and
  33.5% expert wins to -0.005 and 35.0%: inverted both times.

After an adversarial review (verdict: fix first), the second pass. Each item changed the reward rule, a constant or
the expert after the first card. These four are still in the code; their numbers were measured again on the final
code:

- A bound is scored against what the run tested. T (the highest nominal dose the run put on cells) and
  S = min(B, T) are new. Before, a bound above the tested range was paid in full: the review ran the first design
  to 30 uM on seeds 40000 to 40199 and scored 0.444 with honest bounds and 0.678 with every bound raised to "at
  least 100 uM". Now, with the reference design on the same seeds: 0.466 and -0.131. With no run every bound scores
  -0.5, which is why the null row fell from 0.107 to -0.500; the best report with no run is the empty one, at 0.
  T is the nominal top dose whatever the exposure time, so a 24 h token run still counts as having tested 100 uM.
- HEDGE_CAP = 0.5 is new. Before, relabelling a measured CC50 as a lower bound cost nothing (the review: 1,079 of
  1,082 scores unchanged). Now a bound scores at most 0.5 when the run dosed past the true CC50 (L < min(C, T)).
  Every expert answer forced to a bound: 0.720 against 0.825 on seeds 0 to 199, 0.701 against 0.810 on seeds 40000
  to 40199. The cap is a step: a compound whose CC50 sits just under the top dose can honestly look clean there,
  and its bound is capped all the same.
- The expert changed. The review found the first design dominated by cheaper ones. The registered expert is
  1 well per dose and 6 points (it was 2 and 8). The four combinations, on seeds 20000 to 20199, which no card
  uses: 1 well and 6 points 0.820, 1 and 8 0.816, 2 and 6 0.807, 2 and 8 0.795. Quality is 0.853 to 0.857 for all
  four, so the difference is cost. The order is the same on seeds 0 to 199 (0.825, 0.820, 0.816, 0.804) and on
  40000 to 40199 (0.810, 0.802, 0.797, 0.784). Replicates still have no job in this reward: well noise is about 6%
  against a tolerance of 0.5 log units. Cheaper neighbours were tried and not adopted, on the same three seed
  sets: 3 controls of each kind instead of 5 (0.817, 0.824, 0.797), 4 controls (0.818, 0.826, 0.809: within 0.002
  and ahead on one set) and 5 points (0.809, 0.807, 0.796). The expert is the cheapest of the four designs
  compared, not a proven optimum.
- Wording. The edge-wells rule names both mechanisms (fewer cells plated, and water lost), and spec() says which
  plates and wells are charged and how many tips each protocol takes (23 against 18 per compound in the reference
  design).

The second pass also added these, for the unmatched_dmso row. They were removed in the third pass (below):

- A dose range from 100 to 300 uM, where the top well held top dose / 100 percent DMSO (MAX_DMSO_PCT = 3%), with the
  reward's ceiling left at 100 uM (CEILING_UM).
- A design field untreated_controls (cells in medium, no DMSO), at least 3 such wells on every plate, and a gate
  term VEHICLE_HEALTH_GATE = 0.7 on (vehicle mean - full-kill mean) / (untreated mean - full-kill mean). The expert
  got 5 untreated wells per plate, the controls were laid out vehicle, full-kill, untreated in turn, and the random
  policy drew 300 uM top doses and untreated wells.
- Their card: unmatched_dmso still inverted on seeds 0 to 199 (gap -0.004, expert wins 37.5%). The trap is the
  expert with one change, so it still dosed to 100 uM at 1% DMSO, and the health term could only ever fail a
  matched run, because an unmatched run's vehicle wells hold no DMSO.

After the re-review (verdict: fix first), the third pass. Each item is again a change after the first card. All
were made before the cards below were run, and nothing was adjusted afterwards.

- The solvent-matched vehicle is now a QC convention. A plate passes its gate only if its vehicle wells carry the
  DMSO percentage of the top dose on that plate, read from the plate map (the nominal dmso_pct of the vehicle wells
  equals that of the highest-dose wells). A run with match_dmso false leaves its vehicle wells without DMSO, so it
  cannot pass the gate, by rule. This is what a cytotoxicity guideline asks of a solvent control; it is stated here
  as general practice and no clause is quoted. The term has no number to tune. The match-the-dmso rule's "why" in
  rules.toml says the same.
- So the unmatched_dmso row is punished by the gate alone, as no_kill_control is. Its verdict is not a simulated
  effect. Nominal rig unless a rig is named. "With the term switched off" means envs.cytotox_window.vehicle_is_matched
  patched in process to always pass, with Z', controls and everything else as coded:

  | seeds | as coded: expert, trap | gap, expert wins, verdict | With the term switched off: trap | gap, expert wins, verdict |
  |---|---|---|---|---|
  | 0 to 199 (the card) | 0.825, 0.399 | +0.426, 100.0%, visible | 0.830 | -0.005, 41.0%, blind |
  | 40000 to 40199 (fresh card) | 0.810, 0.390 | +0.420, 100.0%, visible | 0.812 | -0.002, 43.0%, blind |
  | 100000 to 103999 | 0.828, 0.400 | +0.428 +- 0.001, 99.95%, visible on 20 of 20 windows of 200 | 0.8315 against 0.8281 | -0.0034 +- 0.0007, 41.4%, blind on 12 windows and inverted on 8 |
  | randomized rigs, 100000 to 100999 | 0.826, 0.400 | +0.426, 100.0%, visible on 5 of 5 | 0.8315 against 0.8260 | -0.0055 +- 0.0016, 40.4%, blind on 2 windows and inverted on 3 |
  | harsh rig, 100000 to 100999 | 0.768, 0.373 | +0.395, 99.9%, visible on 5 of 5 | 0.7772 against 0.7682 | -0.0090 +- 0.0023, 38.1%, inverted (blind on 2 windows, inverted on 3) |

  As coded, the trap came out ahead on 2 of the 4,000 nominal seeds and on 1 of the 1,000 harsh ones: seeds where
  the expert's own plate failed Z', so both runs were halved.
- Why the simulated data cannot carry this row. No design can put more than 1% DMSO in a well. Seeds 0 to 199, with
  the term switched off: the top well of an inert compound reads 0.89 of a DMSO-free vehicle (1.01 of a matched
  one), 14 of 383 clean soluble compounds are called a measured CC50 against 2 when matched, the dip offsets the
  rig's adsorption bias (measured CC50s above 20 uM read 0.035 log units high matched and 0.017 unmatched), and the
  route takes 30 fewer tips (cost 126.8 against 129.8). Quality of the trap 0.862 against 0.857 for the expert
  (0.8635 against 0.8609 on the 4,000 seeds). Across 2,000 random designs each run both ways (seeds 110000 to
  111999), matched minus unmatched is +0.088 +- 0.003 as coded (+0.156 and ahead on 85.6% among the 1,135 designs
  with at least 3 controls of each kind) and -0.005 +- 0.001 with the term switched off (ahead on 22.6%).
- What that means. The reward now teaches "match the DMSO" as a rule to follow, stated in spec(), not as an effect
  to be found in the data. An agent that reads only luminescence has no way to discover it. Read "visible" on this
  row as "the convention fires". A reward that scored a loss of 10 to 30% at the top dose could see the mismatch
  itself; none was built or measured.
- The vehicle-health term, the untreated wells and the 100 to 300 uM range were removed, because with the
  matched-vehicle term in place they were measured to do no job. The measurements were made on the second-pass
  code, each episode run once and its gate evaluated again in process with and without each term:
  - On the card. For the expert and each of the seven traps, the health term alone failed 0 of 200 episodes, on
    seeds 0 to 199 and on 40000 to 40199. Taking it out changed no row but random (0.108 to 0.126, and 0.094 to
    0.105).
  - At or under 100 uM it read evaporation, not DMSO. Random designs with at least 3 wells of every control (seeds
    52000 to 54999): of 901 matched designs that passed Z' the health term failed 134, and 130 of those had no lid
    and outer wells in use. It also failed 62 of 922 unmatched designs, all of them with no lid and outer wells.
  - Above 100 uM nothing beats 100 uM with or without it. Reference design, seeds 50000 to 50199, every bound
    capped at 100 uM in the report; the 100 uM expert scores 0.837. Matched at 150, 200 and 300 uM: 0.805, 0.805
    and 0.737 without the health term (the 100 uM expert wins 69.0%, 82.5% and 84.5% of seeds), 0.484, 0.385 and
    0.357 with it. Unmatched: 0.379, 0.305 and 0.291 either way, because it fails the matched-vehicle term. So
    matched is ahead of unmatched at every dose without the health term too, and the range has no upside because
    the ceiling is 100 uM: the term only turned a small loss into a large one where no design has a reason to go.
  - The re-review also showed the term was leaky: at 10,000 cells per well or more both the vehicle and the
    untreated wells reach confluence and a halved vehicle reads healthy. Its note that the text here held only at
    4,000 cells applies to sentences that went with the term.
- What the removal changed. VEHICLE_HEALTH_GATE, MAX_DMSO_PCT, the untreated_controls field and the requirement of
  3 untreated wells are gone. TOP_LIMIT_UM is 100 uM again, and STANDARD_DMSO_PCT, STANDARD_TOP_UM and CEILING_UM
  are back under their first-card names MATCHED_DMSO_PCT and TOP_LIMIT_UM. The simulator at or under 100 uM was not
  touched: the final code reproduces the second-pass code run with no untreated wells exactly (360 episodes checked
  seed by seed, largest difference 0). What follows from it:
  - The expert changed again: its 5 untreated wells per plate are gone, so its controls are those of the first
    card, 5 vehicle and 5 full-kill. Four compounds now share the first plate and two the second (three and three
    before). Cost 129.8 against 139.8; mean reward 0.825 against 0.826 on seeds 0 to 199.
  - The control layout is vehicle and full-kill alternating again, and the random policy draws the first card's
    choices again (top dose 10, 30 or 100 uM, no untreated wells).
  - edge_wells moved back from visible to weak: 91.0% of seeds in the second pass, 78.5% now (82.5% on the fresh
    seeds, 77.0% on randomized rigs). The second-pass verdict came from the three-kind control layout, which put
    the trap's vehicle wells partly on the outer ring and failed Z' on 62% of its runs. Now its gate passes on
    100% of runs and what is left is the edge physics: quality 0.792 against 0.857. On the harsh rig it is visible
    (99.0%, gate passes 8%).

What the cards show after the third pass (200 seeds, nominal rig, seeds 0 to 199; seeds 40000 to 40199 in brackets):

- Expert 0.825 (0.810). Null -0.500 (-0.500). Random 0.130 (0.100).
- Visible: unmatched_dmso 100.0% (100.0%) and no_kill_control 100.0% (100.0%), both by the gate alone;
  nominal_above_solubility 93.2% (93.5%), mostly by the ceiling convention (below); exact_not_bound 100.0%
  (100.0%), scored against the truth.
- Weak: short_exposure 80.2% (85.0%), edge_wells 78.5% (82.5%), no_interference_control 67.0% (68.5%). These are
  the three mistakes that lose through simulated data, and all three are punished on average, not per episode.
- No row is blind or inverted. With the matched-vehicle term switched off unmatched_dmso would be blind on both
  cards and inverted on 8 of 20 further windows.
- Ten further windows of 200 seeds (100000 to 101999): the four visible rows are visible on all ten
  (nominal_above_solubility 91.2 to 95.5%, near the line), short_exposure and edge_wells are weak on all ten, and
  no_interference_control is weak on nine and blind on one (58.0 to 69.5%).
- Randomized rigs, seeds 0 to 199: expert 0.812, the same four rows visible (unmatched_dmso 99.0%, no_kill_control
  99.5%, nominal_above_solubility 93.8%, exact_not_bound 100.0%), the same three weak. Harsh rig: expert 0.770,
  edge_wells visible as well.

Unchanged conventions, with numbers measured again on seeds 0 to 199 (seeds 40000 to 40199 in brackets):

- A bound is scored against a ceiling, the lower of 100 uM and the solubility limit. That is a scoring convention
  for "what could honestly be claimed", not a measured quantity. Most of what the nominal_above_solubility trap
  loses comes from it: of its 437 (443) bounds above the ceiling, 38 (53) were false.
- The no_kill_control trap is punished by the QC gate alone (quality 0.852 against 0.857). With a background of
  100 counts against about 28,000, leaving out the zero reference barely moves the normalised viability.
- The 0.5 log tolerance, the -0.5 floor, the 0.5 cap and equal weight for a false alarm and a missed liability are
  judgments. A safety team may well weigh a missed liability more.
- Not on the card, because they are not registered traps: 20,000 cells per well instead of 4,000 loses 0.22 and
  the expert wins 93.5% (93.8%) of seeds, because confluent wells hide toxicity; no lid on inner wells loses 0.10
  (0.09) and the expert wins 63.0% (60.5%), part of it through Z', which passes 77% (82%) of those runs.
- A token run that only passes the gate, followed by "at least 100 uM" for every compound, scores 0.229 (0.203),
  against -0.5 for the same answer with no run and 0.825 (0.810) for the expert.
