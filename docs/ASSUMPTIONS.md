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
