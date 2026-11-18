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
- Costs are round numbers that only set the relative price of wells, tips, samples, and synthesis.
- The four original environments share their first hidden draws on a given seed (the review
  environments do not). Do not reuse a seed across environments where one's revealed truth could
  inform another.

## What the rewards do and do not claim

- Truth is known by construction, so every reward is exactly verifiable. That is the whole value of
  simulation here, and also its limit: a policy that scores well has learned this simulator.
- "Expert" means a scripted policy that follows standard bench practice. It is a reference point,
  not an upper bound. The experts of the four original environments are not tuned. The review
  judges are: their thresholds and their list of serious causes were set on seeds 0 to 399, so
  they are reported on other seeds.
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
