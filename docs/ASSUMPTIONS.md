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
