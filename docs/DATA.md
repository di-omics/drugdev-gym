# What data makes these environments informative

The question: which data types, and which ways of generating data, does a post-training RL
environment need before a reward for potency or stability says something true about drug
discovery? Working answer below. It is a plan, not a result.

## An environment needs three things from data

1. **Truth to score against.** A potency or half-life whose value and uncertainty are known.
2. **An observation model.** What the agent sees when it acts: raw wells and raw peak areas,
   with the noise and artifacts a real run has. Not summary IC50s.
3. **A validated reward.** Evidence that the reward separates good practice from bad (envqc).

Simulation gives 1 and 3 for free and fakes 2. Real data is needed mostly to fix 2.

## Data types, by brick

**Potency (dose_response, assay_window)**

- Raw well-level reads with the plate map: both channels, not just the ratio.
- Controls and the on-plate standard curve for every plate.
- The same reference compounds repeated across days, plates, and operators. This gives the
  assay's real reproducibility (minimum significant ratio), which sets the reward tolerance.
- Uniformity plates (all negative or all positive controls) for edge, row, and drift maps.
- The liquid-handler log for the run: volumes, tip use, mix steps, timestamps.

**Stability (stability)**

- Raw peak areas per sample with the internal standard area, not percent remaining.
- Time zero handled as its own measurement, with replicates.
- A positive control compound with a known half-life in every run, and the matrix lot.
- Injection order and time on the autosampler.

**Both, and the part almost nobody publishes**

- **Failure-mode runs.** The same protocol run wrong on purpose: no mix, reused tips, open
  edges, a starved ladder, no internal standard. These are what calibrate the traps, and the
  traps are what make a reward worth training on.

## Ways to generate it, cheapest first

| level | modality | gives | cost |
|---|---|---|---|
| 0 | Mechanistic simulation (this repo) | unlimited seeded tasks, exact truth | none |
| 1 | Instrument QC on the real rig: dye or gravimetric volume checks, uniformity plates | the RigProfile numbers | about a day of robot time, no compounds |
| 2 | Reference replays: public tool compounds with literature potency on a commercial kit, run with the expert protocol and each trap | the sim-to-real gap per brick (sim2real.gap), tolerance from real reproducibility | a few plates |
| 3 | Public data for the truth distribution: potency and stability value ranges, matched pairs for how rough SAR is | realistic hidden landscapes | none |
| 4 | Real-in-the-loop episodes: replay recorded runs (offline), or queue a real plate per episode | a held-out real test set | slow, expensive, small |

The useful ratio is a lot of level 0, a little of levels 1 and 2 chosen to calibrate it, and
level 4 kept as the exam rather than the training set. RL needs far more attempts than a lab
can run live (a plate takes hours to days), so the lab's job is to calibrate and audit the
simulator, not to be the reward loop.

## What "informative" should mean

- The reward tolerance matches real assay reproducibility. A reward tighter than the assay
  trains agents to chase noise.
- Every reward ships with its envqc card, including the mistakes it is blind to.
- A policy trained in simulation is only called transferable after it is scored on level 2 or 4 data.

## Not done yet

Levels 1, 2, and 4 need hardware and reagents. Level 3 needs a pass over public sources
(verify licences and what raw data each one actually exposes before relying on it).
