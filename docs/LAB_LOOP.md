# Lab in the loop: what the lab should measure

The question: when a real lab is added to these environments, which data should it generate
first, so that a reward a model trains on here is true? Working answer below. It is a plan
built from this repository's own tables, not a result. Nothing here has touched hardware.

[DATA.md](DATA.md) lists the data types. This page says which measurements produce them
most cheaply. The table behind it is generated: `drugdev-gym campaign` writes
[results/campaign.md](../results/campaign.md).

The kinds below are drafts. A first draft of this page went through an automated review,
with one reviewer prompted to read as a bench scientist. It moved 14 of the 62 pairs, 10 on
bench grounds and 4 on what the trap actually sets, and it replaced the first sizing. What it
changed is noted where it applies. The expert has not read any of it yet.

## 3 loops on 3 clocks

1. **A model climbs an environment.** Thousands of attempts. The simulator answers in
   milliseconds; a model's attempt costs seconds and from a few thousand to tens of thousands
   of tokens. The lab cannot sit in this loop: a plate takes hours.
2. **An environment climbs its QC card.** Each rule gets a trap, and the card says whether
   the reward can teach it. Minutes on a laptop.
3. **The simulator climbs the lab.** Real measurements replace assumed numbers and show
   which rules matter on a real rig. Days to weeks.

The lab's place is loop 3. Its output is small and has to be chosen so that every
measurement changes what loops 1 and 2 can be trusted for.

## What a mistake changes decides what data prices it

Each trap is the expert protocol with one mistake. Sorted by the kind of choice the mistake
changes, the 62 pairs that count fall into 5 groups:

| group | the mistake changes | what prices it | pairs |
|---|---|---|---|
| from complete records of correct runs | how the same raw data is read, or which data is on the run (standards, time points, a plate, the internal standard, the time zero) | a complete raw record of a correct run, read again or with the same data left out | 15 |
| from grids | a setting of the assay (cells, time, stimulant, exposure), or which concentrations are tested | a grid of conditions repeated over days; a dense grid of concentrations where every well is dosed on its own | 14 |
| on the rig | how the liquid is moved or the sample treated, or where on the plate the wells sit | the physical quantity it changes, measured directly; uniformity plates; then a few same-day paired runs | 14 |
| only when something is wrong | run QC: controls, replicates, an activity control | the same item left out of a run that carries the fault it exists for | 6 |
| not a plate question | a reporting convention, an assay mistake made inside a campaign, campaign strategy, a judgment on finished data | a rule, the assay-level size of the error, a whole campaign, or review cases | 13 |

What each group depends on:

- **Records, 15.** 6 hold on a plain complete record: the raw-ratio, standard-curve and
  curve-top mistakes in the potency environments, and leaving the off-target plates out.
  9 need the record to keep something specific. A time course needs one well per time point
  or a much larger incubation, with the real quench time of every sample, because a record
  denser than one design takes more volume than one incubation well holds. The
  internal-standard pairs hold for one injection order only. The time-zero pair needs both
  zeros on the run. The limit of quantification needs blanks and low calibrators. Solubility
  needs a turbidity read or its own measurement. The enzyme-decay correction needs the
  activity of each fluid lot followed over the incubation, a series of its own; the simulator
  hands that number over as a constant.
- **Grids, 14.** 10 are a setting of the assay. A setting chosen for one day's cells is the
  answer for that day, and the price of not titrating is how far the best setting moves
  between days and cell batches. So a condition grid is repeated, not measured once. 4 are a
  choice of concentrations. 2 of those are the top well of the expert's own ladder and hold
  on a plain record today. The other 2 are a different range, which on a serial ladder is a
  second ladder. It becomes a lookup only where every well is dosed on its own, as on an
  acoustic rig.
- **On the rig, 14.** 10 are handling. The first draft said only a run done wrong on purpose
  prices these. The review's point: most of them change a physical quantity that can be
  measured directly and more sharply, with no cells. A dye ladder with and without mixing,
  and at a small and a large transfer volume. For reused tips, the ladder wells of a sticky
  reference compound read by LC-MS, since a dye does not stick to a tip. Plate weights
  before and after an incubation. Each gastric sample split into the acid and the base
  quench, and injected again after a hold. A control substrate in fresh and aged enzyme. A
  DMSO tolerance curve on the cells. The cost to a potency or a half-life then follows by
  refitting, and a paired run is the confirmation. 4 are position, an effect inside a plate:
  uniformity plates and a split layout price it for the plate-reader assays, and a stable
  marker compound in inner and outer wells of the incubation plate for the stability one.
- **Only when something is wrong, 6.** A correct run has nothing for a control to catch.
  Leaving controls out of a good plate only widens the estimate of its window. The value of
  run QC is how often a fault happens times what an uncaught fault costs, so each control is
  priced on a run that carries the fault it exists for. A replicate catches an unmixed spike
  and an activity control catches aged enzyme, and the handling runs above carry both. Plate
  controls catch a cell, reader or reagent fault. A ladder mistake does not touch the control
  wells, so those need a run with such a fault planted, and the lab's routine run log for
  how often it happens.
- **Not a plate question, 13.** 3 are conventions about what may be claimed where there is no
  data, such as a bound past the top dose. 4 are an assay mistake made inside a campaign. Its
  size comes from the same mistake in its own assay environment. What it costs the campaign
  depends on whether the error changes the ranking of the analogs, which only a campaign, or
  a simulation fed the measured error, can show. 2 are campaign strategy. 4 are judgment,
  priced by review cases, and every run done wrong on purpose is one: a real run with a known
  planted cause.

Where the simulated reward is weakest today: it teaches 1 of the 10 condition pairs outright
and none of the 4 position pairs. The handling pairs look healthy, 8 of 10. But 7 of those 8
rest on physics that was assumed, and 1, unmatched DMSO, on a rule about the plate map. So a
real measurement there checks the size of an effect the simulator already has. Condition and
position are where it may be missing the effect altogether.

What the simulator checks, by 4 rules. A reading or a convention asks the rig for exactly
what the expert asks for; the check cannot tell those 2 apart, so that split is the author's
call. A judgment is made where nothing is run. An upstream trap changes the same assay field,
to the same value, as the assay trap it names. Any other trap changes a design field of its
own kind, or drops a run of that kind, with every other field it changes named and
explained. That is a check on the simulator's actions. Whether the named data would measure
the cost on a bench is a judgment, and `NEEDS` in `campaign.py` records what the review said
each pair depends on.

## What to measure first

A draft order. Step 0 needs no lab. Steps 1 to 5 follow the automated bench review.

0. **Run current models on the environments.** No model has been run here. 2 numbers per
   environment decide what is worth a plate: the mean reward (a pass rate, once a pass rule
   is set), and for each rule how often a model makes that mistake. A mistake no model makes,
   perhaps because the help text of the design field gives the right choice away, is not
   worth pricing. The harness is built (`drugdev-gym baseline`) and checked with the scripted
   policies. It reads 51 of the 62 mistakes from what a model does. Running it on a real
   model has not been done: it spends money and needs a model key.
1. **Liquid-handler metrology.** No cells. The dye ladder and a weight check price the mixing
   and transfer-volume mistakes and fill the rig profile, in about a day. Tip reuse needs
   the sticky reference compound and LC-MS on top.
2. **The plate uniformity study** that validation asks for, extended with a split inner and
   outer layout and plate weights, over 3 days. It prices the position pairs of the
   plate-reader assays and gives the window and the edge and drift maps.
3. **One designed stability record** on LC-MS: one well per time point, real quench times,
   both time zeros, a time zero with and without the mix, the internal standard, a control
   substrate in each fluid, an aged enzyme lot with its age and storage stated, each gastric
   sample split into the acid and the base quench, the quenched plate injected again after a
   hold, blanks and low calibrators, a fast and a slow reference, and a marker compound in
   inner and outer wells. The activity of each fluid lot over the incubation is a second,
   smaller series. Together they price most of the stability pairs. The injection count is
   not worked out here.
4. **The expert protocol with reference compounds on 6 or more days,** with 2 runs done
   right on the same day at least once. The scatter between days sets the reward tolerance.
   The scatter between 2 same-day runs sets the size of every paired comparison. The same
   records price the reading and data mistakes of the potency environments.
5. **A condition grid on 3 or more days,** with a reference curve at each point. The spread
   of the best setting between days is the price of not titrating.
6. **Then runs done wrong on purpose,** paired with the run done right on the same day, for
   what steps 1 to 5 leave open, and a few runs with a planted cell, reader or reagent fault
   for the plate controls. All of them are also review cases.

If each of the 10 handling mistakes were also confirmed by paired runs, the assumed numbers
give about 8 same-day pairs per arm: 120 runs with the 5 reference arms, 4,832 plate wells
and 5,248 LC-MS samples. That is a starting size, not a bound. The 8 is for a shift of 0.25
log units in the mean when one run scatters by 0.15 within a day, by a paired t-test. In
`dose_response` as simulated, one expert potency scatters by about 0.13 within a day at the
median. Reusing tips shifts the mean by about 0.1, which 8 pairs would not resolve, and a
starved ladder returns no fitted value for most compounds, so a shift in the mean is not
defined for it. 10 comparisons read at 0.05 each will also show a false difference somewhere
more often than that, so the reading rule has to be set before the runs. An independent pair
means a different day, cell passage and fresh dilution; several plates from one cell
preparation are one run.

## The unit of data is a run record

A summary value (an EC50, a half-life) can price nothing: it has already made every choice
these environments are about. What one returned run has to hold:

- **The plan as declared.** The design, field for field, as an agent would state it.
- **The run as executed.** The liquid-handler log or transfer report (source, destination,
  volume, tip, mix, time, exceptions), incubation times and temperatures, incubator and
  stack position, the lid, the reader method, minutes from detection reagent to read, and
  the order and real time of every quench and every LC-MS injection, with the temperature
  of the quenched plate and the autosampler tray.
- **The raw reads.** Every read of every well, with the plate map and each well's role. For
  LC-MS, the peak area of the analyte and of the internal standard for every injection, the
  blanks and calibrators, and the vendor raw files.
- **What everything was.** Cell line, passage, viability and density at harvest. Lots of the
  kit, enzymes, matrix, serum and plates. Each compound stock's lot, age, freeze-thaw count
  and solvent, and the DMSO percent in each well. The instrument methods. Operator and date.
  Without these a shift between runs cannot be told from a mistake.
- **The anchors.** Reference compounds on every run. Their known answer is that lab's own
  consensus over several days, with the published value as a sanity check, because a
  published potency depends on the cell line and the kit. The reference set has to carry the
  liabilities the rules are about: a sticky compound, a fast and a slow one, one that
  interferes with the readout, one that is poorly soluble.
- **Any deviation,** planted or noticed, kept apart from the data. The log itself shows a
  missing mix step, so a blind review needs a defined view of the record with the log held
  back.
- **Its use,** assigned before it is run: calibration, or exam. The split is by day and by
  compound, not by plate, because plates from one day share their day.
- **Later, the expert's call** on it, made blind: accept or repeat, the cause, and why.

## Who can run what

- **Complete records and grids are close to work a lab already does,** but the record is a
  custom deliverable. A service lab returns a potency or a percent remaining by default;
  per-transfer logs and raw files have to be agreed, and one pilot record should be checked
  against the list above before a campaign is ordered.
- **Everything measured on the rig belongs to that rig.** The dye ladder, the edge and drift
  maps and the runs done wrong on purpose are specific to one handler, its tips, the plate
  type, the incubator and the reader. Measured at an outside lab, they calibrate that lab's
  rig, not the one a model will later drive. The runs done wrong on purpose also need a lab
  that will follow a deviation sheet.
- **No project compounds should be needed,** provided public compounds with these
  liabilities exist for the target. The question is what a mistake does to a known answer.

## What assay validation already asks for

The NIH Assay Guidance Manual asks a new assay for 2 studies before it is used: a plate
uniformity study over 3 days (2 for a transferred assay), and a replicate-experiment study of
20 to 30 compounds in each of 2 runs. It accepts an assay whose minimum significant ratio
is below 3 and whose plates reach a signal window of at least 2 or a Z-prime of at least 0.4
([HTS Assay Validation](https://www.ncbi.nlm.nih.gov/books/NBK83783/)).

So step 2 extends a study a validating lab runs anyway, and step 4 goes beyond the manual's
2-run study. By the manual's preliminary estimate, MSR = 10^(2 x sqrt(2) x S) with S the
scatter of a log potency between runs, an assay accepted at the limit of 3 scatters by up to
0.17 log units between runs. A 2-run estimate is a floor for the scatter over a campaign of
weeks, which should be measured on the same reference compound over 6 or more runs. That
0.17 is between runs; the 0.15 assumed for the paired runs is within a day, and neither is
measured here. Validation does not include runs done wrong on purpose.

These environments gate a plate at a Z-prime of 0.5. That is stricter than the manual's 0.4
and was chosen by hand.

## What the data becomes

1. **Dials.** Instrument QC and the measured size of each mistake replace assumed defaults
   in [ASSUMPTIONS.md](ASSUMPTIONS.md) (`RigProfile.from_qc` is the first joint). The fit
   should be a spread of rig profiles, with its width from day-to-day scatter, not one point.
2. **A tolerance.** The scatter of the reference compounds across days should set the reward
   tolerance in place of 0.5 chosen by hand.
3. **An exam, not a training set.** A grid or a record is a finite set of answers, and a
   learner trained on it can memorise it. Training stays on the simulator. The measurements
   held back from the fit are the exam, and an exam needs rules of its own: coded compound
   names, truth taken from that rig's consensus and not from the literature, compounds and
   days that were never in the fit, records kept private, and more than one grid so the exam
   has a variance.
4. **Review cases from real runs.** Each run done wrong on purpose is a `plate_review` or
   `timecourse_review` case with real data and a known cause, and the expert's blind calls
   on them are the human reference.
5. **A check that can fail.** Rank agreement between simulator and rig over the scripted
   policies is not one. The nominal and randomized simulated rigs already agree at a rank
   correlation of 0.90 to 1.00 in every environment, because most rungs are easy to tell
   apart. The check that can fail compares, pair by pair, the gap measured on the rig with
   the simulated gap, with an interval, on the pairs that are near ties, and it has to be
   shown to reject a simulator known to be wrong.
6. **A place for what an optimizer finds.** Every arm here is the expert with one scripted
   change, but the failures that held environments back were combinations a learner found. A
   few rig runs should be kept for designs that beat the expert in simulation.

## What rollouts at scale ask of an environment

Post-training on an environment is mostly inference: the model acts, the environment answers,
the reward is logged, over many tasks and many attempts per task. What that asks for, and
where this repository stands:

- **Cheap and reproducible.** Tasks are seeded and plain JSON, and the same task id and
  actions always give the same reward. 50 expert episodes in each of the 10 environments
  take about 16 seconds on a laptop. One instance is one episode; `reset()` does not rewind.
- **Truth out of reach.** Not yet. The package is open and nothing in a task is secret, so a
  search over candidate seeds, about a millisecond each in `dose_response`, recovers the
  hidden truth from a single observation. A real exam needs a private salt held by the
  serving side, with the simulator and its verifier running where the agent cannot read them.
- **A stated suite.** Not yet. There is no manifest of training and exam seed ranges,
  simulator version, reward range or pass rule, and the noise of a task is frozen to its
  seed, so a second attempt at a task sees the same draw.
- **Measured difficulty.** Not yet. `drugdev-gym baseline` runs a model through every
  environment and reports its reward and its mistakes, but no model has been run.
- **A package.** Not built. [OpenEnv](https://github.com/meta-pytorch/OpenEnv) is reset,
  step and state calls with one container image per environment.
  [Harbor](https://huggingface.co/docs/trl/en/harbor) is one directory per task: an
  instruction, an environment image, and tests that write the reward. The environments take
  reset and step with JSON in and out, so the interface is close, but in either packaging
  the simulator has to sit outside the agent's sandbox.

## Not covered

- **Product truth.** A plate also yields a measured value for each compound. How predictive
  models are trained on those values is a different problem and is out of scope here. This
  page is about process truth: how a run was done and how it was read.
- **The acoustic rig.** The kinds are assigned for the 96-well, tip-based environments on
  main. On a 1536 acoustic rig the handling list changes (tips drop out; backfill,
  source-plate handling and the spin come in) and the range pairs become a lookup.
- **Sizes and money.** Only the paired runs are counted, in wells, samples and runs. Steps 1
  to 5 are not sized here, and nothing is priced.
- **One question for the expert.** Whether the `stability` arm with no lid on the outer
  wells is still wanted, now that the lid alone is ruled moot.
