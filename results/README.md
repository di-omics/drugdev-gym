# results

Every other file in this folder is written by a command. None of them is edited by hand: to
change a number, change the code and run the command again. All of it is simulated.

Run the commands from the repository root, in this order. `knowledge` reads `envqc.json` and
`learning.json`, and `campaign` reads `envqc.json`.

| Step | Command | Writes |
|---|---|---|
| 1 | `drugdev-gym qc` | `cards/<env>.md` for all 10 environments, and `envqc.json` |
| 2 | `drugdev-gym qc --env plate_review --env timecourse_review --first-seed 10000` | `cards/plate_review.md` and `cards/timecourse_review.md` again, and their 2 entries in `envqc.json` |
| 3 | `drugdev-gym qc --rig randomized` | `cards/<env>_randomized.md` for all 10 environments, and `envqc_randomized.json` |
| 4 | `drugdev-gym qc --rig randomized --env plate_review --env timecourse_review --first-seed 10000` | `cards/plate_review_randomized.md` and `cards/timecourse_review_randomized.md` again, and their 2 entries in `envqc_randomized.json` |
| 5 | `drugdev-gym learn` | `learning.md`, `learning.json` |
| 6 | `drugdev-gym knowledge` | `knowledge.md`, `knowledge.json` |
| 7 | `drugdev-gym campaign` | `campaign.md`, `campaign.json` |

Steps 2 and 4 are there because the scripted reviewers' thresholds were set on seeds 0 to
399. Their cards are reported on seeds 10000 to 10199 instead. A run with `--env` rewrites
its own cards and keeps the other entries of the JSON file.

## The 28 generated files

| File | What it is | Step |
|---|---|---|
| `cards/<env>.md` (10 files) | one QC card per environment on the nominal rig: every scripted policy of that environment (the expert, each trap, and a random or null policy where it has one) on the same 200 seeds | 1, 2 |
| `cards/<env>_randomized.md` (10 files) | the same cards on randomized rigs | 3, 4 |
| `envqc.json` | every nominal card as data | 1, 2 |
| `envqc_randomized.json` | every randomized card as data | 3, 4 |
| `learning.md`, `learning.json` | the cross-entropy learner on `dose_response`: its learning curve, its 2 designs and the transfer table | 5 |
| `knowledge.md`, `knowledge.json` | each expert rule next to what the reward teaches about it | 6 |
| `campaign.md`, `campaign.json` | each rule-trap pair sorted by the real data that would price it, and the sizing of paired runs | 7 |

The 10 environments are `assay_window`, `dose_response`, `stability`, `lead_opt`,
`plate_review`, `timecourse_review`, `gi_fluid_stability`, `gi_camp_potency`,
`selectivity_panel` and `cytotox_window`.

## To check them

Give each command `--out DIR`, and `campaign` also `--cards DIR/envqc.json`, to write into
an empty folder instead of this one. On the machine that wrote these files, the 7 steps
reproduce all 28 byte for byte. Another machine may differ in the last digits of a number.

## The defaults behind those commands

- `qc`: 200 seeds per policy starting at seed 0, the nominal rig, output in `results/`.
- `learn`: 12 iterations, learner seed 0, 100 unseen test seeds per rig.
- `campaign`: 3 seeds to compare each trap with the expert, scatter 0.15 and shift 0.25 in
  log10 units, alpha 0.05, power 0.8. With any other number it refuses to write here and
  asks for `--out`.

## Not committed

- `results/human/` is where `drugdev-gym score` writes a human expert's scores. Git ignores
  it. No human reference is committed.
- `results/baseline/` is where `drugdev-gym baseline` writes a model run: a report and the
  full transcripts. No model run is committed.
