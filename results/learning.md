# Learning and transfer: dose_response

Cross-entropy method over protocol choices, reward only. Test seeds are unseen by the learner.

## Learning curve (population mean reward)

| iteration | trained on nominal | trained on randomized |
|---|---|---|
| 1 | 0.096 | 0.089 |
| 2 | 0.130 | 0.091 |
| 3 | 0.191 | 0.178 |
| 4 | 0.238 | 0.245 |
| 5 | 0.316 | 0.266 |
| 6 | 0.407 | 0.357 |
| 7 | 0.394 | 0.403 |
| 8 | 0.556 | 0.374 |
| 9 | 0.471 | 0.448 |
| 10 | 0.594 | 0.495 |
| 11 | 0.530 | 0.441 |
| 12 | 0.587 | 0.527 |

## Transfer (mean reward over 100 unseen seeds per rig)

| design | nominal rig | randomized rigs | harsh rig (held out) |
|---|---|---|---|
| learned_on_nominal | 0.723 | 0.710 | 0.652 |
| learned_on_randomized | 0.750 | 0.729 | 0.673 |
| expert | 0.699 | 0.676 | 0.520 |

## Learned designs

- learned_on_nominal: `{"top_conc_nM": 10000.0, "dilution_factor": 10.0, "n_points": 8, "replicates": 1, "transfer_ul": 25.0, "mix": true, "new_tips": true, "avoid_edges": false, "lid": false, "standards": true, "controls_per_side": 8}`
- learned_on_randomized: `{"top_conc_nM": 30000.0, "dilution_factor": 10.0, "n_points": 8, "replicates": 1, "transfer_ul": 25.0, "mix": true, "new_tips": true, "avoid_edges": false, "lid": false, "standards": true, "controls_per_side": 8}`
