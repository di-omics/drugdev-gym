# Learning and transfer: dose_response

Cross-entropy method over protocol choices, reward only. Test seeds are unseen by the learner.

## Learning curve (population mean reward)

| iteration | trained on nominal | trained on randomized |
|---|---|---|
| 1 | 0.116 | 0.114 |
| 2 | 0.140 | 0.118 |
| 3 | 0.223 | 0.202 |
| 4 | 0.257 | 0.283 |
| 5 | 0.288 | 0.335 |
| 6 | 0.378 | 0.403 |
| 7 | 0.423 | 0.410 |
| 8 | 0.514 | 0.455 |
| 9 | 0.482 | 0.456 |
| 10 | 0.580 | 0.497 |
| 11 | 0.475 | 0.535 |
| 12 | 0.593 | 0.563 |

## Transfer (mean reward over 100 unseen seeds per rig)

| design | nominal rig | randomized rigs | harsh rig (held out) |
|---|---|---|---|
| learned_on_nominal | 0.723 | 0.710 | 0.652 |
| learned_on_randomized | 0.706 | 0.674 | 0.531 |
| expert | 0.699 | 0.676 | 0.530 |

## Learned designs

- learned_on_nominal: `{"top_conc_nM": 10000.0, "dilution_factor": 10.0, "n_points": 8, "replicates": 1, "transfer_ul": 25.0, "mix": true, "new_tips": true, "avoid_edges": false, "lid": false, "standards": true, "controls_per_side": 8}`
- learned_on_randomized: `{"top_conc_nM": 10000.0, "dilution_factor": 10.0, "n_points": 8, "replicates": 1, "transfer_ul": 25.0, "mix": true, "new_tips": true, "avoid_edges": true, "lid": false, "standards": true, "controls_per_side": 4}`
