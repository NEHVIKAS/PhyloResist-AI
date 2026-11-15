# Geographic Holdout Validation

Exploratory leave-one-state-out validation was run for the two fungicides. Each outer test fold contains all isolates from one collection state; inner alpha selection also holds out state groups.

## Group coverage

- LA: n=1
- MI: n=1
- MN: n=35
- MT: n=8
- ND: n=33
- NY: n=73
- VA: n=1

Source population labels: {'NA1': 147, 'unassigned': 4, 'NA2': 1}

Population labels are highly imbalanced (the panel is predominantly NA1), so geographic state is used as the only viable documented group field here. Geographic holdout is an internal transfer scenario, not an independent population validation. Tiny held-out state groups make individual fold metrics unstable.

Both traits are analyzed as supplied: least-squares mean EC50 in mM, without log transformation. See each fungicide subfolder for state-fold metrics and isolate predictions.

## Pooled results

- propiconazole: Ridge RMSE 0.015103 mM; training-mean baseline RMSE 0.014513 mM; Ridge R² -0.089.
- tebuconazole: Ridge RMSE 0.007850 mM; training-mean baseline RMSE 0.007300 mM; Ridge R² -0.158.

In this run, pooled Ridge RMSE is slightly higher than the training-mean baseline for both fungicides, so the genotype model did not improve pooled prediction error. The selected alpha was 100 in every outer fold (the largest value tested); this suggests regularization tuning reached the current grid boundary and should be revisited only if a later analysis is justified. These exploratory results do not support a surveillance performance claim.
