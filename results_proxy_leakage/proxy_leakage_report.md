# Proxy Leakage and Bias Audit

## Direct Formula Features

The current heuristic targets are explicitly constructed from these input features. High model performance against proxy targets is therefore expected and should not be interpreted as biological validation.

- chemical: grantham, blosum62, delta_hydropathy, delta_volume, delta_charge, delta_polarity_flag
- conservation: conservation_wt_freq, position_entropy
- structural: af2_neighbor_count_8A, af2_is_buried_proxy, is_globular_domain_region, af2_plddt, delta_volume

## Circularity Assessment

The `mutation_severity_score`, `destabilization_proxy_score`, and derived binary labels are circular with the features listed above. Models trained to predict these targets measure how well an algorithm can reconstruct the engineered formula, not whether it predicts experimental pathogenicity or protein stability.

## Ablation Summary

```csv
feature_group,target_type,feature_count,accuracy_mean,accuracy_std,f1_mean,f1_std,roc_auc_mean,roc_auc_std,average_precision_mean,average_precision_std,mae_mean,mae_std,rmse_mean,rmse_std,r2_mean,r2_std
biochemical_only,high_severity_classification,17,0.787110633727175,0.013644205807261394,0.6149432884368858,0.016684729976095446,0.8642037455977885,0.012903122766977528,0.6677406491297891,0.02673163646345324,,,,,,
biochemical_only,severity_score_regression,17,,,,,,,,,4.508489499331377,0.43592535895795914,6.192270468005727,0.8483261810741104,0.523836036287514,0.06120406042059927
conservation_only,high_severity_classification,5,0.3720730397422126,0.027431782814083597,0.40703588327058055,0.01941888550930457,0.5382031640952748,0.02509060240513182,0.2696258908492501,0.01583997079774438,,,,,,
conservation_only,severity_score_regression,5,,,,,,,,,7.1114652736795545,0.4330299478615455,8.850362480444002,0.6211473233667554,0.020715696320533095,0.05351531376728613
structural_only,high_severity_classification,17,0.704189044038668,0.01391547560962546,0.5476565406647321,0.019829672806758583,0.8038975718320339,0.010425709637423142,0.5656525626392075,0.05206564113177345,,,,,,
structural_only,severity_score_regression,17,,,,,,,,,6.155030453118426,0.35187479519275827,7.7748737634976,0.7142478683428839,0.24392085127259605,0.08004344937953933
combined_all_numeric,high_severity_classification,41,0.9445757250268528,0.008391847469226492,0.8869310406334456,0.01933592480268888,0.9878561667679133,0.002543003938098521,0.9670821993418558,0.0069092691788755385,,,,,,
combined_all_numeric,severity_score_regression,41,,,,,,,,,1.2272097104806492,0.2672484903300999,2.148600447808624,0.9198473093684458,0.9383879820436256,0.045434749717914294
combined_no_formula_features,high_severity_classification,29,0.8781954887218045,0.010135446712855394,0.7489702457350013,0.013698559524745411,0.9329359633804607,0.007653769327530125,0.8313331021216935,0.018370389649910134,,,,,,
combined_no_formula_features,severity_score_regression,29,,,,,,,,,2.867420281436412,0.4910343197443481,4.190752505819739,0.9253015772217599,0.7801304804563125,0.06807885360302757
```


## Interpretation

- If a single feature group performs very well, the proxy target is dominated by that group.
- If `combined_no_formula_features` drops sharply, the original performance was mostly formula leakage.
- The next scientifically grounded step is to replace these targets with external ΔΔG or curated variant labels.
