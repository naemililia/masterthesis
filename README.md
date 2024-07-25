# masterthesis
code corresponding to my master thesis.

In general the code is designed to work for both datasets. However, due to not wanting to adapt paths every time and to keep the code from both datasets, i created duplicates which included minimal changes to adapt for the two datasources EuroSAT and Sen1floods11. Sources for the code used is in the notebooks itself. 

Explanation of what each file does:

- pipeline_labels.ipynb : create the labels (BUILT & LC) for the entire Sen1floods11 weakly labeled dataset
- pipeline_labels_eurosat.ipynb : create the labels (BUILT & LC) for the EuroSAT dataset
- label_creation.ipynb : a step by step creation of a single label (was used to develop the 'method', good to understand what is done when)
- alignment_helpers.py : helper functions for label creation
- Baseline_Sen1floods11.ipynb : train the baseline semantic segmentation model on senfloods hand labeled data
- EuroSAT_baseline.ipynb : train the baseline scene classification model on EuroSAT 'baseline dataset' (the 1/3 that I use)
- train_contrastive_learning_prep.ipynb : pretrain the model via contrastive learning on senfloods dataset
- train_contrastive_learning_eurosat.ipynb : pretrain the model via contrastive learning on eurosat dataset
- train_pretext_1_label.ipynb : train pretext task with one label type on the senfloods dataset
- train_pretext_1_label_eurosat.ipynb : train pretext task with one label type on eurosat dataset
- finetune.ipynb : finetune the pretext and contrastive learning pretrained models on the senfloods hand labeled data. also includes evaluation part (which also exists separate)
- finetune_and_eval_eurosat.ipynb : finetne and evaluate the pretext and contrastive learning pretrained models on the eurosat data
- evaluate.ipynb : evaluate finetuned models on senfloods hand labeled data (you have to reload your test data into memory, makes sense if you want to evaluate multiple models after another, does not make sense if it is just for one model you just finetuned, then use the evaluation part included in the finetune model)


