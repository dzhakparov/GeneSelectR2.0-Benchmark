#!/usr/bin/env python3
"""Stabl stability ranking on one saved outer-training division."""
import argparse
import csv
import time
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
from stabl.stabl import Stabl, group_bootstrap

p = argparse.ArgumentParser()
for key in ('x','y','groups','output'):
    p.add_argument('--'+key, required=True)
p.add_argument('--seed', type=int, required=True)
a = p.parse_args()
x = pd.read_csv(a.x)
y = pd.read_csv(a.y).iloc[:,0].to_numpy(dtype=int)
groups = None if a.groups == 'none' else pd.read_csv(a.groups).iloc[:,0].to_numpy()
assert x.shape[1] == 2000 and len(y) == len(x) and len(np.unique(y)) == 2
assert np.isfinite(x.to_numpy()).all()
params = dict(n_bootstraps=150,
              lambda_grid={'C':np.linspace(0.01,1.0,10)},
              artificial_type='random_permutation',
              artificial_proportion=0.5,
              fdr_threshold_range=np.arange(0.1,1.0,0.01),
              sample_fraction=0.5, replace=False,
              n_jobs=1, random_state=a.seed, verbose=0)
if groups is not None:
    params['bootstrap_func'] = group_bootstrap
start = time.monotonic()
with threadpool_limits(limits=1):
    model = Stabl(**params).fit(x,y,groups=groups)
score = model.stabl_scores_.max(axis=1)
rank = np.argsort(-score, kind='stable')
assert len(rank) == 2000 and len(np.unique(rank)) == 2000
with open(a.output,'w',newline='') as f:
    w=csv.writer(f); w.writerow(['gene','max_selection_frequency'])
    w.writerows((x.columns[i],float(score[i])) for i in rank)
with open(a.output.replace('.csv','_meta.csv'),'w',newline='') as f:
    w=csv.writer(f)
    w.writerow(['seed','native_selected','positive_scores','seconds','min_fdr','fdr_threshold'])
    w.writerow([a.seed,int(model.get_support().sum()),int(np.sum(score>0)),
                round(time.monotonic()-start,3),float(model.min_fdr_),
                float(model.fdr_min_threshold_)])
