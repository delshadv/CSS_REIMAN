#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 19 22:49:45 2022

@author: dv01
"""
#import argparse
#import importlib
#from copy import deepcopy
#from logging import warning

import mne
import h5io
import numpy as np
import pandas as pd
from sklearn.linear_model import RidgeCV
from sklearn.dummy import DummyRegressor
from sklearn.pipeline import make_pipeline
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler, FunctionTransformer
from sklearn.model_selection import KFold, GridSearchCV, cross_validate
from sklearn.metrics import make_scorer, r2_score, mean_absolute_error
import coffeine

from scipy import io
import scipy.misc 
from scipy import io

frequency_bands = {
            "delta": (2, 4),
            "theta": (4.0, 8.0),
            "alpha": (8.0, 12.0),
            "beta_low": (12.0, 30.0),
            "beta_mid": (30.0, 48.0),
            "beta_high": (52.0, 87)
        }


y = io.loadmat("/imaging/henson/users/dv01/y.mat")
y = np.array(y['y'])
rank = 65

data = io.loadmat("/imaging/henson/users/dv01/MEGPLANAR.mat")

covs = np.array(data['covariance_rei_wa_fb'])

A0 = covs[0,0]
A1 = covs[0,1]
A2 = covs[0,2]
A3 = covs[0,3]
A4 = covs[0,4]
A5 = covs[0,5]    

A = np.full((307,204,204,6), 0)

for ii in range(307):
    A[ii,:,:,0] = A0[:,:,ii]

for ii in range(307):
    A[ii,:,:,1] = A1[:,:,ii]

for ii in range(307):
    A[ii,:,:,2] = A2[:,:,ii]

for ii in range(307):
    A[ii,:,:,3] = A3[:,:,ii]
    
for ii in range(307):
    A[ii,:,:,4] = A4[:,:,ii]
    
for ii in range(307):
    A[ii,:,:,5] = A5[:,:,ii]

X = pd.DataFrame(
    {band: list(A[:, ii]) for ii, band in
            enumerate(frequency_bands)})

filter_bank_transformer = coffeine.make_filter_bank_transformer(
    names=list(frequency_bands),
    method='riemann',
    projection_params=dict(scale='auto', n_compo=rank)
    )
model = make_pipeline(
    filter_bank_transformer, StandardScaler())

#cv_params = dict(n_splits=10, shuffle=True, random_state=42)
#cv = KFold(**cv_params)


print("Running cross validation ...")
scores = cross_validate(
    model, X, y, cv=10, scoring = 'balanced_accuracy',
    n_jobs=12)  # XXX too big for joblib
print("... done.")

scores['test_score']

