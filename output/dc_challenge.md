---
title: presentatie data challenge
author: trond
---

# Imputation 1:
1. Linear model with region- and year fixed effects, estimated with OLS
2. Define outliers (Cook's distance)
3. Re-estimate model without outliers (training set)
4. Predict outliers (test data) using model

# Model testing 1: cross validating a time-series model
* Goal: out-of-sample validation of an i-period ahead forecast
* Training data set: series up to t
* Test data set: t+i

# Model testing 2: example 
* 4 different State Space models
* Calculate fitness measure
* Choose model based on out-of-sample fit
