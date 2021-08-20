---
title: "LOCF vs. multiple imputation for addressing sparse measurements"
author: "Ehsan Karim and Belal Hossain"
date: "19 August, 2021"
output:
  html_document:
    keep_md: yes
---



# impute Sparse

Analysis codes, saved results and figures associated with the following article

**Can LOCF be a reasonable approach for addressing sparse measurement issues? An illustration from per-protocol analysis in pragmatic trials**

by Mohammad Ehsanul Karim and Md. Belal Hossain

## Simulation scenarios

- K = 60 # max. follow-up
- MI.m = 10 # number of imputed data


Table: Simulation Scenarios considered (scenario number will match figure numbers).

| scenario|    n| sigma|U          |Role          |Effect |
|--------:|----:|-----:|:----------|:-------------|:------|
|       51| 1000|   2.0|unmeasured |Confounder    |RD     |
|       52| 1000|   5.0|unmeasured |Confounder    |RD     |
|       53| 1000|   0.5|unmeasured |Confounder    |RD     |
|       54| 1000|   3.0|unmeasured |Confounder    |RD     |
|       61| 1000|   2.0|unmeasured |risk factor Y |RD     |
|       71| 1000|   2.0|unmeasured |Confounder    |OR     |
|       81| 1000|   2.0|measured   |Confounder    |RD     |
|       91| 2000|   2.0|unmeasured |Confounder    |RD     |
|       92|  250|   2.0|unmeasured |Confounder    |RD     |

Scenario 51 had additional 2 versions for different MI methods

- Confounder = $\beta1_0=0, \beta2_0=-5, \beta1_1=6, \beta2_1=3$							
- Risk factor for Y = $\beta1_0=2.95, \beta2_0=-5, \beta1_1=0, \beta2_1=0$							
