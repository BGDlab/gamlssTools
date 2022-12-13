

```
Error in if (dv > olddv && itn >= 2 && auto == TRUE) { :
missing value where TRUE/FALSE needed
```

Full context:
```
> gamGMV_fix <-gamlss(formula = GMV ~ fp(log_age, 3) + sex + fs_version + study,
+                 sigma.formula = GMV ~ fp(log_age, 2) + sex + study,
+                 nu.formula = GMV ~ 1,
+                 control = gamlss.control(n.cyc = 200),
+                 family = GG, data=na.omit(cn_pheno_sub), trace = FALSE)
GAMLSS-RS iteration 1: Global Deviance = 1513682
GAMLSS-RS iteration 2: Global Deviance = 1513487
Error in if (dv > olddv && itn >= 2 && auto == TRUE) { :
  missing value where TRUE/FALSE needed
```

- The code runs if the nu and sigma formulas are removed.
