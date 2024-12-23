# EFPBSA-Beta-Prior
Bayesian Penalized Cox-PH using Beta Prior

## To demonstrate the simulation scenarios. Kindly use the following files: 
* Simulation_Gamma_GL_MDPI_1 (.R): This is for the choice Group Lasso
* Simulation_Gamma_EN_MDPI_1 (.R): This is for the choice of Elastic Net
- Each can implement both Gamma or Beta processes prior for the cumulative hazard ratio.

## To demonstrate the real-life datasets. Kindly use the following files: 
* REAL_life_nki70 (.R) 
* REAL_life_ DutchBC(.R)
* Real_Life_DLBCL_chop (.R)
- Here, you can also choose your choice of prior for  the cumulative hazard ratio and select your choice on penalty function (Group Lasso and Elastic Net) 
## Major references
*Lee KH, Chakraborty S, Reeder H, Sun J (2024). _psbcGroup: Penalized Parametric and Semiparametric Bayesian Survival Models with Shrinkage and
  Grouping Priors_. R package version 1.7,
  <https://CRAN.R-project.org/package=psbcGroup>.
* Lee, Kyu Ha, Sounak Chakraborty, and Jianguo Sun. 2015. “Survival Prediction and Variable Selection with Simultaneous Shrinkage and Grouping Priors.” Stat. Anal. Data Min. 8 (2): 114–27. https://doi.org/10.1002/sam.11266.

* 2021. “psbcGroup: Penalized Parametric and Semiparametric Bayesian Survival Models with Shrinkage and Grouping Priors.” R Package Version 1.5. https://CRAN.R-project.org/package=psbcGroup (23 January 2023, date last accessed).

* Zhi Zhao, John Zobolas, Manuela Zucknick, Tero Aittokallio, Tutorial on survival modeling with applications to omics data, Bioinformatics, Volume 40, Issue 3, March 2024, btae132, https://doi.org/10.1093/bioinformatics/btae132

## Suplementary references
* Harrel Jr, F. E. and Lee, K. L. and Mark, D. B. (1996) "Tutorial in biostatistics: multivariable prognostic models: issues in developing models, evaluating assumptions and adequacy, and measuring and reducing error", Statistics in Medicine, 15, pages 361–387.

*Pencina, M. J. and D'Agostino, R. B. (2004) "Overall C as a measure of discrimination in survival analysis: model specific population value and confidence interval estimation", Statistics in Medicine, 23, pages 2109–2123, 2004.
* Heagerty, P.J., Lumley, T., Pepe, M. S. (2000) Time-dependent ROC Curves for Censored Survival Data and a Diagnostic Marker Biometrics, 56, 337 – 344
