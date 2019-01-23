

The code was written in 2016 and first published in my PhD thesis (2017) and is available in McGill University, e-thesis repository.
"Social, genetic and behavioral risk factors of neck cancer - elucidating causal pathways"
http://digitool.Library.McGill.CA:80/R/-?func=dbin-jump-full&object_id=148325&silo_library=GEN01


# 4-Way-decomposition-
Stata command for 4 way decomposition , binary x binary scenario

                           Stata codes for mediation and 4-way decomposition analysis under
                             counterfactual causal framework for case-control study design

                                    Thekke Purakkal, Akhil Soman*, Kaufman, Jay S**.

               *Division of Oral Health and Society, Faculty of Dentistry, McGill University, Montreal Quebec, 
               **Department of Epidemiology, Biostatistics and Occupation Health, McGill University, Montreal, Quebec

                                  Release Date: August 1, 2016; Current Version 1

Stata macros (e.g., PARAMED) for estimating mediation effects in the presence of exposure mediator
interaction under counterfactual causal framework already exist (1). For conducting 4-way decomposition
analysis, although the mathematical equations and SAS codes have been provided by VanderWeele 2014
(2), Stata codes have not been written. The below Stata codes were exclusively written for conducting
mediation analysis under exposure-mediator interaction (alternative method using mathematical equations),
as well as 4-way decomposition analysis for this thesis work. Although the codes given here is specific for
binary outcome, mediator and exposure variables, the code can be easily extended to the case where the
exposure, mediator and outcome are continuous, categorical, binary or their combinations. The codes have
been written using the mathematical formulas for total effect, mediation effects, and 4-way decomposition
effects, as well as for calculating various proportions provided by VanderWeele 2015 and 2016 (3, 4). The
user is encouraged to cross check these codes with the formulas in these references. The codes for
bootstrapping procedure given at the end of the codes can be used to derive the confidence limits for these
estimates.

Let Y be a binary outcome. In a case-control study, Y=1 may represent cases and Y=0 may represent
controls or non-cases. Let A be a binary exposure, and M a binary mediator. Let C1, C2, be continues
covariates, and C3, C4 be binary or categorical. If there are more or fewer covariates, one can add or remove
scalars under “//Covariates”, “//Assigning levels of covariates” and “//calculating bcc” in the below given
code. For a case-control study with rare disease outcome, the line of code for mediator model may be fit
only among controls. Alternatively, or, if the outcome is not rare, one can weight the mediator model using
sampling weights as suggested by VanderWeele and Vansteelandt 2010 (5).

References

1. VanderWeele TJ. Mediation: Introduction and regresion -based approaches. Explanation in Causal
inference: Methods for mediation and Interaction. USA: Oxford University Press; 2015. p. 40-1.
2. VanderWeele TJ. A unification of mediation and interaction: a 4-way decomposition. Epidemiology
(Cambridge, Mass). 2014;25(5):749-61.
3. VanderWeele TJ. A unification of mediation and interaction. Explanation in Causal inference:
Methods for mediation and interaction. USA: Oxford University Press; 2015. p. 371-96.
4. Erratum: A Unification of Mediation and Interaction: A 4-Way Decomposition. Epidemiology
(Cambridge, Mass). 2016;27(5):e36.
5. Vanderweele TJ, Vansteelandt S. Odds ratios for mediation analysis for a dichotomous outcome.
American journal of epidemiology. 2010;172(12):1339-48.
