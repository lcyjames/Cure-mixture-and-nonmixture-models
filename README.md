# Weibull-Cure-mixture-and-nonmixture-models
This package performs parameter estimation for the cure mixture and nonmixture models, with Weibull survival models for the latency component.
The cure mixture model assumes a logistic model for the incidence and Weibull proportional hazards model for the latency part, and an EM algorithm is used for the estimation.
The non-mixture model assumes the complementrary log-log model for the incidence, and the latency part is purely a Weibull distribution, and a regular likelihood maximization is used for the estimation.
