#this module tests the pystra module
import pystra as ra
import numpy as np

#define limit state eqn
limit_state = ra.LimitState(lambda Ms, Mw, Md, My: 1- ((Ms+Mw+Md)/My))

#initialize stochastic model
stochastic_model = ra.StochasticModel()

#define random variables

## Ms is a constant value determined from SDI analysis in GHS
stochastic_model.addVariable(ra.Constant("Ms", 3006))

## Mw is a Gumbel distributed random variable with a mean/nominal ratio of 1 and a COV of 0.15, where nominal value is 27975 from LR rules
Mw_r = 1
Mw_cov = 0.15
Mw_nom = 27975
stochastic_model.addVariable(ra.Gumbel("Mw", Mw_nom*Mw_r, Mw_cov*Mw_nom*Mw_r))

## Md is a Gumbel distributed random variable with a mean/nominal ratio of 1 and a COV of 0.25, where nominal value is 13903 from LR rules
Md_r = 1
Md_cov = 0.25
Md_nom = 13903
stochastic_model.addVariable(ra.Gumbel("Md", Md_nom*Md_r, Md_cov*Md_nom*Md_r))

## My is a Lognormal distributed random variable with a mean/nominal ratio of 1 and a COV of 0.15, where nominal value is provided by design optimization
My_r = 1
My_cov = 0.15
My_nom = 69423.2348766629
stochastic_model.addVariable(ra.Lognormal("My", My_nom*My_r, My_cov*My_nom*My_r))

#initialize reliability analysis
options = ra.AnalysisOptions()
options.setPrintOutput(True)

Analysis = ra.Form(
    analysis_options=options,
    stochastic_model=stochastic_model,
    limit_state=limit_state
)

Analysis.run()

beta = Analysis.getBeta()
failure = Analysis.getFailure()

print(f"Beta is {beta}, corresponding to a failure probability of {failure}")