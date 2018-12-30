The code provided here allows for the stochastic simulation of a nested model of spread of an infectious disease. Two strains of the pathogen can spread: a _wild-type_ sensitive to treatment strain and a _mutated_ resistant one.

The class _NestedSimulator_ can be used to simulate an epidemic and record the number of transmissions of each type, as well as the disease burden.

The class _WithinHostParameterSensitivityAnalysis_ can be used to perform simulations of the within-host model and record the recovery rate or the probability of emergence of the resistance using custom parameters.
