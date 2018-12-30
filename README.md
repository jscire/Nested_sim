The code provided here allows for the stochastic simulation of a nested model of spread of an infectious disease. Two strains of the pathogen can spread: a _wild-type_ sensitive to treatment strain and a _mutated_ resistant one. This code was used to produce results shown in the article 
_Aggressive or moderate drug therapy for infectious diseases? Trade-offs between different treatment goals at the individual and population levels_ by Jérémie Scire, Nathanaël Hozé, Hildegard Uecker (_https://www.biorxiv.org/content/early/2018/04/24/307223_).

The class _NestedSimulator_ can be used to simulate an epidemic and record the number of transmissions of each type, as well as the disease burden.

The class _WithinHostParameterSensitivityAnalysis_ can be used to perform simulations of the within-host model and record the recovery rate or the probability of emergence of the resistance using custom parameters.
