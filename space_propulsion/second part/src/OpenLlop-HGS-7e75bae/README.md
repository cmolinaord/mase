HGS V1.3
========

Code files
----------

| Function | Description
|----------|------------
| hgseq | Finds chemical equilibrium composition of a mixture at given T P
| hgsisentropic | Isentropic expansion, frozen or shifting eq
| hgsprop | Mixture properties
| hgssingle | Properties of a single species
| hgsTp | Adiabatic flame temperature (with dissociation)

| Example | Description
|---------|------------
| 00 | Properties of single elements and mixtures
| 01 | CH4 adiabatic temperature combustion without dissociation
| 02 | Adiabatic H2O2 decomposition
| 02b | V2 Turbopump example (liquid H2O2 adiabatic decomposition)
| 03 | CO2 dissociation equilibrium solved with K
| 04 | H2O equilibrium dissociation for different values of temperature (T)
| 04b | Adiabatic flame temperature with dissociation this code is an example of how hgsTp works
| 05 | Adiabatic H2 / O2 reaction
| 05b | Adiabatic H2 / O2 reaction using hgsTp
| 05c | Adiabatic C3H8 reaction
| 06 | LOX - LH2 combustion
| 07 | isentropic comparison with RPA software
| 08 | ISP vs OF ratio
| 09 | Comparison between hgs / rpa using RP-1
| 10 | Solvers

Internal files:
* BurcatDB.mat: database
* hgsDB
* hgsid
* hgsmix
* hgsfzero
* hgssolve
