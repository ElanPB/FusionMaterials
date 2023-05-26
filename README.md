# FusionMaterials

## Project Description
This project will focus on the procurement of relevant cross sections from ENDF and EXFOR for the effects of high energy neutrons (14 MeV all the way down) on essential materials for Magnetic Confinement Fusion with a focus on materials found in HTS magnets.


The current list of important elements include (all elements found in REBCO tape):
- REBCO: Y, Ba, Cu, O
- Hastelloy: Ni, Cr, Mo
- Substrate: Ag
- Buffer layer: ???

The relevant cross sections for all isotopes of the elements listed above with an abundance greater than 10\% will be looked at

Relevant cross sections include:
- For gas production (any neutron induced reaction with an ejectile of either p or alpha)
  - (n,alpha)
  - (n,n+p)
  - (n,n+alpha)
  - (n,p+alpha)
- For displacement
  - Elastic
  - Inelastic
  - (n,gamma)

## Product
- Plots of data for all the isotopes listed above on all cross sections listed above
- Chi-Squared analysis between ENDF and EXFOR data determining "goodness" of fit for the evaluation

## Files
__RawData__: All raw data from ENDF and EXFOR for all the isotopes looked at in this project<\br>
__ProcessedData__: All plots and info generated from the raw data. Subfolders for each isotope<\br>
__Code__: Code base for this project to create plots and analyze raw data
