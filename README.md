# Bifurcations
AUTO scripts to generate bifurcation diagrams for various ODE systems. AUTO can be downloaded from

    http://indy.cs.concordia.ca/auto/

and the documentation can be found at 

    http://www.dam.brown.edu/people/sandsted/auto/auto07p.pdf

The four models implemented so far are the Goldbeter, Gonzalez-Fernandez & Ermentrout, Meyer & Stryer, and NVU based SMC/EC models. 

To run
------

AUTO must be installed before running. Once everything is set up, simply run

    auto <<script>>
    
in the terminal, where script has the .auto file extension. The models are defined in the .f90 files, the relevant AUTO parameters in the c.* files, and the additional parameters in the autorc files. Having a look at the AUTO documentation linked above is highly recommended.
