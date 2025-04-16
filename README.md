# NKBB

NKBB is a model for the thermal spectrum of thin accretion disks in parametric black hole spacetimes. It employs the Novikov-Thorne model for the description of the accretion disk and the transfer function approach proposed by Cunningham for storing all the relativistic effects of the spacetime metric. Please follow the instructions below carefully in order to ensure a properly working version of the model. This model is designed and tested to work within the X-ray spectral fitting software XSPEC.   
    
If you are using the NKBB model in your work please cite the following papers:    
* Zhou et al., _XSPEC model for testing the Kerr black hole hypothesis using the continuum-fitting method_, [Phys. Rev. D 99, 104031 (2019)](https://doi.org/10.1103/PhysRevD.99.104031)
* Zhou et al., _Thermal spectra of thin accretion discs of finite thickness around Kerr black holes_, [MNRAS 496, 497–503 (2020)](https://doi.org/10.1093/mnras/staa1591)

## Getting started

1. Installing the model:    
    * First, open the file "nkbb.h". You need to change the variable 
    "TR_TABLE_PATH" (line 22) to your current working directory.
    * To create the model, there is the script compile_NKBB.sh. The command lines are
         
            chmod u+r compile_NKBB.sh
            ./compile_NKBB.sh

2. Loading the model in XSPEC:
    * From the directory where the model has been installed, the model can
    simply be loaded inside XSPEC by executing:
  
             lmod nkbb .

3. Calling the model in XSPEC:
    * Execute following in XSPEC to call this model:
     
            model nkbb
    * Please note that NKBB is a model of thermal component spectrum, so the
    norm should be fixed to 1.
    
4. If you want to use the NKBB with variable disk thickness in Kerr spacetime, 
    please set the parameter "defpar_type" to be "0" and change the lower-limit 
    of the parameter "defpar_value" to be "0.0" as well. The scaled "defpar_value"
    is proportional to the mass accretion rate in unit of the Eddington limit.    
    * defpar_value = 0.0 => mass_accretion_rate =  0.0% Eddington limit    
    * defpar_value = 1.0 => mass_accretion_rate = 30.0% Eddington limit   
   
5. For questions and bug reports, please contact <relxill_nk@fudan.edu.cn>

## Usage instructions

To ensure the optimal performance of the model, specific FITS files are mandatory. These files encompass the transfer functions and disk temperature profile of a specific non-Kerr metric and/or accretion disk geometry, which can currently be obtained upon request. A download link will be made available in the future.

The eta parameter regulates the location of the inner edge of the accretion disk and its value cannot be negative. It is defined by the relation $R_{\rm in} = \left( 1 + {\rm eta} \right) R_{\rm ISCO}$, where $R_{\rm in}$ is the radial coordinate of the inner edge of the accretion disk and $R_{\rm ISCO}$ is the radial coordinate of the ISCO.

The def_par_type parameter determines the accretion disk geometry and spacetime metric. It is crucial that this parameter remains fixed during data analysis. The following table outlines the def_par_type values and the corresponding model configurations:

The defpar_type parameter can be set to the following values:

1. defpar_type = 0 corresponds to the Kerr metric with a finite thickness accretion disk, as proposed by Taylor and Reynolds. In this case, the defpar_value parameter is used to specify the mass accretion rate in Eddington units, and the user should fix a lower limit to 0. The range of defpar_value is scaled between 0.0 and 1.0, which corresponds to a mass accretion rate of 0% to 30% of the Eddington limit.    
 
2. defpar_type = 1 corresponds to the Johanssen metric with a non-zero value of α13 specified by the defpar_value parameter.    
 
3. defpar_type = 2 corresponds to the Johanssen metric with a non-zero value of α22 specified by the defpar_value parameter.    
 
4. defpar_type = 3 corresponds to the Johanssen metric with a non-zero value of \epsilon_3 specified by the defpar_value parameter.    
    
*Important Note:* For cases 2, 3, and 4, the defpar_value parameter is scaled within the range of [-1, 1]. The unscale.py and unscale_batch.py [python scripts](https://github.com/ABHModels/relxill_nk/tree/main/unscale_v1.2.4), provided in the relxill_nk package, can be utilized to obtain the actual (unscaled) values of the deformation parameters following the fitting process.
