# avalanche2

## About this file
This is the readme file for my project hosted at [danieltuzes/avalanhe2](https://github.com/danieltuzes/avalanche2). This file is manually created by me and contains some description of this project.
## About this repo
Here you can find my simulations programs and evaluation programs to analyize the data you get. I wrote these during my PhD years.
### names: no avalanche1?
This project is called avalanche2, but there is no avalanche1 or avalanche, it has been removed long ago. avalanche2 is the CA model written to study dislocation avalanches under mean field approximation and avalanche3 is the CA model written to study dislocation patterns.
## so many projects, folders and main files
Here you can find all the folders listed with my programs therein. Each folder contains a program for specific purposes.

1. arithmetic: to make some basic operations on block_data objects stored in files, output is stored in files.
2. **avalanche2**: The dislocation avalanche simulation program. It is based on the continuum model of dislocations in mean field approximation.
3. correlation_integral: to calculate the correlation integrals of the initial point of dislocation avalanches
4. statistics: calculates the
	* av_p:	avalanche size distrubion for events with positive deformation value
	* av_s_av:	average avalanche size in each external stress intervals
	* iDG:	the size of the deformation of the first ith avalanche
	* iDG_av:	the sum and STD of the deformation of the first ith avalanche over different realizations; the ith line gives the property of the ith avalanche
	* ies:	the external stress at the ith avalanche. The results are stored both in incremental corresponding seed_start value and both incremental external stress values
	* ies_av:	the sum and STD of the deformation of the first ith avalanche over different realizations; the ith line gives the property of the ith avalanche
	* log:	beside the stdout, processing info goes there too. Contains the parameters the evaluation program used, the read in files.
	* ssc_DsG:	stress strain curve averaged at strain values (delta sum gamma)
	* ssc_DsGdev:	standard deviation of the stress values at DsG values over different realizations
	* ssc_tau:	stress strain curve averaged at external stress values
	* tau_g_sampl_DsG:	statistics of the lines of tau_g files sampled at strain values
	* tau_g_sampl_tau:	statistics of the lines of tau_g files sampled at external stress values