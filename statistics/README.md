# statistics
This program produces the following files with the described content:
 * av_p:	avalanche size distrubion for events with positive deformation value
 * av_n:	avalanche size distrubion for events with negative deformation value
 * av_poz_part:	avalanche size distrubion where in the size only the positive events are counted
 * av_neg_part:	avalanche size distrubion where in the size only the positive events are counted
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