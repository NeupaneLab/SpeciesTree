# master_seed  = 997758
master_seed  = 995555
run_on_prior = False

simulate_morpho_data = True
simulate_means_from_root_mean = False # simulate means for the tips under Brownian process, if chosen False, a values in the list 'morpho_means' will be chosen

data_size = 100 # number of simulation replicates for each taxa
define_sets = [['E_pulverulenta'],['E_myrtilloides'], ['E_myrtoidea'], ['E_rubra'], ['E_bifida'], ['E_resinosa']] # This list can be changed with corresponding changes on morpho_means and morpho_sigmas  

sim_root_brow = [21.0]*6 #simulating means under brownian process
morpho_means = [6.0, 9.0, 16.0, 12.0, 11.0, 4.0] # use this if means are not simulated under brownian process by specifying 'simulate_root_mean = False' option
morpho_sigma = [2.]*6 #sigmas to generate data from means generated from sim_root_brow or user assigned morphs means

n_gen = 1000000
save_every =  1000

###priors
mean_prior_alpha = 0.0 # mean of the normal distribution prior for mean
mean_prior_sd = 50. # mean of the normal distribution prior for mean


mean_sigma =2. # mean of the lognormal distribution prior for sigma
sig_sd = 1. #sd of the lognormal distribution prior for sigma


brown_mean_prior_alpha =0. # mean of the brownian process 
brown_sig_sd = 50. #sd brownian process 


##starting values 
start_means = [2.]*6
start_sigmas = [0.6]*6


start_root_brow = [4.0] #
rate_multiplier =1. #rate_multiplier for vcv, higher value will increase variance 


start_set = [[18], [37], [60], [73], [99], [114]]
random_start_tree =  "(((E_pulverulenta:1.2,E_myrtilloides:1.2):1.15,(E_myrtoidea:0.05,E_rubra:0.05):2.3):1.25,(E_bifida:0.2,E_resinosa:0.2):3.4);"
#data_file = pd.read_csv("morphometrics2.csv")
