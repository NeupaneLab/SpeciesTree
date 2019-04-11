# master_seed  = 997758
master_seed  = 99554543
run_on_prior = True
simulate_morpho_data = False


data_size = 100 # number of simulation replicates for each taxa
define_sets = [['E_pulverulenta'],['E_myrtilloides'], ['E_myrtoidea'], ['E_rubra'], ['E_bifida'], ['E_resinosa']] # This list cane changed with corresponding changes on morpho_means and morpho_sigmas  

morpho_means = [6.0, 9.0, 16.0, 12.0, 11.0, 4.0] #2 sets defined in this case with 2 means, This list should correspond to the number of sets within define_sets 
morpho_sigma = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1] #2 sets defined in this case with 2 sigmas, This list should correspond to the number of sets within define_sets 

n_gen = 1000000
save_every =  2000

mean_prior_alpha = 50.0 # mean of the normal distribution prior for mean
mean_prior_sd = 10. # mean of the normal distribution prior for mean


mean_sigma =2. # mean of the lognormal distribution prior for sigma
sig_sd = 0.5 #sd of the lognormal distribution prior for sigma



start_means = [5.]*6
start_sigmas = [0.6]*6


brown_mean_prior_alpha =17. # means of the borwn
brown_sig_sd = 10. #sd borwn

sim_root_brow = [21.0]*6


start_root_brow = [4.0] #
rate_multiplier = 1. #rate_multiplier for vcv



start_set = [[18], [37], [60], [73], [99], [114]]
random_start_tree =  "(((E_pulverulenta:1.2,E_myrtilloides:1.2):1.15,(E_myrtoidea:0.05,E_rubra:0.05):2.3):1.25,(E_bifida:0.2,E_resinosa:0.2):3.4);"
data_file = pd.read_csv("morphometrics2.csv")
