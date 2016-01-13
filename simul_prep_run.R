#########################
#
#
# simul_prep_run for simulation work
#
# 1/23/2016
#
# written by Mason Fidino



# source the functions to run everything
source("simul_prep_functions.R")


# load packages, download if necessary
package_load(c("plyr", "reshape2", "rjags",
               "runjags", "mcmcplots"))


# simulate the data, putting the means and standard
# deviations of each parameter group

sim_df <- make_sim_df(nspec = 3, nsite = 80, nyear = 6, nrep = 28,
                      gam = .4, phi = .7, p = 0.15, gam_sd = 1, beta = .5,
                      beta_sd = 1, beta2 = .5, beta_sd2 = .5,  inxs = 0, 
                      inxs_sd = 2, phi_sd = 1, p_sd = 0.5, add_NA = FALSE, 
                      percent_to_NA = 0.2,
                      row_replicate = 400)


# this function will simulate all the needed data to fit both models
# and output a list of lists. Each list contains another list with 4 elements:

# [[1]] sim_list: Contains known values to try to estimate with model
# [[2]] mats: All the different matrices. z, j_list, y_list, and zinit (intial z values)
# [[3]] data_list: data for JAGs
# [[4]] comparison_list: not currently used. Set up for weekly detection instead
#                        of daily detection.

#Note: takes ~0.5 seconds to run for each row of sim_df
#Note: known_inxs fills col wise, give number of inxs equal
#      to size of the entire inxs matrix (we write over the diagonal)

# If test = TRUE, function will only simulate two random rows of sim_df
# and return them.
all_sim <- simulate_from_sim_df(sim_df, test = FALSE,
                                known_inxs = c(1,-0.5,1,-1,1,1,-1.5,-2,1))

# parameters to follow.
params <- c("gam",
            "phi",
            "psi_in", 
            "beta",
            "beta2",
            "p",
            "z",
            "sigma_int")

# load glm module for block updating of parameters (leads to less auto-correlation)

load.module("glm")

# function to batch analyze simulations. Uses runjags under the hood
# so look to ?run.jags for many of the arguments listed in this function.
# for those not in run.jags... Change the number of chains if you have
# less cores!

# model = name / location of interaction model to fit to data.

# comparison_model = name / location of different model you want to fit to
#                    simulated data

# m_name = what you want to have the posterior named. posterior for model will be
#         named "post_matrix_m_name_(simulation number)". Known values will be named
#         "known_m_name_(simulation number)". Posterior for comparison model will be
#         named "post_comp_matrix_m_name_(simulation number)".

# save_path = where you want to save all of these outputs. do not put a "/"
#             at the end of the folder you want to place them in. For example,
#             "C:/simulations" would work, but "C:/simulations/" would not. 

### Note: 400 simulations, with comparions, took roughly 6 days to run.

model <- "interaction_model.R"
comparison_model <- "beta_model.R"


# for predators with competing prey
start_time <- Sys.time()
batch_analyze(all_sim = all_sim, params = params, n_chains = 7,adapt_steps = 3000,
              burn_in = 10000, sample_steps = 10000, thin_steps = 10,
              make_comparisons = TRUE, model = model,
              comparison_model = comparison_model, m_name = "pcp",
              save_path = "C:/simulations")
end_time <- Sys.time()

pcp_time <- end_time - start_time

# for horizontal competition
all_sim <- simulate_from_sim_df(sim_df, test = FALSE,
                                known_inxs = c(1,-1,-2,-0.5,1,-1,0, -0.7,1))

start_time <- Sys.time()
batch_analyze(all_sim = all_sim, params = params, n_chains = 7,adapt_steps = 3000,
              burn_in = 10000, sample_steps = 10000, thin_steps = 10,
              make_comparisons = TRUE, model = model,m_name = "hc",
              save_path = "C:/simulations")
end_time <- Sys.time()

hc_time <- end_time - start_time

# for no interactions
all_sim <- simulate_from_sim_df(sim_df, test = FALSE,
                                known_inxs = c(1,0,0,0,1,0,0,0,1))

start_time <- Sys.time()
batch_analyze(all_sim = all_sim, params = params, n_chains = 7,adapt_steps = 3000,
              burn_in = 10000, sample_steps = 10000, thin_steps = 10,
              make_comparisons = TRUE, model = model,
              comparison_model = comparison_model, m_name = "no_inxs",
              save_path = "C:/simulations")
end_time <- Sys.time()

no_inxs_time <- end_time - start_time
