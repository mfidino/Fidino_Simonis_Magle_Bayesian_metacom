
#bring in the data

library(foreach)
library(ROCR)
library(doParallel)

# bring in the all_sim stuff

known_loc <- "C:/simulations/dynamic_occupancy/knowns/"

#all_sim <- readRDS(paste(known_loc,"all_sim_hc1.txt", sep = ""))

model_loc <- "C:/simulations/dynamic_occupancy/mcmc_matrix"

pull_auc <- function(model_loc = NULL, all_sim_loc = NULL, module = NULL, n_sim = NULL ){
# bring in the all of the simulations to get ROC stuff

# copy over the files onto the hard drive to speed up the read in
to_iter <- list.files(model_loc, pattern = module, full.names = TRUE)

order_iter <- paste( module,1:n_sim, "\\.", sep = "")
# order them 1:n_sim
to_order <- vector("list", length = n_sim)
for(i in 1:length(to_order)){
  to_order[[i]] <- grep(order_iter[i], to_iter)
}

# now they go 1, 2, 3, 4, etc.
to_iter <- to_iter[unlist(to_order)]

# make a cluster
cl <- makeCluster(7)

# register the cluster for parallel computing
registerDoParallel(cl)

known_matrix <- matrix(0, ncol = 24, nrow = length(to_iter))

(strt <- Sys.time())
# loop to go through each simulation and compute the auc
rmse_list <-coverage_list <- precision_list <- known_list <- z_auc_list <- est_list <- vector("list", length = length(to_iter))
for(its in 1:length(to_iter)){
  
file_num <- rep(1:400, each = 2) 
if(its %% 2 != 0){  
all_sim <- readRDS(paste(all_sim_loc, "all_sim_", module,file_num[its], ".txt", sep = ""))
}
  print(its)

# collect the knowns
  
# known_z = simulated z matrix, which is the actual true state  
  
known_z <- matrix(all_sim$mats$z, 
                  nrow = 1, ncol = length(all_sim$mats$z))

# knowns for inxs

inxs_known <- with(all_sim$sim_list,{
  c(l_gam, as.numeric(l_phi_mat),rep(0, nspec), l_beta, l_beta2, p)
})

known_matrix[its,] <- inxs_known
# knowns for no inxs

beta_known <- with(all_sim$sim_list,{
  c(l_gam, l_phi, rep(0, nspec), l_beta, l_beta2, p)
})


# read in table
f_mod <- read.table(to_iter[its], sep = "\t",
                    header = TRUE)


# get number of iterations
iters <- nrow(f_mod)

# pull out the estimated z's from the model
the_z <- f_mod[,grep("z", colnames(f_mod))]

# pull out every other non "z" thing
m_coefs <- f_mod[,2:(min(grep("z", colnames(f_mod)))-1)]

est_list[[its]] <- apply(m_coefs, 2, quantile, probs = c(0.025, 0.5, 0.975))
# make a list for z matrix stuff

ls <- vector('list', length = iters)

# parallel computing for the auc and rmse of the z matrix
ls <- foreach(i = 1:iters, .packages = 'ROCR') %dopar% {
  psi.vals <- as.numeric(the_z[i,])
  z.vals <- known_z
  pred <- prediction(psi.vals, factor(z.vals, levels = c("0", "1")))
  perf <- performance(pred, "auc")
  AUC1 <- perf@y.values[[1]]
  perf <- performance(pred, "tpr", "fpr")
  rmse_z <- sqrt(mean((z.vals-psi.vals)^2))
   c(AUC1, rmse_z)
  
}


# calc rmse, coverage, and precision for parameters




if(its %% 2 == 0){ # even numbers inxs models, odds beta models
  names(inxs_known) <- colnames(m_coefs)
  rmse <- sqrt(apply(sweep(m_coefs,2,inxs_known)^2,2, mean))
  cred_i <- apply(m_coefs, 2, quantile, probs = c(0.025, 0.975))
  precision <- apply(cred_i, 2, dist)
  coverage <- cred_i[1,] <inxs_known & cred_i[2,] > inxs_known
  rmse_list[[its]] <- data.frame(as.list(rmse))
  coverage_list[[its]] <- data.frame(as.list(coverage))
  precision_list[[its]] <- data.frame(as.list(precision))
  known_list[[its]] <- data.frame(as.list(inxs_known))
}
if(its %% 2 != 0){
  names(beta_known) <- colnames(m_coefs)
  rmse <- sqrt(apply(sweep(m_coefs,2,beta_known)^2,2, mean))
  cred_i <- apply(m_coefs, 2, quantile, probs = c(0.025, 0.975))
  precision <- apply(cred_i, 2, dist)
  coverage <- cred_i[1,] <beta_known & cred_i[2,] > beta_known
  rmse_list[[its]] <- data.frame(as.list(rmse))
  coverage_list[[its]] <- data.frame(as.list(coverage))
  precision_list[[its]] <- data.frame(as.list(precision))
  known_list[[its]] <- data.frame(as.list(beta_known))
}


# change ls into a matrix, ls is for the z matrix
my_ans <- matrix(unlist(ls), ncol = 2, nrow = iters, byrow = TRUE)
# name of columns of the matrix
colnames(my_ans) <- c("auc", "rmse_z")


z_auc_list[[its]] <- data.frame(as.list(apply(my_ans, 2, quantile, probs = 0.5)))


}
print(Sys.time() - strt)

stopCluster(cl)

to_return <- list(auc_rmse_z = z_auc_list, 
                  rmse = rmse_list,
                  cov = coverage_list,
                  pres = precision_list,
                  coefs = known_list, 
                  est_list = est_list,
                  known_matrix = known_matrix)
return(to_return)

}



sum_stat_hc <- pull_auc(model_loc, known_loc, module = "hc", n_sim = 400)
sum_stat_no_inxs <- pull_auc(model_loc, known_loc, module = "no_inxs", n_sim = 400)
sum_stat_pcp <- pull_auc(model_loc, known_loc, module = "pcp", n_sim = 400)

saveRDS(sum_stat_hc, "C:/simulations/dynamic_occupancy/summary/sum_stat_hc.txt" )
saveRDS(sum_stat_no_inxs, "C:/simulations/dynamic_occupancy/summary/sum_stat_no_inxs.txt" )
saveRDS(sum_stat_pcp, "C:/simulations/dynamic_occupancy/summary/sum_stat_pcp.txt")

# compile the data in a relevant way to summarize over all the models

known <- plyr::rbind.fill(sum_stat$coefs)
rmse <- plyr::rbind.fill(sum_stat_no_inxs$rmse)
coverage <- plyr::rbind.fill(sum_stat$cov)
pres <- plyr::rbind.fill(sum_stat$pres)
z_rmse <- plyr::rbind.fill(sum_stat$auc_rmse_z)

# remove the psi stuff and merge the phi information, currently
#the beta model has phi.n instead of a matrix. # note this only
# works for this and only this set up.

move_phi_chuck_psi <- function(x){
  x[seq(1, nrow(x),2),c(19, 23, 27)] <- x[seq(1, nrow(x),2),c(4:6)]
  x <- x[,-grep("psi", colnames(x))]
  x <- x[,-c(4:6)]
  return(x)
}

known <- move_phi_chuck_psi(known)
rmse <- move_phi_chuck_psi(rmse)
coverage <- move_phi_chuck_psi(coverage)
pres <- move_phi_chuck_psi(pres)



                   




# now, we need to make comparisons between the beta model(odd rows)
# and the inxs model (even rows)

compare_results <- function(x){
  beta <- apply(x[seq(1,nrow(x),2),], 2, quantile, probs = c(0.025, 0.5, 0.975),
                na.rm = TRUE)
  inxs <- apply(x[seq(2,nrow(x),2),], 2, quantile, probs = c(0.025, 0.5, 0.975),
                na.rm = TRUE)
  to_return <- rbind.fill(data.frame(beta), data.frame(inxs))
  to_return$model <- rep(c("beta", "inxs"), each = 3)
  return(to_return)
}

sum_cov <- function(x){
  beta <- quantile(apply(x[seq(1,nrow(x),2),], 1, sum, na.rm = TRUE,
                na.rm = TRUE)/15, probs = c(0.025, 0.5, 0.975))
  inxs <- quantile(apply(x[seq(2,nrow(x),2),], 1, sum, na.rm = TRUE,
                na.rm = TRUE)/21, probs = c(0.025, 0.5, 0.975))
  to_return <- data.frame(coverage = c(beta, inxs))
  to_return$model <- rep(c("beta", "inxs"), each = 3)
  to_return$probs <- rep(c(0.025, 0.5, 0.975), 2)
  return(to_return)
}

sum_cov_ind <- function(x){
  beta <- apply(x[seq(1,nrow(x),2),], 2, mean, na.rm = TRUE)
  inxs <- apply(x[seq(2,nrow(x),2),], 2, mean, na.rm = TRUE)
  to_return <- data.frame(rbind(beta, inxs))
  to_return$model <- c("beta", "inxs")
  return(to_return)
                        
}

sum_z_rmse <- function(x){
  
}

rmse_ans <- compare_results(rmse)
coverage_ans <- sum_cov(coverage)
coverage_ind <- sum_cov_ind(coverage)
pres_ans <- compare_results(pres)
z_rmse_ans <- compare_results(z_rmse)


test <- data.frame(phi = rmse$phi.1.1., 
                   p1 = known$p.1.,
                   p2 = known$p.2.,
                   p3 = known$p.3.,
                   model = factor(rep(c("beta", "inxs"), 400)))

x <- rmse
beta <- x[seq(1,nrow(x),2),]
inxs <- x[seq(2,nrow(x),2),]

# collect the baseline phi

betap <- beta[,c(13,17,21)]
inxsp <- inxs[,c(13,17,21)]



bpm <- apply(betap, 2, mean)
ipm <- apply(inxsp, 2, mean)

psd <- apply(rbind(betap, inxsp), 2, sd)

(bpm - ipm) / psd
diffs <- betap - inxsp

apply(diffs, 2, mean)

ggplot(test, aes(x = p1, y = phi))+
  stat_smooth(size = 3, aes(color = model))+
  geom_point(aes(color = model))+
  xlab("species 1 detection probability")+
  ylab("confidence interval width, sp 1 persistence")


ls <- vector('list', length = iters)

# parallel computing for the auc and rmse of the z matrix
ls <- foreach(i = 1:400, .packages = 'ROCR') %dopar% {
  psi.vals <- as.numeric(the_z[i,])
  z.vals <- known_z
  pred <- prediction(psi.vals, factor(z.vals, levels = c("0", "1")))
  perf <- performance(pred, "auc")
  AUC1 <- perf@y.values[[1]]
  perf <- performance(pred, "tpr", "fpr")
  rmse_z <- sqrt(mean((z.vals-psi.vals)^2))
  c(AUC1, rmse_z)
  
}
