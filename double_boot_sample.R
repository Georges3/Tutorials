 
rm(list=ls())
setwd("C:/Users/GeorgesB/Downloads/nonparNet-master/nonparNet-master/code/Run_Examples")

source("network_construction.R")
source("random_network.R")
source("resam_patches_functions.R")

net <- random_network(n=20, distrib="ztpois", param = 3, degree = NULL, take.p = 0.05)

#we randomly select 5 seeds for the original sample
set.seed(42)
sam<-patch_sampling(net,n.seed=5,n.wave=2) 
var<-paste("subnet",1:5,sep="_")
library(igraph)
for(i in 1:5){
  assign(var[i],isubgraph(sam[[i]]$subg))
}
par(mfrow=c(3,2),mar=c(1,1,0,0))
for(i in 1:5){
  plot(get(var[i]),arrow_mode="–",label.cex=5)
}

degree(get(var[1]))
degree(subnet_1)

####----Simple bootstrap--------------##
# B_1: Number of Bootstrap samples to take

boot_simple<- function(var,seed,B_1){
  b_var<- array(NA,dim = c(1,seed,B_1))
  
  for (i in 1:B_1) {
    b_var[,,i]<- sample(var,seed,replace = T)
  }
  b_var
}
boot_sample_net <- boot_simple(var,seed=5,B_1=10)


par(mfrow=c(3,2),mar=c(1,1,0,0))
for(i in 1:5){
  
  plot(get(boot_sample_net[,i,1]),arrow_mode="–",label.cex=5)
}

##--------------Double bootstrap-----------

# B_2:Number of bootstrap sample to consider form each bootstrap i: i=1,..., B_1
double_boot <- function(sam_boot,seed,B_1,B_2){
  db_var<- vector("list",B_1)
  for (i in 1:B_1) {
    b_var_new <- array(NA,dim = c(1,seed,B_2))
    for(j in 1:B_2){
      b_var_new[,,j] <- sample(sam_boot[,,j],seed,replace = T)
    }
    db_var[[i]] <- b_var_new
  }
  db_var
}

db_boot_sam_net <- double_boot(sam_boot=boot_sample_net,seed=5,B_1=10,B_2=3)
#----------------------------------------------------------------------------------------------

##---------Mean degree in Original sample---------------------#####

net_degre<- cbind(1:net$n,net$degree)   # The degree sequence of the network
#sam[[1]]$vwaves       # all vertices of the first sampled patch 
#sam[[1]]$vwaves[[1]]  # the seed of the first sampled patch 
nseed=5
#sam[[1]]$vwaves[[2]]   # the first waves of the first sampled patch
sample_mean_degree <- sample_sum_graph<-n_ver<- rep(0,nseed)  

for(i in 1:nseed){ 
sample_seed_waves<-  unlist(c(sam[[i]]$vwaves[[1]],sam[[i]]$vwaves[[2]])) # combine sampled seeds and first waves in each sampled patches
sample_mean_degree[i] <- mean(net_degre[sample_seed_waves,2])       # mean degree in each sampled patch
sample_sum_graph[i] <- sum(net_degre[sample_seed_waves,2])  # sum degree in each sampled patch
n_ver[i] <- length(sample_seed_waves)  # number of vertices in combination of sampled  seeds and first waves
                                        # in each sampled patch
}

######-------------Mean degree in Bootstrap-----------------------##################

all_degrees <- as.numeric(substring(boot_sample_net,8,8)) ## all degrees in the sample
AA <- matrix(all_degrees, ncol = nseed,byrow = T)    # each row gives patches resampled in every bootstrap
boot_mean_degree <- matrix(sample_mean_degree[AA],ncol=nseed,byrow = T) ## mean bootstrap


# Se puede comparar mean of means or sum(sum_gr)/sum(n_ver)

# Instead of sampling patches, why not just resampling the means?

boot_simple_mean<- function(var,seed,B_1){
  b_var_mean<- matrix(NA,nrow = B_1,ncol = seed)
  
  for (i in 1:B_1) {
    b_var_mean[i,]<- sample(var,seed,replace = T)
  }
  b_var_mean
}
boot_sample_net_mean <- boot_simple_mean(sample_mean_degree,seed=nseed,B_1=10)

###------------Mean degree en double bootstrap-----------------------###

 # As db_boot_sam_net is a list
## take first each db_boot_sam_net[[i]] which is an array

dble_bt_mn_dgr <- function(sam_double_boot_net,seed,B_1,B_2){   # sam_double_boot= db_boot_sam_net
Db_mean <- array(NA,dim = c(B_2,seed,B_1))
for(i in 1:B_1){
aa_bouble <- as.numeric(substring(sam_double_boot_net[[i]],8,8))
AA_db <- matrix(aa_bouble, ncol = nseed,byrow = T)    # each row gives patches resampled in every bootstrap
Db_boot_mean_degree <- matrix(sample_mean_degree[AA_db],ncol=nseed,byrow = T) ## mean bootstrap
Db_mean[,,i] <-Db_boot_mean_degree
}
Db_mean
}
double_boot_mean_degree <- dble_bt_mn_dgr(sam_double_boot_net =db_boot_sam_net,seed=5,B_1=10,B_2=3)

corrected_mean_degree <- 3*mean(sample_mean_degree)-3*mean(boot_sample_net_mean)+mean(double_boot_mean_degree)

## Bias corrected Sample Mean degree 
corrected_mean_degree
