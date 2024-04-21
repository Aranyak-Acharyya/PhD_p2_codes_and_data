####CODE__FOR__POWER__CONVERGENCE#####


library(parallel)
library(iterators)
no_of_cores<-detectCores()
clust<-parallel::makeCluster(no_of_cores)
library(foreach)
library(doParallel)
registerDoParallel()


library(Matrix)
library(MASS)
library(MVA)
library(irlba)
library(igraph)
library(smacof)
library(pracma)
library(lsbclust)





e <- new.env()
e$libs <- c("Matrix","MASS","lsbclust",
            "MVA","irlba","igraph","pracma",
            "smacof",.libPaths())
clusterExport(clust, "libs", envir=e)
clusterEvalQ(clust, .libPaths(libs))


RNGkind("L'Ecuyer-CMRG")
set.seed(100)


#setting universal parameters
d<-2
s<-5
D<-d^2

angle<-1
a<-(2^0.5)/sin(angle)
b<-(2^0.5)/cos(angle)




alpha<-2.0
beta<-5.0
sig_ep<-0.1


thres<-qf(0.95,1,s-2)





k_vec<-seq(1,20,1)

rho_func<-function(n)
{
  res<-1/log(n)
  return(res)
}

lambda_func<-function(k)
{
  lambda<-0.95*(0.99)^(k-1)
  return(lambda)
}

n_func<-function(k)
{
  n<-16+4*(k-1)
  return(n)
}

N_func<-function(k)
{
  N<-12+(k-1)
  return(N)
}








clusterExport(clust,list("d","s","a","b",
                         "alpha","beta","sig_ep"))

out_mat<-matrix(,ncol=7)



for(k in k_vec)
{
  
  lambda<-lambda_func(k)
  n<-n_func(k)
  N<-N_func(k)
  N_star<-as.integer(N^0.85)
  
  rho<-rho_func(n)
  
  
  #creating community membership matrix
  Z<-cbind(c(rep(1,n/2),rep(0,n/2)),
           c(rep(0,n/2),rep(1,n/2)))
 
  
  clusterExport(clust,list("lambda","n","rho",
                           "N","Z"))
  
  R<-foreach(i=1:100,.combine = 'rbind') %dopar%
    {
      
      #selecting timepoints
      ts_vec<-runif(s,min=0.25,max=1)
      tm_vec<-runif(N-s,min=0.25,max=1)
      t_vec<-c(ts_vec,tm_vec)
      
      
      #generating responses
      ys<-alpha+beta*ts_vec+rnorm(s,mean=0,
                                    sd=sig_ep)
     
      
      
      
      
      
      
      
      L<-foreach(t=t_vec,.combine = 'cbind') %do%
        {
          
          #defining score matrix
          B<-matrix(c(t/a,t/b,t/b,t/a),nrow=d)
          
          
          
          #defining probability matrix
          P<-Z%*%B%*%t(Z)
          print(P)
          
          #generating adjacency matrix
          pvec<-c(P[lower.tri(P,diag=TRUE)])
          avec<-rbinom(length(pvec),1,pvec)
          A<-matrix(nrow=nrow(P),ncol=ncol(P))
          A[lower.tri(A,diag=TRUE)]<-avec
          A[upper.tri(A)]<-t(A)[upper.tri(A)]
          
          A
          
        }
      
      
      
      M<-foreach(t=t_vec,.combine = 'cbind') %do%
        {
          
          r<-which(t_vec==t)
          A<-L[,((r-1)*n+1):(r*n)]
          
          
          
          
          #obtaining left singular matrix
          A_irlba<-irlba(A,d)
          Vi_hat<-A_irlba$u
          
          Vi_hat
          
        }
      
      
      
      M_irlba<-irlba(M,d)
      V_hat<-M_irlba$u
      
      S<-foreach(t=t_vec,.combine = 'c') %do%
        {
          
          r<-which(t_vec==t)
          A<-L[,((r-1)*n+1):(r*n)]
          
          
          
          
          #estimating sparsity factor individually
          rho_hat_k<-mean(A[upper.tri(A)])
          
          
        }
      
      #estimating sparsity factor globally
      rho_hat<-mean(S)
      
      
      
      Q<-foreach(t=t_vec,.combine = 'rbind') %do%
        {
          #extracting adjacency matrix
          r<-which(t_vec==t)
          A<-L[,((r-1)*n+1):(r*n)]
          
          #estimating score matrix
          R_hat<-t(V_hat)%*%A%*%V_hat
          Q_hat<-(1.0/n)*R_hat
          
          q_hat<-c(Q_hat)
          
          t(q_hat)
          
        }
      
      
      
      
      #declaring matrix to store ASE of the special node
      ASE_s_mat<-Q[1:N_star,]
      
      
      
      
      #creating distance matrix
      B0<-as.matrix(dist(ASE_s_mat,method = "euclidean",
                         diag=TRUE,upper=TRUE))
      
      
      
      
      
      #forming matrix of localization graph
      BB<-ifelse(B0<lambda,B0,0)
      
      
      colnames(BB)<-as.character(seq(1,nrow(B0),1))
      
      
      
      g<-graph_from_adjacency_matrix(BB,
                                     mode="undirected",
                                     weighted=TRUE,
                                     diag=FALSE,
                                     add.colnames = NULL,
                                     add.rownames = NA)
      
      
      
      
      #matrix of shortest path distances
      #D<-shortest.paths(g,v=V(g)[1:s],to=V(g)[1:s])
      D<-distances(g,v=V(g)[1:s],to=V(g)[1:s],mode = "all",
                   weights = NULL,
                   algorithm = c("automatic", "unweighted", 
                                 "dijkstra", "bellman-ford", 
                                 "johnson"))
      Ds<-as.matrix(D)
      
      
      
      
      
      #raw-stress minimization
      MM<-mds(Ds,ndim = 1,type = "ratio",
              weightmat = NULL,
              init = "torgerson")
      
      #finding the scaling factor
      dl<-as.vector(MM$delta)
      dh<-as.vector(MM$dhat)
      fac<-mean(dh/dl)
      
      
      #raw-stress embeddings
      zz<-as.vector(MM$conf)
      z<-zz/fac
      zs<-z[1:s]
      
      
      
      
      
      
      #prediction by true regressors
      df<-as.data.frame(cbind(ts_vec,ys))
      modt<-lm(ys~ts_vec,data = df)
      y_hat_true<-predict(modt)
      
      
      
      
      
      #prediction by embeddings
      df<-as.data.frame(cbind(zs,ys))
      modz<-lm(ys~zs,data=df)
      y_hat_sub<-predict(modz)
      
      
      
      
      
      
      num_true<-sum((y_hat_true-mean(ys))^2)
      den_true<-sum((ys-y_hat_true)^2)
      F_true<-(s-2)*(num_true/den_true)
      
      
      num_sub<-sum((y_hat_sub-mean(ys))^2)
      den_sub<-sum((ys-y_hat_sub)^2)
      F_sub<-(s-2)*(num_sub/den_sub)
      
      
      res<-c(F_true,F_sub)
      
    }
  
  power_true<-length(which(R[,1]>thres))/(nrow(R))
  power_sub<-length(which(R[,2]>thres))/(nrow(R))
  
  power_diff<-abs(power_true-power_sub)
  
  
  
  dec<-c(k,nrow(R),n,N,power_true,power_sub,
         power_diff)
  print(dec)
  
  out_mat<-rbind(out_mat,dec)
  
  
}

stopCluster(clust)



out_mat<-out_mat[-1,]

df<-data.frame(out_mat)
save(df,file = "CopyofP_BalBlock_new14_sparse_est.RData")



load("CopyofP_BalBlock_new14_sparse_est.RData")
out_mat<-as.matrix(df)

rownames(out_mat)<-NULL
colnames(out_mat)<-NULL

print(out_mat)

k_vec<-out_mat[,1]
n_vec<-out_mat[,3]
N_vec<-out_mat[,4]
power_true_vec<-out_mat[,5]
power_sub_vec<-out_mat[,6]
power_diff_vec<-out_mat[,7]






library(ggplot2)
library(reshape2)
library(latex2exp)
library(ggbreak)



df1<-data.frame(k_vec,n_vec,N_vec,power_diff_vec)
ggplot(df1, 
           aes(x=k_vec,y=power_diff_vec)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks=seq(0,20,1)) +
  ylab(TeX("$|\\hat{\\pi}-\\pi^*|$")) +
  xlab(TeX("$K")) 
  

ggsave(file="P2copyofplot14.eps", 
       width = 5, height = 2.5,
       units = "in", dpi = 500)




col1_ASE_vec<-comm2_node_ASE_mat[,1]
col2_ASE_vec<-comm2_node_ASE_mat[,2]

df2<-data.frame(col1_ASE_vec,col2_ASE_vec)
g2<-ggplot(df2, aes(x=col1_ASE_vec,y=col2_ASE_vec)) +
  geom_point() +
  geom_line() +
  xlab(TeX("dim1 of estimated latent positions")) +
  ylab(TeX("dim2 of estimated latent positions"))

g2
  




df<-data.frame(n_vec,risk1_vec,risk2_vec)
dfm<-melt(df, id.vars = 'n_vec')
print(dfm)
ggplot(dfm, aes(x=n_vec, y=value, 
                colour = variable)) +
  geom_point() +
  geom_line() +
  ylab(TeX("sample MSEs")) +
  xlab(TeX("number of nodes(n)")) +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = c("red","orange"),
                      labels=unname(TeX(c(
                        "sample MSE of  $\\hat{\\theta}_{true}$",
                        "sample MSE of $\\hat{\\theta}_{sub}$"))))



ggsave(file="P1plot6.eps", width = 5, height = 2,
       units = "in", dpi = 300)


