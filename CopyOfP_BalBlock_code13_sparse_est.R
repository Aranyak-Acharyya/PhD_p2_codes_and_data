####CODE__FOR__CONVERGENCE__OF__PREDICTED__RESPONSES######

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

RNGkind("L'Ecuyer-CMRG")
set.seed(1234)




e <- new.env()
e$libs <- c("Matrix","MASS","lsbclust",
            "MVA","irlba","igraph","pracma",
            "smacof",.libPaths())
clusterExport(clust, "libs", envir=e)
clusterEvalQ(clust,.libPaths(libs))






#setting universal parameters
d<-2
s<-5
ss<-(s+1)


angle<-1.0
a<-(2^0.5)/sin(angle)
b<-(2^0.5)/cos(angle)



alpha<-2.0
beta<-5.0
sig_ep<-0.01








k_vec<-seq(1,12,1)

rho_func<-function(n)
{
  res<-1/log(n)
  return(res)
}



lambda_func<-function(k)
{
  lambda<-2*(0.99)^(k-1)
  return(lambda)
}

n_func<-function(k)
{
  n<-500+100*(k-1)
  return(n)
}

N_func<-function(k)
{
  N<-15+(k-1)
  return(N)
}











clusterExport(clust,list("d","s","ss","rho_func",
                         "a","b",
                         "alpha","beta","sig_ep"))

out_mat<-matrix(,ncol=6)



for(k in k_vec)
{
  
  lambda<-lambda_func(k)
  n<-n_func(k)
  N<-N_func(k)
  N_star<-as.integer(N^0.75)
  
  rho<-rho_func(n)
  
  
  #creating community membership matrix
  Z<-cbind(c(rep(1,n/2),rep(0,n/2)),
           c(rep(0,n/2),rep(1,n/2)))
 
  
  clusterExport(clust,list("lambda","n","s","rho",
                           "ss","N","Z"))
  
  R<-foreach(i=1:100,.combine = 'c') %dopar%
    {
      
      #selecting timepoints
      ts_vec<-runif(s,min=0.25,max=1)
      tm_vec<-runif(N-s,min=0.25,max=1)
      t_vec<-c(ts_vec,tm_vec)
      
      
      #generating responses
      tss_vec<-t_vec[1:ss]
      yss<-alpha+beta*tss_vec+rnorm(ss,mean=0,
                                    sd=sig_ep)
      ys<-yss[1:s]
      
      y0<-yss[ss]
      
      
      
      
      
      L<-foreach(t=t_vec,.combine = 'cbind') %do%
        {
          
          #defining score matrix
          B<-matrix(c(t/a,t/b,t/b,t/a),nrow=d)
          
          
          
          #defining probability matrix
          P<-(Z%*%B%*%t(Z))*(rho)
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
          Q_hat<-(1.0/rho_hat)*(1.0/n)*R_hat
          
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
      #D<-shortest.paths(g,v=V(g)[1:ss],to=V(g)[1:ss])
      D<-distances(g,v=V(g)[1:ss],to=V(g)[1:ss],mode = "all",
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
      z0<-z[ss]
      zss<-z[1:ss]
      
      
      
      
      
      #prediction by true regressors
      df<-as.data.frame(cbind(tss_vec,yss))
      df_train<-df[1:s,]
      df_test<-df[ss,]
      
      
      
      modt<-lm(yss~tss_vec,data = df_train)
      y_true<-predict(modt,df_test)
      
      
      #prediction by embeddings
      df<-as.data.frame(cbind(zss,yss))
      df_train<-df[1:s,]
      df_test<-df[ss,]
      
      modz<-lm(yss~zss,data=df_train)
      y_sub<-predict(modz,df_test)
      
      
      store<-(y_true-y_sub)^2 
      store
      
    }
  
  se<-sd(R)/10
  
  dec<-c(k,length(R),n,N,mean(R),se)
  print(dec)
  
  out_mat<-rbind(out_mat,dec)
  
  
}

stopCluster(clust)



out_mat<-out_mat[-1,]

df<-data.frame(out_mat)
save(df,file = "CopyofP_BalBlock_new13_sparse_est.RData")


load("CopyofP_BalBlock_new13_sparse_est.RData")
out_mat<-as.matrix(df)

rownames(out_mat)<-NULL
colnames(out_mat)<-NULL

print(out_mat)

k_vec<-out_mat[,1]
n_vec<-out_mat[,3]
N_vec<-out_mat[,4]
risk_vec<-out_mat[,5]
se_vec<-out_mat[,6]





library(ggplot2)
library(reshape2)
library(latex2exp)
library(ggbreak)



df1<-data.frame(k_vec,n_vec,N_vec,risk_vec,se_vec)
ggplot(df1, 
           aes(x=k_vec,y=risk_vec)) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin=risk_vec-se_vec, 
                    ymax=risk_vec+se_vec), 
                alpha=0.4,
                position=position_dodge(0.05)) +
  scale_x_continuous(breaks=seq(0,11,1)) +
  ylab(TeX("$E(\\hat{y}_{sub}-\\hat{y}_{true})^2$")) +
  xlab(TeX("$K")) 
  

ggsave(file="P2plot13_sparse_est.pdf", 
       width = 5, height = 2,
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


