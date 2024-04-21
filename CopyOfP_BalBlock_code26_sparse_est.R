####CODE__FOR__REAL__DATA__ANALYSIS########


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




load("missingEdge-tsg-forAranyak-2.RData")
print(data)
df1<-data[,-2]
print(df1)
col.order<-c("model","score","tg","missing_edge")
df1[,col.order]

load("missingEdge-tsg-forAranyak.RData")
print(data)
df2<-data
data<-rbind(df1,df2)
k_vec<-seq(1,nrow(data),1)
n<-140
d<-3

print(data)




#specifying model-type
mspec<-vector()
for(i in 1:143)
{
  if(i<=130)
  {
    if(i%%10!=0)
    {
      q<-i%/%10
      mspec[i]<-(q+1)
    }
    if(i%%10==0)
    {
      q<-i%/%10
      mspec[i]<-q
    }
  }
  if(i>130)
  {
    mspec[i]<-(i-130)
  }
}



#for(i in 1:143)
#{
#  if(i%%13!=0)
#  {
#    mspec[i]<-i%%13
#  }
#  if(i%%13==0)
#  {
#    mspec[i]<-13
#  }
#}

data$model<-mspec

length(mspec)






#choose threshold for censoring adjacency matrix
p<-0.25




clusterExport(clust,list("d","n","p","k_vec"))

#the quantity 'tim' denotes the position from which a graph from
#each TSG is to be selected, so for our paper tim=40,17,30
tim<-40
  
clusterExport(clust,"tim")


L<-foreach(k=k_vec,.combine='cbind') %do%
{
  
  g<-data$tg[[k]][[tim]]
  
  
  
  g<-as.undirected(
    g,
    mode =  "each", 
    edge.attr.comb = igraph_opt("edge.attr.comb")
  )
  
   
  weight<-E(g)$weight
  
  
  
  
  A<-as_adjacency_matrix(
    g,
    type = "both",
    attr = "weight",
    names = FALSE
  )
  
 
  vec_modA<-as.vector(abs(A))
  vnz_modA<-vec_modA[!vec_modA %in% c(0)]
  
  
  qq<-quantile(vnz_modA,probs=p)
  
  
  
  AA<-matrix(ifelse(as.vector(abs(A))>qq,1,0),
             nrow=n,byrow=FALSE)
  
  AA
}








M<-foreach(k=k_vec,.combine = 'cbind') %do%
  {
    
    r<-which(k_vec==k)
    A<-L[,((r-1)*n+1):(r*n)]
    
    
    
    
    
    
    #obtaining left singular matrix
    A_irlba<-irlba(A,d)
    Vi_hat<-A_irlba$u
    
    Vi_hat
    
  }



M_irlba<-irlba(M,d)
V_hat<-M_irlba$u

S<-foreach(k=k_vec,.combine = 'c') %do%
  {
    
    #extracting the adjacency matrix
    r<-which(k_vec==k)
    A<-L[,((r-1)*n+1):(r*n)]
    
    
    rho_hat_k<-mean(A[upper.tri(A)])
    
  }

rho_hat<-mean(S)

clusterExport(clust,"rho_hat")


Q<-foreach(k=k_vec,.combine = 'rbind') %do%
  {
    #extracting adjacency matrix
    r<-which(k_vec==k)
    A<-L[,((r-1)*n+1):(r*n)]
    
    #estimating score matrix
    R_hat<-t(V_hat)%*%A%*%V_hat
    Q_hat<-(1/rho_hat)*(1.0/n)*R_hat
    
    
    q_hat<-Q_hat[upper.tri(Q_hat,diag=TRUE)]
    
    t(q_hat)
    
  }

stopCluster(clust)

dQ<-as.data.frame(Q)
print(dQ)
#pairs(dQ)


Chat_corr_coeff_mat<-xicor(dQ)
print(Chat_corr_coeff_mat)


library(ggplot2)
library(GGally)






#declaring matrix to store ASE of the special node
ASE_s_mat<-Q
#ASE_s_mat<-Q[,c(1,2,3,5,6,9)]









gr1<-ggpairs(as.data.frame(ASE_s_mat)) +
theme(axis.line=element_blank(),
      axis.text=element_blank(),
      axis.ticks=element_blank())

gr1

#ggsave(file="P2pairs_all_t21.pdf",
#       plot=gr1,
#       width=3,height=3,
#       units="in",dpi=500)


x1<-ASE_s_mat[,1]
x2<-ASE_s_mat[,4]
dfr<-data.frame(x1,x2,mspec)

gr2<-ggplot(dfr, 
            aes(x=x1,y=x2)) +
  geom_point() +
  #geom_line() +
  scale_x_continuous(breaks=scales::pretty_breaks(n=5)) +
  scale_y_continuous(breaks=scales::pretty_breaks(n=5)) +
  xlab(TeX("$(1,1)$-th entry")) +
  ylab(TeX("$(1,3)$-th entry")) +
  theme(legend.position = "none")


gr2




ggsave(file="P2_scatter_relevant.pdf",
       plot = gr2,
       width = 3, height = 3,
       units = "in", dpi = 500)





#creating distance matrix
B0<-as.matrix(dist(ASE_s_mat,method = "euclidean",
                   diag=TRUE,upper=TRUE))


lambda<-0.005


#forming matrix of localization graph
BB<-ifelse(B0<lambda,B0,0)


colnames(BB)<-as.character(seq(1,nrow(B0),1))



g<-graph_from_adjacency_matrix(BB,
                               mode="undirected",
                               weighted=TRUE,
                               diag=FALSE,
                               add.colnames = NULL,
                               add.rownames = NA)

a<-is.connected(g)


while(a==FALSE)
{
  lambda<-lambda+0.005
  
  #forming matrix of localization graph
  BB<-ifelse(B0<lambda,B0,0)
  
  
  colnames(BB)<-as.character(seq(1,nrow(B0),1))
  
  
  
  g<-graph_from_adjacency_matrix(BB,
                                 mode="undirected",
                                 weighted=TRUE,
                                 diag=FALSE,
                                 add.colnames = NULL,
                                 add.rownames = NA)
  
  a<-is.connected(g)
  
}



#matrix of shortest path distances
D<-distances(g,v=V(g),to=V(g),mode = "all",
             weights = NULL,
             algorithm = c("automatic", "unweighted", 
                           "dijkstra", "bellman-ford", 
                           "johnson")
)
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

print(z)


y_eff<-data$score[k_vec]


df<-cbind(k_vec,y_eff,z,as.factor(mspec))
print(df)

colnames(df)<-c("k","y","z","mspec")

df<-as.data.frame(df)
print(df)

y<-df$y
z<-df$z

mod<-lm(y~z,data=df)

summary(mod)



p<-summary(mod)$coefficients[2,4]

rh<-cor(y,z)

res<-c(tim,p,rh,lambda)

print(res)








library(ggplot2)
library(reshape2)
library(latex2exp)
library(ggpmisc)
library(ggpubr)
library(broom)


lb<-c("one","two","three","four",
      "five","six","seven","eight",
      "nine","ten","eleven","twelve",
      "thirteen")

#gr4<-ggplot(df,aes(x=z,y=y,label=lb)) + 
#  geom_point() +
#  geom_text(label=lb) +
#  ylab(TeX("$y_i$")) +
#  xlab(TeX("$\\hat{z}_i$"))+
#  stat_smooth(method = "lm", 
#              formula = y ~ x, 
#              geom = "smooth") 

length(y)







gr4<-ggplot(df,aes(x=z,y=y)) + 
  geom_point() +
  ylab(TeX("$y_i$")) +
  xlab(TeX("$\\hat{z}_i$"))+
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth") +
  stat_correlation(mapping = use_label(c("R", "n","p","R2")),
                   label.x = "left",
                   label.y = "top") 
  #stat_poly_eq(aes(label =  paste(after_stat(eq.label), "*\" with \"*", 
  #                                after_stat(rr.label), "*\", \"*",
  #                               after_stat(p.value.label), "*\".\"",
  #                                sep = "")),
  #             formula = y~x, size = 3)




gr4









ggsave(file="P2pairs_all_t40.pdf",
       plot=gr1,
       width=5,height=5,
       units="in",dpi=500)


ggsave(file="P2res_vs_embedding_all_t40.pdf",
       plot=gr4,
       width=5,height=4,
       units="in",dpi=500)




data$model<-as.factor(data$model)

gr6<-ggplot(data,aes(x=model,y=score,fill=model)) +
  geom_violin() +
  geom_boxplot(width=0.1) +
  labs(x="model",y="learning score") +
  theme(legend.position="none")

gr6


ggsave(file="P2model_vs_score.pdf",
       plot=gr6,
       width=5,height=4,
       units="in",dpi=500)






library(gridExtra)
library(grid)
library(ggplotify)

#Gr1<-as.ggplot(~gr1)

gr5<-arrangeGrob(gr2,gr4,ncol=2)

plot(gr5)



ggsave(file="P2intro_plot_all_t40.pdf",gr5,
       width=7,height=3.5,
       units="in",dpi=500)



#####THE__END__FOR__REAL__DATA__SECTION__NEXT__NONPARAMETRIC########



library(ggplot2)
library(reshape2)
library(latex2exp)
library(np)

z_modnp<-z
print(z)
print(z_modnp)
df_modnp<-data.frame(y,z,z_modnp)
print(df_modnp)


modnp<-npreg(y~z_modnp,
             regtype="ll",
             gradients=TRUE,
             data=df_modnp)

summary(modnp)

plot(modnp,gradients=TRUE,plot.errors.method="asymptotic")

#yhatnp<-predict(modnp,newdata=df)



z_modnp<-seq(min(z)-0.005,max(z)+0.005,length.out=400)
print(z_modnp)

yhatnp<-predict(modnp,newdata=as.data.frame(z_modnp))

print(yhatnp)

dfnp2<-data.frame(z_modnp,yhatnp)



#dfnp<-data.frame(z,y,yhatnp)

#dfnpm<-reshape::melt(dfnp,id.vars = 'z')

grnp<-ggplot(df,aes(x=z,y=y)) +
  geom_point() +
  geom_line(data=dfnp2,aes(x=z_modnp,y=yhatnp),color="blue") +
  xlab(TeX("$\\hat{z}_i$")) +
  ylab(TeX("$y_i$"))

grnp

ggsave(file="P2res_vs_embedding_np.pdf",
       grnp,
       width=5,
       height=4,
       units="in")







library(gridExtra)
library(grid)

#grid.arrange(gr2,gr4,ncol=2)









gr5<-arrangeGrob(gr2,gr4,ncol=2)

plot(gr5)



ggsave(file="P2intro_plot.pdf",gr5,
       width=7,height=3.5,
       units="in",dpi=500)







