library(ape)
library(deSolve)
library("strap")
library("stringr")
library("ips")
library("ggtree")
require(treeio)
require(ggplot2)
require(ggtree)
library("reshape2")
library("ggstance")
library("diversitree")
library("phytools")
library("ggimage")
library("ggpubr")
library("dplyr")
library("patchwork")
library(tidyverse)

#############################################
#Code to recreate the result from the paper.
#"Early origin of viviparity and multiple reversions to oviparity
#in squamate reptiles" by R. Alexander Pyron and Frank T. Burbrink

#load the tree
phy <- ape::read.tree(file='data.txt')
#plot the tree using ggtree
ggplot(phy)+geom_tree()+theme_tree2()


#load the tip states 
data1 <- read.csv("data.csv")

#note there are only 4000 tips, but we have more than 8000 species in the excel,
#the phylogenetic tree is pruned, need to match the tips first then use
#BiSEE for phylogenetic analysis.

match1 <- match(phy$tip.label,data1$Species)
num_na <- str_count(match1, "NA")

#There are NA which means not matched, according to Pyron, "We pruned the tree
#to include only those species for which we had parity data (3951 species total)"

#Thus we need to prun the tree again

#tree pruning algorithm
j=0
v=c()
for (i in 1:length(match1)) {
  if(is.na(match1[i])){
    v <- c(v,phy$tip.label[i])
    j=j+1
  }
}

#new phylogenetic tree after pruning
newphy <- drop.tip(phy,v,trim.internal=TRUE)
match2 <- match(newphy$tip.label,data1$Species)
ggplot(newphy)+geom_tree()+theme_tree2()
 

#now need to add states in phylo such that we could apply BiSSE
#notice that there are 3962 tips, but in the paper they have 3951 in total
#probably because there are some uncertain tips (0/1) in the data and they also 
#get rid of those tips, we'll see...
subdata = subset(data1,select=-c(Reference))
comb_tree <- inner_join(as_tibble(newphy),subdata,by = c('label'='Species'))

id <- grep("0/1",comb_tree$Parity)
length(id)

v2=c()
for (i in 1:length(id)){
  v2 <- c(v2,newphy$tip.label[id[i]])
}

finaltree <- drop.tip(newphy,v2,trim.internal=TRUE)

test <- as_tibble(finaltree)
comb_tree <- left_join(as_tibble(finaltree),subdata,by = c('label'='Species'))

new_tree <- as.phylo(comb_tree)
ggplot(new_tree)+geom_tree()+theme_tree2()

########################################################
#pruning complete
#new_tree --- final phylogenetic tree
#comb_tree$Parity[1:3951] --- Parity mode


#visualization
tree1 <- as.treedata(comb_tree)

true_tree<-ggtree(tree1,right=TRUE,options(ignore.negative.edge=TRUE),branch.length='none')+
  geom_point(aes(color=Parity))+scale_color_manual(values=c("#0072BD","#D95319"))+
  ggtitle("Parity Mode of Squamate Reptiles") + 
  theme(text = element_text(size = 12.0,family = "serif"),plot.title = element_text(size=12))+labs(colour="Parity Mode")
true_tree

#ggsave("Pruned_tree.pdf",width = 15,height = 20) 


#######################################################
#Find parameter estiamtions using BiSSE (ML & Bayesian)
#first force the tree to be ultrametric due to rounding error (otherwise BiSSE won't work)
phy1 <- force.ultrametric(new_tree,method = "nnls")

ggplot(phy1)+geom_tree()+theme_tree2()

states <- setNames(as.numeric(comb_tree$Parity[1:3951]) ,new_tree$tip.label[1:3951])
lik <- make.bisse(phy1,states,sampling.f = c(0.4668068,0.630988))
lik(pars)


#perform ML
p <- starting.point.bisse(phy1)
p
#full model
fit <- find.mle(lik,p)
fit$lnLik
round(coef(fit),5)

#different hypothesis settings 

#no reversal 
lik.r <- constrain(lik, q10 ~ 0)
fit.r <- find.mle(lik.r, p[argnames(lik.r)])
fit.r$lnLik
round(coef(fit.r),5)

#equal lambda
lik.l <- constrain(lik, lambda1 ~ lambda0)
fit.l <- find.mle(lik.l, p[argnames(lik.l)])
fit.l$lnLik
round(coef(fit.l),5)

#equal mu
lik.m <- constrain(lik, mu1 ~ mu0)
fit.m <- find.mle(lik.m, p[argnames(lik.m)])
fit.m$lnLik
round(coef(fit.m),5)

#equal q
lik.q <- constrain(lik, q10 ~ q01)
fit.q <- find.mle(lik.q, p[argnames(lik.q)])
fit.q$lnLik
round(coef(fit.q),5)

#MK2
lik.mk2 <- constrain(lik, lambda1 ~ lambda0, mu1 ~ mu0)
fit.mk2  <- find.mle(lik.mk2 , p[argnames(lik.mk2 )])
fit.mk2 $lnLik
round(coef(fit.mk2 ),5)

#MK1
lik.mk1 <- constrain(lik, lambda1 ~ lambda0, mu1 ~ mu0, q01~ q10)
fit.mk1  <- find.mle(lik.mk1 , p[argnames(lik.mk1 )])
fit.mk1 $lnLik
round(coef(fit.mk1 ),5)

#mk1 no reversal
lik.mk1r <- constrain(lik, lambda1 ~ lambda0, mu1 ~ mu0, q10~0)
fit.mk1r  <- find.mle(lik.mk1r , p[argnames(lik.mk1r )])
fit.mk1r $lnLik
round(coef(fit.mk1r ),5)


#For the highest AIC score, to account for parameter uncertainty, using MCMC sampling

#######################################################
#ASR (ancestral state reconstruction)


#BiSSE model
st <- asr.marginal(lik,coef(fit))
plot(phy1, main = "Marginal reconstruction of ancestral states under the best-fit 6-parameter BiSSE model", show.tip.label = FALSE, show.node.label = FALSE)
#add to the internal nodes
nodelabels(thermo = t(st), piecol = 1:2, cex = 0.5)


#MK2
likmk <- make.mk2(phy1,states)
fitmk <- find.mle(likmk,pars[5:6],method = "subplex")
round(coef(fitmk),5)
stmk <- asr.marginal(likmk,coef(fitmk))
st <- asr.marginal(lik,coef(fitmk))
plot(phy1, main = "Marginal reconstruction of ancestral states under the best-fit 6-parameter BiSSE model", show.tip.label = FALSE, show.node.label = FALSE)
nodelabels(thermo = t(st), piecol = 1:2, cex = 0.5)


