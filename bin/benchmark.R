#!/usr/bin/env Rscript

library(argparse)
library(corpcor)
library(data.table)
library(ggplot2)
library(ggpubr)
library(MASS)

# parse arguments
parser = ArgumentParser(description='Benchmark shrinkage and imputation')
parser$add_argument('-i1', '--input1', type='character', help="Input sample count data - truth")
parser$add_argument('-i2', '--input2', type='character', help="Input sample count data - bench")
parser$add_argument('-o', '--output', type='character', help="Output mse")
parser$add_argument('--transf', type='character', help='clr or alr')
parser = parser$parse_args()


# read input files ----------------------------------------------------------------------

count_truth = fread(parser$input1)
count_truth = as.matrix(count_truth)
count_bench = fread(parser$input2)
count_bench = as.matrix(count_bench)

print(dim(count_truth))
print(dim(count_bench))

# functions -----------------------------------------------------------------------------
noShrink=function(M){
    Cov=cov(M) 
    Cor=cov2cor(Cov)
    PC=cor2pcor(Cov)
    
    C=list(Cov,Cor,PC)
    names(C)=c('cov', 'corr', 'pcor')
    return(C)
}
aShrink=function(M){
    Cov=cov.shrink(M,verbose=FALSE)  
    Cor=cov2cor(Cov)
    PC=cor2pcor(Cov)
    
    C=list(Cov,Cor,PC)
    names(C)=c('cov.shrink', 'corr.shrink', 'pcor.shrink')
    return(C)
}
bShrink=function(M,intype,outtype="clr"){
    if (intype=="alr"){
        D=ncol(M)+1
    } else{
        D=ncol(M)
    }
    N=nrow(M)
    
    if (intype=="rco"){
        
        print("Warning: estimating compositions from counts via shrinkage")
        P=matrix(0,N,D)
        for (i in 1:nrow(M)){
            P[i,]=shrink(M[i,])
        }
        B=log(P)
    } else if (intype=="alr"){
        
        P=exp(M)/(1+apply(exp(M),1,sum))
        P=cbind(P,1-apply(P,1,sum))
        B=log(P)
        
    } else if (intype=="clr"){
        f=log(apply(exp(M),1,sum)) #backtransform to log P:
        B=M-f
    }
    
    Cb=cov.shrink(B,verbose=FALSE)  
    if (outtype=="alr"){
        F=cbind(diag(rep(1,D-1)),rep(-1,D-1))
        Cov=F%*%Cb%*%t(F)
    } else if (outtype=="clr"){
        G=diag(rep(1,D))-matrix(1/D,D,D)
        Cov=G%*%Cb%*%G
    } else{
        die("outtype unsupported")
    }
    Cor=cov2cor(Cov)
    PC=cor2pcor(Cov)
    
    
    C=list(Cov,Cor,PC)
    names(C)=c('cov.basis', 'corr.basis', 'pcor.basis')
    return(C)
    
}
shrink <- function(n){
        N=sum(n)
        q=n/N       
        D=length(n)
        T=rep(1/D,D)
        
        l=(1-sum(q^2))/((N-1)*sum((T-q)^2))
        if (l<0){
            l=0
        }
        if(l>1){
            l=1
        }
        qs=l*T+(1-l)*q
        
        return(qs)
}
calculate_mse <- function(true, pred){
    l = list()
    for (i in c('cov', 'cov.shrink', 'cov.basis')){
        l[[i]] = mse(true[['cov']], pred[[i]])
    }
    for (i in c('corr', 'corr.shrink', 'corr.basis')){
        l[[i]] = mse(true[['corr']], pred[[i]])
    }
    for (i in c('pcor', 'pcor.shrink', 'pcor.basis')){
        l[[i]] = mse(true[['pcor']], pred[[i]])
    }
    return(l)
}
mse <- function(true, pred){
    true = true[lower.tri(true)]
    pred = pred[lower.tri(pred)]
    return( mean((true-pred)^2) )
}
get_cor_name <- function(name){
    if (name %in% c('cov', 'cov.shrink', 'cov.basis')){
        'covariance'
    }else if (name %in% c('corr', 'corr.shrink', 'corr.basis')){
        'correlation'
    }else if (name %in% c('pcor', 'pcor.shrink', 'pcor.basis')){
        'partial correlation'
    }
}
get_shrink_name <- function(name){
    if (name %in% c('cov', 'corr', 'pcor')){
        'no shrinkage'
    }else if (name %in% c('cov.shrink', 'corr.shrink', 'pcor.shrink')){
        'shrink std'
    }else if (name %in% c('cov.basis', 'corr.basis', 'pcor.basis')){
        'shrink basis'
    }
}
get_clr <- function(count){
    lgm = apply(log(count),1,mean)
    clr = log(count) - lgm
    return(clr)
}
get_alr <- function(count){
    alr = log(count)-log(count[,ncol(count)])
    alr = alr[,-ncol(count)]
    return(alr)
}
get_lr <- function(count, transf=c('clr','alr')){
    if (transf == 'clr'){
        lr = get_clr(count)
    }else if (transf == 'alr'){
        lr = get_alr(count)
    }
    return(lr)
}

# benchmark -----------------------------------------------------------------------------

# calculate log ratios
lr_truth = get_lr(count_truth, transf=parser$transf)
lr_bench = get_lr(count_bench, transf=parser$transf)

# calculate correlations
bcor_truth = noShrink(lr_truth)
bcor_bench = noShrink(lr_bench)
bcor_bench = c(bcor_bench, aShrink(lr_bench))
bcor_bench = c(bcor_bench, bShrink(lr_bench, parser$transf, parser$transf))

# calculate mse
tmp = calculate_mse(bcor_truth, bcor_bench)

# organize and save output data frame
df = data.frame(mse=numeric(), corr=character(), shrink=character())
for (method in names(tmp)){
    row = data.frame(mse=tmp[[method]], corr=get_cor_name(method), shrink=get_shrink_name(method))
    df  = rbind(df, row)
}
write.csv(df, parser$output, row.names=F, quote=F)