#!/usr/bin/env Rscript

library(argparse)
library(corpcor)
library(data.table)
library(ggplot2)
library(ggpubr)
library(MASS)
library(propr)

# parse arguments
parser = ArgumentParser(description='Benchmark shrinkage')
parser$add_argument('-i', '--input', type='character', help="Input count data")
parser$add_argument('-o', '--outdir', type='character', help="Output directory")
parser$add_argument('--ngene', type='integer', help='Number of genes considered')
parser$add_argument('--ncell', type='integer', help='Number of cells considered')
parser$add_argument('--nsamp', type='integer', default=200, help='Number of samples')
parser$add_argument('--simulation', action='store_true', default=FALSE, help='If false (default), only sample cells. If true, simulate cells.')
parser = parser$parse_args()


# read input files ----------------------------------------------------------------------

count = fread(parser$input)
count = as.matrix(count)

print(dim(count))
print(sum(count==0))

if (!file.exists(parser$outdir)){
    dir.create(parser$outdir, recursive = TRUE)
}

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
simulate_lr <- function(count, ncell){
    alr= get_alr(count)
    S  = cov(alr) 
    m  = apply(alr, 2, mean)
    al = mvrnorm(n = ncell, mu=m, Sigma=S)  # alr samples
    P  = exp(al)/(1+apply(exp(al),1,sum))
    P  = cbind(P,1-apply(P,1,sum))          # sampled compositions
    lg = apply(log(P),1,mean)
    cl = log(P)-lg  
    return(list(cl, al))
}
sample_cells <- function(count, ncell){
    csamp = sample(1:nrow(count), ncell)
    count = count[csamp,]
    return(count)
}


# benchmark -----------------------------------------------------------------------------
df = list()
df[['clr']] = data.frame(mse=numeric(), corr=character(), shrink=character())
df[['alr']] = data.frame(mse=numeric(), corr=character(), shrink=character())
for (i in 1:parser$nsamp){
    
    # sample genes
    gsamp  = sample(1:ncol(count), parser$ngene)
    dat    = count[,gsamp]

    # clr/alr ground truth
    lr_truth = list()
    lr_truth[['clr']] = get_clr(dat)
    lr_truth[['alr']] = get_alr(dat)

    # simulate or sample cells, and compute clr and alr
    lr_bench = list()
    if (parser$simulation){
        lr = simulate_lr(dat, parser$ncell)
        lr_bench[['clr']] = lr[[1]]
        lr_bench[['alr']] = lr[[2]]
    }else{
        dat = sample_cells(dat, parser$ncell)
        lr_bench[['clr']] = get_clr(dat)
        lr_bench[['alr']] = get_alr(dat)
    }

    bcor_truth = list(); bcor_bench = list()
    for (transf in c('clr', 'alr')){
        # compute correlation
        bcor_truth[[transf]] = noShrink(lr_truth[[transf]])
        bcor_bench[[transf]] = noShrink(lr_bench[[transf]])
        bcor_bench[[transf]] = c(bcor_bench[[transf]], aShrink(lr_bench[[transf]]))
        bcor_bench[[transf]] = c(bcor_bench[[transf]], bShrink(lr_bench[[transf]], transf, transf))

        # calculate mse
        tmp = calculate_mse(bcor_truth[[transf]], bcor_bench[[transf]])
        for (method in names(tmp)){
            row = data.frame(mse=tmp[[method]], corr=get_cor_name(method), shrink=get_shrink_name(method))
            df[[transf]] = rbind(df[[transf]], row)
        }
    }
}
for (transf in c('clr','alr')){
    out = paste0(parser$outdir, '/mse_', parser$ncell, '_', parser$ngene, '_', parser$nsamp, '_', transf, '.csv')
    write.csv(df[[transf]],out, row.names=F, quote=F)
}