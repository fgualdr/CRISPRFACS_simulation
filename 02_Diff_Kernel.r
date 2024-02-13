#!/usr/bin/env Rscript 

# suppress verbosity when loading libraries:
options(warn=-1)
library(optparse)
library(dplyr)
library(tidyr)
library(GeneralNormalizer)
library(ComplexHeatmap)
library(tidyr)
library(tibble)
library(ggplot2)
library(data.table)
library(MASS)
library(BiocParallel)
library(snow)

# # # Custom functions:
`%!in%` = Negate(`%in%`)
# Function to perform the orthogonal projection of a point onto a line:

orthogonal_projection <- function(p, l1, l2) {
  dir <- l2 - l1
  numerator <- sum((p - l1) * dir)
  denominator <- sum(dir^2)
  lambda <- numerator / denominator
  projection <- l1 + lambda * dir
  transformation <- t(t(projection) - p)
  return(list(projection=projection, transformation=transformation))
}

rows_in_lim <- function(m,ld,lu){
        m[m<ld] = 0
        m[m>lu] = 0
        m[m!=0] = 1
        return(rownames(m)[which(rowSums(m) == ncol(m))])
}
split_tibble <- function(tibble, col = 'col') tibble %>% split(., .[, col])
# For each results we want to retrieve the number of hits per taget when varying the Log2FC. 
# There are 10 potential gRNA (hits) per target listed as "target name"_"n° of the gRNA"
# We want to retrieve the number of hits per target when varying the Log2FC.
# We will use the function "get_hits_per_target" to do so.
# For each L2FC we will any hit with a L2FC > L2FC.
# We will then count the number of hits per target i.e. remove _[0-9]* from the target name and count the number of unique target name.
get_hits_per_target <- function(res,log2fc){
    # res is the DESEQ2 results having the hits as raw names
    # log2fc is the log2fc threshold
    # the res has column log2FoldChange and padj
    # Subset the res:
    df_subset <- res[which(abs(as.numeric(as.character(res$log2FoldChange))) > log2fc),]
    # If there is no hit we return NA and finish:
    if(nrow(df_subset) != 0){   
        # Remove the _[0-9]* from the target name:
        df_subset$target_name <- gsub("_[0-9]*","",rownames(df_subset))
        # we aggregate the df_subset by target name and count the number of hits i.e raws:
        df_subset <- aggregate(df_subset$target_name,by=list(df_subset$target_name),length)
        names(df_subset) <- c("target_name","n_hits")
        # Specify the log2fc:
        df_subset$log2fc <- log2fc
        return(df_subset)
    }
}

# 

random_gen <- function(l){
        # define internally get_hits_per_target
        get_hits_per_target <- function(res,log2fc){
            df_subset <- res[which(abs(as.numeric(as.character(res$log2FoldChange))) > log2fc),]
            # If there is no hit we return NA and finish:
            if(nrow(df_subset) != 0){   
                # Remove the _[0-9]* from the target name:
                df_subset$target_name <- gsub("_[0-9]*","",rownames(df_subset))
                # we aggregate the df_subset by target name and count the number of hits i.e raws:
                df_subset <- aggregate(df_subset$target_name,by=list(df_subset$target_name),length)
                names(df_subset) <- c("target_name","n_hits")
                # Specify the log2fc:
                df_subset$log2fc <- log2fc
                return(df_subset)
            }
        }

        res_up_r = l[["res_up_r"]]
        res_down_r = l[["res_down_r"]]
        max_oligo_per_target = l[["max_oligo_per_target"]]
        l2c_cuts = l[["l2c_cuts"]]

        # how many groups do we need:
        real_names <- c(rownames(res_up_r),rownames(res_down_r))

        d = round(  (length(real_names) / max_oligo_per_target)+0.5,0)
        
        artifical_names = rep(paste0("ScraTarget",1:d),each=max_oligo_per_target)
        artifical_names <- split(artifical_names,artifical_names)
        artifical_names <- lapply(artifical_names,function(x){
            x <- paste0(x,"_",1:length(x))
            return(x)
        })
        artifical_names <- unlist(artifical_names)
        artifical_names <- artifical_names[sample(1:length(artifical_names))]
        artifical_names <- artifical_names[1:length(real_names)]
        # shuffle the real names 
        real_names <- real_names[sample(1:length(real_names))]
        names(artifical_names) <- real_names
        # Shuffle the rows in the res_up_r:
        res_shuffle_up = res_up_r[names(artifical_names[which(names(artifical_names) %in% rownames(res_up_r))]),]
        rownames(res_shuffle_up) <- artifical_names[which(names(artifical_names) %in% rownames(res_up_r))]
        # Shuffle the rows in the res_down_r:
        res_shuffle_down = res_down_r[names(artifical_names[which(names(artifical_names) %in% rownames(res_down_r))]),]
        rownames(res_shuffle_down) <- artifical_names[which(names(artifical_names) %in% rownames(res_down_r))]
        # Get the hits per target:
        res_shuffle_up <- lapply(l2c_cuts,function(x) get_hits_per_target(res_shuffle_up,x))
        res_shuffle_up <- do.call(rbind,res_shuffle_up)
        res_shuffle_down <- lapply(l2c_cuts,function(x) get_hits_per_target(res_shuffle_down,x))
        res_shuffle_down <- do.call(rbind,res_shuffle_down)
        res_shuffle_down$log2fc <- -1*res_shuffle_down$log2fc
        # 
        xy <- rbind(res_shuffle_up,res_shuffle_down)
        rownames(xy) = 1:dim(xy)[1]
        return(xy)
    }

# mass_kde2_parallel is a function that save to file all the MASS::kde2d results with a unique name in the working directory:
# the input is a list with the following elements:
# - xy: the xy data frame
# - n_x: the number of x grid points
# - n_y: the number of y grid points
# - n_v: the number of unique values in the y axis
# the output will be NULL as it save to file instead of adding to the list - this is because we 
# use bplapply and we want to save the results to file and save memory usage:
# it will need an increamental number so the list will include the "id" slot too:
mass_kde2_parallel <- function(l){
        h = c(MASS::bandwidth.nrd(l[["xy"]]$n_hits), MASS::bandwidth.nrd(l[["xy"]]$log2fc))
        # bandwidth.nrd returns the estimated bandwidth of a univariate distribution using a modified version of Silverman's rule of thumb.
        # when dealing with discrete data the bandwidth is 1:
        # as n_hits are discrete entity the badnwith will always be 1 so we set it to 1:
        h[1] = 1
        # on the othet hand if the log2fc has estimated bandwidth of 0 we through an error:
        if(h[2] == 0){
            stop("ERROR: bandwidth for log2fc is 0")
        }
        kde2d <- MASS::kde2d( 
            l[["xy"]]$n_hits, 
            l[["xy"]]$log2fc, 
            n=c( l[["n_x"]]+2, l[["n_y"]] ), 
            h = h,
            lims = c(
                    0,l[["n_x"]]+1, # for the x axis we consider the number of hits
                    c(min(l[["n_v"]])*2, max(l[["n_v"]])*2) # for the y axis we consider the log2fc
                    )
            )
        saveRDS(kde2d,paste0("kde2d_",Sys.Date(),"_",Sys.time(),"_",Sys.getpid(),"_",l[["id"]],".rds"))
        # delete the kde2d object to free up memory:
        rm(kde2d)
        return(NULL)
    }

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # Main script:
# Parse the options with optparse:

setwd("/hpcnfs/data/GN2/fgualdrini/Master_batch_scripts/CRISPR_TOOLS_BANCH/Simulations_KERNEL/")


opt <- list()
opt$file <- "/hpcnfs/data/GN2/fgualdrini/Master_batch_scripts/CRISPR_TOOLS_BANCH/Simulations_KERNEL/groundTruth_reps_ALL.txt"
opt$processes <- 3
opt$skip_outliers <- TRUE
opt$target <- "LOW"
opt$control <- "HIGH"
opt$basename_save <- gsub("\\.txt","",opt$file)
opt$use_nt_scramble <- TRUE
opt$scramble_id <- "scrambled"
opt$grna_id_col <- 1


# Parse the cpus:
cpus=as.numeric(as.character(opt$processes)) - 1
cat("Number of cpus:",cpus,"\n")
if(cpus <= 1){cpus = 2} ## at least 2 cpus
if(!is.null(cpus)){
    param <- BiocParallel::SnowParam(workers=cpus,tasks=0,stop.on.error=TRUE,progressbar = TRUE,type="SOCK")
    param_serial <- BiocParallel::SerialParam(stop.on.error = TRUE,progressbar = TRUE)
}else{
    param <- NULL
    param_serial <- BiocParallel::SerialParam(stop.on.error = TRUE,progressbar = TRUE)
}

# define target and control:
opt$target = gsub("\\.bam$","",opt$target)
opt$control = gsub("\\.bam$","",opt$control)
target=opt$target
control=opt$control
check_levels <- c(target,control)
check_levels <- strsplit(check_levels, "_")
check_levels <- do.call(rbind, check_levels)
# this will be a matrix with two columns:
check_levels <- check_levels[, apply(check_levels, 2, function(x) length(unique(x)) > 1)]
# re-assign to target and control:
target <- check_levels[1]
control <- check_levels[2]

names(check_levels) <- c(opt$target,opt$control)

cat("Target:",target,"\n")
cat("Control:",control,"\n")

# Saving folder will always be the workdirector:
save_folder=paste0(getwd(),paste0("/CrisprRes",target,"_vs_",control,"/"))
dir.create(save_folder)

# We read the input file:
x = read.delim(opt$file,sep="\t",header=TRUE,row.names=NULL,stringsAsFactors=FALSE)
rownames(x) = x[,"elements_rep1"]
# select poisson:
x = x[,grep("nbinom",colnames(x))]
# first two columns are start and end:
samples.vec <- colnames(x)
samples.vec <- samples.vec[!grepl("start|end|gene",tolower(samples.vec))]

# check if replicates are present, they should be ending in R1 or R2 etc.. or r1 or r2 etc..

deg_design = do.call(rbind,strsplit(samples.vec,"_"))
# remove columns with one value:
deg_design = deg_design[,apply(deg_design,2,function(x){length(unique(x)) > 1})]
colnames(deg_design) = c("Sample_Condition","Sample_Replicate")
deg_design = as.data.frame(deg_design,stringsAsFactors=FALSE)
rownames(deg_design) = samples.vec
deg_design$Sample_ID = rownames(deg_design)

x = x[,rownames(deg_design)]
cs = apply(x,2,function(y){!all(is.na(y))})
x = x[,which(cs)]
cs = apply(x,2,function(y){!all(y==0)})
x = x[,which(cs)]
rn = rownames(x)
x = as.data.frame(apply(x,2,function(y){
    return(as.numeric(as.character(y)))
}),stringsAsFactors=FALSE)
rownames(x) = rn

##
if(length(grep("presort",tolower(colnames(x)))) != 0){
    ref = colnames(x)[grep("presort",tolower(colnames(x)))]
}else{
    ref = colnames(x)[1]
}

# We normalize the data:
ori_mat = round(x,0)
Result = RunNorm(ori_mat,deg_design,fix_reference=deg_design$Sample_ID,row_name_index=1,saving_path=save_folder,n_pop=1,n_pop_reference=1,BiocParam=param)
colData = Result$scaling_factors
mat = Result$norm_mat

deg_design = deg_design[order(deg_design$Sample_Replicate),]
# Norm on PLASMID BY Replicate:


# if we have replicates with generate Log2FC from DESEQ2 fitting:
if(length(unique(deg_design$Sample_Replicate)) > 1){

    rn = lapply(split_tibble(colData,"Sample_Condition"),function(x){
        return( rows_in_lim( ori_mat[,rownames(x)],10,2000) )
    })
    rn =unique(unlist(rn))
    mat = round(ori_mat[rn,rownames(colData)],0)
    mat = mat + abs(min(mat))
    mat = mat[which(rowMeans(mat)<2000 & rowMeans(mat)>10),]
    formula = ~ Sample_Replicate + Sample_Condition   
    dds_ex <- DESeq2::DESeqDataSetFromMatrix(countData = ori_mat[rownames(mat),], colData = colData, design = formula)
    dds_ex <- DESeq2::estimateSizeFactors(dds_ex)

    dds_ex <- DESeq2::estimateDispersions(dds_ex,fitType="local",maxit=10000) # parametric local mean
    dds_ex <- DESeq2::nbinomWaldTest(dds_ex,maxit = 10000)       
    vsd <- DESeq2::vst(dds_ex, blind=FALSE)
    cat("Run DESEQ2:\n")
    pdf(paste0(save_folder,"/MAplot.pdf"))
        DESeq2::plotDispEsts(dds_ex)
    dev.off()
    mat <- SummarizedExperiment::assay(vsd)
    mat <- limma::removeBatchEffect(mat, vsd$Batch)
    SummarizedExperiment::assay(vsd) <- mat
    counts_batch_corrected <- SummarizedExperiment::assay(vsd)

    # We save the data by contrast:
    cond = unique(deg_design$Sample_Condition)
    cond = cond[!grepl("PLASMID",cond)]
    cat("All conditions are : ",paste0(cond,collapse=","),"\n")

    save_folder_sub_stat = paste0(save_folder,"/Testing_stat/")
    dir.create(save_folder_sub_stat)
    # Up until here we can free up space. We remove all object but contrasts, dds_ex , save_folder_sub_stat which are define above and use in the loop below:
    # we keep also the defined functions: orthogonal_projection, get_hits_per_target, random_gen, mass_kde2_parallel
    # rm(list=setdiff(ls(),c("contrasts","dds_ex","save_folder_sub_stat","orthogonal_projection","get_hits_per_target","random_gen","mass_kde2_parallel","%!in%","param","target","control","opt")))
    cat("Run DESEQ2 on contrast:",paste0(c(target,control),collapse=","),"\n")
    target_cond = unique(colData[grep(tolower(target),tolower(colData$Sample_Condition)),]$Sample_Condition)
    control_cond = unique(colData[grep(tolower(control),tolower(colData$Sample_Condition)),]$Sample_Condition)
    res = DESeq2::results(dds_ex, contrast=c("Sample_Condition",target_cond,control_cond))
    nm <- paste0(save_folder_sub_stat,"/",target,"_vs_",control,"_RESULTS.txt")
    pdf(paste0(save_folder_sub_stat,"/",target,"_vs_",control,"_RESULTS.pdf"))
            DESeq2::plotMA(res,alpha = 0.01, main = "", xlab = "mean of normalized counts", MLE = FALSE)
    dev.off()

    write.table(res,paste0(save_folder_sub_stat,"/Multivariate_",target,"_vs_",control,"_DESEQ2_RES.txt"),sep="\t",col.names=NA)
    res = as.data.frame(res,stringsAsFactors=FALSE)
}else{
    save_folder_sub_stat = paste0(save_folder,"/Testing_stat/")
    dir.create(save_folder_sub_stat)
    # here we have no replicates so the formula is just the condition:
    formula = ~ Sample_Condition
    dds_ex <- DESeq2::DESeqDataSetFromMatrix(countData = ori_mat, colData = colData, design = formula)
    dds_ex <- DESeq2::estimateSizeFactors(dds_ex)
    SummarizedExperiment::colData(dds_ex)$sizeFactor <- colData$scaling
    # applying the vst() and then you can just calculate the simple difference between the two samples. 
    # This will give you a log2 fold change that is moderated, such that small counts do not produce large LFC (avoiding having to pick a pseudocount). 
    vst_dat <- DESeq2::vst(dds_ex, blind=TRUE)
    # get the vst counts:
    counts_vst <- SummarizedExperiment::assay(vst_dat)
    l2fc <- counts_vst[,names(check_levels)[1]] - counts_vst[,names(check_levels)[2]]
    counts_vst <- cbind(counts_vst,l2fc)
    counts_vst <- as.data.frame(counts_vst,stringsAsFactors=FALSE)
    colnames(counts_vst ) <- c(colnames(counts_vst)[1:ncol(counts_vst)-1],"log2FoldChange")
    res = counts_vst

}

rownames(res) = gsub("_","",rownames(res))
rownames(res) = gsub("gRNA","gRNA_",rownames(res))
# We compute this for a range of cut offs:
l2c_cuts <- seq(min(abs(res$log2FoldChange)),max(abs(res$log2FoldChange)),by=0.001)


# We separate up and down regulated effects:
res_up = res[which(res$log2FoldChange > 0),]
res_down = res[which(res$log2FoldChange < 0),]
res_down$log2FoldChange <- -1*res_down$log2FoldChange

if(opt$use_nt_scramble){
    # we separate the "scaramble|random" targets from the rest:
    res_up_r = res_up[grepl("^rand|^scram|^non_target",rownames(res_up)),]
    res_up = res_up[!grepl("^rand|^scram|^non_target",rownames(res_up)),]
    res_down_r = res_down[grepl("^rand|^scram|^non_target",rownames(res_down)),]
    res_down = res_down[!grepl("^rand|^scram|^non_target",rownames(res_down)),]
}else{
    # if no scrrambled to be used we randomly all the gRNA by grouping in random sets
    res_up_r = res_up
    res_up = res_up[!grepl("^rand|^scram|^non_target",rownames(res_up)),]
    res_down_r = res_down
    res_down = res_down[!grepl("^rand|^scram|^non_target",rownames(res_down)),]
}


# for the target we get the number of hits per target for each l2fc:
res_up <- lapply(l2c_cuts,function(x) get_hits_per_target(res_up,x))
res_up <- do.call(rbind,res_up)
res_down <- lapply(l2c_cuts,function(x) get_hits_per_target(res_down,x))
res_down <- do.call(rbind,res_down)
# as we considered the abs of the res_down their associated l2fc are negative:
res_down$log2fc <- -1*res_down$log2fc

# Ok given the res_up and res_down we compute the bivariate distribution of the number of hits per target for each pair of log2fc:
# We aggregate by n_hits and log2fc and count the number of elements:
xy <- rbind(res_up,res_down)
rownames(xy) = 1:dim(xy)[1]

# Ok now we do the same for the random hits - but we bootstrap the number of hits per target:
# basically we perform 1000 iteractions where we group randomly the scamble into groups having the same number of targets as the real data:
# and each time we compute the number of hits per target for each l2fc:
# We do this for the up and down regulated target so for the res_up_r and res_down_r:
# the randomization should be done without replacement:
max_oligo_per_target = max(c(res_up$n_hits,res_down$n_hits))

# to make it consistent with the design and data - we can consider that for one target we might have that out of the total oligos some will go up and some will go down:
# so the sampling should go together - so that names are spread across the up and down regulated targets:

# to use bplapply we make a list object containing the res_up_r and res_down_r 1000 times:
# so will be a named list:
l1 = list( 
    "res_up_r"=res_up_r,
    "res_down_r"=res_down_r,
    "max_oligo_per_target"=max_oligo_per_target,
    "l2c_cuts"=l2c_cuts
            )
# repeat this list 1000 times so we have a list of list:
l_rep = rep(list(l1),10)
# select top 10 elelemts:
list_r <-  BiocParallel::bplapply(l_rep,random_gen,BPPARAM = param)

# for each of the 1000 list we compute the MASS::kde2d
# for consistency we need to set "n" and "lim" coherently to the xy of the real data:
n_grid_x <- length(unique(xy$n_hits))
# for the y grid we consider the total number of all oligos in the real data:
n_grid_y <- dim(res)[1] # like saying each oligo has a different log2fc
# we combine the list_r with n_grid_x and n_grid_y so to have a list of list each element with the "xy" of the list_r and the "n_grid_x" and "n_grid_y":
list_r_m <- lapply(1:length(list_r),function(i){
    l = list()
    l[["id"]] = i
    l[["xy"]] = list_r[[i]]
    l[["n_x"]] = n_grid_x
    l[["n_y"]] = n_grid_y
    l[["n_v"]] = unique(xy$log2fc)
    return(l)
})

# to use the multiple scrambled bivariate we will:
# 1) For each of the replicate datasets, apply the MASS::kde2d function to estimate the bivariate kernel density:
setwd(save_folder_sub_stat)
# generate a temporary folder:
dir.create("tmp")
setwd("tmp")
list_r_kde <- lapply(list_r_m,mass_kde2_parallel)

# list_r_kde we list the files generated in tmp:
list_r_kde_f <- list.files(pattern = "^kde2d_")
if(length(list_r_kde_f)==0){
    # show possible errors in list_r_kde
    print(list_r_kde)
    cat(length(list_r_kde_f),"\n")
}

# cat("Get the maximum Machine number")
cat("Get the maximum Machine number \n")
list_r_kde_norm <- lapply(list_r_kde_f,function(x){
        # read the kde2d:
        x = readRDS(x)
        x$z <- x$z + (.Machine$double.xmin)
        xlim <- range(x$x)
        ylim <- range(x$y)
        cell_size <- (diff(xlim) / length(x$x)) * (diff(ylim) / length(x$y))
        # Compute the normalization function:
        norm <- sum(x$z) * cell_size  # normalizing constant    
        x$z = (x$z*cell_size)/norm 
        return(x)
})

# remove the tmp folder:
setwd(save_folder_sub_stat)
unlink(paste0(save_folder_sub_stat,"/tmp/"), recursive = TRUE, force = TRUE)

cat("Calculate the mean and standard deviation of the kernel density estimates\n")
# 2) Calculate the mean and standard deviation of the kernel density estimates.
# the "z" is a matrix we want to preserve the dimensionality of the matrix so we use sapply. sapply will get the 
mean_density <- Reduce("+", lapply(list_r_kde_norm, "[[", "z")) / length(list_r_kde_norm)
sd_density <- Reduce("+", lapply(list_r_kde_norm, function(x) (x[["z"]] - mean_density)^2) ) / (length(list_r_kde_norm) - 1)

# 3) Calculate the z-score for each combination of n_hits and log2fc
# We use the kde2d function from the kde2d package to compute the bivariate distribution:
kde2d_mean <- list()
kde2d_mean[["x"]] <- list_r_kde_norm[[1]][["x"]]
kde2d_mean[["y"]] <- list_r_kde_norm[[1]][["y"]]
kde2d_mean[["z"]] <- mean_density
xlim <- range(kde2d_mean$x)
ylim <- range(kde2d_mean$y)
cell_size <- (diff(xlim) / length(kde2d_mean$x)) * (diff(ylim) / length(kde2d_mean$y))
# Compute the normalization function:
norm <- sum(kde2d_mean$z) * cell_size  # normalizing constant    
kde2d_mean$z = (kde2d_mean$z*cell_size)/norm 
# Plot the bivariate distribution:
cat("Plot the bivariate distribution\n")
pdf(paste0(save_folder_sub_stat,"/Multivariate_",target,"_vs_",control,"_KDE_RANDOM.pdf"))
    par( mfrow= c(2,2) )
    filled.contour( 
                    kde2d_mean,
                    color.palette=colorRampPalette(c('white','blue','yellow','red','darkred')),
                    ylim = range(c(-0.5,0.5), finite = TRUE),
                    zlim = range(kde2d_mean$z, finite = TRUE)
                    )
    persp(kde2d_mean, phi = 40, theta = 70, d = 5,col="#B8EBF2",border=NULL)
dev.off()

png(paste0(save_folder_sub_stat,"/Multivariate_",target,"_vs_",control,"_KDE_RANDOM.png"))
    par( mfrow= c(2,2) )
    filled.contour( 
                    kde2d_mean,
                    color.palette=colorRampPalette(c('white','blue','yellow','red','darkred')),
                    ylim = range(c(-0.5,0.5), finite = TRUE),
                    zlim = range(kde2d_mean$z, finite = TRUE)
                    )
    #persp(kde2d_mean, phi = 40, theta = 70, d = 5,col="#B8EBF2",border=NULL)
dev.off()


# Define the function such as that given a pair of hits and log2fc it returns the probability after having interpolated the density function:
kde2dAkima = kde2d_mean$z
colnames(kde2dAkima) = kde2d_mean$y
rownames(kde2dAkima) = kde2d_mean$x
kde2dAkima = reshape2::melt(kde2dAkima)
colnames(kde2dAkima) = c("n_hits","log2fc","density")

# For evey point in the xy data frame we compute the probability of the point given the density function:
# The probability is reported as 1-cdf and in particular we compute the cdf as the sum of the density from the point to the upper/lower limit

Res = xy

# up
Res_up = Res[Res$log2fc >0,]
Res_up = Res_up[order(Res_up$log2fc,decreasing=TRUE),]
Res_up = Res_up %>% distinct(target_name,n_hits, .keep_all = T)
prob = lapply(1:nrow(Res_up),function(i){
    if (Res_up$log2fc[i] < 0){
        # Get the x_values and y_values:
        x_idx <- which(kde2dAkima$n_hits == Res_up$n_hits[i])
        # as the l2fc is less than 0 we need to sum the density from the point to the lower limit meanin any point below the one selected
        y_idx = which(kde2dAkima$log2fc <= Res_up$log2fc[i])
        # Subset the density values:
        summed_p = sum(kde2dAkima$density[intersect(x_idx,y_idx)])
        return(summed_p)
    } else {
        # HEre we get cdf as the sum of the density from the point to the upper limit
        x_idx <- which(kde2dAkima$n_hits == Res_up$n_hits[i])
        y_idx = which(kde2dAkima$log2fc >= Res_up$log2fc[i])
        # Su up all the values 
        summed_p = sum(kde2dAkima$density[intersect(x_idx,y_idx)])
        return(summed_p)
    }
})
Res_up$Bivariate_prob = unlist(prob)

# down
Res_down = Res[Res$log2fc <0,]
Res_down = Res_down[order(Res_down$log2fc,decreasing=FALSE),]
Res_down = Res_down %>% distinct(target_name,n_hits, .keep_all = T)
prob = lapply(1:nrow(Res_down),function(i){
    if (Res_down$log2fc[i] < 0){
        # Get the x_values and y_values:
        x_idx <- which(kde2dAkima$n_hits == Res_down$n_hits[i])
        # as the l2fc is less than 0 we need to sum the density from the point to the lower limit meanin any point below the one selected
        y_idx = which(kde2dAkima$log2fc <= Res_down$log2fc[i])
        # Subset the density values:
        summed_p = sum(kde2dAkima$density[intersect(x_idx,y_idx)])
        return(summed_p)
    } else {
        # HEre we get cdf as the sum of the density from the point to the upper limit
        x_idx <- which(kde2dAkima$n_hits == Res_down$n_hits[i])
        y_idx = which(kde2dAkima$log2fc >= Res_down$log2fc[i])
        # Su up all the values 
        summed_p = sum(kde2dAkima$density[intersect(x_idx,y_idx)])
        return(summed_p)
    }
})
Res_down$Bivariate_prob = unlist(prob)
Res = rbind(Res_up,Res_down)
Res$adjp_fdr_bivariate = p.adjust(Res$Bivariate_prob, method = "fdr")

write.table(Res,paste0(save_folder_sub_stat,"/Multivariate_",target,"_vs_",control,"_prob.txt"),sep="\t",col.names=NA)


##################################
# Plot Global result:
##################################

Res_by_target = split(Res,Res$target_name)
v_hits = as.character(1:10)

# UP
up_fc = lapply(Res_by_target,function(x){
    y = x[x$log2fc > 0,"log2fc"]
    names(y) = as.character(x[x$log2fc > 0,"n_hits"])
    y = y[v_hits]
    names(y) = v_hits
    return(y)
})
up_fc = as.data.frame(do.call(cbind,up_fc),stringsAsFactors=FALSE)
up_fc[is.na(up_fc)] = 0

up_pval = lapply(Res_by_target,function(x){
    y = -log10(x[x$log2fc > 0,"adjp_fdr_bivariate"])
    names(y) = as.character(x[x$log2fc > 0,"n_hits"])
    y = y[v_hits]
    names(y) = v_hits
    return(y)
})
up_pval = as.data.frame(do.call(cbind,up_pval),stringsAsFactors=FALSE)
up_pval[is.na(up_pval)] = 0

# DOWN
down_fc = lapply(Res_by_target,function(x){
    y = x[x$log2fc < 0,"log2fc"]
    names(y) = as.character(x[x$log2fc < 0,"n_hits"])
    y = y[v_hits]
    names(y) = v_hits
    return(y)
})
down_fc = as.data.frame(do.call(cbind,down_fc),stringsAsFactors=FALSE)
down_fc[is.na(down_fc)] = 0
down_pval = lapply(Res_by_target,function(x){
    y = -log10(x[x$log2fc < 0,"adjp_fdr_bivariate"])
    names(y) = as.character(x[x$log2fc < 0,"n_hits"])
    y = y[v_hits]
    names(y) = v_hits
    return(y)

})
down_pval = as.data.frame(do.call(cbind,down_pval),stringsAsFactors=FALSE)
down_pval[is.na(down_pval)] = 0


# Criteria at least 2 oligos with log2fc >= log2(1.1) and pval <= 0.01

w_tuf = apply(up_fc,2,function(x){
    return(length(which(abs(x) >= log2(1.20))))
    })
w_tuf = w_tuf[which(w_tuf>=2)]
w_tup = apply(up_pval[2:nrow(up_pval),],2,function(x){max(abs(x))})
w_tup = w_tup[which(w_tup>0.01)]
w_tu = intersect(names(w_tuf),names(w_tup))


w_tdf = apply(down_fc,2,function(x){
    return(length(which(abs(x) >= log2(1.20))))
    })
w_tdf = w_tdf[which(w_tdf>=2)]
w_tdp = apply(down_pval[2:nrow(down_pval),],2,function(x){max(abs(x))})
w_tdp = w_tdp[which(w_tdp>0.01)]
w_td = intersect(names(w_tdf),names(w_tdp))


w_a = union(w_tu,w_td)

# Ranks:

ord = Res
ord = ord[ord$target_name %in% w_a,]
ord = ord[ord$n_hits >=2,]
ord = ord[order(ord$log2fc),]
ord$pos = 1:dim(ord)[1]
ord = split(ord,ord$target_name)

ord = unlist(lapply(ord,function(x){
    return(sum(x$pos))
}))

ord = ord[order(ord)]

down_fc = down_fc[,names(ord)]
down_pval = down_pval[,names(ord)]
up_fc = up_fc[,names(ord)]
up_pval = up_pval[,names(ord)]

m = rbind(down_fc,up_fc)
ord = colSums(m)
names(ord) = colnames(down_fc)
ord = ord[order(ord)]

down_fc = down_fc[,names(ord)]
down_pval = down_pval[,names(ord)]
up_fc = up_fc[,names(ord)]
up_pval = up_pval[,names(ord)]

# PLOT TIME
# DOWN
# Scale to plateaux at 10^-10
down_pval[down_pval>=10] = 10
down_pval = (down_pval-min(down_pval))/max((down_pval-min(down_pval)))
down_pval = sqrt(down_pval/pi)
down_pval = (down_pval-min(down_pval))/max((down_pval-min(down_pval)))

col_fun_down = circlize::colorRamp2(c( -1 , -0.1  ), c( "blue" , "white" ) ) 
dds =  20/max(nrow(down_fc),ncol(down_fc))
HRdown = Heatmap(          down_fc,
                        width = ncol(down_fc)*unit(dds, "cm"), 
                        height = nrow(down_fc)*unit(dds, "cm"),
                        rect_gp = gpar(type = "none"), 
                        cell_fun = function(j, i, x, y, width, height, fill) {
                            grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = NA, fill = NA))
                            grid.circle(x = x, y = y, r = abs(down_pval[i, j])/2.0 * unit(dds, "cm"), 
                            gp = gpar(fill = col_fun_down(down_fc[i, j]), col = NA))
                        }, 
                        show_row_names = TRUE,
                        show_column_names = TRUE,
                        column_title_gp = gpar(fontsize = 5),
                        row_title_gp = gpar(fontsize = 5),
                        row_names_gp = gpar(fontsize = 5),
                        heatmap_legend_param = list(title = "Calls",direction = "horizontal"),
                        column_names_gp = gpar(fontsize = 5),
                        row_dend_width = unit(2, "mm"),
                        border=TRUE,
                        col = col_fun_down ,
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE,
                        row_dend_reorder = FALSE,
                        column_dend_reorder=FALSE
                        )

# UP
# Scale to plateaux at 10^-10
up_pval[up_pval>=10] = 10
up_pval = (up_pval-min(up_pval))/max((up_pval-min(up_pval)))
up_pval = sqrt(up_pval/pi)
up_pval = (up_pval-min(up_pval))/max((up_pval-min(up_pval)))

col_fun_up = circlize::colorRamp2(c( 0.1, 1 ), c( "white" , "red" ) ) 
dds =  20/max(nrow(up_fc),ncol(up_fc))
HRup = Heatmap(          up_fc,
                        width = ncol(up_fc)*unit(dds, "cm"), 
                        height = nrow(up_fc)*unit(dds, "cm"),
                        rect_gp = gpar(type = "none"), 
                        cell_fun = function(j, i, x, y, width, height, fill) {
                            grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = NA, fill = NA))
                            grid.circle(x = x, y = y, r = abs(up_pval[i, j])/2.0 * unit(dds, "cm"), 
                            gp = gpar(fill = col_fun_up(up_fc[i, j]), col = NA))
                        }, 
                        show_row_names = TRUE,
                        show_column_names = TRUE,
                        column_title_gp = gpar(fontsize = 5),
                        row_title_gp = gpar(fontsize = 5),
                        row_names_gp = gpar(fontsize = 5),
                        heatmap_legend_param = list(title = "Calls",direction = "horizontal"),
                        column_names_gp = gpar(fontsize = 5),
                        row_dend_width = unit(2, "mm"),
                        border=TRUE,
                        col = col_fun_up ,
                        cluster_rows = FALSE, 
                        cluster_columns = FALSE,
                        row_dend_reorder = FALSE,
                        column_dend_reorder=FALSE
                        )

ht_list = HRup %v% HRdown
pdf(paste0(save_folder_sub_stat,"/Multivariate_",target,"_vs_",control,"_final_results_graph.pdf"),useDingbats = FALSE , width=10   , height=10)
        draw(ht_list)
dev.off()  



## from the above code we are assessing the probability to observe a determinied number of gRNA "scoring" at a specific log2fc cutoff per target.
## To get an FDR we need to compare the above probability to the probability of observing the same number of gRNA "scoring" at a specific log2fc cutoff per target in the 
## random data - this is already done as the random is used to generate the bi-variate frequency distribution.
# To get a global score per target and therefore rank them providing:

