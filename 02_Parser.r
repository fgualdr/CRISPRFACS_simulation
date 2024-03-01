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
source("/hpcnfs/data/GN2/fgualdrini/Master_batch_scripts/CRISPR_TOOLS_BANCH/Simulations_KERNEL/accessory_utiles_kernel.r")
setwd("/hpcnfs/data/GN2/fgualdrini/Master_batch_scripts/CRISPR_TOOLS_BANCH/Simulations_KERNEL/")

#### params:

use_nt_scramble = TRUE

# 1) CPUS:
cpus=4
if(cpus <= 1){cpus = 2} ## at least 2 cpus
if(!is.null(cpus)){
    param <- BiocParallel::SnowParam(workers=cpus,tasks=0,stop.on.error=TRUE,progressbar = TRUE,type="SOCK")
    param2 <- BiocParallel::SnowParam(workers=2,tasks=0,stop.on.error=TRUE,progressbar = TRUE,type="SOCK")
    param_serial <- BiocParallel::SerialParam(stop.on.error = TRUE,progressbar = TRUE)
}else{
    param <- NULL
    param_serial <- BiocParallel::SerialParam(stop.on.error = TRUE,progressbar = TRUE)
    param2 <- BiocParallel::SnowParam(workers=2,tasks=0,stop.on.error=TRUE,progressbar = TRUE,type="SOCK")
}

# 2) Input file:
all_simulations <- list.files("/hpcnfs/data/GN2/fgualdrini/Master_batch_scripts/CRISPR_TOOLS_BANCH/Simulations_KERNEL/Simulated_data_counts/",pattern=".txt$",full.names=TRUE,recursive=FALSE)
names(all_simulations) <- gsub(".txt","",basename(all_simulations))
simulation_path <- unique(dirname(all_simulations))
meta_test <- read.delim("/hpcnfs/data/GN2/fgualdrini/Master_batch_scripts/CRISPR_TOOLS_BANCH/Simulations_KERNEL/meta_data.txt",sep="\t",header=TRUE,row.names=NULL,stringsAsFactors=FALSE)

# loop through the meta data - for each row we will run the analysis:
# we sequence from the last row to the first:
for(i in nrow(meta_test):1){
    cat(i,"\n")
    # get the row:
    simulation_path_test <- paste0(simulation_path,"/",meta_test[i,"file_id"],"/")
    # check if the path does not exist:
    if(!dir.exists(simulation_path_test)){
        cat("Processing:",meta_test[i,"file_id"],"\n")
        file = all_simulations[meta_test[i,"file_id"]]
        cat("File:",file,"\n")
        cat("Simulation path:",simulation_path_test,"\n")

        dir.create(simulation_path_test)

        # get the input file:
        x = read.delim(file,sep="\t",header=TRUE,row.names=NULL,stringsAsFactors=FALSE)
        rownames(x) = x$sgRNA
        # We keep the counts only:
        samples.vec <- colnames(x)[grep("_r[0-9]$|_rep[0-9]$",tolower(colnames(x)))]
        # samples should be within c("plasmid","presort","low","high"):
        samples.vec <- samples.vec[grepl("^plasmid|^presort|^low|^high",tolower(samples.vec))]
        x <- x[,samples.vec]

        # the analysis will be executed pair-wise:
        list_comparisons <- list(
                high_vs_low = c("high","low")
                #high_vs_presort = c("high","presort"),
                #low_vs_presort = c("low","presort")
        )

        for(comparison in list_comparisons){
            cat("Comparison:",comparison[1],"vs",comparison[2],"\n")
            save_folder <- paste0(simulation_path_test,paste0("/CrisprRes",comparison[1],"_vs_",comparison[2],"/"))
            unlink(save_folder, recursive = TRUE)
            dir.create(save_folder)
            # 3) taget / control:
            target=comparison[1]
            control=comparison[2]
            check_levels <- c(target,control)
            check_levels <- strsplit(check_levels, "_")
            check_levels <- do.call(rbind, check_levels)
            # this will be a matrix with two columns:
            check_levels <- check_levels[, apply(check_levels, 2, function(x) length(unique(x)) > 1)]
            # re-assign to target and control:
            target <- check_levels[1]
            control <- check_levels[2]
            names(check_levels) <- c("target","control")
            cat("Target:",target,"\n")
            cat("Control:",control,"\n")

            # 4) Make the design:
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

            # 5) Normalize the data:
            ori_mat = round(x,0)
            Result = RunNorm( ori_mat,
                            deg_design,
                            fix_reference=deg_design$Sample_ID,
                            row_name_index=1,
                            saving_path=save_folder,
                            n_pop=1,
                            n_pop_reference=1,
                            BiocParam=param,
                            Save_results=FALSE)
            colData = Result$scaling_factors
            mat = Result$norm_mat
            deg_design = deg_design[order(deg_design$Sample_Replicate),]

            # 6) Compute DESEQ2
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

            # 7) Compute the bivariate distribution:
            # We compute this for a range of cut offs:
            l2c_cuts <- seq(min(abs(res$log2FoldChange)),max(abs(res$log2FoldChange)),by=log2(1.01)) 
            # 0.01 means 2^0.01 = 1.007 i.e. 0.7% change less then 1% change - so we check increasing the l2fc by 0.01 already quite granular
            # checnging this will increase matrix size and computation time - and could exceed memory - also in the Kde2d function we limit the bandwidth to and minimum of 0.01
            # else huge memory usage and long computation time

            # We separate up and down regulated effects:
            res_up = res[which(res$log2FoldChange > 0),]
            res_down = res[which(res$log2FoldChange < 0),]
            res_down$log2FoldChange <- -1*res_down$log2FoldChange

            if(use_nt_scramble){
                # we separate the "scaramble|random" targets from the rest:
                res_up_r = res_up[grepl("^rand|^scram|^nontargeting",rownames(res_up)),]
                res_up = res_up[!grepl("^rand|^scram|^nontargeting",rownames(res_up)),]
                res_down_r = res_down[grepl("^rand|^scram|^nontargeting",rownames(res_down)),]
                res_down = res_down[!grepl("^rand|^scram|^nontargeting",rownames(res_down)),]
            }else{
                # if no scrrambled to be used we randomly all the gRNA by grouping in random sets
                res_up_r = res_up
                res_up = res_up[!grepl("^rand|^scram|^nontargeting",rownames(res_up)),]
                res_down_r = res_down
                res_down = res_down[!grepl("^rand|^scram|^nontargeting",rownames(res_down)),]
            }
            # for the target we get the number of hits per target for each l2fc:
            res_up <- BiocParallel::bplapply(l2c_cuts,function(x) get_hits_per_target(res_up,x),BPPARAM = param)
            res_up <- as.data.frame(do.call(rbind,res_up),stringsAsFactors=FALSE)
            res_down <- BiocParallel::bplapply(l2c_cuts,function(x) get_hits_per_target(res_down,x),BPPARAM = param)
            res_down <- as.data.frame(do.call(rbind,res_down),stringsAsFactors=FALSE)
            # as we considered the abs of the res_down their associated l2fc are negative:
            res_down$log2fc <- -1*res_down$log2fc
            # clear the memory with gc() keeping anythin to be used downstream:
            gc() 

            # Ok given the res_up and res_down we compute the bivariate distribution of the number of hits per target for each pair of log2fc:
            # We aggregate by n_hits and log2fc and count the number of elements:
            xy <- as.data.frame(rbind(res_up,res_down),stringsAsFactors=FALSE)
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
            gc()

            # for each of the 1000 list we compute the MASS::kde2d
            # for consistency we need to set "n" and "lim" coherently to the xy of the real data:
            n_grid_x <- length(unique(xy$n_hits))
            # for the y grid we consider the total number of all oligos in the real data:
            n_grid_y <- dim(res)[1] # like saying each oligo has a different log2fc
            # we combine the list_r with n_grid_x and n_grid_y so to have a list of list each element with the "xy" of the list_r and the "n_grid_x" and "n_grid_y":
            list_r_m <- lapply(1:length(list_r),function(i){
                l = list()
                l[["id"]] = i
                l[["xy"]] = data.table(list_r[[i]])
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
            gc()

            list_r_kde <- BiocParallel::bplapply(list_r_m, function(x) BiocParallel::bptry(mass_kde2_parallel(x)), BPPARAM = param_serial)
            list_r_kde_f <- list.files(pattern="rds$",full.names=TRUE,recursive=TRUE)
            # We should account for 
            cat("Get the maximum Machine number \n")
            list_r_kde_norm <- lapply(list_r_kde_f,function(x){
                    # read the kde2d:
                    x = readRDS(x)
                    # is better to add only to the "0"
                    x$z[x$z == 0] <- .Machine$double.xmin
                    # x$z <- x$z + (.Machine$double.xmin)
                    xlim <- range(x$x)
                    ylim <- range(x$y)
                    cell_size <- (diff(xlim) / length(x$x)) * (diff(ylim) / length(x$y))
                    # Compute the normalization function:
                    norm <- sum(x$z) * cell_size  # normalizing constant    
                    x$z = (x$z*cell_size)/norm 
                    x$z_log = log(x$z)
                    return(x)
            })

            # remove the tmp folder:
            setwd(save_folder_sub_stat)
            unlink(paste0(save_folder_sub_stat,"/tmp/"), recursive = TRUE, force = TRUE)

            cat("Calculate the mean and standard deviation of the kernel density estimates\n")
            # 2) Calculate the mean and standard deviation of the kernel density estimates.
            # the "z" is a matrix we want to preserve the dimensionality of the matrix so we use sapply. sapply will get the 
            mean_density <- Reduce("+", lapply(list_r_kde_norm, "[[", "z")) / length(list_r_kde_norm)
            # to do the same on the z_log the average of logged values should be lp1 = log(p1); lp2 = log(p2); ...; lpn = log(pn) and then the average of the logs is log((p1+p2+...+pn)/n)

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
            y_lim <- summary(kde2d_mean$y)
            # consider 1st and 3rd quartile:
            y_lim <- quantile(kde2d_mean$y,probs=c(0.3,0.7))
            pdf(paste0(save_folder_sub_stat,"/Multivariate_",target,"_vs_",control,"_KDE_RANDOM.pdf"))
                filled.contour( 
                                kde2d_mean,
                                color.palette=colorRampPalette(c('white','blue','yellow','red','darkred')),
                                ylim = range(y_lim, finite = TRUE),
                                zlim = range(kde2d_mean$z, finite = TRUE)
                                )
            dev.off()

            # Define the function such as that given a pair of hits and log2fc it returns the probability after having interpolated the density function:
            kde2dAkima = kde2d_mean$z
            colnames(kde2dAkima) = kde2d_mean$y
            rownames(kde2dAkima) = kde2d_mean$x
            kde2dAkima = reshape2::melt(kde2dAkima)
            colnames(kde2dAkima) = c("n_hits","log2fc","density")
            kde2dAkima = kde2dAkima[which(kde2dAkima$n_hits > 0 & kde2dAkima$n_hits <= max_oligo_per_target),]
            # scale the density to sum to 1:
            kde2dAkima$density = kde2dAkima$density/sum(kde2dAkima$density)

            # For evey point in the xy data frame we compute the probability of the point given the density function:
            # The probability is reported as 1-cdf and in particular we compute the cdf as the sum of the density from the point to the upper/lower limit

            Res = xy
            # for each element in Res$target_name we keep those rows having highest n_hits and highest abs(log2fc)
            Res = Res[order(Res$n_hits,abs(Res$log2fc),decreasing=TRUE),]
            # for each pair of target_name, n_hits we keep the highest log2fc 
            Res_up = Res[Res$log2fc >0,]
            # use goup by and select the single element s with highest log2fc
            Res_up = Res_up %>% group_by(target_name,n_hits) %>% filter(log2fc == max(log2fc))
            # convert back to data frame:
            Res_up = as.data.frame(Res_up,stringsAsFactors=FALSE)
            
            # or the lowest log2fc for the down
            Res_down = Res[Res$log2fc <0,]
            Res_down = Res_down %>% group_by(target_name,n_hits) %>% filter(log2fc == min(log2fc))
            Res_down = as.data.frame(Res_down,stringsAsFactors=FALSE)

            # up
            prob = lapply(1:nrow(Res_up),function(i){
                # display a progress bar:
                if(i %% 1000 == 0){
                    cat(paste0("Progress: ",i,"/",nrow(Res_up),"\n"))
                }
                # Here we get cdf as the sum of the density from the point to the upper limit
                x_idx <- which(kde2dAkima$n_hits == Res_up$n_hits[i])
                # here we need to sum the density from the point to the upper limit so any point above the one selected
                # like saying the probability to see such a change as extreme as the one observed
                y_idx = which(kde2dAkima$log2fc >= Res_up$log2fc[i])
                # Su up all the values 
                summed_p = sum(kde2dAkima$density[intersect(x_idx,y_idx)])
                return(summed_p)
            })
            cat("\n")
            Res_up$Bivariate_prob = unlist(prob)
            # down
            prob = lapply(1:nrow(Res_down),function(i){
                if(i %% 1000 == 0){
                    cat(paste0("Progress: ",i,"/",nrow(Res_up),"\n"))
                }
                # Get the x_values and y_values:
                x_idx <- which(kde2dAkima$n_hits == Res_down$n_hits[i])
                # as the l2fc is less than 0 we need to sum the density from the point to the lower limit i.e. any point below the one selected
                # like saying the probability to see such a change as extreme as the one observed
                y_idx = which(kde2dAkima$log2fc <= Res_down$log2fc[i])
                # Subset the density values:
                summed_p = sum(kde2dAkima$density[intersect(x_idx,y_idx)])
                return(summed_p)
            })
            cat("\n")
            Res_down$Bivariate_prob = unlist(prob)
            Res = rbind(Res_up,Res_down)
            Res$adjp_fdr_bivariate = p.adjust(Res$Bivariate_prob, method = "fdr")
            write.table(Res,paste0(save_folder_sub_stat,"/Multivariate_",target,"_vs_",control,"_prob.txt"),sep="\t",col.names=NA)

            # the above will give multiple results per target - we want to integrate across all hits:
            # we get the list of targets only make a one column data frame and then we aggregate by target_name and n_hits
            Res_t = Res[,c("target_name"),drop=FALSE]
            # get unique rows only:
            Res_t = unique(Res_t)
            Res_t$Cumulative_p_up <- 1
            Res_t$Cumulative_p_down <- 1
            # for each direction:
            up_p <- lapply(1:nrow(Res_t),function(x){
                if(x %% 100 == 0){
                    cat(paste0("Progress: ",x,"/",nrow(Res_t),"\n"))
                }
                if(length(which(Res_up$target_name %in% Res_t$target_name[x])) >0){
                    # we integrate only for the log2fc at the various n_hits:
                    y = Res_up[which(Res_up$target_name %in% Res_t$target_name[x]),]
                    rownames(y) = paste0(y$target_name,"_",y$n_hits)
                    # we make a data frame empty with target_name, n_hits, log2fc and Bivariate_prob
                    # target_name the same as in y, n_hits 1:max(kde2dAkima$n_hits)
                    y_d <- data.frame( target_name=rep(Res_t$target_name[x],length(1:max(kde2dAkima$n_hits))),
                                        n_hits=1:max(kde2dAkima$n_hits),
                                        log2fc=rep(NA,length(1:max(kde2dAkima$n_hits))),
                                        Bivariate_prob=rep(NA,length(1:max(kde2dAkima$n_hits))))
                    rownames(y_d) = paste0(y_d$target_name,"_",y_d$n_hits)
                    y_d$Bivariate_prob = y[rownames(y_d),"Bivariate_prob"]
                    y_d$log2fc = y[rownames(y_d),"log2fc"]
                    # now we might have a situation where we need to fill in the NAs starting from the higher n_hits:
                    # if the higher n_hits has log2fc == NA we need to give all the kde2dAkima$density having n_hits == with the same n_hits and lof2fc >=0 so we sum all of them 
                    # if subsequent n_hits have log2fc == NA we need to give all the kde2dAkima$density having n_hits == with the same n_hits and lof2fc >= log2fc of the previous (i.e. the one with higer n_hits)
                    # we do this for all the n_hits:
                    c = max(kde2dAkima$n_hits)
                    for(i in length(1:max(kde2dAkima$n_hits)):1){
                        if(is.na(y_d$log2fc[i]) & i == c){
                            # get the x_values and y_values:
                            x_idx <- which(kde2dAkima$n_hits == y_d$n_hits[i])
                            y_idx = which(kde2dAkima$log2fc >= 0)
                            # Su up all the values
                            summed_p = sum(kde2dAkima$density[intersect(x_idx,y_idx)])
                            y_d$Bivariate_prob[i] = summed_p
                        }else{
                            if(is.na(y_d$log2fc[i])){
                                # get the x index:
                                x_idx <- which(kde2dAkima$n_hits == y_d$n_hits[i])
                                # we get the lfc of the previous n_hits : ie.e the higher one:
                                lfc = y_d$log2fc[i+1]
                                y_idx = which(kde2dAkima$log2fc >= lfc)
                                # Su up all the values
                                summed_p = sum(kde2dAkima$density[intersect(x_idx,y_idx)])
                                y_d$Bivariate_prob[i] = summed_p
                                y_d$log2fc[i] = lfc
                            }
                        }
                    }
                    return(sum(y_d$Bivariate_prob))
                }else{
                    return(1)
                }
            })
            Res_t$Cumulative_p_up = unlist(up_p)
            # for the down 
            down_p <- lapply(1:nrow(Res_t),function(x){
                if(x %% 100 == 0){
                    cat(paste0("Progress: ",x,"/",nrow(Res_t),"\n"))
                }
                if(length(which(Res_down$target_name %in% Res_t$target_name[x])) >0){
                    # we integrate only for the log2fc at the various n_hits:
                    y = Res_down[ which(Res_down$target_name %in% Res_t$target_name[x]),]
                    rownames(y) = paste0(y$target_name,"_",y$n_hits)
                    # we make a data frame empty with target_name, n_hits, log2fc and Bivariate_prob
                    # target_name the same as in y, n_hits 1:max(kde2dAkima$n_hits)
                    y_d <- data.frame( target_name=rep(Res_t$target_name[x],length(1:max(kde2dAkima$n_hits))),
                                        n_hits=1:max(kde2dAkima$n_hits),
                                        log2fc=rep(NA,length(1:max(kde2dAkima$n_hits))),
                                        Bivariate_prob=rep(NA,length(1:max(kde2dAkima$n_hits))))
                    rownames(y_d) = paste0(y_d$target_name,"_",y_d$n_hits)
                    y_d$Bivariate_prob = y[rownames(y_d),"Bivariate_prob"]
                    y_d$log2fc = y[rownames(y_d),"log2fc"]
                    # now we might have a situation where we need to fill in the NAs starting from the higher n_hits:
                    # if the higher n_hits has log2fc == NA we need to give all the kde2dAkima$density having n_hits == with the same n_hits and lof2fc >=0 so we sum all of them 
                    # if subsequent n_hits have log2fc == NA we need to give all the kde2dAkima$density having n_hits == with the same n_hits and lof2fc >= log2fc of the previous (i.e. the one with higer n_hits)
                    # we do this for all the n_hits:
                    c = max(kde2dAkima$n_hits)
                    for(i in length(1:max(kde2dAkima$n_hits)):1){
                        if(is.na(y_d$log2fc[i]) & i == c){
                            # get the x_values and y_values:
                            x_idx <- which(kde2dAkima$n_hits == y_d$n_hits[i])
                            y_idx = which(kde2dAkima$log2fc <= 0)
                            # Su up all the values
                            summed_p = sum(kde2dAkima$density[intersect(x_idx,y_idx)])
                            y_d$Bivariate_prob[i] = summed_p
                        }else{
                            if(is.na(y_d$log2fc[i])){
                                # get the x index:
                                x_idx <- which(kde2dAkima$n_hits == y_d$n_hits[i])
                                # we get the lfc of the previous n_hits : ie.e the higher one:
                                lfc = y_d$log2fc[i+1]
                                y_idx = which(kde2dAkima$log2fc <= lfc)
                                # Su up all the values
                                summed_p = sum(kde2dAkima$density[intersect(x_idx,y_idx)])
                                y_d$Bivariate_prob[i] = summed_p
                                y_d$log2fc[i] = lfc
                            }
                        }
                    }
                    return(sum(y_d$Bivariate_prob))
                }else{
                    return(1)
                }
            })
            Res_t$Cumulative_p_down = unlist(down_p)

            Res_t$adjp_fdr_bivariate_up = p.adjust(Res_t$Cumulative_p_up, method = "fdr")
            Res_t$adjp_fdr_bivariate_down = p.adjust(Res_t$Cumulative_p_down, method = "fdr")
            write.table(Res_t,paste0(save_folder_sub_stat,"/Multivariate_",target,"_vs_",control,"_prob_per_TARGET.txt"),sep="\t",col.names=NA)
        }
        # save the kde2dAkima, Res , kde2d_mean and xy
        save(kde2dAkima,Res_t,Res_up,Res_down,kde2d_mean,xy,file=paste0(save_folder_sub_stat,"/Multivariate_",target,"_vs_",control,"_prob_per_TARGET.RData"))
    }
}
