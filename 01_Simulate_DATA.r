library(ggplot2)
library(reshape)
library(MAUDE)
library(truncnorm)
library(ggplot2)
library(data.table)
library(dplyr)
library(DESeq2)
set.seed(76484377)

# qsub -I -q nocg_interact -l select=1:ncpus=5
# singularity shell -B /hpcnfs docker://fgualdr/docker_maude


setwd("/hpcnfs/data/GN2/fgualdrini/Master_batch_scripts/CRISPR_TOOLS_BANCH/Simulations_KERNEL")

# add a folder to keep the tests:
simulated_folder <- paste0(getwd(),"/Simulated_data/")
unlink(simulated_folder, recursive = TRUE)
dir.create(simulated_folder)

# considering the above code we need to modify so that we can consider some other scenarios and conditions to test.
# We want to Simulate again the data so that we can generate the "expression" distribution which will follow a normal distribution
# this will be the global distribution of expression associated with all the cells in the experiment
# we need to generate cells harbouring TRUE regulators in both positive and negative directions
# so modelling as a normal distribution we can consider the global having mean 0 and sd 1 (later we can generate different distributions with different sd)
# we can simulate to sort the top 15% and bottom 15% of the distribution.
# When assembling the Simulation we will control for:
# number of genes: n_genes (these will include TRUE regulators and non-regulators)
# n_pos_regulators: n_pos_regulators
# n_neg_regulators: n_neg_regulators
# number of non-targeting gRNA: n_nt_gRNA
# number of gRNA per gene: n_guide
# so in total we will have tot_gRNA = (n_genes * n_guide) + n_nt_gRNA
# Then we will consider the distribution of the expression signal for the TRUE regulators which will be modelled as skew normal distribution having:
# mu_pos = 0.5 and sd_pos = 0.1
# mu_neg = -0.5 and sd_neg = 0.1 .. these values can be tested to evaluate the impact of asymmetric disposition and number of pos and neg
# so we will assume that the sum of the skewed distribution with the global will generate still a normal distribution..
# basically the skewed pos and neg are inscribed within the global distribution of expression.
# So we simulate that within the bins we will get both true and non-true regulators.
# By summing the 3 distribution we need to make them at the end normally distributed.
# we will then simulate read counts per cells sorted.

# Define the composition of the library, the proportion between positive and negative regulators, the number of non-targeting gRNA and the number of gRNA per gene
total_genes <- c(700,1500) # n° of target genes
n_no_effect <- c(0.9,0.7)  # percentage of genes not having an effect
pos_neg_ratio <- c(0.5,0.9,0.1)  # ratio of positive and negative regulators
n_gRNA <- c(10,4)  # number of gRNA per target gene 
n_non_targeting_gRNA <- 1000  # number of non-targeting gRNA
# Define the distribution of the effect of the gRNA: 
# The non essential effect is modelled as a normal distribution with mean 0 and variance 0.5
non_essential_mean <- 0       # mean of non-essential genes
non_essential_variance <- 0  # variance of non-essential genes
# Positive and Negative effects are defined as truncated normal distributions
# a) positive effect: truncated normal distribution with mean 0.1, variance 2, truncated between 0.1 and 4
pos_effect_lower_bound <- 0.1  # lower bound of the positive effect # this can be reduced to check sensitivity of the Tool in detecting small effects
pos_effect_upper_bound <- 3   # upper bound of the positive effect
pos_effect_variance <- c(0.5,1,2.5,5)
mean_pos_effect <- rev(c(0.1,0.25,0.5,1))
# b) negative effect: truncated normal distribution with mean -0.5, variance 2, truncated between -4 and -0.1
neg_effect_lower_bound <- -3   # lower bound of the negative effect - cAN BE USED TO SIMULATE OUTLIERS
neg_effect_upper_bound <- -0.1  # upper bound of the negative effect
neg_effect_variance <- c(0.5,1,2.5,5)
mean_neg_effect <- rev(c(0.1,0.25,0.5,1) )
# We deine the variability of the gRNA efficiency
guide_efficiency_mean = 0.8
guide_efficiency_variance = 0.05  # variance of the gRNA efficiency (so if a group of gRNA targeting the same gene will have the gene specific effect we should add a variability to the gRNA efficiency)
# Define the total number of cells per bins
# n° of replicate experiments
n_reps = 3  # number of replica experiments
# n° of read counts per gRNA (this is the average count per gRNA)
depth_factors = c(100,200) # read per gRNA
interceptSD <- 0.8 # this is the spead of expression levels in log2 scale
# add parameters to include variability in registering the read-out by the FACS machine
noise_mean = 0
noise_variance = c(0.05)
# sorting parameters:
tailP=0.00000001

fraction_sorted = c(0.1,0.15,0.20,0.25)  # fraction of the distribution to be sorted
# this parameter is done by each replicate! so we can later compare..

# nbinom dispersion:
dispersion <- 0.05
# we can loop through a grid of parameters to generate different scenarios

Combo_params <- expand.grid(total_genes = total_genes, 
                            n_no_effect = n_no_effect, 
                            pos_neg_ratio = pos_neg_ratio, 
                            n_gRNA = n_gRNA, 
                            pos_effect_variance = pos_effect_variance, 
                            mean_pos_effect = mean_pos_effect, 
                            neg_effect_variance = neg_effect_variance, 
                            mean_neg_effect = mean_neg_effect, 
                            depth_factors = depth_factors, 
                            noise_variance = noise_variance)

# add a test_ID column which will report the column name and the value of the parameter all collapse by "_" so that we can use it as file names when saving the data
Combo_params$test_ID = apply(Combo_params, 1, function(x) paste0(names(x),"_",x, collapse = "_"))

for(test in 1:nrow(Combo_params)){

  print(paste0("Test: ", test))
  print(Combo_params$test_ID[test])

  # generate the saving subfolder
  subfolder = paste0(simulated_folder,Combo_params$test_ID[test],"/")
  dir.create(subfolder)

  # we can then extract the parameters for the test
  total_genes = Combo_params$total_genes[test]
  n_no_effect = Combo_params$n_no_effect[test]
  pos_neg_ratio = Combo_params$pos_neg_ratio[test]
  n_gRNA = Combo_params$n_gRNA[test]
  pos_effect_variance = Combo_params$pos_effect_variance[test]
  mean_pos_effect = Combo_params$mean_pos_effect[test]
  neg_effect_variance = Combo_params$neg_effect_variance[test] 
  mean_neg_effect = Combo_params$mean_neg_effect[test] * -1
  depth_factors = Combo_params$depth_factors[test]
  noise_variance = Combo_params$noise_variance[test]


  n_true_targets <- total_genes*(1-n_no_effect)  # number of genes having a TRUE effect

  # test the below with a single set of parameters set:


  # Define the effect distribution for the non-essential genes, the positive and negative regulators and the non-targeting gRNA
  # a) non-essential genes: normal distribution with mean 0 and variance 0.5
  non_essential_effect <- rnorm(  n = round(total_genes*n_no_effect)  , 
                                  mean = non_essential_mean , 
                                  sd = non_essential_variance
                                )
  # b) non-targeting gRNA: normal distribution with mean 0 and variance 0.5
  scrambled_effect <- rnorm(  n = n_non_targeting_gRNA  , 
                              mean = non_essential_mean , 
                              sd = non_essential_variance
                            )
  # c) positive regulators: truncated normal distribution with mean 0.1, variance 2, truncated between 0.1 and 4
  positive_effect <- rtruncnorm(  n = round(n_true_targets*pos_neg_ratio) , 
                                  a = pos_effect_lower_bound, 
                                  b = pos_effect_upper_bound, 
                                  mean = mean_pos_effect, 
                                  sd = pos_effect_variance
                                )
  # d) negative regulators: truncated normal distribution with mean -0.5, variance 2, truncated between -4 and -0.1
  negative_effect <- rtruncnorm(  n = round(n_true_targets*(1-pos_neg_ratio)) , 
                                  a = neg_effect_lower_bound, 
                                  b = neg_effect_upper_bound, 
                                  mean = mean_neg_effect, 
                                  sd = neg_effect_variance
                                )

  # The above are the effects - so random considering "TRUE" positve and "TRUE" negative regulators and non-essential genes + non-targeting gRNA
  # We assemble the ground truth for the simulation - which will be the ground truth related to the effect of each TARGET gene
  # Assembled in data.frames
  groundTruth_non_essential = data.frame(gid = paste0("non_essential_",1:(round(total_genes*n_no_effect))), 
                                        meanEffect = non_essential_effect,
                                        type = "Nonessential"
                                        )
  groundTruth_non_targeting = data.frame( gid = paste0("non_targeting_",1:n_non_targeting_gRNA), 
                                          meanEffect = scrambled_effect,
                                          type = "Nontargeting"
                                        )                               
  groundTruth_pos = data.frame(   gid = paste0("truepositive_",1:round(n_true_targets*pos_neg_ratio)), 
                                  meanEffect = positive_effect,
                                  type = "positive"
                                )
  groundTruth_neg = data.frame(   gid = paste0("truenegative_",1:round(n_true_targets*(1-pos_neg_ratio))), 
                                  meanEffect = negative_effect,
                                  type = "negative"
                                )

  # Combine in a single data.frame

  groundTruth = rbind(groundTruth_non_essential, groundTruth_pos, groundTruth_neg, groundTruth_non_targeting)
  groundTruth$n_gRNA = n_gRNA
  # the non-targeting have only 1 gRNA
  groundTruth$n_gRNA[groundTruth$type == "non_targeting"] = 1

  p <- ggplot(groundTruth, aes(x=meanEffect)) + 
    geom_density(adjust = 1/5) +
    facet_wrap(~type) + 
    theme_classic() + 
    xlab("gRNA effect") + 
    scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand=c(0,0)) + 
    coord_cartesian(xlim=c(min(groundTruth$meanEffect), max(groundTruth$meanEffect)));
  ggsave(paste0(subfolder,"01_ground_truth_per_type.png"), p, width=8, height=6, units="in", dpi=300);

  p <- ggplot(groundTruth, aes(x=meanEffect)) + 
    geom_density(adjust = 1/5) +
    theme_classic() + 
    xlab("gRNA effect") + 
    scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand=c(0,0)) + 
    coord_cartesian(xlim=c(min(groundTruth$meanEffect), max(groundTruth$meanEffect)));
  ggsave(paste0(subfolder,"01_ground_truth_total.png"), p, width=8, height=6, units="in", dpi=300);

  # expand groundTruth each row by the number of gRNA
  rownames(groundTruth) = groundTruth$gid
  groundTruth = groundTruth[rep(row.names(groundTruth), groundTruth$n_gRNA),]
  rownames(groundTruth) = NULL
  # Add the elements column as the gid + _ + number of gRNA
  groundTruth$elements = paste0(groundTruth$gid, "_gRNA_", 1:groundTruth$n_gRNA)
  # the mean TRUE effect associated with perturbing the target gene will be scaled by the gRNA efficiency
  # guide effieincy is defined randomly with a mean of 0.8 and a variance of 0.1 and capped at 1
  # so we can model it as a truncated normal distribution with mean 0.8, variance 0.1, lower bound 0 and upper bound 1
  guide_efficiency = rtruncnorm(  n = nrow(groundTruth),
                                  a = 0, 
                                  b = 1, 
                                  mean = guide_efficiency_mean, 
                                  sd = guide_efficiency_variance
                                )

  # we can then scale the mean effect by the guide efficiency
  groundTruth$meanEffect_pergRNA = groundTruth$meanEffect * guide_efficiency # as % of the effect

  ### OK now we need to generate replicates including each time the noise associated with the read-out of the FACS - so given the gorund truth we can simulate the read counts
  ll_reps = list()

  n = nrow(groundTruth)
  average_counts <- depth_factors # average count per gene
  interceptMean <- log2(average_counts)
  BaseMean = rnorm(n,interceptMean,interceptSD)
  m = 2 # two conditions only
  dispMeanRel=function(x) 4/x + .01 # here instead of doing this we can vary the dispersion parameter

  for(i in 1:n_reps){

    cat("Replica: ", i, "\n")

    groundTruth_rep = groundTruth

    # this is the FACS read-out noise
    groundTruth_rep$meanEffect_pergRNA_noise = groundTruth_rep$meanEffect_pergRNA + rnorm(nrow(groundTruth_rep), noise_mean, noise_variance)
    # here we apply centering again
    groundTruth_rep$meanEffect_pergRNA_noise = groundTruth_rep$meanEffect_pergRNA_noise - mean(groundTruth_rep$meanEffect_pergRNA_noise)

    p <- ggplot(groundTruth_rep, aes(x=meanEffect_pergRNA_noise)) + 
      geom_density(adjust = 1/5) +
      facet_wrap(~type) + 
      theme_classic() + 
      xlab("gRNA effect") + 
      scale_x_continuous(expand=c(0,0)) + 
      scale_y_continuous(expand=c(0,0)) + 
      coord_cartesian(xlim=c(min(groundTruth_rep$meanEffect_pergRNA_noise), max(groundTruth_rep$meanEffect_pergRNA_noise)));
    ggsave(paste0(subfolder,"02_ground_truth_per_type_rep",i,".png"), p, width=8, height=6, units="in", dpi=300);

    p <- ggplot(groundTruth_rep, aes(x=meanEffect_pergRNA_noise)) + 
      geom_density(adjust = 1/5) +
      theme_classic() + 
      xlab("gRNA effect") + 
      scale_x_continuous(expand=c(0,0)) + 
      scale_y_continuous(expand=c(0,0)) + 
      coord_cartesian(xlim=c(-1,1));
    ggsave(paste0(subfolder,"02_ground_truth_total_rep",i,".png"), p, width=8, height=6, units="in", dpi=300);
    groundTruth_repORI = groundTruth_rep
    for(sorting_frac in fraction_sorted){
      groundTruth_rep = groundTruth_repORI
      cat("Sorting fraction: ", sorting_frac, "\n")
      subfolder_frac = paste0(subfolder,"sorted_",sorting_frac,"/")
      # we can then simulate the number of reads per gRNA per sorted bin:
      binBounds = data.frame(Bin=c("LOW","HIGH"), fraction=rep(sorting_frac,2))
      # "LOW" BIN will be the bottom 15% of the distribution
      # "HIGH" BIN will be the top 15% of the distribution
      # so in theory from -Inf to 0.15 and from 0.85 to Inf
      binBounds$binStartQ = c(tailP, 1-sorting_frac); # start at : LOW = 0.15, HIGH = 0.85
      binBounds$binEndQ = c(sorting_frac, 1-tailP); # end at : LOW = 0.30, HIGH = 1
      # binBounds$binStartZ and binBounds$binEndZ are going to be based on the quantiles of the actual distribution fit_dist
      binBounds$binStartZ = quantile(groundTruth_rep$meanEffect_pergRNA_noise, probs = binBounds$binStartQ);
      binBounds$binEndZ = quantile(groundTruth_rep$meanEffect_pergRNA_noise, probs = binBounds$binEndQ);

      pos_sel = groundTruth_rep[grepl("positive_", groundTruth_rep$elements),]
      neg_sel = groundTruth_rep[grepl("negative_", groundTruth_rep$elements),]
      non_essential_non_targeting = groundTruth_rep[grepl("non_essential_", groundTruth_rep$elements) | grepl("non_targeting_", groundTruth_rep$elements),]

      #plot the select examples and show the bin structure
      p = ggplot(groundTruth_rep, aes(x=meanEffect_pergRNA_noise,y = after_stat(count)))+
        geom_density(alpha=0.2,adjust = 1/5) + 
        geom_vline(xintercept = sort(unique(c(binBounds$binStartZ,binBounds$binEndZ))),colour="gray") + 
        theme_classic() + 
        scale_fill_manual(values=c("red","darkgray")) + 
        xlab("Target expression") + 
        scale_x_continuous(expand=c(0,0)) + 
        scale_y_continuous(expand=c(0,0)) + 
        coord_cartesian(xlim=c(-1,1)) +
        geom_segment(data=binBounds, aes(x=binStartZ, xend=binEndZ, colour=Bin, y=0, yend=0), size=5, inherit.aes = FALSE); 

      p <- p + geom_density(data=pos_sel, aes(x=meanEffect_pergRNA_noise,y = after_stat(count)), alpha=0.2, fill="red",adjust = 1/5) + 
        geom_density(data=neg_sel, aes(x=meanEffect_pergRNA_noise,y = after_stat(count)), alpha=0.2, fill="blue",adjust = 1/5) + 
        geom_density(data=non_essential_non_targeting, aes(x=meanEffect_pergRNA_noise,y = after_stat(count)), alpha=0.2, fill="darkgray",adjust = 1/5) +
        theme_classic() + 
        scale_fill_manual(values=c("red","darkgray")) + 
        xlab("Target expression") + 
        scale_x_continuous(expand=c(0,0)) + 
        scale_y_continuous(expand=c(0,0)) + 
        coord_cartesian(xlim=c(-1,1)) +
        geom_segment(data=binBounds, aes(x=binStartZ, xend=binEndZ, colour=Bin, y=0, yend=0), size=5, inherit.aes = FALSE);
      
      ggsave(paste0(subfolder_frac,"03_ground_truth_per_type_rep",i,"_bins.png"), p, width=8, height=6, units="in", dpi=300);

      # Assign gRNA to BINs
      groundTruth_rep$BIN = "NS"
      groundTruth_rep$BIN[groundTruth_rep$meanEffect_pergRNA_noise < binBounds$binEndZ[1]] = "LOW"
      groundTruth_rep$BIN[groundTruth_rep$meanEffect_pergRNA_noise > binBounds$binStartZ[2]] = "HIGH"
      rownames(groundTruth_rep) = groundTruth_rep$elements

      cat("LOW: ", sum(groundTruth_rep$BIN == "LOW"), "HIGH: ", sum(groundTruth_rep$BIN == "HIGH"), "\n")

      # get high counts:
      beta <- as.data.frame(cbind(  BaseMean, # baseMean is the mean expression in log2 scale
                      betaSD = groundTruth_rep$meanEffect_pergRNA_noise - binBounds[binBounds$Bin == "HIGH",]$binStartZ # this is the mean effect by shifting the mean effect to the lower bound of the HIGH bin we simulate the sorting
                      ))
      # sumbratracting to the lower bound will change the betaSD so the recorded effect but practically this is the same ... like we are not changin the 
      beta$dispersion <- dispMeanRel(2^(beta[,1])) # dispersion is a function of the mean
      
      # HIGH
      colData <- DataFrame(condition=factor(rep(c("PRESORT","HIGH"),times=c(ceiling(1),floor(1))), levels=c("PRESORT", "HIGH"))) 
      x <- stats::model.matrix.default(~ colData$condition)
      mu <- t(2^(x %*% t(beta[,1:2])) ) # the %*% is the matrix multiplication between the model matrix and the beta
      countData <- matrix(rnbinom(m*n, mu=mu, size=1/beta$dispersion), ncol=m)
      mode(countData) <- "integer"
      # PRESORT and HIGH
      pre_sort <- countData[,1:(m/2)]
      high_sort <- countData[,(m/2+1):m]
      colnames(mu) <- c("mu_presort","mu_h")
      beta_h <- as.data.frame(cbind(beta[,1:3],mu))

      # model low counts - we keep the same random generated BaseMean counts as in HIGH
      beta[,2] <- (groundTruth_rep$meanEffect_pergRNA_noise*-1) - (binBounds[binBounds$Bin == "LOW",]$binEndZ*-1)
      # LOW
      colData <- DataFrame(condition=factor(rep(c("PRESORT","LOW"),times=c(ceiling(1),floor(1))), levels=c("PRESORT", "LOW")))
      x <- if (m > 1) {
        stats::model.matrix.default(~ colData$condition)
      } else {
        cbind(rep(1,m),rep(0,m))
      }
      mu <- t(2^(x %*% t(beta[,1:2])) ) # here mu is the mean of the negative binomial distribution
      countData <- matrix(rnbinom(m*n, mu=mu, size=1/dispersion), ncol=m) # the size is the dispersion parameter
      mode(countData) <- "integer"
      # PRESORT and LOW
      pre_sort <- round(rowMeans(cbind(pre_sort, (countData[,1:(m/2)]))),0)
      low_sort <- round((countData[,(m/2+1):m]),0)
      colnames(mu) <- c("mu_presort","mu_low")
      beta_l <- as.data.frame(cbind(beta[,1:3],mu))

      ###################################
      ###################################
      ###################################
      # merge the beta_h and beta_l so that nonredundant columns are kept
      betaS <- cbind(beta_h[,1],beta_h[,2:3],beta_l[,2:3],beta_h[,4:5],beta_l[,4:5])
      colnames(betaS) <- c("BaseMean","betaSD_H","dispersion_H","betaSD_L","dispersion_L","mu_presort_H","mu_h","mu_presort_L","mu_l")
      groundTruth_rep <- as.data.frame(cbind(groundTruth_rep, betaS))

      groundTruth_rep$PRESORT = pre_sort
      groundTruth_rep$HIGH = high_sort
      groundTruth_rep$LOW = low_sort

      write.table(groundTruth_rep, file=paste0(subfolder_frac,"groundTruth_rep_",i,".txt"), sep="\t", quote=FALSE, row.names=FALSE)
      # separate by fraction:

      ll_reps[[paste0("sorting_fraction_",sorting_frac)]][[i]] = groundTruth_rep

      # plot a pair wise scatter : groundTruth_rep$PRESORT vs groundTruth_rep$HIGH; groundTruth_rep$PRESORT vs groundTruth_rep$LOW; groundTruth_rep$HIGH vs groundTruth_rep$LOW
      # colour by type scale x y as log2
      groundTruth_rep_plot = groundTruth_rep
      # make the positive and negative in type be on top of other points:
      groundTruth_rep_plot$type = factor(groundTruth_rep_plot$type, levels = rev(c("positive","negative","Nonessential","Nontargeting")))
      # order by levels:
      # make the colour scheme so that are all grey except for negative and positive that are blue and red respectively
      groundTruth_rep_plot$colour = "grey"
      groundTruth_rep_plot$colour[groundTruth_rep_plot$type == "positive"] = "red"
      groundTruth_rep_plot$colour[groundTruth_rep_plot$type == "negative"] = "blue"
      col_palette = c("grey","red","blue","grey")
      names(col_palette) = c("Nonessential","positive","negative","Nontargeting")
      
      # to ensure that Nonessential and Nontargeting are on the bottom we need to set the levels of the factor type
      groundTruth_rep_plot$type = factor(groundTruth_rep_plot$type, levels = c("Nonessential","Nontargeting","positive","negative"))

      data_special <- groundTruth_rep_plot[!(groundTruth_rep_plot$type %in% c("Nonessential", "Nontargeting")), ][,c("PRESORT","HIGH","LOW","type")]
      data_other <- groundTruth_rep_plot[groundTruth_rep_plot$type %in% c("Nonessential", "Nontargeting"), ][,c("PRESORT","HIGH","LOW","type")]
      
      p <- ggplot() +
        geom_point(data=data_other, aes(x=log2(PRESORT), y=log2(HIGH), colour=type), alpha=0.5) +
        geom_point(data=data_special, aes(x=log2(PRESORT), y=log2(HIGH), colour=type), alpha=0.5) +
        geom_point(alpha=0.5) + 
        theme_classic() + 
        xlab("PRESORT") + 
        ylab("HIGH") + 
        scale_x_continuous(expand=c(0,0)) + 
        scale_y_continuous(expand=c(0,0)) +
        scale_colour_manual(values=col_palette)
      ggsave(paste0(subfolder_frac,"04_groundTruth_rep_",i,"_PRESORT_vs_HIGH.png"), p, width=8, height=6, units="in", dpi=300);

      p <- ggplot() +
        geom_point(data=data_other, aes(x=log2(PRESORT), y=log2(LOW), colour=type), alpha=0.5) +
        geom_point(data=data_special, aes(x=log2(PRESORT), y=log2(LOW), colour=type), alpha=0.5) +
        geom_point(alpha=0.5) + 
        theme_classic() + 
        xlab("PRESORT") + 
        ylab("LOW") + 
        scale_x_continuous(expand=c(0,0)) + 
        scale_y_continuous(expand=c(0,0)) +
        scale_colour_manual(values=col_palette)
      ggsave(paste0(subfolder_frac,"04_groundTruth_rep_",i,"_PRESORT_vs_LOW.png"), p, width=8, height=6, units="in", dpi=300);

      p <- ggplot() +
        geom_point(data=data_other, aes(x=log2(HIGH), y=log2(LOW), colour=type), alpha=0.5) +
        geom_point(data=data_special, aes(x=log2(HIGH), y=log2(LOW), colour=type), alpha=0.5) +
        geom_point(alpha=0.5) + 
        theme_classic() + 
        xlab("HIGH") + 
        ylab("LOW") + 
        scale_x_continuous(expand=c(0,0)) + 
        scale_y_continuous(expand=c(0,0)) +
        scale_colour_manual(values=col_palette)
      ggsave(paste0(subfolder_frac,"04_groundTruth_rep_",i,"_HIGH_vs_LOW.png"), p, width=8, height=6, units="in", dpi=300);
    }
  }

  for(sorting_frac in names(ll_reps)){
    subfolder_frac = paste0(subfolder,"sorted_",gsub("sorting_fraction_","",sorting_frac),"/")
    ll_reps_sel = ll_reps[[sorting_frac]]
    Combo <- lapply(1:length(ll_reps_sel), function(x) {
        y = ll_reps_sel[[x]]
        colnames(y) = paste0(colnames(y), "_rep", x)
        return(y)
      })
      Combo = do.call(cbind, Combo)
      write.table(Combo, file=paste0(subfolder_frac,"groundTruth_rep_all.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  }
}


