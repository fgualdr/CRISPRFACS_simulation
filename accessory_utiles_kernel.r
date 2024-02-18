
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
        df_subset$target_name <- gsub("_.*","",rownames(df_subset))
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
                df_subset$target_name <- gsub("_.*","",rownames(df_subset))
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
            # exit to any apply or lapply or bplapply
            return(NULL) 
            }else{
                if(h[2] <= log2(1.05)){
                    h[2] = log2(1.05) # as we binned the L2FC in 0.01 we set the bandwidth to 0.01 in case the estimate exceeds 0.01 - this control vectors dimensions
                }
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
        gc()
        return(NULL)
    }
