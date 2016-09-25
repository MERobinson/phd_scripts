# Script: perm_test.R  
# Author: Mark Robinson
# Date: Jan 2016
# Description: Function to calc probability of observing average param change accross genes

# General function for quoting
q <- function(...) {
  sapply(match.call()[-1], deparse)
}
if(!suppressMessages(require(doMC, warn.conflicts=F))) {
  install.packages('doMC')
  if(!suppressMessages(require(doMC, warn.conflicts=F))) {
    stop('Loading package doMC failed!')
  }
}
if(!suppressMessages(require(utils, warn.conflicts=F))) {
  install.packages('utils')
  if(!suppressMessages(require(utils, warn.conflicts=F))) {
    stop('Loading package utils failed!')
  }
}
if(!suppressMessages(require(plyr, warn.conflicts=F))) {
  install.packages('plyr')
  if(!suppressMessages(require(plyr, warn.conflicts=F))) {
    stop('Loading package plyr failed!')
  }
}
# Register doMC with CPU number.
registerDoMC()


# General function for loading data
load_data <- function(regex, sep="\t") {
  
  files <- list.files(path=getwd(), pattern=regex, full.names=T, recursive=FALSE)
  results <- list()
  for (i in 1:length(files)) {
    basename=basename(files[i])
    condition <- sub("^([A-Za-z]{3,4})_.+$","\\1",basename)
    data <- read.delim(files[i],header=T,sep=sep)
    results[[condition]] = data
  }
  return(results)
}

# p val for perm function
prob <- function(row_idx, param_idx, data, permuted, upper, min_n, perm) { 

  nuc_n <- data[row_idx,paste0(param_idx,"_n")] # number of nuc in current gene
  if (nuc_n < min_n) { # if less than minimum number of nuc, set as NA
    p = NA
  } else if (upper==T) { # otherwise calc p
    p = sum( permuted[,nuc_n] > data[row_idx,param_idx] )/( perm+1 )
  } else if (upper==F) {
    p = sum( permuted[,nuc_n] < data[row_idx,param_idx] )/( perm+1 )
  }
  
  # change zero's to lowest possible pval
  p <- ifelse (p == 0, 1/(perm+1), p)
}

# Permutation test function
perm_test <- function(perm, min_n, regex, genic_only=F) {
  
  print(paste0("Performing permuation test with ",perm," iterations"))
  print(paste0("Min nucleosome number = ",min_n))
  
  # list all files matching regex
  files <- list.files(path=getwd(), pattern=regex, full.names=T, recursive=FALSE)
  
  print("Files for testing:")
  print(files)
  
  # combine data by condition
  canonical <- data.frame()
  combined <- list(ChdA=data.frame(), ChdB=data.frame(), ChdC=data.frame(), Isw=data.frame())
  for (i in 1:length(files)) {
    basename=basename(files[i])
    condition <- sub("^([A-Za-z]{2,4})_.+$","\\1",basename)
    data <- read.table(files[[i]], header=T, sep='\t', 
                       colClasses=R.utils::colClasses("cnnnnnnnnncf---------------------"))
    data$cond <- rep(condition)
    if (condition=="WT") {
      exp <- sub("^.+_(Exp[A-Z])_.+$","\\1",basename)
      if (exp=="ExpA") {
        combined[["ChdC"]] <- rbind(combined[["ChdC"]], data)
      } else if (exp=="ExpB") {
        combined[["ChdA"]] <- rbind(combined[["ChdA"]], data)
        combined[["ChdB"]] <- rbind(combined[["ChdB"]], data)
      } else if (exp=="ExpC") {
        combined[["Isw"]] <- rbind(combined[["Isw"]], data)
      }
      canonical <- rbind(canonical, data)
    } else {
      combined[[condition]] <- rbind(combined[[condition]], data)
    }
  }
  
  # calculate average canonical nucleosome params
  canonical <- ddply(canonical, .(chr, canonical), summarize, occ=mean(occ,na.rm=T),
                     dist=mean(dist,na.rm=T), size=mean(size,na.rm=T) )
   
  results <- list()
  for (i in 1:length(combined)) {
    
    condition <- names(combined)[i]
    data <- combined[[i]]
    
    # make pos change absolute since direction unmeaningful
    data$pos <- abs(data$pos)
    
    # check if nucleosome changes have been calculated already 
    if (file.exists(paste0(condition,"_parameter_changes.txt"))) {
      data <- read.delim(paste0(condition,"_parameter_changes.txt"), sep="\t")
    } else {
      # convert other param values to change compared with canonical
      data <- foreach(j=1:nrow(data), .combine='rbind', .multicombine=T) %dopar% {
        row <- data[j,]
        matched <- canonical[canonical$chr==row$chr & canonical$canonical==row$canonical,]
        row[,7:9] <- row[,7:9]-matched[,3:5]
        row
      }
      write.table(data, paste0(condition,"_parameter_changes.txt"),sep="\t",row.names=F,quote=F)
    }
    
    # limit analysis to just genic nucleosomes
    if (genic_only==TRUE) {
      data <- subset(data, data$nuc_category > 1 & data$nuc_category < 9)
    }
    
    # seperate WT and test conditions
    genes <- ddply(data, .(gene_id,cond), summarise,
                   pos_n=sum(!is.na(pos)), occ_n=sum(!is.na(occ)),
                   dist_n=sum(!is.na(dist)), size_n=sum(!is.na(size)),
                   pos=mean(pos,na.rm=T), occ=mean(occ,na.rm=T),
                   dist=mean(dist,na.rm=T), size=mean(size,na.rm=T))
    write.table(genes, paste0(condition,"_gene_changes_summary.txt"), sep="\t", row.names=F, quote=F)
    ref <- genes[genes$cond=="WT",]
    test <- genes[genes$cond!="WT",]
    
    # for each param...
    for (param in q(pos, occ, dist, size)) {
      
      print(paste0("Permuting ",condition," ",param," values"))
      
      all_obs <- na.omit(data[,param]) # get all possible non-na values for current param 
      
      max_n <- max(test[,paste0(param,"_n")]) # get max number of nucleosomes
      
      permuted <- matrix(nrow = perm, ncol = max_n) # initialise matrix to store permuted averages
      
      # calc perm matrix
      permuted <- foreach(i=1:ncol(permuted), .combine='cbind', .multicombine=T) %dopar% {
        perm_res <- sapply(1:nrow(permuted), function(x) { mean(sample(all_obs, i, replace=T)) })
      }
                                
      # calc p values (no. of permutations with values >/< gene value divided by no. of perm)
      if(param=="pos") {
        test$pos_upper <- sapply(1:nrow(test), function(x) { prob(x, param, data=test, permuted=permuted, upper=T, min_n=min_n, perm=perm) })
      } else {
        test[,paste0(param,"_upper")] <- sapply(1:nrow(test), function(x) { 
          prob(x, param, data=test, permuted=permuted, upper=T, min_n=min_n, perm=perm) })
        test[,paste0(param,"_lower")] <- sapply(1:nrow(test), function(x) { 
          prob(x, param, data=test, permuted=permuted, upper=F, min_n=min_n, perm=perm) })
      }
    }
    
    # adjust pval for multiple comparisons
    test[,11:ncol(test)] <- sapply(test[,11:ncol(test)],function(x) { p.adjust(x,method="BH") })
    
    # write tables
    date <- Sys.Date()
    filename <- paste0(condition,"_permutation_test_results_",date,".txt")
    write.table(test,filename,sep="\t",row.names=F,col.names=T,quote=F)  
    
    # Add results to list
    results[[ condition ]] <- test
  } 
  
  return(results)
}

# Individual permutation test function
indi_perm_test <- function(condition, perm, min_n, regex, genic_only=F) {
  
  print(paste0("Performing permuation test with ",perm," iterations"))
  print(paste0("Min nucleosome number = ",min_n))
  
  # load data
  files <- list.files(path=getwd(), pattern=regex, full.names=T, recursive=FALSE)
  print("Files for testing:"); print(files)
  combined <- list()
  combined[[condition]] <- data.frame()
  canonical <- data.frame()
  for (i in 1:length(files)) {
    basename=basename(files[i])
    cond <- sub("^([A-Za-z]{2,5})_.+$","\\1",basename)
    data <- read.table(files[[i]], header=T, sep='\t', 
                       colClasses=R.utils::colClasses("cnnnnnnnnncf---------------------"))
    data$cond <- rep(cond)
    if (cond == condition) {
      combined[[cond]] <- rbind(combined[[cond]], data)
    } else if (cond == "WT") {
      exp <- sub("^.+_(Exp[A-Z])_.+$","\\1",basename)
      if (exp=="ExpC") {
        combined[[condition]] <- rbind(combined[[condition]], data)
      }
      canonical <- rbind(canonical, data)
    }
  }
  
  # calculate average canonical nucleosome params
  canonical <- ddply(canonical, .(chr, canonical), summarize, occ=mean(occ,na.rm=T),
                     dist=mean(dist,na.rm=T), size=mean(size,na.rm=T) )
  
  results <- list()
  for (i in 1:length(combined)) {
    
    cond <- names(combined)[i]
    data <- combined[[i]]
    
    # make pos change absolute since direction unmeaningful
    data$pos <- abs(data$pos)
    
    # check if nucleosome changes have been calculated already 
    if (file.exists(paste0(condition,"_parameter_changes.txt"))) {
      data <- read.delim(paste0(condition,"_parameter_changes.txt"), sep="\t")
    } else {
      # convert other param values to change compared with canonical
      data <- foreach(j=1:nrow(data), .combine='rbind', .multicombine=T) %dopar% {
        row <- data[j,]
        matched <- canonical[canonical$chr==row$chr & canonical$canonical==row$canonical,]
        row[,7:9] <- row[,7:9]-matched[,3:5]
        row
      }
      write.table(data, paste0(condition,"_parameter_changes.txt"),sep="\t",row.names=F,quote=F)
    }
    
    # limit analysis to just genic nucleosomes
    if (genic_only==TRUE) {
      data <- subset(data, data$nuc_category > 1 & data$nuc_category < 9)
    }
    
    # seperate WT and test conditions
    genes <- ddply(data, .(gene_id,cond), summarise,
                   pos_n=sum(!is.na(pos)), occ_n=sum(!is.na(occ)),
                   dist_n=sum(!is.na(dist)), size_n=sum(!is.na(size)),
                   pos=mean(pos,na.rm=T), occ=mean(occ,na.rm=T),
                   dist=mean(dist,na.rm=T), size=mean(size,na.rm=T))
    write.table(genes, paste0(condition,"_gene_changes_summary.txt"), sep="\t", row.names=F, quote=F)
    ref <- genes[genes$cond=="WT",]
    test <- genes[genes$cond!="WT",]
    
    # for each param...
    for (param in q(pos, occ, dist, size)) {
      
      print(paste0("Permuting ",condition," ",param," values"))
      
      all_obs <- na.omit(data[,param]) # get all possible non-na values for current param 
      
      max_n <- max(test[,paste0(param,"_n")]) # get max number of nucleosomes
      
      permuted <- matrix(nrow = perm, ncol = max_n) # initialise matrix to store permuted averages
      
      # calc perm matrix
      permuted <- foreach(i=1:ncol(permuted), .combine='cbind', .multicombine=T) %dopar% {
        perm_res <- sapply(1:nrow(permuted), function(x) { mean(sample(all_obs, i, replace=T)) })
      }
      
      # calc p values (no. of permutations with values >/< gene value divided by no. of perm)
      if(param=="pos") {
        test$pos_upper <- sapply(1:nrow(test), function(x) { prob(x, param, data=test, permuted=permuted, upper=T, min_n=min_n, perm=perm) })
      } else {
        test[,paste0(param,"_upper")] <- sapply(1:nrow(test), function(x) { 
          prob(x, param, data=test, permuted=permuted, upper=T, min_n=min_n, perm=perm) })
        test[,paste0(param,"_lower")] <- sapply(1:nrow(test), function(x) { 
          prob(x, param, data=test, permuted=permuted, upper=F, min_n=min_n, perm=perm) })
      }
    }
    
    # adjust pval for multiple comparisons
    test[,11:ncol(test)] <- sapply(test[,11:ncol(test)],function(x) { p.adjust(x,method="BH") })
    
    # write tables
    date <- Sys.Date()
    filename <- paste0(condition,"_permutation_test_results_",date,".txt")
    write.table(test,filename,sep="\t",row.names=F,col.names=T,quote=F)  
    
    # Add results to list
    results[[ condition ]] <- test
  } 
  
  return(results)
}