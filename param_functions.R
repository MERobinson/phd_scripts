# param_functions.R

require(reshape2)
require(ggplot2)
require(plyr)
require(gridExtra)

###### General function for quoting
q <- function(...) {
  sapply(match.call()[-1], deparse)
}

###### General function for loading data
load_data <- function( regex, sep="\t", header=T, path=getwd(), check.names=T ) {
  
  files <- list.files(path=path, pattern=regex, full.names=T, recursive=FALSE)
  results <- list()
  for (i in 1:length(files)) {
    basename=basename(files[i])
    condition <- sub(regex,"\\1",basename)
    data <- read.delim(files[i],header=header,sep=sep,check.names=check.names)
    results[[condition]] = data
  }
  return(results)
}

###### General functions for assigning multiple variables

# Generic form
'%=%' = function(l, r, ...) UseMethod('%=%')

# Binary Operator
'%=%.lbunch' = function(l, r, ...) {
  Envir = as.environment(-1)
  
  if (length(r) > length(l))
    warning("RHS has more args than LHS. Only first", length(l), "used.")
  
  if (length(l) > length(r))  {
    warning("LHS has more args than RHS. RHS will be repeated.")
    r <- extendToMatch(r, l)
  }
  
  for (II in 1:length(l)) {
    do.call('<-', list(l[[II]], r[[II]]), envir=Envir)
  }
}

# Used if LHS is larger than RHS
extendToMatch <- function(source, destin) {
  s <- length(source)
  d <- length(destin)
  
  # Assume that destin is a length when it is a single number and source is not
  if(d==1 && s>1 && !is.null(as.numeric(destin)))
    d <- destin
  
  dif <- d - s
  if (dif > 0) {
    source <- rep(source, ceiling(d/s))[1:d]
  }
  return (source)
}

# Grouping the left hand side
g = function(...) {
  List = as.list(substitute(list(...)))[-1L]
  class(List) = 'lbunch'
  return(List)
}

##### Combine data by condition
comb_cond <- function( data, ref="WT" ) {

  names <- names(data)
  samples <- list()
  results <- list()
  
  # initiate data frame for each condition and record experiment (for reference data)
  for (i in 1:length(data)) {
    g(cond,rep,exp) %=% unlist(strsplit(names[i],"_")) # split string and assign to seperate variables
    if ( cond != ref ) {
      if ( exp %in% names(samples) ) {
        samples[[ exp ]] = c( samples[[exp]], cond )
      } else {
        samples[[ exp ]] = cond
      }
      
      results[[ cond ]] = data.frame()
    } 
  }
  samples <- lapply(samples,unique) # remove duplicates
  
  # merge data by 
  for (i in 1:length(data)) {
    g(cond,repl,exp) %=% unlist(strsplit(names[i],"_"))
    data_copy <- cbind( data[[i]], cond=rep(cond,nrow(data[[i]])),
                        rep=rep(repl,nrow(data[[i]])),
                        exp=rep(exp,nrow(data[[i]])) )
    # if wt, add to all cond from that exp
    if ( cond == ref) {
      for ( j in 1:length(samples[[exp]]) ) {
        cond = samples[[exp]][j]
        results[[ cond ]] = rbind( results[[cond]], data_copy  )
      }
    } else {
      results[[ cond ]] = rbind( results[[cond]], data_copy )
    }
  }
  # add element for faceting (to group WT with samples)
  for (i in 1:length(results)) { 
    results[[i]] = cbind( results[[i]], set=rep(names(results)[i],nrow(results[[i]])) )
  }
  # melt
  for (i in 1:length(results)) {
    results[[i]]$pos <- abs(results[[i]]$pos)
    results[[i]] <- melt(results[[i]],
                         id.vars=c("chr","start","end","dyad","canonical","gene_id","nuc_category","cond","rep","exp","set"),
                         measure.vars=q(pos,occ,dist,size),
                         variable.name="param")
  } 
  return(results)
}  

##### summarise parameters per condition
require(effsize)
summary <- function(data,abs=F) {
  
  results <- list()
  for (i in 1:length(data)) {
    temp_data <- data[[i]]
    if (abs==TRUE) {
      temp_data$value <- abs(temp_data$value)
    }
    cond_summary <- ddply( temp_data, .(cond,param), summarize, 
                           N=sum(!is.na(value)), Mean=round(mean(value,na.rm=T),3),
                           SD=round(sd(value,na.rm=T),3), SE=SD/sqrt(N), CV = SD/Mean)
    test <- ddply( temp_data, .(param), summarise, Cohens_d=cohen.d(d=value,f=cond,na.rm=T)$estimate,
                   effsize=cohen.d(d=value,f=cond,na.rm=T)$magnitude)
    cond_summary$Cohens_d = c(test$Cohens_d,rep(NA,4))
    cond_summary$effsize = c(test$effsize,rep(NA,4))
    results[[ names(data)[i] ]] = cond_summary
  }
  return(results)
} 

##### plot density
require(dplyr)
plot_density <- function(data,palette,levels) {
  
  combined <- data.frame()
  for (cond in names(data)) {
    temp <- data[[cond]]
    limits <- data.frame(param = levels(temp$param),
                         min = c(0,0,0,120),
                         max = c(25,600,15,180))
    temp <- inner_join(temp, limits) %>% filter(value > min, value < max)
    temp$colour <- paste(temp$cond,temp$rep,sep=" ")
    combined <- rbind(combined,temp)
  }
  combined$colour <- factor(combined$colour, levels=levels)
  combined$set2 <- paste(combined$set,combined$param,sep=" ")
  levels <- list()
  for (set in unique(combined$set)) {
    temp <- list( paste0(set," pos"), paste0(set," occ"), paste0(set," dist"), paste0(set," size"))
    levels <- c(levels,temp)
  }
  combined$set2 <- factor(as.factor(combined$set2), levels=levels)
  plot <- ggplot(data=combined, aes(x=value,fill=colour))+geom_density(adjust=1,alpha=.5,show.legend=T)+facet_wrap(~set2,scales="free",ncol=4)
  plot <- plot+theme_bw()+scale_fill_manual(values=palette)
  return(plot)
}

##### plot scatter
plot_scatter <- function(data,ref,high="white",low="blue",levels) {
  combined <- data.frame()
  for (cond in names(data)) {
    temp <- data[[cond]]
    ref_sub <- subset(temp,temp$Condition=="WT")
    test_sub <- subset(temp,temp$Condition==cond)
    merged <- merge(ref_sub,test_sub,by=q(Chrn,Canonical_Position,Replicate,Experiment,Set,Parameter))
    merged$colour <- paste(merged$Set,merged$Replicate,sep=" ")
    combined <- rbind(combined,merged)
  }
  combined$colour <- factor(combined$colour, levels=levels)
  combined$set2 <- paste(combined$Set,combined$Parameter,sep=" ")
  levels <- list()
  for (set in unique(combined$Set)) {
    temp <- list( paste0(set," Position"), paste0(set," Occupancy"), paste0(set," Distribution"), paste0(set," Size"))
    levels <- c(levels,temp)
  }
  combined$set2 <- factor(as.factor(combined$set2), levels=levels)
  plot <- ggplot(data=combined, aes(x=value.x,y=value.y))+facet_wrap(~set2,scales="free",ncol=4)+geom_smooth( method="lm", se=F)
  plot <- plot+theme_bw()+stat_binhex(na.rm=T)+scale_fill_gradient(trans="log",low="goldenrod4",high="goldenrod1")
  return(plot)
}

plot_scatter2 <- function(data,ref,palette) {
  plots <- list()
  for (i in 1:length(data)) {
    temp <- data[[i]]
    name <- names(data)[i]
    ref_sub <- subset(temp,temp$Condition=="WT")
    test_sub <- subset(temp,temp$Condition==name)
    merged <- merge(ref_sub,test_sub,by=q(Chrn,Canonical_Position,Replicate,Experiment,Set,Parameter))
    merged$colour <- paste(merged$Set,merged$Replicate,sep=" ")
    combined <- merged
    combined$set2 <- paste(combined$Set,combined$Parameter,sep=" ")
    levels <- list( paste0(name," Position"), paste0(name," Occupancy"), paste0(name," Distribution"), paste0(name," Size") )
    combined$set2 <- factor(as.factor(combined$set2), levels=levels)
    high = palette[(i*2)]
    low = palette[(i*2)-1]
    plot <- ggplot(data=combined, aes(x=value.x,y=value.y))+facet_wrap(~set2,scales="free",ncol=4)
    plot <- plot+theme_bw()+stat_binhex(na.rm=T,bins=30)+scale_fill_gradient(trans="log",low=low,high=high)
    plots[[name]] = plot
  }
  plots[["nrow"]] = length(data)
  do.call(grid.arrange, plots)
}

##### peak-to-peak/NRL
repeat_length <- function(data) {
  results <- data.frame()
  for (cond in names(data)) {
    temp <- data[[cond]]
    temp <- temp[grep("^chr[0-9]+$", temp$Chrn), ]
    temp$pos <- temp$Canonical_Position+temp$Position
    temp$nrl <- c( NA, tail(temp$pos,-1)-head(temp$pos,-1) )
    print(head(temp))
    print(tail(temp))
    print(range(temp$nrl))
    res <- data.frame( Condition=cond, NRL=mean(temp$nrl,na.rm=T), SD=round(sd(temp$nrl,na.rm=T),2), SE=sd(temp$nrl,na.rm=T)/sqrt(nrow(temp)-1) )
    results <- rbind( results, res )
    #results[[cond]] = ddply(data[[cond]], .(Condition,Replicate), summarise,  NRL=mean(head(Canonical_Position,-1)-tail(value,-1),na.omit=T) )
  }
  return(results)
}

##### finding genic overlaps
suppressMessages(library(ChIPpeakAnno))
