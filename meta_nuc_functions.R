require(reshape2)
require(plyr)
require(gridExtra)
require(ggplot2)

q <- function(...) {
  sapply(match.call()[-1], deparse)
}

meta_plot <- function(df, ref="WT", test, palette, subset=NA) {
  
  # if only want to include subset of data
  if(!is.na(subset[[1]])) {
    df <- subset(df, df$gene_id %in% subset)
  }
  
  df$cond2 <- paste(df$cond,df$rep,sep="_")
  
  ## plot average read coverages for nuc categories
  data <- df[,-c(1:11,(ncol(df)-4):(ncol(df)-1))]
  data <- melt(data, id.vars=q(cond2,nuc_category) )
  data <- na.omit(data)
  data <- ddply(data, .(cond2,nuc_category,variable), summarize,  
                N=length2(value,na.rm=T), 
                mean=round(mean(value,na.rm=T),3),  
                sd=round(sd(value,na.rm=T),3), 
                se=round(sd/sqrt(N),3), 
                ci=round(se * qt(.95/2 + .5, N-1),3) )
  data$variable <- as.numeric(as.character(data$variable))
  min = min(data$variable); max = max(data$variable)
  nuc_plot <- ggplot(data=data, aes(x=variable, y=mean, group=cond2, colour=cond2))+geom_line()+facet_wrap(~nuc_category,nrow=1)
  nuc_plot <- nuc_plot+scale_color_manual(values=palette)+geom_ribbon(aes(ymin=mean-ci, ymax=mean+ci, fill=cond2,color=NA),alpha=0.2)
  nuc_plot <- nuc_plot+scale_fill_manual(values=palette)+theme_bw()+scale_x_continuous(breaks=c(min,0,max))
  
  data <- melt(df,id.vars=q(cond,rep,nuc_category),measure.vars=q(pos,occ,dist,size))
  stats <- data.frame()
  for (nuc_cat in unique(data$nuc_category)) {
    for (param in q(pos,occ,dist,size)) {
      for (rep in unique(data$rep)) {
        population <- data[data$variable==param & data$nuc_category==nuc_cat,]
        population <- population[population$rep==rep,]
        ref_data = population[population$cond==ref,]
        test_data = population[population$cond==test,]
        mean_diff = mean(test_data$value,na.rm=T)-mean(ref_data$value,na.rm=T)
        sd2 = (sd(population$value, na.rm=T))^2 
        n1 = length(na.omit(ref_data$value))
        n2 = length(na.omit(test_data$value))
        sed=sqrt( (sd2/n1)+(sd2/n2) )
        temp <- data.frame(variable=param, nuc_category=nuc_cat, rep=rep, mean=mean_diff, sed=sed)
        stats <- rbind( stats, temp )
      }
    }
  }
  stats
  
  p1 <- ggplot(data=subset(stats,stats$variable=="pos"),aes(x=nuc_category,y=mean, group=rep, colour=rep))+ylab("Pos. Change (bp)")+expand_limits(y=0)+ scale_y_continuous(expand=c(0.15,0),labels=fmt())+scale_colour_manual(values=palette)+xlab("")
  p1 <- p1+theme_bw()+geom_linerange(aes(min=mean-sed,max=mean+sed),size=2,position=position_dodge(.25))+geom_point(size=4,position=position_dodge(.25))+geom_point(col="dimgrey",position=position_dodge(.25),size=2)+geom_hline(yintercept=0,linetype=2,colour="red")
  p1 <- p1+theme(axis.text.x=element_blank(),axis.ticks=element_blank(),plot.margin = rep(unit(0,"null"),4),panel.margin = unit(0,"null"))
  p2 <- ggplot(data=subset(stats,stats$variable=="occ"),aes(x=nuc_category,y=mean, group=rep, colour=rep))+ylab("Occ. Change")+expand_limits(y=0)+ scale_y_continuous(expand=c(0.15,0),labels=fmt())+scale_colour_manual(values=palette)+xlab("")
  p2 <- p2+theme_bw()+geom_linerange(aes(min=mean-sed,max=mean+sed),size=2,position=position_dodge(.25))+geom_point(size=4,position=position_dodge(.25))+geom_point(col="dimgrey",position=position_dodge(.25),size=2)+geom_hline(yintercept=0,linetype=2,colour="red")
  p2 <- p2+theme(axis.text.x=element_blank(),axis.ticks=element_blank(),plot.margin = rep(unit(0,"null"),4),panel.margin = unit(0,"null"))
  p3 <- ggplot(data=subset(stats,stats$variable=="dist"),aes(x=nuc_category,y=mean, group=rep, colour=rep))+ylab("Dist. Change")+expand_limits(y=0)+ scale_y_continuous(expand=c(0.15,0),labels=fmt())+scale_colour_manual(values=palette)+xlab("")
  p3 <- p3+theme_bw()+geom_linerange(aes(min=mean-sed,max=mean+sed),size=2,position=position_dodge(.25))+geom_point(size=4,position=position_dodge(.25))+geom_point(col="dimgrey",position=position_dodge(.25),size=2)+geom_hline(yintercept=0,linetype=2,colour="red")
  p3 <- p3+theme(axis.text.x=element_blank(),axis.ticks=element_blank(),plot.margin = rep(unit(0,"null"),4),panel.margin = unit(0,"null"))
  p4 <- ggplot(data=subset(stats,stats$variable=="size"),aes(x=nuc_category,y=mean, group=rep, colour=rep))+ylab("Size Change (bp)")+expand_limits(y=0)+ scale_y_continuous(expand=c(0.15,0),labels=fmt())+scale_colour_manual(values=palette)+xlab("")
  p4 <- p4+theme_bw()+geom_linerange(aes(min=mean-sed,max=mean+sed),size=2,position=position_dodge(.25))+geom_point(size=4,position=position_dodge(.25))+geom_point(col="dimgrey",position=position_dodge(.25),size=2)+geom_hline(yintercept=0,linetype=2,colour="red")
  p4 <- p4+theme(axis.text.x=element_blank(),axis.ticks=element_blank(),plot.margin = rep(unit(0,"null"),4),panel.margin = unit(0,"null"))
  
  ## to align axes properly...
  plots <- list(p1,p2,p3,p4)
  grobs <- list()
  widths <- list()
  for (i in 1:length(plots)){
    grobs[[i]] <- ggplotGrob(plots[[i]])
    widths[[i]] <- grobs[[i]]$widths[2:5]
  }
  maxwidth <- do.call(grid::unit.pmax, widths)
  for (i in 1:length(grobs)){
    grobs[[i]]$widths[2:5] <- as.list(maxwidth)
  }
  params <- do.call("grid.arrange", c(grobs, ncol = 1))
  
  # final
  grid.arrange(nuc_plot, params, ncol=1)
}

fmt <- function(){
  f <- function(x) as.character(round(x,2))
  f
}

length2 <- function (x, na.rm=FALSE) {
  if (na.rm) sum(!is.na(x))
  else       length(x)
}

# General function for loading data
load_data <- function(regex, sep, path=getwd(), header=T, col.names=col.names, check.names=T) {
  files <- list.files(path=path, pattern=regex, full.names=T, recursive=FALSE)
  results <- list()
  for (i in 1:length(files)) {
    basename=basename(files[i])
    condition <- sub(regex,"\\1",basename)
    data <- read.delim(files[i],header=header,sep=sep,col.names=col.names,check.names=check.names,row.names=NULL)
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
    if ( cond == ref ) {
      for ( j in 1:length(samples[[exp]]) ) {
        cond = samples[[exp]][j]
        if (!is.null(cond)) results[[ cond ]] = rbind(results[[cond]], data_copy)
      }
    } else if ( cond != ref ) {
      results[[ cond ]] = rbind( results[[cond]], data_copy )
    }
  }
  # add element for faceting (to group WT with samples)
  for (i in 1:length(results)) { 
    results[[i]] = cbind( results[[i]], set=rep(names(results)[i],nrow(results[[i]])) )
  }
  return(results)
} 


##### metplot of just average nuc profile split by gene sets
meta_plot_sets <- function(df, ref="WT", test, palette, set) {
  
  df$cond2 <- paste(df$cond,df$rep,sep="_")
  df$set <- rep(NA)
  
  for (i in 1:length(levels(set$set))) {
    name <- levels(set$set)[i] 
    df[df$gene_id %in% set[set$set==name,]$gene_id,]$set <- name  
  }
  df$set <- as.factor(df$set)
  
  ## plot average read coverages for nuc categories
  data <- df[,-c(1:11,(ncol(df)-4):(ncol(df)-2))]
  data <- melt(data, id.vars=q(cond2,nuc_category,set) )
  data <- na.omit(data)
  data <- ddply(data, .(cond2,nuc_category,variable,set), summarize,  
                N=length2(value,na.rm=T), 
                mean=round(mean(value,na.rm=T),3),  
                sd=round(sd(value,na.rm=T),3), 
                se=round(sd/sqrt(N),3), 
                ci=round(se * qt(.95/2 + .5, N-1),3) )
  print(head(data))
  data$variable <- as.numeric(as.character(data$variable))
  min = min(data$variable); max = max(data$variable)
  nuc_plot <- ggplot(data=data, aes(x=variable, y=mean, group=cond2, colour=cond2))+geom_line()+facet_wrap(~set+nuc_category,ncol=11)
  nuc_plot <- nuc_plot+scale_color_manual(values=palette)+geom_ribbon(aes(ymin=mean-ci, ymax=mean+ci, fill=cond2),alpha=0.2)
  nuc_plot <- nuc_plot+scale_fill_manual(values=palette)+theme_bw()+scale_x_continuous(breaks=c(min,0,max))
  return(nuc_plot)
}


# ## plot of quartiles
# meta_plot <- function(df, palette) {
#   
#   ## plot average read coverages for nuc categories
#   data <- df[,-c(1:11,35:36)]
#   data <- melt(data, id.vars=q(condition,nuc_category,set) )
#   data <- na.omit(data)
#   data <- data[data$set %in% c(1,4),]
#   data <- ddply(data, .(condition,nuc_category,set,variable), summarize,  
#                 N=length2(value,na.rm=T), 
#                 mean=round(mean(value,na.rm=T),3),  
#                 sd=round(sd(value,na.rm=T),3), 
#                 se=round(sd/sqrt(N),3), 
#                 ci=round(se * qt(.95/2 + .5, N-1),3) )
#   data$variable <- as.numeric(as.character(data$variable))
#   min = min(data$variable); max = max(data$variable)
#   nuc_plot <- ggplot(data=data, aes(x=variable, y=mean, group=interaction(condition, set), colour = interaction(condition, set)))
#   (nuc_plot <- nuc_plot+geom_line()+facet_wrap(~nuc_category,nrow=1))
#   #nuc_plot <- nuc_plot+scale_color_manual(values=palette)+geom_ribbon(aes(ymin=mean-ci, ymax=mean+ci, fill=cond2,color=NA),alpha=0.2)
#   #nuc_plot <- nuc_plot+scale_fill_manual(values=palette)+theme_bw()+scale_x_continuous(breaks=c(min,0,max))
#   
#   data <- melt(df,id.vars=q(cond,rep,nuc_category),measure.vars=q(pos,occ,dist,size))
#   stats <- data.frame()
#   for (nuc_cat in unique(data$nuc_category)) {
#     for (param in q(pos,occ,dist,size)) {
#       for (rep in unique(data$rep)) {
#         population <- data[data$variable==param & data$nuc_category==nuc_cat,]
#         population <- population[population$rep==rep,]
#         ref_data = population[population$cond==ref,]
#         test_data = population[population$cond==test,]
#         mean_diff = mean(test_data$value,na.rm=T)-mean(ref_data$value,na.rm=T)
#         sd2 = (sd(population$value, na.rm=T))^2 
#         n1 = length(na.omit(ref_data$value))
#         n2 = length(na.omit(test_data$value))
#         sed=sqrt( (sd2/n1)+(sd2/n2) )
#         temp <- data.frame(variable=param, nuc_category=nuc_cat, rep=rep, mean=mean_diff, sed=sed)
#         stats <- rbind( stats, temp )
#       }
#     }
#   }
#   stats
#   
#   p1 <- ggplot(data=subset(stats,stats$variable=="pos"),aes(x=nuc_category,y=mean, group=rep, colour=rep))+ylab("Pos. Change (bp)")+expand_limits(y=0)+ scale_y_continuous(expand=c(0.15,0),labels=fmt())+scale_colour_manual(values=palette)+xlab("")
#   p1 <- p1+theme_bw()+geom_linerange(aes(min=mean-sed,max=mean+sed),size=2,position=position_dodge(.25))+geom_point(size=4,position=position_dodge(.25))+geom_point(col="dimgrey",position=position_dodge(.25),size=2)+geom_hline(yintercept=0,linetype=2,colour="red")
#   p1 <- p1+theme(axis.text.x=element_blank(),axis.ticks=element_blank(),plot.margin = rep(unit(0,"null"),4),panel.margin = unit(0,"null"))
#   p2 <- ggplot(data=subset(stats,stats$variable=="occ"),aes(x=nuc_category,y=mean, group=rep, colour=rep))+ylab("Occ. Change")+expand_limits(y=0)+ scale_y_continuous(expand=c(0.15,0),labels=fmt())+scale_colour_manual(values=palette)+xlab("")
#   p2 <- p2+theme_bw()+geom_linerange(aes(min=mean-sed,max=mean+sed),size=2,position=position_dodge(.25))+geom_point(size=4,position=position_dodge(.25))+geom_point(col="dimgrey",position=position_dodge(.25),size=2)+geom_hline(yintercept=0,linetype=2,colour="red")
#   p2 <- p2+theme(axis.text.x=element_blank(),axis.ticks=element_blank(),plot.margin = rep(unit(0,"null"),4),panel.margin = unit(0,"null"))
#   p3 <- ggplot(data=subset(stats,stats$variable=="dist"),aes(x=nuc_category,y=mean, group=rep, colour=rep))+ylab("Dist. Change")+expand_limits(y=0)+ scale_y_continuous(expand=c(0.15,0),labels=fmt())+scale_colour_manual(values=palette)+xlab("")
#   p3 <- p3+theme_bw()+geom_linerange(aes(min=mean-sed,max=mean+sed),size=2,position=position_dodge(.25))+geom_point(size=4,position=position_dodge(.25))+geom_point(col="dimgrey",position=position_dodge(.25),size=2)+geom_hline(yintercept=0,linetype=2,colour="red")
#   p3 <- p3+theme(axis.text.x=element_blank(),axis.ticks=element_blank(),plot.margin = rep(unit(0,"null"),4),panel.margin = unit(0,"null"))
#   p4 <- ggplot(data=subset(stats,stats$variable=="size"),aes(x=nuc_category,y=mean, group=rep, colour=rep))+ylab("Size Change (bp)")+expand_limits(y=0)+ scale_y_continuous(expand=c(0.15,0),labels=fmt())+scale_colour_manual(values=palette)+xlab("")
#   p4 <- p4+theme_bw()+geom_linerange(aes(min=mean-sed,max=mean+sed),size=2,position=position_dodge(.25))+geom_point(size=4,position=position_dodge(.25))+geom_point(col="dimgrey",position=position_dodge(.25),size=2)+geom_hline(yintercept=0,linetype=2,colour="red")
#   p4 <- p4+theme(axis.text.x=element_blank(),axis.ticks=element_blank(),plot.margin = rep(unit(0,"null"),4),panel.margin = unit(0,"null"))
#   
#   ## to align axes properly...
#   plots <- list(p1,p2,p3,p4)
#   grobs <- list()
#   widths <- list()
#   for (i in 1:length(plots)){
#     grobs[[i]] <- ggplotGrob(plots[[i]])
#     widths[[i]] <- grobs[[i]]$widths[2:5]
#   }
#   maxwidth <- do.call(grid::unit.pmax, widths)
#   for (i in 1:length(grobs)){
#     grobs[[i]]$widths[2:5] <- as.list(maxwidth)
#   }
#   params <- do.call("grid.arrange", c(grobs, ncol = 1))
#   
#   # final
#   grid.arrange(nuc_plot, params, ncol=1)
# }

fmt <- function(){
  f <- function(x) as.character(round(x,2))
  f
}