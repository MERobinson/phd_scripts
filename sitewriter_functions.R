## Functions for sitewriter plotting comparing multiple sites

q <- function(...) {
  sapply(match.call()[-1], deparse)
}

# sitewriter calling function - generic
sitewriter <- function(ref, test, sites, out) {
  system(paste(c("perl ~/Scripts/Perl/site_writer.pl -o",out,paste(c(ref,test),collapse=","),
                 paste(sites,collapse=",")),collapse=" "))
}

extract_sites <- function(data, sites, cond, param) {
  
  if (param == "expression") {
    up <- subset(data, data$padj < 0.05 & data$log2FoldChange > 1)$gene_id
    do <- subset(data, data$padj < 0.05 & data$log2FoldChange < -1)$gene_id
    mis <- sites[sites$gene_id %in% up,]
    not <- sites[sites$gene_id %in% do,]
    print(paste0("Number of upregulated genes: ",nrow(mis)))
    print(paste0("Number of downregulated genes: ",nrow(not)))
  } else {
    if (param == "all") {
      genes <- data[(data[,"pos_upper"]<0.05 & !is.na(data[,"pos_upper"]) |
                       data[,"occ_upper"]<0.05 & !is.na(data[,"occ_upper"]) |
                       data[,"occ_lower"]<0.05 & !is.na(data[,"occ_lower"]) |
                       data[,"dist_upper"]<0.05 & !is.na(data[,"dist_upper"]) |
                       data[,"dist_lower"]<0.05 & !is.na(data[,"dist_lower"]) |
                       data[,"size_upper"]<0.05 & !is.na(data[,"size_upper"]) |
                       data[,"size_lower"]<0.05 & !is.na(data[,"size_lower"]) ),]$gene_id
    } else {
      genes <- data[data[,param]<0.05 & !is.na(data[,param]),]$gene_id
    }
    mis <- sites[sites$gene_id %in% genes,]
    not <- sites[!sites$gene_id %in% genes,]
    print(paste0("Number of mismodeled genes: ",nrow(mis)))
    print(paste0("Number of non-mismodeled genes: ",nrow(not)))
  }
  dir.create(paste0(cond,"_",param),showWarnings=F)
  write.table(mis, paste0(cond,"_",param,"/par_mismodeled_sites.txt"), sep="\t", row.names=F, quote=F, col.names=F)
  write.table(not, paste0(cond,"_",param,"/non_mismodeled_sites.txt"), sep="\t", row.names=F, quote=F, col.names=F)
}

# list .sgr files
list_sgr <- function(conditions, path) {
  base=path
  dirs <- list.dirs(path = base, full.names = TRUE, recursive = FALSE)
  names(dirs) <- sapply(dirs, function(x) { condition <- sub("^.+/([A-Za-z0-9]{2,5})$", "\\1", x) })
  dirs <- dirs[names(dirs) %in% conditions]
  sgr <- lapply(dirs, function(dir) {
    condition <- sub("^.+/([A-Za-z0-9]{2,5})$", "\\1", dir)
    files <- list.files(path=dir, pattern="*_b5_nuc_map.sgr", full.names=T, recursive=F)
    if (condition=="WT") {
      files <- list("WT1"=c("rep1"=files[1], "rep2"=files[2]),
                    "WT2"=c("rep1"=files[3], "rep2"=files[5]),
                    "WT3"=c("rep1"=files[6], "rep2"=files[7]))
    } else if (condition=="Isw" | condition=="ChdAB") {
      files <- files[2:3]
      names(files) <- paste0("rep",seq_along(files))
    } else {
      names(files) <- paste0("rep",seq_along(files))
    }
    files
  })
  return(sgr)
}

# sitewriter callin function - misodeled
call_sitewriter <- function(sgr, set, out) {
  condition <- sub("^([A-Za-z]{2,5})_.+$", "\\1", set)
  sites <- list.files(path=set, pattern="*_sites.txt", full.names=T, recursive=F)
  names(sites) <- sapply(sites, function(x) { tmp <- sub("^.+(mis|par)_modeled_sites.txt$", "\\1", x) })
  if (condition=="ChdC") wt <- sgr[["WT"]][["WT1"]]
  if (condition=="ChdA" | condition=="ChdB") wt <- sgr[["WT"]][["WT2"]]
  if (condition=="Isw" | condition=="ChdAB") wt <- sgr[["WT"]][["WT3"]]
  cond <- sgr[[condition]]
  system(paste(c("perl ~/Scripts/Perl/site_writer.pl -o",set,paste(c(cond,wt),collapse=","),
                 paste(sites,collapse=",")),collapse=" "))
}

sitewriter_data <- function(set, split="sites") {
  files <- list.files(path=set, pattern="*_sitewriter.txt", full.names=T, recursive=F)
  data=list()
  for (i in seq_along(files)) {
    file <- files[[i]]
    name <- sub("^.+/([A-Za-z]{2,5}_(0H|rep[0-9])_Exp[A-Z])_.+.txt$", "\\1", file)
    cond <- sub("^([A-Za-z]{2,5})_.+$", "\\1", name)
    site <- sub("^.+/.+_((par|non)_mismodeled)_.+.txt$", "\\2", file)
    temp <- read.delim(file, sep="\t")
    if (i==1) {
      data[[1]] <- as.data.frame(temp[,1])
      colnames(data[[1]]) <- "bin"
      if (split!="none") {
        data[[2]] <- as.data.frame(temp[,1])
        colnames(data[[2]]) <- "bin" 
      }
    }
    if (split=="sites") {
      if (site=="par") {
        data[[1]][,name] <- temp[,2]
      } else if (site=="non") {
        data[[2]][,name] <- temp[,2]
      }
    } else if (split=="cond") {
      if (cond=="WT") {
        data[[1]][,name] <- temp[,2]
      } else {
        data[[2]][,name] <- temp[,2]
      }
    } else if (split=="none") {
      data[[1]][,name] <- temp[,2]
    }
  }
  return(data)
}

plot_sitewriter <- function(data, colors) {
  require(ggplot2)
  require(reshape2)
  require(gridExtra)
  parameters <- list(geom_rect(data=NULL,aes(xmin=-Inf,xmax=0,ymin=-Inf,ymax=Inf),fill="grey80",color="grey40"),
                     geom_rect(data=NULL,aes(xmin=0,xmax=Inf,ymin=-Inf,ymax=Inf),fill="white",colour="grey40"),
                     scale_y_continuous(name="Nucleosome Score", limits=c(0,10)), xlab("Distance to Gene Start (bp)"),
                     geom_line(size=1.3), scale_color_manual(values=colors),
                     theme(text = element_text(size=16, color="grey60"), axis.title=element_text(size=18, color="grey30")) )
  
  data1 <- melt(data[[1]], id.vars="bin")
  (p1 <-ggplot(data1, aes(x=bin,y=value,color=variable))+parameters)
  
  data2 <- melt(data[[2]], id.vars="bin")
  p2 <- ggplot(data2, aes(x=bin,y=value,color=variable))+parameters
  
  grid.arrange(p1,p2)
}

plot_sitewriter_indi <- function(data, colors) {
  require(ggplot2)
  require(reshape2)
  parameters <- list(geom_rect(data=NULL,aes(xmin=-Inf,xmax=0,ymin=-Inf,ymax=Inf),fill="grey80",color="grey40"),
                     geom_rect(data=NULL,aes(xmin=0,xmax=Inf,ymin=-Inf,ymax=Inf),fill="white",colour="grey40"),
                     scale_y_continuous(name="Nucleosome Score", limits=c(0,10)), xlab("Distance to Gene Start (bp)"),
                     geom_line(size=1.3), scale_color_manual(values=colors),
                     theme(text = element_text(size=16, color="grey60"), axis.title=element_text(size=18, color="grey30")) )
  
  data1 <- melt(data[[1]], id.vars="bin")
  (p1 <-ggplot(data1, aes(x=bin,y=value,color=variable))+parameters)
  return(p1)
}