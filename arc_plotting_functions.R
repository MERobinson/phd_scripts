## Functions for Arc-Plotting
#require(devtools)
#install_github('arcdiagram',  username='gastonstat')
require(arcdiagram)
require(devtools)
require(reshape2)
require(gridExtra)

comb = function(n, x) {
  return(factorial(n) / (factorial(x) * factorial(n-x)))
}

arc_params <- function (df,groups,palette) {
  
  # initialise variables
  count=0
  vlabels <- character(ncol(df))
  vgroup <- numeric(ncol(df))
  vfill <- character(ncol(df))
  vborders <- character(ncol(df))
  edgelist <- matrix(nrow=comb(ncol(df),2),ncol=2)
  values <- numeric(comb(ncol(df),2))
  degrees <- numeric(ncol(df))
  pval <- numeric(comb(ncol(df),2))
  line_group <- numeric(comb(ncol(df),2))
  table <- matrix(nrow=comb(ncol(df),2),ncol=4)
  
  # loop through reference columns
  for (i in 1:ncol(df)) {
    
    # get just relevent bits of colname as label
    name1 <- colnames(df)[i]
    vlabels[i] <- sub("^([pos|occ|dist|size])_.+$","\\U\\1",name1,perl=TRUE) ## took this out to manually add labels
    
    # set colours depending on sample
    vfill[i] = palette[i]
    vborders[i] = palette[i]
    vgroup[i] = groups[i]
    degrees[i] <- sum(df[,i])/20
    
    # loop through test columns
    for (j in (i+1):(ncol(df)) ) {
      
      if (j>ncol(df)) break
      
      name2 <- colnames(df)[j]
      count = count+1
      summed = df[,i] + df[,j]
      edgelist[count,1] = i
      edgelist[count,2] = j
      table[count,1] = name1
      table[count,2] = name2
      population = sum(df[,i])+sum(df[,j])
      values[count] = (sum(summed==2)/population)*100
      table[count,3] = sum(summed==2)
      line_group[count] = vgroup[i]
      
      # calc hypergeometric
      success = sum(summed==2)
      if (success==0) {
        p = 1
      } else {
        white_balls = sum(df[,i])
        black_balls = nrow(df)-white_balls
        draws = sum(df[,j])
        p = 1-phyper(success,white_balls,black_balls,draws) 
      }
      pval[count] = p
      table[count,4] = p
    }
  }
  graph_param <- list("vlabels"=vlabels,"vgroup"=vgroup,"edgelist"=edgelist,"values"=values,"degrees"=degrees,"vfill"=vfill,
                      "vborders"=vborders,"pval"=pval,"table"=table,"linegroup"=line_group) 
  return(graph_param)
}

arc_plotting <- function(param) {
  
  x = data.frame(param$vgroup, param$degrees, param$vlabels, ind=1:length(param$vlabels))
  layout(matrix(c(1,2),2,1))
  
  # plot arc diagram of relative overlaps
  arcplot(param$edgelist, ordering=x$ind, labels=param$vlabels, cex.labels=1,
          show.nodes=TRUE, col.nodes=param$vborders, bg.nodes=param$vfill,
          cex.nodes = (log(param$degrees)+1.5)/2, pch.nodes=21,
          lwd.nodes = 5, line=-0.5,
          col.arcs = hsv(0, 0, 0.2, 0.25), lwd.arcs = 1 * param$values)
  
  # adjust p values and plot significance
  padj <- p.adjust(param$pval,method="BH")
  col <- rep(0.2,length(padj))
  col <- replace(col,(padj<0.05),0.8)
  arcplot(param$edgelist, ordering=x$ind, labels=param$vlabels, cex.labels=1,
          show.nodes=TRUE, col.nodes=param$vborders, bg.nodes=param$vfill,
          cex.nodes = (log(param$degrees)+1.5)/2, pch.nodes=21,
          lwd.nodes = 5, line=-0.5,
          col.arcs = hsv(0, col, col, 0.25), lwd.arcs = 1 * -log10(padj+1e-20))
}