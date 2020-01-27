# Functions for plotting networks with respect to binary phenotype

# Function to create igraph object for case/control
create_cc_igraph  = function(w.adj0, w.adj1, col, lay=NULL, ...) {
  require(igraph)
  
  palette = categorical_pal(max(col))
  
  control_graph = graph.adjacency(w.adj0, mode="undirected", weighted=TRUE, diag = FALSE)
  E(control_graph)$color[E(control_graph)$weight>0] = "forestgreen"
  E(control_graph)$color[E(control_graph)$weight<=0] = "orangered"
  V(control_graph)$color = palette[col]
  # Will use layout for case network
  
  case_graph = graph.adjacency(w.adj1, mode="undirected", weighted=TRUE, diag = FALSE)
  E(case_graph)$color[E(case_graph)$weight>0] = "forestgreen"
  E(case_graph)$color[E(case_graph)$weight<=0] = "orangered"
  V(case_graph)$color = palette[col]
  if (is.function(lay)) {
    l_case = lay(case_graph, ...)
  } else {
    l_case = lay
  }
  
  return(list(case=case_graph, control=control_graph, lay=l_case, pal=palette))
}


plotNetworks = function(counts, prec0, prec1, labs=NULL, scale_line_width=30, lay=layout.star, 
                        lab0="Group 0", lab1="Group 1", ci0=NULL, ci1=NULL, vertex.size0=50, vertex.size1=50,
                        col=NULL, seed=NULL, reorder=NULL, vertex.label.cex=NULL, cutoff=NULL, phyla_names=NULL,
                        prec=TRUE,cor=FALSE,cex.main=1,side.legend=TRUE, ...) {
  require(igraph)
  
  n.spec = NCOL(counts)
  
  if (is.null(labs)) {
    labs = substr(colnames(counts)[-n.spec], 4, nchar(colnames(counts)[-n.spec]))
  }
  if (is.null(col)) {
    col = rainbow(n.spec-1)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  sgn <- ifelse(prec, -1, 1)
  
  if (!is.null(ci0) & !is.null(ci1)) {
    sig.links0 = ci_to_adj(ci0[[1]],
                           ci0[[2]])
    sig.links1 = ci_to_adj(ci1[[1]],
                           ci1[[2]])
    partial_cor0 = sgn*cov2cor(prec0)*sig.links0
    partial_cor1 = sgn*cov2cor(prec1)*sig.links1
  } else {
    partial_cor0 = sgn*cov2cor(prec0)
    partial_cor1 = sgn*cov2cor(prec1)
  }
  
  if (!is.null(cutoff)){
    partial_cor0[abs(partial_cor0)<cutoff] = 0
    partial_cor1[abs(partial_cor1)<cutoff] = 0
  }
  
  diag(partial_cor0) = 0
  diag(partial_cor1) = 0
  
  palette = categorical_pal(max(col))
  
  control_graph = graph.adjacency(partial_cor0, mode="undirected", weighted=TRUE, diag = FALSE)
  E(control_graph)$color[E(control_graph)$weight>0] = "forestgreen"
  E(control_graph)$color[E(control_graph)$weight<=0] = "orangered"
  V(control_graph)$color = palette[col]
  if (is.function(lay)) {
    l_control = lay(control_graph)
  } else {
    l_control = lay
  }
  
  case_graph = graph.adjacency(partial_cor1, mode="undirected", weighted=TRUE, diag = FALSE)
  E(case_graph)$color[E(case_graph)$weight>0] = "forestgreen"
  E(case_graph)$color[E(case_graph)$weight<=0] = "orangered"
  V(case_graph)$color = palette[col]
  if (is.function(lay)) {
    l_case = lay(case_graph, ...)
  } else {
    l_case = lay
  }
  
  unique_col = palette[col[!duplicated((col))]]
  
  layout(matrix(c(1,2,3,1,2,3,4,4,4), ncol=3, byrow=TRUE), heights=c(8, 0.5), widths=c(6,6,2.5,2))
  par(mai=rep(0.2, 4))
  plot(control_graph, layout=l_case, 
              edge.width=abs(E(control_graph)$weight)*scale_line_width,
              vertex.size=vertex.size0, vertex.shape="circle",
              vertex.label=labs, vertex.label.color="black", vertex.label.cex=vertex.label.cex)
  title(lab0, cex.main=cex.main)
  par(mai=rep(0.2, 4))
  plot(case_graph, layout=l_case, 
              edge.width=abs(E(case_graph)$weight)*scale_line_width,
              vertex.size=vertex.size1, vertex.shape="circle",
              vertex.label=labs, vertex.label.color="black", vertex.label.cex=vertex.label.cex)
  title(lab1, cex.main=cex.main)
  plot.new()
  par(mai=c(0,0,0,0))
  if (side.legend) {
  legend("center", legend=phyla_names[!duplicated(phyla_names)], pch = 21, pt.cex=6, 
         col="black", pt.bg=unique_col, title = "Phylum", cex=1, box.col="white", xjust=1, 
         y.intersp = 4, x.intersp = 3)
  }
  plot.new()
  legend("center", legend=c("Positive assoc.", "Negative assoc."), 
         lty=c(1,1), col=c("forestgreen","orangered"), ncol=2, cex=2, lwd = 6, box.col="white")
  
  return(l_case)
}

radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}

