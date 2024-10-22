### Plot interaction networks (downloaded from stringd)

library("igraph")
library("plotrix")

## Function for plotting an elliptical node
myellipse <- function(coords, v=NULL, params) {
  vertex.color = params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color = vertex.color[v]
  }
  vertex.size = 1/30 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size = vertex.size[v]
  }
  
  draw.ellipse(x=coords[,1], y=coords[,2],
               a = vertex.size, b=vertex.size/2, col=vertex.color)
}

## Register the shape with igraph
add_shape("ellipse", clip=shapes("circle")$clip, plot=myellipse)


### Load gene list
imp.gene.symbol = as.character(read.table("data/clinical/sev_important_genes_symbol_exp.txt")$V1)

### Load STRING interactions
string.network = read.table("data/clinical/STRING/string_interactions.tsv", sep="\t")
string.network.edges = cbind(as.character(string.network$V1), as.character(string.network$V2))
  
string.network.nodes = data.frame(gene=unique(as.character(c(string.network.edges[,1], string.network.edges[,2]))))
string.network.nodes$label = rep(0, length(string.network.nodes))

string.network.nodes$label[string.network.nodes$gene %in% imp.genes] = 1
string.network.nodes$gene = as.character(string.network.nodes$gene)

graph = graph.data.frame(string.network.edges, vertices = string.network.nodes, directed = F)
V(graph)$color = ifelse(V(graph)$label==1,"red","grey")
V(graph)$label.cex = 0.5

pdf("Figures/network.pdf", width=10, height=10)
plot(graph, vertex.label=V(graph)$name, vertex.shape="ellipse", vertex.size=2)
dev.off()
