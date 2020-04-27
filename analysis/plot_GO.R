### Plot p-value of GO analysis (downloaded from DAVID)

plot_pval <- function(GO.table, count=10){
  count = min(count, nrow(GO.table))
  pval = -log(GO.table$PValue)
  names(pval) = GO.table$Term
  barplot(pval[1:count], names.arg=NULL, horiz=T,las=2)
}



KEGG.table = read.table("data/clinical/sev_important_genes_KEGG.txt", header=T, sep="\t")

pdf("Figures/pval.KEGG.pdf", width=8, height=6)
plot_pval(KEGG.table)
dev.off()


BP.table = read.table("data/clinical/sev_important_genes_BP.txt", header=T, sep="\t")

pdf("Figures/pval.BP.pdf", width=8, height=6)
plot_pval(BP.table)
dev.off()
