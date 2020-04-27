#cat top_mapped_filtered.txt | cut -f1 | sort | uniq | sed '/GENBANK_ACCESSION/d'  > top_genes_merged.txt

while read line
do grep -w "${line}" ../GPL6244.annot | awk '{print $1}' | sort | uniq >> top_genes_probe.txt
 done < top_genes_merged.txt

cat top_genes_probe.txt | sort | uniq > top_genes_merged_probe.txt
