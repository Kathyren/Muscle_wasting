#!/bin/bash

USECASE=$1
FILE=$2
CSV_FILE="network_count.csv"

if [[ -z "$USECASE" || -z "$FILE" || ! -f "$FILE" ]]; then
    echo "Usage: $0 usecase file.cyj"
    exit 1
fi

echo "Analyzing file: $FILE (Use case: $USECASE)"
echo "--------------------------------"

# Total nodes
num_nodes=$(jq '.elements.nodes | length' "$FILE")
# Total edges
num_edges=$(jq '.elements.edges | length' "$FILE")
# miRNAs
num_mirnas=$(jq '[.elements.nodes[] | select(.type=="mirna" or .data.type=="mirna")] | length' "$FILE")
# Genes
num_genes=$(jq '[.elements.nodes[] | select(.type=="gene" or .data.type=="gene")] | length' "$FILE")
# DEGs
num_deg=$(jq '[.elements.nodes[] | select(
  (.data | to_entries[]? 
   | select(.key|startswith("dds_original_")) 
   | .value|tonumber 
   | if . > 2 or . < -2 then true else false end)
  )] | length' "$FILE")
# Proportion of DEGs
deg_ratio=$(awk -v d="$num_deg" -v n="$num_nodes" 'BEGIN { if (n==0) print 0; else printf "%.4f", d/n }')

# Pathways (top 3, comma-separated)
top_pathways=$(jq -r '.elements.nodes[].data.pathways[]? // empty' "$FILE" | \
  sort | uniq -c | sort -nr | head -n 3 | awk '{$1=""; sub(/^ /, ""); print}' | paste -sd ',' -)

# Density: 2E / N(N-1)
density=$(awk -v e="$num_edges" -v n="$num_nodes" 'BEGIN {
  if (n<=1) print 0; else printf "%.5f", (2*e)/(n*(n-1))
}')

# Print summary
echo "Number of nodes: $num_nodes"
echo "Number of edges: $num_edges"
echo "Number of miRNAs: $num_mirnas"
echo "Number of genes: $num_genes"
echo "Number of DEGs: $num_deg"
echo "Proportion of DEGs: $deg_ratio"
echo "Density: $density"
echo "Top pathways: $top_pathways"

# Add CSV header if it doesn't exist
if [[ ! -f "$CSV_FILE" ]]; then
    echo -e "Usecase\tFile\tNetwork Final size\tmiRNA amount\tGenes amount\tEdges\tDensity\tDEGs\tProportion of DEG\tMost common pathways" > "$CSV_FILE"
fi

# Append data
echo -e "$USECASE\t$(basename "$FILE")\t$num_nodes\t$num_mirnas\t$num_genes\t$num_edges\t$density\t$num_deg\t$deg_ratio\t$top_pathways" >> "$CSV_FILE"
