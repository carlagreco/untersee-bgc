import json
import csv
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description='Extract biosynthetic genes from antiSMASH JSON output.')
parser.add_argument('input_file', type=str, help='Path to the antiSMASH JSON input file')
parser.add_argument('output_file', type=str, help='Path to save the biosynthetic genes TSV output file')
args = parser.parse_args()

# Load the antiSMASH JSON file
input_file = args.input_file
output_file = args.output_file

# Open and parse the JSON file
with open(input_file, 'r') as file:
    data = json.load(file)

# Prepare to collect biosynthetic gene information
biosynthetic_genes = []

# Create a dictionary to map gene IDs to their PFAM domains
gene_to_pfams = {}

# First, collect all PFAM domains and associate them with gene IDs
for record in data.get('records', []):
    for feature in record.get('features', []):
        if feature.get('type') == 'PFAM_domain':
            qualifiers = feature.get('qualifiers', {})
            gene_id = qualifiers.get('label', ['<unknown>'])[0]
            pfam_id = qualifiers.get('db_xref', ['<unknown>'])[0]
            
            if gene_id not in gene_to_pfams:
                gene_to_pfams[gene_id] = []
            
            gene_to_pfams[gene_id].append(pfam_id)

# Iterate through records to find relevant regions and genes
for record in data.get('records', []):
    cluster_category = None
    cluster_genes = []

    for feature in record.get('features', []):
        if feature.get('type') == 'region':
            # Extract cluster category from the "product" field
            cluster_category = ', '.join(feature.get('qualifiers', {}).get('product', ['Unknown']))

        if feature.get('type') == 'CDS':
            qualifiers = feature.get('qualifiers', {})
            gene_id = qualifiers.get('ID', ['<unknown>'])[0]
            gene_kind = qualifiers.get('gene_kind', ['Unknown'])[0].lower()

            if gene_kind in ['biosynthetic', 'biosynthetic-additional']:
                # Get PFAM domains for this gene
                pfam_domains = gene_to_pfams.get(gene_id, [])
                
                cluster_genes.append({
                    'gene_id': gene_id,
                    'gene_kind': gene_kind,
                    'cluster_category': cluster_category,
                    'pfam_domains': ';'.join(pfam_domains)
                })

    # Only keep records with a protocluster and relevant genes
    if cluster_category and cluster_genes:
        biosynthetic_genes.extend(cluster_genes)

# Write the output to a TSV file
with open(output_file, 'w', newline='') as tsvfile:
    fieldnames = ['gene_id', 'gene_kind', 'cluster_category', 'pfam_domains']
    writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter='\t')

    writer.writeheader()
    writer.writerows(biosynthetic_genes)

print(f"Biosynthetic gene information has been saved to {output_file}")