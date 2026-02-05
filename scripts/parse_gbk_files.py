import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

def extract_and_calculate_features(record, feature_type):
    """ Extract features and calculate CDS count and GC content specific to the feature type """
    features_list = []
    
    # Determine the correct qualifier key based on the feature type
    if feature_type == 'region':
        name_key = 'region_number'
    elif feature_type == 'candidate_cluster':
        name_key = 'candidate_cluster_number'
    elif feature_type == 'protocluster':
        name_key = 'protocluster_number'
    
    for feature in record.features:
        if feature.type == feature_type:
            start, end = feature.location.start, feature.location.end
            length = end - start + 1
            name = feature.qualifiers[name_key][0]  # Directly get the name
            num_cds = sum(1 for f in record.features if f.type == 'CDS' and start <= f.location.start <= f.location.end <= end)
            gc_content = gc_fraction(record.seq[start:end])
            
            # Get all products and format as a list
            products = feature.qualifiers.get('product', [])
            product_list = ";".join(products) if products else "None"
            
            features_list.append({
                'contig_id': record.id, 'start': start, 'end': end, 
                'length': length, 'num_cds': num_cds, 'gc_content': gc_content,
                'product': product_list, 
                'on_contig_edge': feature.qualifiers.get('contig_edge', ['False'])[0] == 'True', 
                'name': name
            })
    return features_list

def process_gbk_files(directory):
    """ Process all GenBank files in a directory and extract feature data """
    regions, cand_clusters, protoclusters = [], [], []
    for filename in os.listdir(directory):
        if filename.endswith('.gbk') and 'region' in filename:
            for record in SeqIO.parse(os.path.join(directory, filename), "genbank"):
                regions.extend(extract_and_calculate_features(record, 'region'))
                cand_clusters.extend(extract_and_calculate_features(record, 'candidate_cluster'))
                protoclusters.extend(extract_and_calculate_features(record, 'protocluster'))
    
    return pd.DataFrame(regions), pd.DataFrame(cand_clusters), pd.DataFrame(protoclusters)

# Directory with GenBank files
gbk_directory = "LU_1kb"

# Process files and save to CSV
regions_df, cand_clusters_df, protoclusters_df = process_gbk_files(gbk_directory)
regions_df.to_csv("regions_summary.csv", index=False)
cand_clusters_df.to_csv("cand_clusters_summary.csv", index=False)
protoclusters_df.to_csv("protoclusters_summary.csv", index=False)

print("Dataframes saved to 'regions_summary.csv', 'cand_clusters_summary.csv', and 'protoclusters_summary.csv'.")
