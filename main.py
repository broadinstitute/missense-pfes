
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
import argparse

# Figure plotting
plt.rcParams['figure.dpi'] = 500
plt.rcParams['font.family'] = 'Arial'

def plot_heatmap(df_score, data_folder, title = 'Protein Feature Enrichment Score', max = None):
    # plot the df_score heatmap for the same protein
    plt.figure(figsize=(23,6))
    pos_list = [df_score.RefAA[x] + str(x+1) for x in df_score.index]
    df_score_plot = df_score.set_index('RefAA')
    df_score_plot.index = pos_list
    df_score_plot.index.name = 'RefAA'
    if 'average' in df_score_plot.columns:
        df_score_plot = df_score_plot.drop(columns=['average'])
    
    if max == None:
        max = np.max([abs(df_score_plot.max().max()), abs(df_score_plot.min().min())])
    ax = sns.heatmap(df_score_plot.T.iloc[:, 0:], cmap='coolwarm', vmax=max, vmin=-max)
    
    # xlabel 
    ax.set_xlabel('Reference Amino Acid', fontsize=15)
    ax.set_ylabel('Mutated Amino Acid', fontsize=15)
    ax.set_title(title, fontsize=20)

    # # Set x-axis labels
    every = int(df_score_plot.shape[0]/15)
    tics = np.arange(df_score_plot.shape[0])[:-every:every].tolist() + [df_score_plot.shape[0]-1]
    ax.set_xticks(tics)
    ax.set_xticklabels(df_score_plot.iloc[tics].index)

    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)

    # set colorbar font 
    cbar = ax.collections[0].colorbar
    cbar.set_label('Score', fontsize=15)
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=15)
    
    plt.tight_layout()
    plt.savefig(f'{data_folder}/heatmap_{gene_name}_{uid}.png')
    # plt.show()
    
def plot_heatmap_1D(df_score, title = 'Protein Feature Enrichment Score', max = None):
    # plot the df_score heatmap for the same protein
    plt.figure(figsize=(24,0.5))
    pos_list = [df_score.RefAA[x] + str(x+1) for x in df_score.index]
    df_score_plot = df_score.set_index('RefAA')
    df_score_plot.index = pos_list
    df_score_plot.index.name = 'RefAA'
    if 'average' in df_score_plot.columns:
        df_score_plot = df_score_plot.drop(columns=['average'])
    
    if max == None:
        max = np.max([abs(df_score_plot.max().max()), abs(df_score_plot.min().min())])
    ax = sns.heatmap(df_score_plot.T.iloc[:, 0:], cmap='coolwarm', vmax=max, vmin=-max)
    
    # # xlabel 
    ax.set_xlabel('Reference Amino Acid', fontsize=15)
    ax.set_ylabel('  ', fontsize=15)
    ax.set_title(title, fontsize=20)

    # # # Set x-axis labels
    every = int(df_score_plot.shape[0]/15)
    tics = np.arange(df_score_plot.shape[0])[:-every:every].tolist() + [df_score_plot.shape[0]-1]
    ax.set_xticks(tics)
    ax.set_xticklabels(df_score_plot.iloc[tics].index)
    ax.set_yticks([])
    
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    
    # set colorbar font 
    cbar = ax.collections[0].colorbar
    cbar.set_label('Score', fontsize=15)
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=15)
    
    # remove colorbar
    # ax.collections[0].colorbar.remove()
    plt.show()
    
# Data file 
bucket = 'missense-score/PFES_NOV2024'
data_url = f'https://storage.googleapis.com/{bucket}'
unimeta = pd.read_csv(f'{data_url}/uniprot_metadata.tsv',delimiter='\t')
unimeta.index = unimeta.UniprotKB_Entry.tolist()

genemeta = pd.read_csv(f'{data_url}/gene_metadata.tsv',delimiter='\t')
genemeta.index = genemeta.HGNC_symbol.tolist()


# Initialize the parser
parser = argparse.ArgumentParser(description="Process either gene name or UniProt ID, along with mutation.")

# Create a mutually exclusive group for gene_name and uniprot
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-gene", help="The gene name to process", default='None')
group.add_argument("-uniprot", help="The UniProt ID to process", default='None')

# Add mutation argument
# Create a mutually exclusive group for gene_name and uniprot
group_m = parser.add_mutually_exclusive_group()
group_m.add_argument("-mutation", help="Single mutation (one letter code, e.g. M1A)")
group_m.add_argument("-mutations", help="The list of mutations to process")

# Parse the arguments - gene, uniprot, mutation
args = parser.parse_args()

gene_name = args.gene
uid = args.uniprot
mutation = args.mutation
mutations = args.mutations

if gene_name != 'None':
    uid = genemeta.loc[gene_name].UniprotKB_Entry
if uid != 'None':
    gene_name = unimeta.loc[uid].HGNC_symbol
    
print(f'Gene Name: {gene_name}')
print(f'UniProt ID: {uid}')

if mutation: 
    print(f'Mutation: {mutation}')
    
if mutations:
    print(f'Mutations: read from {mutations}')
    # Read the mutations from the file
    with open(mutations, 'r') as f:
        mutations = f.readlines()
    mutations = [m.strip() for m in mutations]
    print(mutations)
    
import requests
def url_exists(url):
    response = requests.head(url)
    return response.status_code == 200

exist_in_gcs = url_exists(f'{data_url}/PFES/scores_class_{uid}.txt')
if exist_in_gcs:
    # download file from GCS
    data_folder = f'data_{gene_name}_{uid}'
    if not os.path.exists(data_folder):
        os.system(f"mkdir {data_folder}")

    os.system(f"gsutil -q cp gs://{bucket}/PFES/scores_class*_{uid}.txt {data_folder}")

    if os.path.exists(f"{data_folder}/scores_class_{uid}.txt"):
        print(f"Reading scores for {gene_name}/{uid}")
        
        df_score = pd.read_csv(f"{data_folder}/scores_class_{uid}.txt",sep='\t')
        df_score_list = pd.read_csv(f"{data_folder}/scores_class_list_{uid}.txt",sep='\t')
        plot_heatmap(df_score, data_folder)
        total_max = df_score.drop(columns=['RefAA']).max().max()
    
else:
    print(f"No scores found for {gene_name}/{uid}")
    exit(0)
