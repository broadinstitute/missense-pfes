import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
import argparse

# Figure plotting
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'Arial'

def plot_heatmap(df_score, data_folder, key, title = 'Protein Feature Enrichment Score', max = None):
    # plot the df_score heatmap for the same protein
    plt.figure(figsize=(15,6))
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
    plt.savefig(f'{data_folder}/heatmap_{key}_{gene_name}_{uid}.png', bbox_inches='tight')
    print ('heatmap generated:', f'{data_folder}/heatmap_{key}_{gene_name}_{uid}.png')
    
def plot_heatmap_1D(df_score, data_folder, key, title = 'Protein Feature Enrichment Score', max = None):
    # plot the df_score heatmap for the same protein
    plt.figure(figsize=(20,0.5))
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

    # plt.subplots_adjust(top=0.6, bottom=0.4)
    # plt.tight_layout()
    plt.savefig(f'{data_folder}/heatmap_{key}_{gene_name}_{uid}.png', bbox_inches='tight')
    print ('heatmap generated:', f'{data_folder}/heatmap_{key}_{gene_name}_{uid}.png')
    
# Data file 
bucket = 'missense-score/PFES_NOV2024'
data_url = f'https://storage.googleapis.com/{bucket}'

unimeta = pd.read_csv(f'../files/uniprot_metadata_2024_04.tsv',delimiter='\t')
unimeta.index = unimeta.UniprotKB_Entry.tolist()

genemeta = pd.read_csv(f'../files/gene_metadata_2024_04.tsv',delimiter='\t')
genemeta.index = genemeta.HGNC_symbol.tolist()

df_plp_like = pd.read_csv("../files/plp_likelihood.tsv",delimiter='\t')

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


import requests
def url_exists(url):
    response = requests.head(url)
    return response.status_code == 200

exist_in_gcs = url_exists(f'{data_url}/PFES/scores_class_{uid}.txt')
if exist_in_gcs:
    # download file from GCS
    data_folder = f'../data_{gene_name}_{uid}'
    if not os.path.exists(data_folder):
        os.system(f"mkdir {data_folder}")

    os.system(f"gsutil -q cp gs://{bucket}/PFES/scores_class*_{uid}.txt {data_folder}")
    os.system(f"gsutil -q cp gs://{bucket}/PFES/scores_decomp_class*_{uid}.txt {data_folder}")

    if os.path.exists(f"{data_folder}/scores_class_{uid}.txt"):
        print(f"Reading scores for {gene_name}/{uid}")
        
        df_score = pd.read_csv(f"{data_folder}/scores_class_{uid}.txt",sep='\t')
        df_score_list = pd.read_csv(f"{data_folder}/scores_class_list_{uid}.txt",sep='\t')
        
        # Generate Heatmap 
        plot_heatmap(df_score, data_folder, 'all')
        total_max = df_score.drop(columns=['RefAA']).max().max()
        
        # Generate Heatmap - decomposed score
        mops = ['Physicochemical', 'Function', 'Domain','Structure', 'Modification', 'PPI']
        df_score_decomp = pd.read_csv(f"{data_folder}/scores_decomp_class_{uid}.txt",sep='\t')
        
        mop = 'Physicochemical'
        df_score_plot = pd.DataFrame()
        df_score_plot['RefAA'] = df_score_decomp['RefAA'] 
        aas = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K','M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
        for alt in aas: 
            df_score_plot[alt] = df_score_decomp[alt].apply(lambda x: eval(x)[mop])
        plot_heatmap(df_score_plot, data_folder, mop, title = f'Protein Feature Enrichment Score: {mop}', max=total_max/2.0)

        df_score_plot = pd.DataFrame()
        df_score_plot['RefAA'] = df_score_decomp['RefAA']
        for mop in ['Function', 'Domain','Structure', 'Modification', 'PPI']:
            df_score_plot[f'{mop}'] = df_score_decomp['A'].apply(lambda x: eval(x)[mop]).tolist()
            plot_heatmap_1D(df_score_plot[['RefAA',f'{mop}']], data_folder, mop, title = f'Protein Feature Enrichment Score: {mop}', max=total_max/2.0)

        # Pathogenicity likelihood
        # smooth line (2nd polynomial) fit to the data points
        z = np.polyfit(df_plp_like.Center, df_plp_like['PLP dominance'], 2)
        p_plp = np.poly1d(z)

        from scipy.optimize import curve_fit
        # Define the sigmoid function
        def sigmoid(x, x0, k):
            return 100 / (1 + np.exp(-k * (x - x0)))
        # Extract data from your DataFrame
        x_data = df_plp_like.Center
        y_data = df_plp_like['PLP dominance']
        # Fit the sigmoid curve to the data
        popt, pcov = curve_fit(sigmoid, x_data, y_data, method='dogbox', bounds=([-np.inf, 0], [np.inf, 1]))
        # popt contains the optimal values for x0 and k
        x0, k = popt
        # Create a callable sigmoid function similar to `p_plp`
        p_plp = lambda x: sigmoid(x, x0, k)

        print ('')
        if mutations:
            print(f'Mutations: read from {mutations}')
            # Read the mutations from the file
            with open(mutations, 'r') as f:
                mutations = f.readlines()
            mutations = [m.strip() for m in mutations]
            print(mutations)

        if mutation: 
            print(f'Mutation: {mutation}')
            mutations = [mutation]
        
        # Mutations
        df_mut_score = pd.DataFrame(columns=['Mutation'])
        df_mut_score['Mutation'] = mutations
        
        # Check if valid
        for mutation in mutations:
            
            refaa = mutation[0]
            resid = int(mutation[1:-1])
            mutaa = mutation[-1]

            pos_list = [df_score.RefAA[x] + str(x+1) for x in df_score.index]

            if mutation[0:-1] in pos_list:
                try: 
                    this_mut_score = df_score.iloc[resid-1][mutaa]
                    this_mut_list = df_score_list.iloc[resid-1][mutaa].split(",")
                    this_mut_decomp = df_score_decomp.iloc[resid-1][mutaa]
                    
                    # print(f'Mutation {mutation} is a valid mutation')
                    # print (f'   Total PFE Score: {this_mut_score}')
                    
                    closest_row = df_plp_like.iloc[(df_plp_like['Center'] - this_mut_score).abs().argmin()]#['Center']
                    # print (f"     :correspoding to {p_plp(this_mut_score):.1f}% likelihood of pathogenicity \n")
                        
                    # print (f'   Decomposed PFE Score\n     :{eval(this_mut_decomp)} \n')
                    # print (f'   List of features;Odd ratio\n     : {this_mut_list} \n')

                    df_mut_score.loc[df_mut_score.Mutation == mutation, 'Total PFE Score'] = this_mut_score
                    df_mut_score.loc[df_mut_score.Mutation == mutation, 'Pathogenicity Likelihood'] = p_plp(this_mut_score)
                    for key in eval(this_mut_decomp).keys():
                        df_mut_score.loc[df_mut_score.Mutation == mutation, f'{key} PFE Score'] = eval(this_mut_decomp)[key]
                        
                    df_mut_score.loc[df_mut_score.Mutation == mutation, 'List of features;Odd ratio'] = "|".join(this_mut_list)
                    
                except:
                    print(f'Error: Mutation {mutation} not a valid mutation')

                df_mut_score.to_csv(f"{data_folder}/mutations_PFES_{gene_name}_{uid}.csv", index=False)
                print (f'CSV file for mutations generated: {data_folder}/mutations_PFES_{gene_name}_{uid}.csv')
                
            else:
                print(f'Error: Mutation {mutation} not found in the score file, reference AA is {df_score.iloc[resid-1].RefAA}{resid}')


else:
    print(f"No scores found for {gene_name}/{uid}")
    exit(0)
    
    

