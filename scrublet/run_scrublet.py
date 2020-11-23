import scrublet as scr
import scipy.io
import numpy as np
import pandas as pd
import os


# Run this script from within the seurat output folder /outs/
# will create a folder /scrublet/ that contains the predictions for each barcode.



print(os.getcwd())

input_dir = 'filtered_feature_bc_matrix/'

print(os.path.exists(input_dir))

counts_matrix = scipy.io.mmread(input_dir + 'matrix.mtx.gz').T.tocsc()

print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))

scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.15)

doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
                                                          

scores = pd.DataFrame({"doublet_scores":doublet_scores, "predicted_doublets":predicted_doublets})

output_dir = "scrublet/"

if(os.path.exists("scrublet/") == False):
  os.mkdir(output_dir)

scores.to_csv(output_dir + "predicted_doublets.csv")

print(len(doublet_scores))

