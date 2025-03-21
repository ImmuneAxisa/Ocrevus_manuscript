##########################################################
################ Command line arguments ##################
##########################################################

# To do: add options for PHATE parameters

import argparse

# Initiate the parser
parser = argparse.ArgumentParser(description="this script runs PHATE on a PCA matrix. It returns PHATE embeddings.")
 
# Add long and short argument
parser.add_argument("--input", "-i", help="path to PCA matrix in csv format, first column must be cell IDs", required = True)
parser.add_argument("--output", "-o", help="path to PHATE embeddings output table", required = True)
parser.add_argument("--prefix", "-p", help="prefix for PHATE embeddings output table", default = "")

# Read arguments from the command line
args = parser.parse_args()

# Print input and output path
print("Set input file is %s" % args.input)

if args.prefix == "":
    out = args.output + "/" + 'data_phate.csv'
    print("Set output is %s" % out)
else:
    out = args.output + "/" + args.prefix + '_data_phate.csv'
    print("Set output is %s" % out)


##########################################################
####################### Set up ###########################
##########################################################

import pandas as pd
import numpy as np
import graphtools as gt
import phate
import magic
import scprep
import meld
import cmocean
import sklearn
import scipy
import seaborn as sns

#import gzip
#import pickle5 as pickle

import matplotlib.pyplot as plt
import diffxpy.api as de


import datetime
import os, tempfile
import types


########### PHATE options ###############


knn=10
decay=10
n_jobs=-1
n_components = 3
mds_solver = "smacof"


##########################################################
####################### PHATE ############################
##########################################################

# Read input PC matrix 
# NB: load_csv will store first column as index by default, where cell IDs should be
data_pca = scprep.io.load_csv(args.input)


# Run phate
phate_op = phate.PHATE(knn=knn, decay=decay, n_jobs=n_jobs, n_components = n_components, mds_solver = mds_solver)
data_phate = phate_op.fit_transform(data_pca)

# kmeans clustering
clusters = phate.cluster.kmeans(phate_op)

# export phate embeddings and clusters

#create header
header = ["PHATE" + str(dim) for dim in range(1,n_components+1)]
#convert to dataframe
data_phate_df = pd.DataFrame(data=data_phate, index=data_pca.index, columns=header)
#Add clusters
data_phate_df["PHATE_clusters"] = clusters
#save file
data_phate_df.index.name = "index"
data_phate_df.to_csv(out)



