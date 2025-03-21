# This script MELD and VFC
# Returns MELD densities and likelihoods, and VFC clusters
# Modified from: https://nbviewer.jupyter.org/github/KrishnaswamyLab/MELD/blob/master/notebooks/Wagner2018_Chordin_Cas9_Mutagenesis.ipynb

##########################################################
################ Command line arguments ##################
##########################################################


import argparse

# Initiate the parser
parser = argparse.ArgumentParser(description="this script runs MELD and  VFC on a PCA matrix. It returns MELD densities and likelihoods, and VFC clusters.")

# Add long and short argument
parser.add_argument("--input", "-i", help="path to PCA matrix in csv format, first column must be cell IDs", required = True)
parser.add_argument("--meta", "-m", help="path to metadata table in csv format, first column must be cell IDs", required = True)
parser.add_argument("--labels", "-l", help="name of metadata column to use as label for MELD (categorical)", required = True)
parser.add_argument("--cluster", "-c", help="name of cluster metadata column to use for subsetting before VFC is run (categorical)")
parser.add_argument("--benchmarking", "-b", help="path to benchmarking output table in csv format")
parser.add_argument("--output", "-o", help="path to output table", required = True)
parser.add_argument("--prefix", "-p", help="prefix for output table", default = "")

# Read arguments from the command line
args = parser.parse_args()

# Print input and output path
print("Set input file is %s" % args.input)
print("Set metadata file is %s" % args.meta)
print("Sample label to test is %s" % args.labels)

if args.prefix == "":
    out = args.output + "/"
    print("Set output stem is %s" % out)
else:
    out = args.output + "/" + args.prefix + '_'
    print("Set output stem is %s" % out)


##########################################################
######################## Set up ##########################
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

import os, tempfile
import types
from datetime import datetime
#import gzip
#import pickle5 as pickle

# setting defaults for matplotlib font sizes
import matplotlib.pyplot as plt
plt.rc('font', size=14)

# making sure plots & clusters are reproducible
np.random.seed(42)

#%load_ext autoreload
#%autoreload 2

import diffxpy.api as de



##########################################################
######################### Options ########################
##########################################################


################ inputs #################

data_pca = scprep.io.load_csv(args.input)

metadata = scprep.io.load_csv(args.meta)

sample_labels = metadata[args.labels]


########### MELD options ################



########### VFC options ################

if args.cluster:
    # if cluster input reassign cluster variable to clusterID
    print("Setting cluster column to be used by VFC")
    metadata['clusterID'] = metadata[args.cluster]
    
else:
    print("No cluster provided. VFC will be run on the whole dataset")
    # Else create an articial cluster that is constant across all the dataset
    metadata['clusterID'] = 0
    


##########################################################
########################## MELD ##########################
##########################################################


if args.benchmarking:
    # get results:
    results = scprep.io.load_csv(args.benchmarking)
    
    # We want to take the average of each set of random seeds for each combination of beta and knn values
    results_wide = results.groupby(['beta', 'knn']).mean().sort_values(by='MSE').reset_index()
    
    # get best result
    top_result = results_wide.sort_values('MSE').iloc[0]
    
else:
    print("No benchmarking. Using default beta and knn values")
    top_result = pd.Series([60,5], index = ['beta','knn'])


print(top_result)
############## run MELD ###############


G = gt.Graph(data_pca, knn=int(top_result['knn']), use_pygsp=True)


meld_op = meld.MELD(beta=top_result['beta'])

sample_densities = meld_op.fit_transform(G, sample_labels=sample_labels)


# add cell ids to index
sample_densities.index = data_pca.index

sample_densities.index.name = "index"

# Function from MELD manual, wrapper for sklearn.preprocessing.normalize
sample_likelihoods = meld.utils.normalize_densities(sample_densities)


# add cell ids to index
sample_likelihoods.index = data_pca.index

sample_likelihoods.index.name = "index"


# Save outputs

sample_densities.to_csv(out + 'MELD_densities.csv')

sample_likelihoods.to_csv(out + 'MELD_LLH.csv')


##########################################################
###################### VFC clustering ####################
##########################################################

print("finished running MELD, now running VFC")

np.random.seed(0)

## -------------------------------------------------------------------------------------------------------
# partial eigen decomposition for compute_fourier_basis.
# See: https://github.com/epfl-lts2/pygsp/issues/27
# Implementation for partial decomposition written but not released
# Pasting the updated function here and replacing the instance method with this new version

def compute_fourier_basis_new(self, n_eigenvectors=100, recompute=False):
        r"""Compute the (partial) Fourier basis of the graph (cached).
        The result is cached and accessible by the :attr:`U`, :attr:`e`,
        :attr:`lmax`, and :attr:`mu` properties.
        Parameters
        ----------
        n_eigenvectors : int or `None`
            Number of eigenvectors to compute. If `None`, all eigenvectors
            are computed. (default: None)
        recompute: bool
            Force to recompute the Fourier basis if already existing.
        Notes
        -----
        'G.compute_fourier_basis()' computes a full eigendecomposition of
        the graph Laplacian :math:`L` such that:
        .. math:: L = U \Lambda U^*,
        or a partial eigendecomposition of the graph Laplacian :math:`L`
        such that:
        .. math:: L \approx U \Lambda U^*,
        where :math:`\Lambda` is a diagonal matrix of eigenvalues and the
        columns of :math:`U` are the eigenvectors.
        *G.e* is a vector of length `n_eigenvectors` :math:`\le` *G.N*
        containing the Laplacian eigenvalues. The largest eigenvalue is stored
        in *G.lmax*. The eigenvectors are stored as column vectors of *G.U* in
        the same order that the eigenvalues. Finally, the coherence of the
        Fourier basis is found in *G.mu*.
        References
        ----------
        See :cite:`chung1997spectral`.
        Examples
        --------
        >>> G = graphs.Torus()
        >>> G.compute_fourier_basis(n_eigenvectors=64)
        >>> G.U.shape
        (256, 64)
        >>> G.e.shape
        (64,)
        >>> G.compute_fourier_basis()
        >>> G.U.shape
        (256, 256)
        >>> G.e.shape
        (256,)
        >>> G.lmax == G.e[-1]
        True
        >>> G.mu < 1
        True
        """
        if n_eigenvectors is None:
            n_eigenvectors = self.N

        if (hasattr(self, '_e') and hasattr(self, '_U') and not recompute
                and n_eigenvectors <= len(self.e)):
            return

        assert self.L.shape == (self.N, self.N)
        if self.N**2 * n_eigenvectors > 3000**3:
            self.logger.warning(
                'Computing the {0} eigendecomposition of a large matrix ({1} x'
                ' {1}) is expensive. Consider decreasing n_eigenvectors '
                'or, if using the Fourier basis to filter, using a '
                'polynomial filter instead.'.format(
                    'full' if n_eigenvectors == self.N else 'partial',
                    self.N))

        # TODO: handle non-symmetric Laplacians. Test lap_type?
        if n_eigenvectors == self.N:
            self._e, self._U = np.linalg.eigh(self.L.toarray())
        else:
            # fast partial eigendecomposition of hermitian matrices
            self._e, self._U = scipy.sparse.linalg.eigsh(self.L,
                                                   n_eigenvectors,
                                                   which='SM')
        # Columns are eigenvectors. Sorted in ascending eigenvalue order.

        # Smallest eigenvalue should be zero: correct numerical errors.
        # Eigensolver might sometimes return small negative values, which
        # filter's implementations may not anticipate. Better for plotting
        # too.
        assert -1e-12 < self._e[0] < 1e-12
        self._e[0] = 0

        if self.lap_type == 'normalized':
            # Spectrum bounded by [0, 2].
            assert self._e[-1] <= 2

        assert np.max(self._e) == self._e[-1]
        self._lmax = self._e[-1]
        self._mu = np.max(np.abs(self._U))


## -------------------------------------------------------------------------------------------------------
# Initialize list to store VFC results on each cluster (if no cluster it'll run once on the whole dataset)
vfc_results = list()

  
for cluster in np.unique(metadata['clusterID']):
    
    now = datetime.now()
    print("Start processing cluster " + str(cluster) + " at " + now.strftime("%d/%m/%Y %H:%M:%S"))
    
    pca = data_pca.loc[metadata['clusterID'] == cluster]
    
    if pca.shape[0] > 200 and pca.shape[0] < 40000:
    
        ## -----------------------------------------------------------------------------
        # If cluster to subset on have been provided rebuild graph on the subset. 
        # If running VFC on full set we already had the graph when running MELD
        if args.cluster:
            print("building cluster specific graph")
            G = gt.Graph(pca, knn=int(top_result['knn']), use_pygsp=True)
    
        ## -----------------------------------------------------------------------------
        # Modify default fourier method to partil method and compute fourier basis
        G.compute_fourier_basis = types.MethodType(compute_fourier_basis_new, G)
        
        G.compute_fourier_basis()
        
        now = datetime.now()
        print("finished compute_fourier_basis at " + now.strftime("%d/%m/%Y %H:%M:%S"))
        
        ## -----------------------------------------------------------------------------
        # Set VFC operator
        vfc = meld.VertexFrequencyCluster(n_clusters = 5)
        
        
        ## -----------------------------------------------------------------------------
        # Binarise the sample labels using the first category in the likelihood output from MELD
        indicator = sample_labels.loc[metadata['clusterID'] == cluster] == sample_likelihoods.columns[0]
        
        # subset likelihoods
        likelihoods_subset = sample_likelihoods.loc[metadata['clusterID'] == cluster]
        
        
        ## -----------------------------------------------------------------------------
        ################################## Run VFC #####################################
        vfc.fit_transform(G, indicator, likelihoods_subset.iloc[:, 1])
        
        now = datetime.now()
        print("finished fit_transform at " + now.strftime("%d/%m/%Y %H:%M:%S"))
        
        
        ## -----------------------------------------------------------------------------
        # test a few Kmeans values for clustering
        clusters_by_n = {}
        for n in [2,3,4,5,6,7,8,9,10]:
            clusters_by_n["KMeans_" + str(n)] = vfc.predict(n)
            now = datetime.now()
            print("finished KMeans_" + str(n) + " at " + now.strftime("%d/%m/%Y %H:%M:%S"))
        
        
        ## -----------------------------------------------------------------------------
        # Save outputs 
        
        # Create data frame of results with the different K means
        clusters_by_n_df = pd.DataFrame(clusters_by_n, index = data_pca.loc[metadata['clusterID'] == cluster].index)
        
        # Add index name
        clusters_by_n_df.index.name = "index"
        
        # Add cluster iteration 
        clusters_by_n_df['clusterID'] = cluster
        
        ## -----------------------------------------------------------------------------
        ################################# PHATE ########################################
        
        
        ## -----------------------------------------------------------------------------
        # PHATE options
        
        knn=10
        decay=10
        n_jobs=-1
        n_components = 3
        mds_solver = "smacof"
        
        
        ## -----------------------------------------------------------------------------
        # Run phate
        phate_op = phate.PHATE(knn=knn, decay=decay, n_jobs=n_jobs, n_components = n_components, mds_solver = mds_solver)
        data_phate = phate_op.fit_transform(pca)
        
        ## -----------------------------------------------------------------------------
        ############## export phate embeddings and clusters ############################
        
        #create header
        header = ["VCF_PHATE" + str(dim) for dim in range(1,n_components+1)]
        #convert to dataframe
        data_phate_df = pd.DataFrame(data=data_phate, index=pca.index, columns=header)
        
        data_phate_df.index.name = "index"
        
        # Concat VCF cluster calls and PHATE embeddings run on the cluster
        merged_df=pd.concat([clusters_by_n_df, data_phate_df], axis=1)
        
        # Append data frame for the current iteration to the result list
        vfc_results.append(merged_df)
        
        






    
# Concatenate list of dataframes into one dataframe and save to file
clusters_by_n_df = pd.concat(vfc_results)
    
clusters_by_n_df.to_csv(out + 'VCF_clusters.csv')
