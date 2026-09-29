import sys
from multiprocessing import Pool, TimeoutError, Array
from contextlib import closing
import csv
import gzip
import os
import scipy.io
from scipy.stats import boxcox
import numpy as np
import numpy.linalg as LA
from sklearn.decomposition import PCA
import pandas as pd
from datetime import datetime
import math
from scipy.stats import boxcox
import random
# Set the seed for reproducibility
random.seed(42)
np.random.seed(42)

print("Start time:")
print(datetime.now())

#define arguments values
workspace_path = sys.argv[1] # 
filename_fitness = sys.argv[2] # "lst_Scer_strains_Fitness.csv" # 
filename_mtx_pres_abs = sys.argv[3] # "mtx_Scer_strains_Gene_Presence_Absence.csv" # 
filename_mtx_expr = sys.argv[4] # "mtx_Scer_strains_Expression.csv" # 
filename_mtx_phylo_dist = sys.argv[5] #16 species recovered in YPD out of 21 # "mtx_Scer_strains_Phylo_dist_Kernel.csv" # 
dataset_name = sys.argv[6] # 
nb_subsamples = int(sys.argv[7]) # 
nb_cpus = int(sys.argv[8]) # 
sample_size = int(sys.argv[9]) # 

# Original EM REML n x n formulation
def reml_em(Kernel,X,y,sig_estimate=None,verbose=True,n_iter=100):
    # Calculating new y
    Q,R = LA.qr(X)
    M = lambda O : O - Q @ (Q.T @ O)
    resid = M(y)
    y_new = (resid - resid.mean()) / resid.std(ddof=1)

    n = Kernel[0].shape[0]
    p = len(Kernel)
    if sig_estimate is None:
        sig_estimate = np.zeros(p+1, dtype=np.float64)
        sig_estimate[:] = np.var(y_new,ddof=1) / (p+1)
        #print(sig_estimate)
    iteration_estimate = []
    # Convergence
    exit_code = 1

    for i in range(n_iter):
        prev_sig = sig_estimate.copy()

        V = sig_estimate[p] * np.eye(n, dtype=np.float64) # for I
        for j in range(p):
            V += sig_estimate[j] * Kernel[j]
        
        V_inv = LA.inv(V)
        R = V_inv - V_inv.dot(X).dot(LA.inv(X.T.dot(V_inv).dot(X))).dot(X.T).dot(V_inv)

        trace = np.zeros(p+1, dtype=np.float64)
        quad = np.zeros(p+1, dtype=np.float64)
        for j in range(p):
            quad[j] = y_new.T.dot(R).dot(Kernel[j]).dot(R).dot(y_new)
            trace[j] = np.trace( R.dot(Kernel[j]) )

        quad[p] = y_new.T.dot(R).dot(R).dot(y_new)
        trace[p] = np.trace( R )

        sig_estimate = prev_sig - ( (prev_sig**2) * (trace - quad) ) / n
        iteration_estimate.append(sig_estimate)

        print("Quad:", quad)
        print("Trace (RK):", trace)

        if verbose:
            print_str = "\t EM REML vanilla round " + str(i) + ": "
            for j in range(p + 1):
                print_str += "  {:.4f}".format(sig_estimate[j])
            print(print_str, flush=True)
        
        diff = np.max( np.abs(sig_estimate - prev_sig) )
        if diff < 1e-4 :
            exit_code = 0
            if verbose:
                print("\t Estimates converged, exitting - diff: %.6f" % (diff) )
            break

    iteration_estimate = np.array(iteration_estimate)
    return sig_estimate,iteration_estimate[0,:],iteration_estimate,exit_code

v_comp_a = []
v_comp_b = []
v_comp_c = []
v_comp_d = []
v_comp_e = []
v_comp_f = []
v_comp_g = []
v_comp_residuals = []
v_comp_Gene_content = []
v_comp_Expression = []
v_comp_Phylo_dist_ypd = []

#import data
Fitness_ypd = pd.read_csv("{0}/{1}".format(workspace_path,filename_fitness),sep=",",header=None).to_numpy().squeeze()
Fitness_ypd = (Fitness_ypd - np.mean(Fitness_ypd))/np.std(Fitness_ypd)
org_Fitness_ypd = Fitness_ypd

mtx_phylo_dist_ypd = pd.read_csv("{0}/{1}".format(workspace_path,filename_mtx_phylo_dist),sep=",",header=None).to_numpy()
#mtx_phylo_dist_ypd = (pd.DataFrame((mtx_phylo_dist_ypd+1.0)[:,np.where(np.sum(mtx_phylo_dist_ypd,axis=0)!=0)[0]] ).apply(lambda x: boxcox(x)[0])).to_numpy() # (mtx_phylo_dist_ypd - np.mean(mtx_phylo_dist_ypd,axis=0))/np.std(mtx_phylo_dist_ypd,axis=0) # 
org_mtx_phylo_dist_ypd = mtx_phylo_dist_ypd
#print(np.shape(mtx_phylo_dist_ypd))

E = pd.read_csv("{0}/{1}".format(workspace_path,filename_mtx_expr),sep=",",header=None).to_numpy()
#print(np.shape(E))
#arr_ind_genes_never_Expressed = np.where(np.sum(E,axis=0)==0)[0]
#print(arr_ind_genes_never_Expressed)
#arr_ind_genes_Expressed = np.where(np.sum(E,axis=0)!=0)[0]
#E = E[:,arr_ind_genes_Expressed]
#E = (pd.DataFrame((E+1.0)[:,np.where(np.sum(E,axis=0)!=0)[0]] ).apply(lambda x: boxcox(x)[0])).to_numpy() # (E - E.mean(axis=0))/np.std(E,axis=0) # 
#print(np.shape(E))
#print(E[0:3,0:3])
#print(np.std(E,axis=0))
#print(np.any(np.isnan(E)))
org_E = E

Gene_content = pd.read_csv("{0}/{1}".format(workspace_path,filename_mtx_pres_abs),sep=",",header=None).to_numpy()
#print(np.shape(Gene_content))
#arr_ind_genes_never_present = np.where(np.sum(Gene_content,axis=0)==0)[0]
#print(arr_ind_genes_never_present)
#print(len(arr_ind_genes_never_present))
#arr_ind_genes_present = np.where(np.sum(Gene_content,axis=0)!=0)[0]
#Gene_content = Gene_content[:,arr_ind_genes_present]
#(Gene_content[:,np.where(np.nansum(Gene_content,axis=0)!=0)[0]]- np.nanmean(Gene_content[:,np.where(np.nansum(Gene_content,axis=0)!=0)[0]],axis=0))
#print(np.shape(Gene_content))
#print(Gene_content[0:3,0:3])
#print(np.any(np.isnan(Gene_content)))
org_Gene_content = Gene_content

for id_subsample in np.arange(nb_subsamples):
    lst_ids_current_subsample = random.sample(np.arange(np.shape(org_mtx_phylo_dist_ypd)[0]).tolist(), sample_size)
    Fitness_ypd = org_Fitness_ypd[lst_ids_current_subsample]
    mtx_phylo_dist_ypd = org_mtx_phylo_dist_ypd[lst_ids_current_subsample,:]
    E= org_E[lst_ids_current_subsample,:]
    Gene_content = org_Gene_content[lst_ids_current_subsample,:]
    
    #USE the phylogenetic VCV matrix as phylogenetic feature instead of the phylogenetic distances covariance
    mtx_cov_phylo_dist = np.cov(mtx_phylo_dist_ypd) #pd.read_csv("{0}/Phylo_vcv_kernel_21_spcs_YPD.csv".format(workspace_path),sep=",",header=None).to_numpy() #

    N = np.shape(mtx_phylo_dist_ypd)[0]

    #Define kernels
    K_Gene_content = pd.DataFrame(Gene_content.T).cov().to_numpy() / (np.shape(Gene_content)[1])#(Gene_content @ Gene_content.T) / (np.shape(Gene_content)[1])
    K_E = np.cov(E) / (np.shape(E)[1])#(E @ E.T) / (np.shape(E)[1])
    K_phylo_dist_ypd = mtx_cov_phylo_dist / (np.shape(mtx_cov_phylo_dist)[1])#(mtx_cov_phylo_dist @ mtx_cov_phylo_dist.T) / (np.shape(mtx_cov_phylo_dist)[1]) #/719591#

    #single
    Kernels_Gene_content = [K_Gene_content]
    Kernels_Expression = [K_E]
    Kernels_phylo_dist_ypd = [K_phylo_dist_ypd]

    #pairs
    Kernels_Gene_content_and_Expression = [K_Gene_content,K_E]
    Kernels_Gene_content_and_phylo_dist_ypd = [K_Gene_content,K_phylo_dist_ypd]
    Kernels_Expression_and_phylo_dist_ypd = [K_E,K_phylo_dist_ypd]

    #full model (ALL EXPLANATORY VARIABLES)
    Kernels_full_model = [K_Gene_content,K_E,K_phylo_dist_ypd]

    #Representation of no fixed effect
    X = np.ones(N).reshape((N,1))
    
    #Test REML VarPart on each component INDIVIDUALLY
    print("##########################Test REML VarPart on each component INDIVIDUALLY################################")
    sigma_sqs_reml_Pheno_vs_Gene_content = reml_em(Kernel=Kernels_Gene_content, X=X, y=Fitness_ypd, verbose=True, n_iter=500)[0]
    sigma_sqs_reml_Pheno_vs_Expression = reml_em(Kernel=Kernels_Expression, X=X, y=Fitness_ypd, verbose=True, n_iter=500)[0]
    sigma_sqs_reml_Pheno_vs_Phylo_dist_ypd = reml_em(Kernel=Kernels_phylo_dist_ypd, X=X, y=Fitness_ypd, verbose=True, n_iter=500)[0]

    #Test REML VarPart on each component PAIRS
    print("##########################Test REML VarPart on each component PAIRS################################")
    sigma_sqs_reml_Pheno_vs_Gene_content_and_Expression = reml_em(Kernel=[K_Gene_content, K_E], X=X, y=Fitness_ypd, verbose=True, n_iter=500)[0]
    sigma_sqs_reml_Pheno_vs_Gene_content_and_Phylo_dist_ypd = reml_em(Kernel=[K_Gene_content,K_phylo_dist_ypd], X=X, y=Fitness_ypd, verbose=True, n_iter=500)[0]
    sigma_sqs_reml_Pheno_vs_Expression_and_Phylo_dist_ypd = reml_em(Kernel=[K_E,K_phylo_dist_ypd], X=X, y=Fitness_ypd, verbose=True, n_iter=500)[0]

    #Test REML VarPart (FULL)
    print("##########################Test REML VarPart (FULL)################################")
    sigma_sqs_reml_Full = reml_em(Kernel=[K_Gene_content, K_E,K_phylo_dist_ypd], X=X, y=Fitness_ypd, verbose=True, n_iter=500)[0]
    
    R2_Pheno_vs_Gene_content = (sigma_sqs_reml_Pheno_vs_Gene_content[0])/np.sum(sigma_sqs_reml_Pheno_vs_Gene_content)
    R2_Pheno_vs_Expression = (sigma_sqs_reml_Pheno_vs_Expression[0])/np.sum(sigma_sqs_reml_Pheno_vs_Expression)
    R2_Pheno_vs_Phylo_dist_ypd = (sigma_sqs_reml_Pheno_vs_Phylo_dist_ypd[0])/np.sum(sigma_sqs_reml_Pheno_vs_Phylo_dist_ypd)

    R2_Pheno_vs_Gene_content_and_Expression = np.sum(sigma_sqs_reml_Pheno_vs_Gene_content_and_Expression[0:(len(sigma_sqs_reml_Pheno_vs_Gene_content_and_Expression)-1)])/np.sum(sigma_sqs_reml_Pheno_vs_Gene_content_and_Expression)
    R2_Pheno_vs_Gene_content_and_Phylo_dist_ypd = np.sum(sigma_sqs_reml_Pheno_vs_Gene_content_and_Phylo_dist_ypd[0:(len(sigma_sqs_reml_Pheno_vs_Gene_content_and_Phylo_dist_ypd)-1)])/np.sum(sigma_sqs_reml_Pheno_vs_Gene_content_and_Phylo_dist_ypd)
    R2_Pheno_vs_Expression_and_Phylo_dist_ypd = np.sum(sigma_sqs_reml_Pheno_vs_Expression_and_Phylo_dist_ypd[0:(len(sigma_sqs_reml_Pheno_vs_Expression_and_Phylo_dist_ypd)-1)])/np.sum(sigma_sqs_reml_Pheno_vs_Expression_and_Phylo_dist_ypd)

    R2_Full = np.sum(sigma_sqs_reml_Full[0:(len(sigma_sqs_reml_Full)-1)])/np.sum(sigma_sqs_reml_Full)
    
    #individual prediction contributions
    print("##########################individual prediction contributions to phenotype variation################################")
    print(R2_Pheno_vs_Gene_content)
    print(R2_Pheno_vs_Expression)
    print(R2_Pheno_vs_Phylo_dist_ypd)

    #pairs contributions
    print("##########################pairs contributions to phenotype variation################################")
    print(R2_Pheno_vs_Gene_content_and_Expression)
    print(R2_Pheno_vs_Gene_content_and_Phylo_dist_ypd)
    print(R2_Pheno_vs_Expression_and_Phylo_dist_ypd)

    #full model contribution to phenotype variation
    print("##########################full model contribution to phenotype variation################################")
    print(R2_Full)
    
    #Varpart solution
        #full model and full predictors contributions to response variable variance
    comp_full = R2_Full
    comp_d_plus_e = R2_Pheno_vs_Gene_content + R2_Pheno_vs_Phylo_dist_ypd - R2_Pheno_vs_Gene_content_and_Phylo_dist_ypd
    comp_b_plus_e = R2_Pheno_vs_Gene_content + R2_Pheno_vs_Expression - R2_Pheno_vs_Gene_content_and_Expression
    comp_e_plus_f = R2_Pheno_vs_Expression + R2_Pheno_vs_Phylo_dist_ypd - R2_Pheno_vs_Expression_and_Phylo_dist_ypd
        #single var contribution only (independently of all others) = Full - pairwise
    comp_a = R2_Full - R2_Pheno_vs_Expression_and_Phylo_dist_ypd
            #adjust negative variances to 0 (due to multicollinearity or bas normalization)
    if (comp_a<0):
        comp_a = 0
    comp_c = R2_Full - R2_Pheno_vs_Gene_content_and_Phylo_dist_ypd
            #adjust negative variances to 0 (due to multicollinearity or bas normalization)
    if (comp_c<0):
        comp_c = 0
    comp_g = R2_Full - R2_Pheno_vs_Gene_content_and_Expression
            #adjust negative variances to 0 (due to multicollinearity or bas normalization)
    if (comp_g<0):
        comp_g = 0
        #all component intersection
    comp_e = comp_a + comp_b_plus_e + comp_d_plus_e - R2_Pheno_vs_Gene_content
            #adjust negative variances to 0 (due to multicollinearity or bas normalization)
    if (comp_e<0):
        comp_e = 0
        #pairs intersection
    comp_b = comp_b_plus_e - comp_e
            #adjust negative variances to 0 (due to multicollinearity or bas normalization)
    if (comp_b<0):
        comp_b = 0
    comp_d = comp_d_plus_e - comp_e
            #adjust negative variances to 0 (due to multicollinearity or bas normalization)
    if (comp_d<0):
        comp_d = 0
    comp_f = comp_e_plus_f - comp_e
            #adjust negative variances to 0 (due to multicollinearity or bas normalization)
    if (comp_f<0):
        comp_f = 0
        #residuals
    residual_comp = 1 - (comp_a+comp_b+comp_c+comp_d+comp_e+comp_f+comp_g)

    #save components explained variance
    v_comp_a.append(comp_a)
    v_comp_b.append(comp_b)
    v_comp_c.append(comp_c)
    v_comp_d.append(comp_d)
    v_comp_e.append(comp_e)
    v_comp_f.append(comp_f)
    v_comp_g.append(comp_g)
    v_comp_residuals.append(residual_comp)
    v_comp_Gene_content.append(R2_Pheno_vs_Gene_content)
    v_comp_Expression.append(R2_Pheno_vs_Expression)
    v_comp_Phylo_dist_ypd.append(R2_Pheno_vs_Phylo_dist_ypd)
    
    #timestamp for the end of the current iteration
    print("End time of the current iteration:")
    print(datetime.now())

#save Phenotype variance explained in table
data = {'id_subsample' : np.arange(nb_subsamples).tolist(),
        'comp_a': v_comp_a,
        'comp_b' : v_comp_b,
        'comp_c': v_comp_c,
        'comp_d': v_comp_d,
        'comp_e': v_comp_e,
        'comp_f': v_comp_f,
        'comp_g': v_comp_g,
        'residual_comp': v_comp_residuals,
        'comp_Gene_content': v_comp_Gene_content,
        'comp_Expression': v_comp_Expression,
        'comp_Phylo_dist_ypd': v_comp_Phylo_dist_ypd}

df_out = pd.DataFrame(data)
print(df_out)
df_out.to_csv("{0}/Table_phenotype_variance_explained_in_the_dataset_{1}.tsv".format(workspace_path,dataset_name), sep='\t',na_rep="NA",header=True,index=False)

#timestamp for the end
print("End time:")
print(datetime.now())
