# Helper Module for Network Analysis Thresholding

import pandas as pd
import numpy as np
import scipy
import matplotlib.pyplot as plt

from numpy import linalg as LA

import empyricalRMT as rmt
from empyricalRMT.construct import generate_eigs
from empyricalRMT.eigenvalues import Eigenvalues

def get_similarity_matrix(df_counts_rel):
    # make covariance matrix (OTUs x OTUs)
    df_counts_rel_corr = df_counts_rel.corr(method="pearson")
    # make similarity data
    S_corr = df_counts_rel_corr.abs()
    return S_corr


def get_thresholded_matrix(S_corr, s_t):
     # QUESTION: should this make -r --> -1
    # apply threshold
    A = S_corr >= s_t
    A = A.astype(int)
    # replace values on diagonals with 0
    np.fill_diagonal(A.values, 0)
    # keep non-zero rows and columns
    A = A.loc[:, (A != 0).any(axis=0)]
    A = A.loc[(A != 0).any(axis=1), :]
    # get a sense of how many OTUs remain
    print("shape", A.shape)
    return A


def get_unfolded_eigs(S_corr, s_t, verbose=False):
    A = get_thresholded_matrix(S_corr, s_t)
    # calculate eigenvalues
    w, _ = LA.eig(A)
    # sort eigenvalues
    # To test NNSD distribution, order the eigenvalue, smallest first
    w_sorted = np.sort(w)
    # use empiricalRMT library
    eigs = Eigenvalues(w_sorted)
    # unfold the eigenvalues by fitting to cubic spline
    unfolded = eigs.unfold(smoother="poly", degree=3)
    if verbose:
        # plot some classic observables and compare to theory
        ensembles = ["poisson", "goe"]  # theoretically expected curves to plot
        unfolded.plot_nnsd(ensembles=ensembles)  # nearest neighbours spacings
    return unfolded.vals


def get_expected_poisson(unfolded):
    p = np.pi
    _spacings = np.diff(unfolded)
    s = np.linspace(_spacings.min(), _spacings.max(), len(unfolded))
    poisson = np.exp(-s)
    return poisson


def get_expected_GOE(unfolded):
    p = np.pi
    _spacings = np.diff(unfolded)
    s = np.linspace(_spacings.min(), _spacings.max(), len(unfolded))
    goe = ((p * s) / 2) * np.exp(-(p / 4) * s * s)
    return goe


def chi_squared(observed, expected):
    d = (observed-expected)**2
    d = d/expected
    d = d.sum()
    return d


def find_RMT_threshold(df_counts_rel, s_tb, alpha=0.05):
    """
    df_counts_rel: pandas DataFrame, Samples x OTUs relative counts
    s_tb: baseline theshold starting point
    alpha: critcal value for chi-square test

    Use Random Matrix Theory to find a threshold to 
    separate random noise (following GOE NNSD) 
    from non-random signal (following Poisson NNSD).
    
    H0: P(d) follows the Poisson distribution
    H1: P(d) does not follow the Poisson distribution
    """

    S_corr = get_similarity_matrix(df_counts_rel)
    
    # set significance threshold
    s_t = s_tb
    
    coarse_step = 0.1
    fine_step = 0.01
    is_finding_finer_theshold = False    
    
    s_t_vals = []
    unfolded_len_vals = []
    Xsq_crit_vals = []
    Xsq_poisson_vals = []
    Xsq_GOE_vals = []
    
    while(s_t < 1.0):
        unfolded = get_unfolded_eigs(S_corr, s_t, verbose=False)
        if len(unfolded) < 2:
            break

        Xsq_crit = scipy.stats.chi2.ppf(1-alpha, df=len(unfolded)-1)
        
        # chi-square for poisson distribution
        expected_poisson = get_expected_poisson(unfolded)
        Xsq_poisson = chi_squared(unfolded, expected_poisson)
        
        # chi-square for gaussian distribution
        expected_GOE = get_expected_GOE(unfolded)
        Xsq_GOE = chi_squared(unfolded, expected_GOE)
        
        # savepoint
        s_t_vals.append(s_t)
        unfolded_len_vals.append(len(unfolded))
        Xsq_crit_vals.append(Xsq_crit)
        Xsq_poisson_vals.append(Xsq_poisson)
        Xsq_GOE_vals.append(Xsq_GOE)        
        print("s_t", s_t)
        print("len(unfolded)", len(unfolded))
        print("Xsq_crit", Xsq_crit)
        print("Xsq_poisson", Xsq_poisson)
        print("Xsq_GOE", Xsq_GOE)
        
        if (is_finding_finer_theshold):
            if (Xsq_poisson > Xsq_crit):
                # reject null hypothesis that NNSD is Poisson
                print("--reject Poisson")
                s_t += fine_step
            
            else:
                # do not reject null hypothesis that NNSD is Poisson
                print("--do not reject Poisson (fine)!")
                break
    
        else:
            if (Xsq_poisson > Xsq_crit):
                # reject null hypothesis that NNSD is Poisson
                print("--reject Poisson")
                s_t += coarse_step

            else:
                # do not reject null hypothesis that NNSD is Poisson
                # find a finer threshold
                print("--do not reject Poisson (coarse)")
                s_t -= coarse_step
                is_finding_finer_theshold = True
    
    print("s_t_vals", s_t_vals)
    print("unfolded_len_vals", unfolded_len_vals)
    print("Xsq_crit_vals", Xsq_crit_vals)
    print("Xsq_poisson_vals", Xsq_poisson_vals)
    print("Xsq_GOE_vals", Xsq_GOE_vals)

    print("final threshold:", s_t)     
    return s_t


def visualize_RMT_threshold(df_counts_rel, threshold, alpha=0.05, show_poisson=False):
    """
    df_counts_rel: pandas DataFrame, Samples x OTUs relative counts
    threshold: theshold on similarity correlation matrix
    alpha: critcal value for chi-square test
    """
    S_corr = get_similarity_matrix(df_counts_rel)
    unfolded = get_unfolded_eigs(S_corr, threshold, verbose=True)
    Xsq_crit = scipy.stats.chi2.ppf(1-alpha, df=len(unfolded)-1)

    # chi-square for poisson distribution
    expected_poisson = get_expected_poisson(unfolded)
    Xsq_poisson = chi_squared(unfolded, expected_poisson)

    print(f"Xsq_poisson {Xsq_poisson}")
    print(f"Xsq_crit {Xsq_crit} at alpha={alpha}")

    if show_poisson:
        plt.figure()
        plt.hist(get_expected_poisson(unfolded))
        plt.title("Expected Unfolded NNSD Poisson Distribution")
    



