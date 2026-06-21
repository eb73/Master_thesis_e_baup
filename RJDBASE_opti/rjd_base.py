"""Functions to run the original RJD-BASE algorithm.
Reference code was given by Artmeis Pados to ensure the use of the same built-in Python functions.
"""
from sklearn.metrics import normalized_mutual_info_score
from sklearn.cluster import KMeans
import numpy as np
from scipy.sparse.linalg import eigsh
import time


def gen_normalized_random_vector(m):
    """ Generate a random vector of size m (entries follow Uniform(0,1)), and normalize it to have sum 1. """
    mu_ = np.random.rand(m)
    mu = mu_/sum(mu_)
    return mu


def rjd_base_iterate(L_list, k):
    """Peforms one trial of RJD-BASE."""
    m = len(L_list)
    #n = L_list[0].shape[0]

    mu = gen_normalized_random_vector(m) # generate mu randomly, and normalize it
    L = sum(mu_i*L_i for mu_i, L_i in zip(mu, L_list)) # compute the weighted sum of L_i's by mu_i's
    w, X = eigsh(L, k=k, which='SA')  # compute the bottom-k eigenvalues and corresponding eigenvectors
    obj = sum(w) # compute sum of bottom-k eigenvalues
    
    # Q, _ = np.linalg.qr(X) # eigsh already returns normalized eigenvectors

    return X, obj

def rjd_base(L_list, k, num_trials=10):
    """Performs num_trials trials of RJD-BASE."""
    X = None
    obj = -np.inf
    for _ in range(num_trials):
        X_trial, obj_trial = rjd_base_iterate(L_list, k)
        if obj_trial > obj:
            obj = obj_trial
            X = X_trial
            
    return X

def rjd_base_timed(L_list, k, num_trials=10):
    """Performs RJD-BASE with num_trials trials and keep track of time per iteration."""
    X = None
    obj = -np.inf
    vec_time = []
    for _ in range(num_trials):
        t0 = time.time()
        X_trial, obj_trial = rjd_base_iterate(L_list, k)
        t1 = time.time()
        vec_time.append(t1-t0)
        if obj_trial > obj:
            obj = obj_trial
            X = X_trial
            
    return X, np.mean(vec_time)

def rjd_base_with_memory(L_list, k, num_trials, Y_true):
    """Performs RJD-BASE with num_trials and keep track of NMI/objective at each iteration."""
    # Outputs:
    nmi_vec = [] #list of NMI from every trial
    #V_list = [] # list of embedding matrix from every trial
    obj_vec = [] # sum of bottom-k eigenvalues from every trial (warning: lambda_0 is 0, so actually take bottom-(k+1) eigenvalues)

    obj_best = -np.inf # highest eigenvalue sum
    X_best = None # embedding matrix corresponding to highest eigenvalue sum
    nmi_best = None # NMI corresponding to highest eigenvalue sum
    for _ in range(num_trials):
        # compute RJD-BASE iterate
        X_trial, obj_trial = rjd_base_iterate(L_list, k)
        # compute predicted clusters at this trial
        labels = KMeans(n_clusters=k, random_state=0).fit_predict(X_trial)
        nmi_trial = normalized_mutual_info_score(Y_true, labels)
        # append values of objective and NMI to results vectors
        obj_vec.append(obj_trial)
        nmi_vec.append(nmi_trial)
        # update best value
        if obj_trial > obj_best:
            obj_best = obj_trial
            X_best = X_trial
            nmi_best = nmi_trial
            
    return X_best, obj_best, nmi_best, obj_vec, nmi_vec
    
