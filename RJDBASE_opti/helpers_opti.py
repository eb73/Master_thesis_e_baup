import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import normalized_mutual_info_score as nmi
from sklearn.metrics import adjusted_mutual_info_score as ami
from sklearn.metrics import rand_score as RI
from sklearn.metrics import adjusted_rand_score as ARI
from tqdm import tqdm
import time

## ============================= OBJECTIVE AND GRADIENT =============================
def softmax(alpha, x):
    """Compute the alpha-softmax function (smooth approximation of max)
    alpha is a scalar that controls the smoothness of the approximation. As alpha -> inf, softmax -> max.
    x is a list or array of values to compute the softmax over.
    """
    top = sum(xi*np.exp(alpha*xi) for xi in x)
    bottom = sum(np.exp(alpha*xi) for xi in x)
    return top/bottom

def quasimax(alpha, x):
    """Compute the alpha-quasimax function (smooth approximation of max)
    alpha is a scalar that controls the smoothness of the approximation. As alpha -> inf, quasimax -> max.
    x is a list or array of values to compute the quasimax over.
    """
    s = sum(np.exp(alpha*xi) for xi in x)
    return (1/alpha)*np.log(s)

def f(L,X):
    return np.trace(X.T @ L @ X)

def d_f(L,X):
    return 2*L @ X

def sG(X, L_list): # TODO change name to sG or infinity norm s_L
    """Compute s_G(X) = ||vec(Tr(X^T L_i X))||_inf, where L_list is the list of matrices L_i."""
    m = len(L_list)
    fX = np.array([f(L_list[i], X) for i in range(m)]) # compute fX = [Tr(X^T L_1 X), ..., Tr(X^T L_m X)]
    return np.linalg.norm(fX, ord=np.inf) # np.max(fX) gives the same result when L are Laplacians (symmetric positive semi-definite)

def weights_softmax(alpha, L_list, X):
    """Compute the alpha-softmax weights for the function f(L,X)
    The softmax approximation can be written g(X) = sum w_i * f(L_i, X).
    """
    m = len(L_list)
    fX = np.array([f(Li,X) for Li in L_list])
    z = alpha*fX
    expz = np.exp(z)
    weights = np.zeros(m)
    for i in range(m):
        weights[i] = expz[i]/sum(expz[j] for j in range(m))
    return weights

def g(alpha, L_list, X, f):
    """Compute the function g(X) = sum w_i * f(L_i, X) = softmax(alpha, vec(Tr(X^TLX))) where w_i are the softmax weights."""
    m = len(L_list)
    weights = weights_softmax(alpha, L_list, X)
    return sum(weights[i]*f(L_list[i], X) for i in range(m))

def d_g(alpha, L_list, X, f, df):
    """Compute the gradient of g(X) = sum w_i * f(L_i, X)."""
    m = len(L_list)
    weights = weights_softmax(alpha, L_list, X)
    fX = np.array([f(L_list[i], X) for i in range(m)])
    gX = g(alpha, L_list, X, f)
    return sum(weights[i]*(1+alpha*(fX[i]-gX))*df(L_list[i], X) for i in range(m))

def h(alpha, L_list, X, f):
    """Compute the function h(X) = quasimax(alpha, vec(Tr(X^TLX)))."""
    m = len(L_list)
    fX = np.array([f(L_list[i], X) for i in range(m)])
    return (1/alpha)*np.log(sum(np.exp(alpha*fX[i]) for i in range(m)))

def d_h(alpha, L_list, X, df):
    """Compute the gradient of h(X) = quasimax(alpha, vec(Tr(X^TLX)))."""
    m = len(L_list)
    weights = weights_softmax(alpha, L_list, X)
    return sum(weights[i]*df(L_list[i], X) for i in range(m))

## ============================= SCORES =============================
def subspace_invariance(X, L_list):
    """ Compute the maximum Frobenius norm of LX-XX^TLX, for L in L_list.
    This measures how much the subspace spanned by the columns of X is invariant under the action of the matrices in L_list."""
    P = X @ X.T
    vals = []
    for L in L_list:
        diff = L@X - P@L@X
        vals.append(np.linalg.norm(diff, 'fro'))
    return max(vals)

def predict_labels(X, k):
    """ Compute the predicted labels of each data point, based on approximated embedding X. k=number of clusters."""
    labels = KMeans(
        n_clusters=k,random_state=0,n_init=1,algorithm="elkan",max_iter=100
        ).fit_predict(X)
    return labels

def compute_scores(X, L_list, k, Y_true):
    """Utility function to compute directly the scores sG(X), subspace invariance, NMI & AMI."""
    labels = predict_labels(X, k)
    return [sG(X, L_list), subspace_invariance(X, L_list), nmi(Y_true, labels), ami(Y_true, labels)] 
# the RI and ARI performs similarly as nmi and ami, so it's not relevant to add RI(Y_true,labels), ARI(Y_true, labels)]

def print_scores(v):
    print(f"sG={v[0]:.3}, subspace invariance={v[1]:.3f}, nmi={v[2]:.3f}, ami={v[2]:.3f}")
    return

## ============================= RIEMANNIAN GRADIENT DESCENT =============================
def sym(A):
    """ Symmetric part of a matrix. """
    return (A + A.T)/2


def riemannian_gd(obj_func, dobj_func, L_list, X0, step_size, max_iter, tol_conv, iter_window=5):
    """Perform riemannian gradient descent to minimize the objective function obj_func with derivative dobj_func."""
    X_old = X0
    X_new = None
    obj_vec = []
    sG_vec = []
    time_vec = []
    for i in tqdm(range(max_iter)):
        t0 = time.time()
        euclidean_grad = dobj_func(L_list, X_old)
        riemann_proj = euclidean_grad - X_old @ (sym(X_old.T @ euclidean_grad)) # project to tangent space
        grad_norm = np.linalg.norm(riemann_proj, 'fro')
        X_new = X_old - step_size*riemann_proj # gradient step
        X_new, _ = np.linalg.qr(X_new) # orthonormalize
        t1 = time.time()
        time_vec.append(t1-t0)
        obj_vec.append(obj_func(L_list, X_new))
        sG_vec.append(sG(X_new, L_list))
        
        if grad_norm < tol_conv: 
            print("Convergence reached at iteration ", i, " (grad norm)")
            break
        
        if i>iter_window and abs(obj_vec[-1] - obj_vec[-iter_window]) / abs(obj_vec[-iter_window]) < tol_conv:
            print("Convergence reached at iteration ", i, " (objective stable)")
            break
        
        X_old = X_new
    i +=1 # because start at 0
    return X_new, obj_vec, sG_vec, i, np.mean(time_vec)