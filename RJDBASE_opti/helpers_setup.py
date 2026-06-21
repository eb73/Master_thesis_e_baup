"""Helpers for loading or generating input data."""

## ============================= IMPORT =============================
import numpy as np
from scipy.spatial.distance import cdist
from pathlib import Path
import numpy as np
from scipy.io import loadmat
from sklearn.metrics import pairwise_distances

## ============================= CREATE LAPLACIANS AND Y_TRUE =============================
def load_example_caltech():
    data_dir = Path("caltech-101")
    lap_path = data_dir / "laplacians.npy"
    y_path = data_dir / "Y_true.npy"
    mat_path = Path(data_dir / "Caltech101-7.mat") # this is the file that contains pre-processed modalities from the data

    if not lap_path.exists() or not y_path.exists():
        _build_caltech_npy(data_dir, mat_path)

    L_caltech = np.load(lap_path)
    Y_true = np.load(y_path)
    k = 7
    return L_caltech, k, Y_true


def make_laplacian(X, k):
    """
    input: X Nxd feature matrix
    output: Laplacian L (self-tuning spectral clustering style)
    """
    N = X.shape[0]

    dists = pairwise_distances(X, metric='euclidean', squared=True)  # (N, N)

    sorted_dists = np.sort(dists, axis=1)
    sigmas = np.sqrt(sorted_dists[:, k])  # (N,)

    sigma_matrix = np.outer(sigmas, sigmas)
    W = np.exp(-dists / sigma_matrix)
    np.fill_diagonal(W, 0)

    D_inv_sqrt = np.diag(1.0 / np.sqrt(W.sum(axis=1) + 1e-10))
    L_sym = np.eye(N) - D_inv_sqrt @ W @ D_inv_sqrt
    return (L_sym + L_sym.T) / 2


def _build_caltech_npy(data_dir: Path, mat_path: Path):
    """
    Builds laplacians.npy and Y_true.npy from the raw .mat for caltech data.
    Only called when the .npy files don't already exist.
    """
    if not mat_path.exists():
        raise FileNotFoundError(
            f"Missing {mat_path}. Download it from "
            "https://github.com/eb73/Master_thesis_e_baup"
        )

    mat = loadmat(mat_path, struct_as_record=False, squeeze_me=True)

    # build laplacians from the feature views in the .mat file 
    X_views = mat["X"] # shape: (num_modalities,), each entry Nxd
    Y_true = mat["Y"] # ground truth labels

    k = 7 # number of clusters = number of selected categories
    laplacians = [make_laplacian(X_views[i], k) for i in range(X_views.shape[0])]
    L = np.stack(laplacians, axis=0)  #shape: (num_modalities, N, N)

    data_dir.mkdir(parents=True, exist_ok=True)
    np.save(data_dir / "laplacians.npy", L)
    np.save(data_dir / "Y_true.npy", Y_true)



def example_weighted_sbm(N=300, k=6, m=4, sigma1=1, sigma2=1e6, alpha=0.9, beta=0.05, gamma=0.06, theta=0.7, delta=0.2, eps=0.005):
    """Create weighted SBM dataset (synthetic graph 2)."""
    rng = np.random.default_rng(42) # to replicate results
    # sample imbalanced cluster sizes via Dirichlet
    props = rng.dirichlet(np.ones(k))
    sizes = np.round(props * N).astype(int)
    sizes[-1] += N - sizes.sum() # fix rounding error

    Y = np.concatenate([np.full(s, c) for c, s in enumerate(sizes)]) # labels=Y
    rng.shuffle(Y) # random permutation of the cluster labels
    
    # generate modalities
    Ls = []
    for moda in range(m):
        # feature vector xi in R^N following iid N(0,1) 
        x = rng.normal(0,1,size=(N,1)) 
        # block matrix B
        B = np.zeros((k,k))
        sigma = sigma1 
        if moda==0: #B1
            B = np.diag([alpha]*3 + [beta]*3)
            B+=eps
        elif moda==1: #B2
            B = np.diag([beta]*3 + [alpha]*3)
            B+=eps
        elif moda==2: #B3
            B = gamma*np.ones((k,k))
            B+=eps
            sigma = sigma2
        else: #B4
            B = theta*np.eye(k) + delta*(np.ones((k,k))-np.eye(k))
        
        #similarity matrix S
        dist2 = cdist(x,x,metric="sqeuclidean")
        S = np.exp(-dist2/(2*sigma**2))
        #mask matrix C
        C = B[Y[:, None], Y[None, :]]
        #weighted adjacency matrix W
        W = S*C # want element-wise product hence * and not @
        np.fill_diagonal(W,0) #force diagonal to be 0 to avoid self-loops
        W = (W+W.T)/2 #ensure symmetry
        
        #laplacian
        L = sym_norm_laplacian(W)
        Ls.append(L)
    
    return Ls, k, Y
            
        
def gen_stochastic_block_graph(sizes_clusters, p_in, p_out, w_in, w_out, eps):
    N = sum(sizes_clusters)
    labels = np.concatenate([np.full(s,k) for k,s in enumerate(sizes_clusters)])

    W = eps*np.ones((N,N))

    for i in range(N):
        W[i,i] = 0 # force diagonal to be 0 to avoid self-loops
        for j in range(i+1,N):

            p = p_in if labels[i]==labels[j] else p_out
            w = w_in if labels[i]==labels[j] else w_out
            if np.random.rand() < p:
                W[i,j]+=w
                W[j,i]+=w

    return W


def load_example_1(pin, pout, N1=50, eps=0, sigma=0.0):
    L1 = sym_norm_laplacian(gen_stochastic_block_graph([N1,N1,N1], pin, pout, 1, 0.1, eps))
    L2 = sym_norm_laplacian(gen_stochastic_block_graph([2*N1,N1],   pin, pout, 1, 0.1, eps))
    L3 = sym_norm_laplacian(gen_stochastic_block_graph([N1,2*N1],   pin, pout, 1, 0.1, eps))

    def add_noise(L, sigma):
        if sigma == 0: return L
        noise = np.random.randn(*L.shape) * sigma
        noise = (noise + noise.T) / 2  # keep symmetric
        return L + noise

    Ls = [add_noise(L, sigma) for L in [L1, L2, L3]]
    Y_true = np.concatenate([np.ones(N1), 2*np.ones(N1), 3*np.ones(N1)])
    return Ls, 3, Y_true



## ============================= AUXILIARIES =============================
def laplacian(W):
    D=np.diag(W.sum(1))
    return D-W


def sym_norm_laplacian(W):
    D_inv_sqrt = np.diag(1/np.sqrt(W.sum(1)))
    return np.eye(W.shape[0]) - D_inv_sqrt @ W @ D_inv_sqrt