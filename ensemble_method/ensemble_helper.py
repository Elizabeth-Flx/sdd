import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

import networkx as nx
import community as community_louvain

# cluster label encoding to ensure labels are consecutive integers starting from 0 and nonans

def encode_cluster_labels(df: pd.DataFrame):
    # Drop missing rows (simplest policy)
    print('Rows with Nan values: ', df.isna().any(axis=1).sum())
    df = df.dropna(axis=0, how="any")
    print('Remaining rows: ', df.shape[0])

    spot_ids = df.index.astype(str).tolist()
    methods = df.columns.astype(str).tolist()

    cols = []
    # Encode each column to integer labels
    # returns lists 
    for col in df.columns:
        codes, _ = pd.factorize(df[col], sort=True)
        cols.append(codes.astype(np.int32))
    baseCls = np.stack(cols, axis=1)  # N x M
    return baseCls, spot_ids, methods

# need to make cluster ids globally unique across methods: since build ONE big membership matrix of all methods
def make_global_cluster_ids(baseCls: np.ndarray):
    N, M = baseCls.shape
    nClsOrig = baseCls.max(axis=0) + 1
    offsets = np.concatenate([[0], np.cumsum(nClsOrig)[:-1]]).astype(np.int32)
    bcs = baseCls + offsets
    C_total = int(nClsOrig.sum())
    return bcs, nClsOrig, C_total


# builing cluster segment membership matrix
def build_cluster_segments(bcs: np.ndarray, C_total: int):
    N, M = bcs.shape
    rows = bcs.reshape(-1)  # row-major: (spot0 all methods), then spot1...
    cols = np.repeat(np.arange(N, dtype=np.int32), M)  # match row-major
    data = np.ones_like(rows, dtype=np.float32)
    return sp.coo_matrix((data, (rows, cols)), shape=(C_total, N)).tocsr()


def jaccard_similarity_rows(A: sp.csr_matrix) -> np.ndarray:
    """
    Jaccard similarity between rows of sparse binary matrix A (C x N).
    Returns dense (C x C) float array.
    """
    A = A.astype(np.float32)
    inter = (A @ A.T).toarray()  # (C x C)
    sizes = np.asarray(A.sum(axis=1)).reshape(-1)
    union = sizes[:, None] + sizes[None, :] - inter
    with np.errstate(divide="ignore", invalid="ignore"):
        S = np.divide(inter, union, out=np.zeros_like(inter, dtype=np.float32), where=(union > 0))
    return S

import numpy as np
import scipy.sparse as sp

def nmi_similarity_rows(A: sp.csr_matrix, eps: float = 1e-12) -> np.ndarray:
    """
    NMI similarity between rows of sparse binary matrix A (C x N).
    Each row is a cluster-segment membership indicator over N spots.
    Returns dense (C x C) in [0, 1] (approximately; numerical eps may cause tiny negatives).
    """
    A = A.astype(np.float32)
    C, N = A.shape
    # n11 = intersections between cluster segments
    n11 = (A @ A.T).toarray().astype(np.float64)  # (C x C)
    # row sums = cluster sizes
    a = np.asarray(A.sum(axis=1)).reshape(-1).astype(np.float64)
    b = a
    # build all 2x2 counts for each pair (i,j)
    n10 = a[:, None] - n11
    n01 = b[None, :] - n11
    n00 = N - (a[:, None] + b[None, :]) + n11
    # to props
    p11 = n11 / N
    p10 = n10 / N
    p01 = n01 / N
    p00 = n00 / N
    # marginals
    px1 = a / N
    px0 = 1.0 - px1
    py1 = px1
    py0 = 1.0 - py1

    # helper: term = p(xy) * log( p(xy) / (p(x)*p(y)) ), with safe masking
    def mi_term(pxy, px, py):
        denom = (px[:, None] * py[None, :]) + eps
        ratio = (pxy + eps) / denom
        return np.where(pxy > 0, pxy * np.log(ratio), 0.0)

    MI = (
        mi_term(p11, px1, py1) +
        mi_term(p10, px1, py0) +
        mi_term(p01, px0, py1) +
        mi_term(p00, px0, py0)
    )

    # entropies H(X) and H(Y) for each row-variable (binary)
    def H_binary(p1):
        p1 = np.clip(p1, eps, 1.0 - eps)
        p0 = 1.0 - p1
        return -(p1 * np.log(p1) + p0 * np.log(p0))

    Hx = H_binary(px1)  # (C,)
    Hy = Hx.copy()

    denom = np.sqrt(Hx[:, None] * Hy[None, :]) + eps
    NMI = MI / denom

    # numeric cleanup
    NMI = np.clip(NMI, 0.0, 1.0)

    # set perfect self-similarity
    np.fill_diagonal(NMI, 1.0)
    return NMI.astype(np.float32)


# random walk enhancements
def compute_pts_ii(S: np.ndarray, t: int = 20) -> np.ndarray:
    """
    SpatialESD computePTS_II implementation for dense S (C x C).
    """
    S = S.astype(np.float64, copy=True)
    S[S < 0] = 0.0
    np.fill_diagonal(S, 0.0)
    
    row_sum = S.sum(axis=1)
    isolated = row_sum <= 1e-12 # MATLAB uses == 0 before division

    P = np.zeros_like(S)
    P[~isolated] = S[~isolated] / row_sum[~isolated, None]
    P[P < 0] = 0.0  # MATLAB: P(P<0)=0

    tmpP = P.copy()
    inProdP = P @ P.T
    for _ in range(t - 1):
        tmpP = tmpP @ P
        inProdP += tmpP @ tmpP.T

    diag = np.diag(inProdP)
    denom = np.sqrt(np.outer(diag, diag))
    newS = np.divide(inProdP, denom, out=np.zeros_like(inProdP), where=(denom > 0))

    # MATLAB isolation check after normalization
    isolated2 = P.sum(axis=1) < 1e-10
    newS[isolated2, :] = 0.0
    newS[:, isolated2] = 0.0

    return newS.astype(np.float32)


# builing similarity*similarity matrix
def build_method_spot_matrices(bcs: np.ndarray, clsSimRW: np.ndarray) -> list[np.ndarray]:
    """
    For each method m, build A^(m) = clsSimRW[ bcs[:,m], bcs[:,m] ] (N x N).
    Returns list of N x N dense float matrices.
    """
    clsSimRW = clsSimRW.copy()
    np.fill_diagonal(clsSimRW, 1.0)

    N, M = bcs.shape
    mats = []
    for m in range(M):
        c = bcs[:, m]
        A = clsSimRW[np.ix_(c, c)].astype(np.float32, copy=False)
        mats.append(A)
    return mats


def frobenius_norm(X: np.ndarray) -> float:
    return float(np.sqrt(np.sum(X * X)))

def ensdd_weighted_consensus(
    mats: list[np.ndarray],
    lambda_: float | None = None,
    prob_quantile: float = 0.5,
    niter: int = 100,
    epsilon: float = 1e-5
) -> tuple[np.ndarray, np.ndarray, dict]:
    """
    Returns:
      H: consensus matrix (N x N)
      w: weights (M,)
      info: dict with loss history etc.
    """
    M = len(mats)
    if M == 0:
        raise ValueError("No method matrices provided.")
    N = mats[0].shape[0]
    for A in mats:
        if A.shape != (N, N):
            raise ValueError("All matrices must have same shape (N,N).")

    w = np.full(M, 1.0 / M, dtype=np.float64)
    H = sum(w[m] * mats[m] for m in range(M)).astype(np.float64)

    losses = []
    lambdas = []
    diffs = []

    prev_loss = None
    for k in range(1, niter + 1):
        d = np.array([frobenius_norm(mats[m].astype(np.float64) - H) for m in range(M)], dtype=np.float64)

        if k == 1 and lambda_ is None:
            lambda_ = float(np.quantile(d, prob_quantile))
            if lambda_ <= 0:
                lambda_ = 1e-8

        z = -d / lambda_
        z -= z.max()
        w = np.exp(z)
        w /= w.sum()

        # update H
        H = sum(w[m] * mats[m] for m in range(M)).astype(np.float64)

        # recompute distances for updated H  
        d_new = np.array([frobenius_norm(mats[m].astype(np.float64) - H) for m in range(M)], dtype=np.float64)

        loss_main = float(np.sum(w * d_new))
        loss_entropy = float(np.sum(w * np.log(np.clip(w, 1e-15, 1.0))))
        loss_all = loss_main + lambda_ * loss_entropy
        losses.append(loss_all)
        lambdas.append(lambda_)
        if prev_loss is None:
            diffs.append(np.nan)
        else:
            diffs.append(abs(loss_all - prev_loss))

        if prev_loss is not None and abs(loss_all - prev_loss) < epsilon:
            break
        prev_loss = loss_all

    #info = {"loss": np.array(losses), "lambda": np.array(lambdas), "diff": np.array(diffs)}
    return H.astype(np.float32), w.astype(np.float32)

def dense_similarity_to_csr_knn(
    H: np.ndarray,
    k: int = 30,
    keep_diagonal: bool = False,
    make_symmetric: bool = True,
    dtype=np.float32,
    tau : float = 0.0,
) -> sp.csr_matrix:
    """
    Build a kNN graph (per row keep top-k similarities).
    Produces a sparse CSR suitable for Leiden.
    """
    H = np.asarray(H)
    N = H.shape[0]

    if make_symmetric:
        H = 0.5 * (H + H.T)

    if not keep_diagonal:
        np.fill_diagonal(H, -np.inf)

    rows = []
    cols = []
    data = []

    for i in range(N):
        row = H[i]
        # indices of top-k (excluding -inf)
        idx = np.argpartition(row, -k)[-k:]
        # filter out invalid
        idx = idx[np.isfinite(row[idx]) & (row[idx] > tau)]
        rows.extend([i] * len(idx))
        cols.extend(idx.tolist())
        data.extend(row[idx].astype(dtype, copy=False).tolist())

    A = sp.coo_matrix((np.array(data, dtype=dtype),
                       (np.array(rows, dtype=np.int32), np.array(cols, dtype=np.int32))),
                      shape=(N, N)).tocsr()

    # symmetrize by max (common for neighbor graphs)
    if make_symmetric:
        A = A.maximum(A.T)

    # set diagonal to 0
    A.setdiag(0.0)
    A.eliminate_zeros()
    return A

# edge based consensus ##

# make dense to sparse csr igraph
def csr_to_igraph(A: sp.csr_matrix):
    import igraph as ig
    A = A.tocsr()
    A.eliminate_zeros()

    # use only upper triangle edges
    A_triu = sp.triu(A, k=1).tocoo()
    edges = list(zip(A_triu.row.tolist(), A_triu.col.tolist()))
    weights = A_triu.data.astype(float).tolist()

    g = ig.Graph(n=A.shape[0], edges=edges, directed=False)
    g.es["weight"] = weights
    return g

# leiden:: RESOLUTION
def leiden_from_csr(A: sp.csr_matrix, resolution: float = 1.0, seed: int = 0,
                    partition_type: str = "rb"):
    """
    partition_type:
      - "rb": RBConfigurationVertexPartition (common, like 'resolution' modularity)
      - "cpm": CPMVertexPartition (often more stable for controlling cluster sizes)
    """
    import leidenalg as la

    # enforce symmetry + clean diagonal
    A = A.maximum(A.T).tocsr()
    A.setdiag(0.0)
    A.eliminate_zeros()

    g = csr_to_igraph(A)

    if partition_type == "rb":
        part = la.find_partition(
            g, la.RBConfigurationVertexPartition,
            weights=g.es["weight"],
            resolution_parameter=float(resolution),
            seed=int(seed),
        )
    elif partition_type == "cpm":
        part = la.find_partition(
            g, la.CPMVertexPartition,
            weights=g.es["weight"],
            resolution_parameter=float(resolution),
            seed=int(seed),
        )
    else:
        raise ValueError("partition_type must be 'rb' or 'cpm'")

    labels = np.asarray(part.membership, dtype=np.int32)
    return labels

# leiden: K-SPECIFIC
def tune_resolution_for_k(A: sp.csr_matrix, K: int,
                          r_min: float = 0.01, r_max: float = 5.0,
                          max_iter: int = 30, tolerance: int = 0,
                          seed: int = 0, partition_type: str = "rb"):
    """
    Returns: labels, resolution_found, history(list of dicts)
    """
    history = []
    best = None  # (abs_diff, resolution, labels)

    lo, hi = float(r_min), float(r_max)
    for it in range(max_iter):
        mid = (lo + hi) / 2.0
        labels = leiden_from_csr(A, resolution=mid, seed=seed, partition_type=partition_type)
        k_mid = int(np.unique(labels).size)
        diff = abs(k_mid - K)

        history.append({"iter": it, "resolution": mid, "k": k_mid, "diff": diff})

        if best is None or diff < best[0]:
            best = (diff, mid, labels)

        # stop if within tolerance
        if k_mid >= K - tolerance and k_mid <= K + tolerance:
            return labels, mid, history

        # adjust bounds
        if k_mid > K + tolerance:
            hi = mid
        else:
            lo = mid

    # fallback: 
    return best[2], best[1], history

# see which to use leiden
def leiden_final(A: sp.csr_matrix,
                 resolution: float | None = None,
                 K: int | None = None,
                 #auto_bic: bool = False,
                 seed: int = 0,
                 partition_type: str = "rb"):
    """
    Choose exactly one mode:
      - resolution provided -> mode 1
      - K provided -> mode 2
      - auto_bic=True -> mode 3
    Returns: labels, info(dict)
    """
    if sum([resolution is not None, K is not None]) != 1:
        raise ValueError("Choose exactly one of: resolution, K")

    if resolution is not None:
        labels = leiden_from_csr(A, resolution=resolution, seed=seed, partition_type=partition_type)
        return labels, {"mode": "resolution", "resolution": float(resolution),
                        "k": int(np.unique(labels).size)}

    if K is not None:
        labels, res_found, hist = tune_resolution_for_k(
            A, K=K, seed=seed, partition_type=partition_type
        )
        return labels, {"mode": "target_k", "resolution": float(res_found),
                        "k": int(np.unique(labels).size), "history": hist}

    # auto bic
    #labels, res_best, table = select_resolution_by_bic(
    #    A, seed=seed, partition_type=partition_type
    #)
    #return labels, {"mode": "auto_bic", "resolution": float(res_best),
    #                "k": int(np.unique(labels).size), "scan": table}
    

def get_cluster_centroids_from_H(H: np.ndarray, labels: np.ndarray):
    """
    Build cluster-cluster similarity matrix S (k0 x k0) from consensus H (N x N)
    by averaging H over pairs of points across clusters.
    """
    labels = np.asarray(labels)
    clus = np.unique(labels)
    k0 = len(clus)
    idx = [np.where(labels == c)[0] for c in clus]

    S = np.zeros((k0, k0), dtype=np.float64)
    for i in range(k0):
        Ii = idx[i]
        for j in range(i, k0):
            Ij = idx[j]
            # mean similarity between clusters i and j
            val = H[np.ix_(Ii, Ij)].mean()
            S[i, j] = val
            S[j, i] = val
    np.fill_diagonal(S, 1.0)
    return S, clus

def merge_to_exact_k_using_H(H: np.ndarray, labels: np.ndarray, K: int, linkage_method="average"):
    """
    Merge Leiden clusters down to exactly K using hierarchical clustering on (1 - similarity).
    """
    S, clus = get_cluster_centroids_from_H(H, labels)  # k0 x k0 similarity
    k0 = S.shape[0]
    if k0 < K:
        raise ValueError(f"Cannot merge up: only {k0} clusters but need K={K}. Overcluster first.")

    # convert similarity to distance
    D = 1.0 - S
    np.fill_diagonal(D, 0.0)

    # condensed distance vector for scipy linkage
    dvec = squareform(D, checks=False)
    Z = linkage(dvec, method=linkage_method)

    # cluster the k0 “Leiden clusters” into exactly K groups
    super_labels = fcluster(Z, t=K, criterion="maxclust")  # values 1..K

    # map each original cell/spot label to merged label
    mapping = {clus[i]: int(super_labels[i] - 1) for i in range(k0)}  # 0..K-1
    merged = np.vectorize(mapping.get)(labels).astype(np.int32)
    return merged

def leiden_overcluster_then_merge(A_csr: sp.csr_matrix, H: np.ndarray, K: int,
                                  seed: int = 0, partition_type: str = "rb",
                                  r_start: float = 0.1, r_max: float = 20.0,
                                  step: float = 1.25):
    """
    1) run Leiden with increasing resolution until k0 >= K
    2) merge clusters down to exactly K using consensus H
    """
    r = r_start
    last = None
    target = 2 * K

    while r <= r_max:
        labels = leiden_from_csr(A_csr, resolution=r, seed=seed, partition_type=partition_type)
        k0 = np.unique(labels).size
        last = (labels, r, k0)
        if k0 >= target:
            break
        r *= step

    if last is None:
        raise RuntimeError("Leiden did not run.")
    labels, r_used, k0 = last

    if k0 < K:
        raise RuntimeError(
            f"Could not reach at least K clusters. Max resolution {r_max} gave k0={k0} < K={K}."
        )

    merged = merge_to_exact_k_using_H(H, labels, K, linkage_method="average")
    return merged, {"leiden_resolution": r_used, "k0_before_merge": k0, "K_final": K}

# hierarchical for small cases
def hc_from_consensus(H: np.ndarray, K: int, linkage_method: str = "average") -> np.ndarray:
    # symmetrize + clean
    H = 0.5 * (H + H.T)
    np.fill_diagonal(H, 1.0)

    # distance = 1 - similarity
    D = 1.0 - H
    np.fill_diagonal(D, 0.0)

    # scipy linkage needs condensed distance vector
    dvec = squareform(D, checks=False)
    Z = linkage(dvec, method=linkage_method)

    # exactly K clusters (labels 0..K-1)
    labels = fcluster(Z, t=K, criterion="maxclust") - 1
    return labels.astype(np.int32)