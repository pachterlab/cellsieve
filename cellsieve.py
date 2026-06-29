#!/usr/bin/env python3
"""
Cell Sieve

Description:
Combinatorially compressed gene panels to distinguish cell types in scRNAseq data.
This script selects optimal gene panels for cell type classification from scRNA-seq data 
using an Elastic Net global prior combined with a Mixed Integer Linear Programming (MILP) solver.

Example Commands:

1. Default Run (runs a single panel size of 10):
   cellsieve --adata-path /path/to/data.h5ad --out-dir ./output --n-workers 8

2. Custom Single Panel Size:
   cellsieve --adata-path /path/to/data.h5ad --out-dir ./output --panel-size 40

3. Custom Titration (runs multiple specific sizes sequentially, e.g., 5, 15, and 30):
   cellsieve --adata-path /path/to/data.h5ad --out-dir ./output --titrate-panels 5 15 30

4. Misannotation Detection & Correction (flags and cleans bad annotations):
   cellsieve --adata-path /path/to/data.h5ad --out-dir ./output --correct-misannotations --panel-size 20

5. Find Redundant Gene Swaps (finds correlated substitutes for selected genes):
   cellsieve --adata-path /path/to/data.h5ad --out-dir ./output --find-redundant-swaps --corr-threshold 0.80
"""

import argparse
import time
import sys
import os
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.linear_model import LogisticRegression, SGDClassifier
from sklearn.metrics import classification_report
from sklearn.model_selection import StratifiedShuffleSplit, StratifiedKFold, cross_val_predict, cross_val_score
from sklearn.preprocessing import StandardScaler, MaxAbsScaler
import warnings
import gurobipy as gp
from gurobipy import GRB


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

def as_dense_array(matrix, dtype=float):
    return np.asarray(matrix.toarray() if hasattr(matrix, "toarray") else matrix, dtype=dtype)

def parse_args():
    parser = argparse.ArgumentParser(description="Unified Elastic Net MILP Gene-Panel Selector.")
    parser.add_argument("--celltype-column", default="celltype", help="Column in adata.obs containing cell type labels.")
    parser.add_argument("--adata-path", required=True, help="Path to your adata .h5ad file")
    parser.add_argument("--out-dir", default="./elastic_sieve_output", help="Directory to save outputs.")
    parser.add_argument("--random-seed", type=int, default=42)

    # Mode Toggles & Alpha
    parser.add_argument("--alpha", type=float, default=0.7, help="Specify the Elastic Net alpha (L1 Ratio) from 0.0 to 1.0. Default 0.7")
    parser.add_argument("--titrate-panels", nargs="+", type=int, default=None, help="Run multiple panel sizes sequentially. Usage: --titrate-panels 5 10 20 40")
    parser.add_argument("--panel-size", type=int, default=20, help="Target panel size if not titrating. Default is 10.")
    parser.add_argument("--max-cells", type=int, default=15000, help="Max cells to use for the Global Prior fitting.")
    parser.add_argument("--correct-misannotations", action="store_true", help="Flag cells with confidently incorrect labels, save a CSV report, and output a cleaned h5ad file.")
    parser.add_argument("--min-correction-margin", type=float, default=0.30, help="Minimum probability margin to override an existing annotation. Default is 0.30.")
    
    # Redundant Swaps
    parser.add_argument("--find-redundant-swaps", action="store_true", help="Find genes highly correlated with the MILP-selected genes to use as substitutes.")
    parser.add_argument("--corr-threshold", type=float, default=0.75, help="Minimum Pearson R to consider a valid swap. Default is 0.75.")

    # MILP Solver Parameters
    parser.add_argument("--gurobi-license-file", default=None, help="Optional path to gurobi.lic. Sets GRB_LICENSE_FILE for this run.")
    parser.add_argument("--mip-gap", type=float, default=0.10)
    parser.add_argument("--time-limit", type=float, default=300)
    parser.add_argument("--n-workers", type=int, default=24)
    return parser.parse_args()


# =============================================================================
# Elastic Net Prior Generation (Simplified)
# =============================================================================

def compute_global_elasticnet(adata, y, args):
    log(f"Using provided Alpha (L1 Ratio): {args.alpha:.4f}")
    
    # 1. Prepare stable sparse matrix representation with Stratified Sampling
    all_cells = np.arange(adata.n_obs)
    if len(all_cells) > args.max_cells:
        log(f"Subsampling {len(all_cells)} cells down to {args.max_cells} using StratifiedShuffleSplit...")
        sss = StratifiedShuffleSplit(n_splits=1, train_size=args.max_cells, random_state=args.random_seed)
        train_idx, _ = next(sss.split(all_cells, y))
    else:
        train_idx = all_cells

    X_scaled = MaxAbsScaler().fit_transform(adata[train_idx].X)
    y_sub = y[train_idx]
    
    sgd_model = SGDClassifier(
        loss="log_loss", penalty="elasticnet", class_weight="balanced", 
        l1_ratio=args.alpha, alpha=0.001, max_iter=1000, tol=1e-3, 
        learning_rate="invscaling", eta0=0.1, random_state=args.random_seed, n_jobs=args.n_workers 
    )

    start_time = time.time()
    log("Fitting Global Elastic Net prior model...")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        W = sgd_model.fit(X_scaled, y_sub).coef_
        
    log(f"Global Elastic Net fit completed in {time.time() - start_time:.1f} seconds.")

    # 2. Extract non-zero coefficients
    all_indices = np.arange(adata.n_vars)
    global_nonzero_mask = np.any(np.abs(W) > 1e-4, axis=0)
    
    pos_indices = all_indices[global_nonzero_mask]
    pos_scores = np.max(np.abs(W[:, global_nonzero_mask]), axis=0)
        
    return pos_indices, pos_scores


def localized_subset_refit(adata, y, candidate_pool_idx, candidate_scores, args, panel_size_val):
    n_top_candidates = min(panel_size_val * 4, len(candidate_pool_idx))
    
    top_local_idx = np.argsort(candidate_scores)[::-1][:n_top_candidates]
    final_candidate_idx = candidate_pool_idx[top_local_idx]
    final_candidate_priors = candidate_scores[top_local_idx] 
    
    log(f"Forwarding top {len(final_candidate_idx)} candidate genes to the MILP matrix.")
    
    X_sub = as_dense_array(adata[:, final_candidate_idx].X)
    sub_scaler = StandardScaler().fit(X_sub)
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sub_ref_model = LogisticRegression(
            penalty="l2", class_weight="balanced", 
            n_jobs=args.n_workers, max_iter=2000, random_state=args.random_seed
        ).fit(sub_scaler.transform(X_sub), y)
    
    return final_candidate_idx, sub_scaler, sub_ref_model, final_candidate_priors
    


# =============================================================================
# Shared Constraint Setup & MILP Functions
# =============================================================================

def select_mrmr_constraints(y_true, decision_scores, n_pairs_target=3000):
    n_cells, n_classes = len(y_true), decision_scores.shape[1]
    
    true_scores = decision_scores[np.arange(n_cells), y_true]
    margins = true_scores[:, np.newaxis] - decision_scores
    margins[np.arange(n_cells), y_true] = np.inf 
    
    potential_idx = []
    for i in range(n_cells):
        hardest_wrong = np.argsort(margins[i, :])[:5]
        for w in hardest_wrong:
            potential_idx.append([i, y_true[i], w, margins[i, w]])
            
    df_pot = pd.DataFrame(potential_idx, columns=['cell_idx', 'true_label', 'wrong_label', 'margin'])
    df_pot['relevance'] = 1.0 / (df_pot['margin'] + 1e-3)
    df_pot['relevance'] = (df_pot['relevance'] - df_pot['relevance'].min()) / (df_pot['relevance'].max() - df_pot['relevance'].min() + 1e-8)
    
    selected_indices = [df_pot['relevance'].idxmax()]
    cell_hit_count = np.zeros(n_cells)
    pair_hit_count = np.zeros((n_classes, n_classes))
    
    row = df_pot.iloc[selected_indices[0]]
    cell_hit_count[int(row.cell_idx)] += 1
    pair_hit_count[int(row.true_label), int(row.wrong_label)] += 1
    
    relevance_arr = df_pot['relevance'].values
    cell_idx_arr = df_pot['cell_idx'].values.astype(int)
    true_arr = df_pot['true_label'].values.astype(int)
    wrong_arr = df_pot['wrong_label'].values.astype(int)
    
    for _ in range(1, n_pairs_target):
        redundancy = (cell_hit_count[cell_idx_arr] * 2.0) + (pair_hit_count[true_arr, wrong_arr])
        score = (relevance_arr) - (redundancy / (redundancy.max() + 1e-8))
        score[selected_indices] = -np.inf
        
        best_idx = np.argmax(score)
        selected_indices.append(best_idx)
        cell_hit_count[cell_idx_arr[best_idx]] += 1
        pair_hit_count[true_arr[best_idx], wrong_arr[best_idx]] += 1

    final_pairs = df_pot.iloc[selected_indices].copy()
    counts = final_pairs.groupby(['true_label', 'wrong_label']).size().reset_index(name='pair_count')
    final_pairs = final_pairs.merge(counts, on=['true_label', 'wrong_label'])
    final_pairs['macro_weight'] = 1.0 / final_pairs['pair_count']
    return final_pairs

def build_milp_matrix(adata, candidate_idx, scaler, model, pair_df):
    unique_cells = np.asarray(sorted(pair_df["cell_idx"].unique()), dtype=int)
    cell_to_row = {cell: i for i, cell in enumerate(unique_cells)}
    X = scaler.transform(as_dense_array(adata[unique_cells, candidate_idx].X, dtype=float))
    W, b = model.coef_, model.intercept_
    A = np.zeros((len(pair_df), len(candidate_idx)), dtype=np.float64)
    bias = np.zeros(len(pair_df), dtype=np.float64)
    
    for s, row in enumerate(pair_df.itertuples(index=False)):
        r = cell_to_row[int(row.cell_idx)]
        t, w = int(row.true_label), int(row.wrong_label)
        A[s, :] = X[r, :] * (W[t] - W[w])
        bias[s] = b[t] - b[w]
    return A, bias


def solve_milp(A, bias, weights, args, panel_size, candidate_priors):
    n_constraints, n_genes = A.shape

    model = gp.Model("unified_gene_panel_milp")
    model.Params.TimeLimit = args.time_limit
    model.Params.MIPGap = args.mip_gap
    model.Params.Threads = args.n_workers
    model.Params.MIPFocus = 1
    model.Params.Method = 3        
    model.Params.Cuts = 0          
    model.Params.Heuristics = 0.95  
    
    z = model.addVars(n_genes, vtype=GRB.BINARY, name="z")
    slack = model.addVars(n_constraints, lb=0.0, vtype=GRB.CONTINUOUS, name="slack")
    
    for s in range(n_constraints):
        expr = gp.LinExpr(float(bias[s]))
        nz = np.flatnonzero(A[s, :])
        if len(nz): expr.addTerms([float(A[s, g]) for g in nz], [z[int(g)] for g in nz])
        expr.addTerms(1.0, slack[s])
        model.addConstr(expr >= 1.0, name=f"margin_{s}")
        
    gene_count = gp.quicksum(z[g] for g in range(n_genes))
    model.addConstr(gene_count == panel_size, name="exact_panel_size")

    weighted_slack = gp.quicksum(float(weights[s]) * slack[s] for s in range(n_constraints))
    # Feed the raw priors directly; no normalization
    prior_reward = gp.quicksum(float(candidate_priors[g]) * z[g] for g in range(n_genes))
    
    # Static weights: The solver balances margin vs. raw prior without alpha bias
    margin_focus = 1.0
    prior_focus = 1.0  # You may need to tune this depending on the raw magnitude of your coefficients

    obj_expr = margin_focus * weighted_slack - prior_focus * prior_reward
    
    model.setObjective(obj_expr)
    model.ModelSense = GRB.MINIMIZE

    model.optimize()
    if model.SolCount == 0: raise RuntimeError("Gurobi produced no solution.")
    return np.array([g for g in range(n_genes) if z[g].X > 0.5], dtype=int)


# =============================================================================
# Redundant Solutions (Gene Swaps)
# =============================================================================

import scipy.sparse as sp
from joblib import Parallel, delayed

def compute_redundant_swaps(adata, selected_genes, out_dir, args):
    log(f"Calculating redundant gene swaps (Pearson R >= {args.corr_threshold}) using {args.n_workers} workers...")
    all_genes = adata.var_names.tolist()
    
    # 1. Isolate the small matrix of selected genes and make ONLY that dense
    valid_selected = [g for g in selected_genes if g in all_genes]
    if not valid_selected:
        log("No valid selected genes found for swap calculation.")
        return
        
    X_sel = adata[:, valid_selected].X
    if sp.issparse(X_sel):
        X_sel = X_sel.toarray()
        
    # Standardize selected genes
    X_sel_std = (X_sel - X_sel.mean(axis=0)) / (X_sel.std(axis=0) + 1e-8)
    
    n_cells, n_all_genes = adata.shape
    correlations = np.zeros((n_all_genes, len(valid_selected)))
    
    # 2. Define the worker function for a single chunk
    def _correlate_chunk(start_idx, end_idx):
        X_chunk = adata[:, start_idx:end_idx].X
        if sp.issparse(X_chunk):
            X_chunk = X_chunk.toarray()
            
        # Standardize the chunk
        X_chunk_std = (X_chunk - X_chunk.mean(axis=0)) / (X_chunk.std(axis=0) + 1e-8)
        
        # Vectorized Pearson correlation for the chunk
        return start_idx, end_idx, np.dot(X_chunk_std.T, X_sel_std) / float(n_cells)

    # 3. Process chunks in parallel using threading (avoids copying adata in memory)
    chunk_size = 2000
    results = Parallel(n_jobs=args.n_workers, backend="threading")(
        delayed(_correlate_chunk)(i, min(i + chunk_size, n_all_genes)) 
        for i in range(0, n_all_genes, chunk_size)
    )
    
    # Reassemble the results into the correlation matrix
    for start_idx, end_idx, chunk_corr in results:
        correlations[start_idx:end_idx, :] = chunk_corr
        
    # 4. Extract the high-correlation pairs
    rows = []
    for j, sel_gene in enumerate(valid_selected):
        high_corr_idx = np.where(correlations[:, j] >= args.corr_threshold)[0]
        
        for idx in high_corr_idx:
            cand_gene = all_genes[idx]
            if cand_gene != sel_gene:
                rows.append({
                    "MILP_Selected_Gene": sel_gene,
                    "Substitute_Gene": cand_gene,
                    "Pearson_Correlation": round(float(correlations[idx, j]), 4)
                })
                
    # 5. Save results
    if rows:
        df = pd.DataFrame(rows)
        df = df.sort_values(by=["MILP_Selected_Gene", "Pearson_Correlation"], ascending=[True, False])
        out_csv = Path(out_dir) / "redundant_gene_swaps.csv"
        df.to_csv(out_csv, index=False)
        log(f"Found {len(df)} redundant gene swaps. Saved to {out_csv.name}")
    else:
        log(f"No redundant genes found above correlation {args.corr_threshold}.")


# =============================================================================
# Validation Outputs
# =============================================================================

def generate_outputs(out_dir, size_str, adata, y_full, selected_idx, unique_types, args):
    sub_out_dir = Path(out_dir) / f"panel_size_{size_str}_alpha_{args.alpha:.2f}"
    sub_out_dir.mkdir(parents=True, exist_ok=True)
    
    selected_names = adata.var_names[selected_idx].tolist()
    X_final_dense = as_dense_array(adata[:, selected_idx].X)
    X_final_scaled = StandardScaler().fit_transform(X_final_dense)
    
    # 1. Generalization Power (Stratified CV Split)
    class_counts = pd.Series(y_full).value_counts()
    n_splits = min(3, class_counts.min()) 
    cv_macro_f1 = np.nan
    if n_splits >= 2:
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=args.random_seed)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cv_scores = cross_val_score(LogisticRegression(class_weight="balanced", max_iter=2000, random_state=args.random_seed), 
                                        X_final_scaled, y_full, cv=cv, scoring="f1_macro", n_jobs=args.n_workers)
            cv_macro_f1 = np.mean(cv_scores)

    # 2. Compression Capacity (Global Matrix Fit)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        final_model = LogisticRegression(class_weight="balanced", max_iter=2000, random_state=args.random_seed).fit(X_final_scaled, y_full)
    pd.DataFrame(final_model.coef_, index=unique_types, columns=selected_names).to_csv(sub_out_dir / "model_coefficients.csv")

    final_pred = final_model.predict(X_final_scaled)
    report = classification_report(y_full, final_pred, target_names=unique_types, digits=4)
    report_dict = classification_report(y_full, final_pred, target_names=unique_types, output_dict=True, zero_division=0)
    global_macro_f1 = report_dict["macro avg"]["f1-score"]
    
    # 3. Write Metrics Summary
    pd.DataFrame([{
        "adata": args.adata_path,
        "num_celltypes": len(unique_types),
        "alpha_complexity": args.alpha,
        "requested_panel_size": size_str,
        "actual_panel_size": len(selected_names),
        "global_macro_f1": global_macro_f1,
        "cv_macro_f1": cv_macro_f1,
        "accuracy": report_dict["accuracy"]
    }]).to_csv(sub_out_dir / "summary_metrics.csv", index=False)
    
    pd.Series(selected_names).to_csv(sub_out_dir / "selected_genes.txt", index=False, header=False)
    (sub_out_dir / "classification_report.txt").write_text(report)
    
    print(f"\n{'='*50}\nFINAL REPORT: {len(selected_names)} GENES (Alpha: {args.alpha:.2f} | Size: {size_str})\n{'='*50}")
    print(f"Number of Cell Types:              {len(unique_types)}")
    print(f"Compression Capacity (Global F1):  {global_macro_f1:.4f}")
    print(f"Generalization Power (CV F1):      {cv_macro_f1:.4f}\n")

    # --- NEW: Redundant Gene Swaps ---
    if args.find_redundant_swaps:
        compute_redundant_swaps(adata, selected_names, sub_out_dir, args)

    # --- UPDATED: Misannotation Detection & Correction ---
    if args.correct_misannotations:
        log("Generating unbiased Out-of-Fold predictions to evaluate labels...")
        
        clean_eval_model = LogisticRegression(class_weight="balanced", max_iter=2000, random_state=args.random_seed)
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            probs = cross_val_predict(
                clean_eval_model, X_final_scaled, y_full, 
                cv=5, method="predict_proba", n_jobs=args.n_workers
            )
        
        annotated_probs = probs[np.arange(len(y_full)), y_full]
        predicted_probs = np.max(probs, axis=1)
        final_oof_pred = np.argmax(probs, axis=1)
        
        margins = predicted_probs - annotated_probs
        
        mis_df = pd.DataFrame({
            "cell_barcode": adata.obs_names,
            "annotated_label": [unique_types[y] for y in y_full],
            "predicted_label": [unique_types[p] for p in final_oof_pred],
            "annotated_prob": annotated_probs,
            "predicted_prob": predicted_probs,
            "margin": margins
        })
        
        suspicious = mis_df[mis_df["annotated_label"] != mis_df["predicted_label"]].copy()
        suspicious = suspicious[suspicious["margin"] > args.min_correction_margin]
        
        # 1. Output the CSV report
        suspicious.to_csv(sub_out_dir / "flagged_misannotations.csv", index=False)
        log(f"Flagged {len(suspicious)} suspicious cells. Saved to flagged_misannotations.csv")
            
        # 2. Update and save the AnnData object
        log("Applying corrected labels to the AnnData object...")
        
        # Create a safe new column name based on the input column
        corrected_col = f"{args.celltype_column}_corrected"
        
        # Default everything to the original labels
        adata.obs[corrected_col] = adata.obs[args.celltype_column].copy()
        
        # Create a dictionary mapping only the suspicious barcodes to their new predicted labels
        correction_map = suspicious.set_index("cell_barcode")["predicted_label"].to_dict()
        
        # Apply the mapping only to cells that exist in the correction_map
        adata.obs[corrected_col] = adata.obs_names.map(
            lambda x: correction_map.get(x, adata.obs.loc[x, args.celltype_column])
        )
        
        # Save the new cleaned h5ad file into the output directory
        out_name = f"CLEANED_{Path(args.adata_path).name}"
        out_path = sub_out_dir / out_name
        adata.write_h5ad(out_path)
        log(f"Saved {len(suspicious)} corrections to new column '{corrected_col}'. File saved to: {out_name}")
    # ------------------------------------

# =============================================================================
# Main Orchestrator
# =============================================================================

def main():
    args = parse_args()

    if args.gurobi_license_file:
        os.environ["GRB_LICENSE_FILE"] = str(Path(args.gurobi_license_file).expanduser())
    
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    log(f"Reading Dataset: {args.adata_path}")
    adata = sc.read_h5ad(args.adata_path)

    celltype_series = pd.Series(adata.obs[args.celltype_column].astype(str).values)
    unique_types = sorted(celltype_series.unique())
    y_full = celltype_series.map({ct: i for i, ct in enumerate(unique_types)}).values.astype(int)

    # Determine runtime loops based on args
    if args.titrate_panels:
        sizes_to_run = args.titrate_panels
        log(f"Titrating across panel sizes: {sizes_to_run}")
    else:
        sizes_to_run = [args.panel_size]
        log(f"Running single panel size: {sizes_to_run[0]}")
        
    # --- Execute Pipeline ---
    # Step 1: Base Feature Extraction / Prior Generation (Elastic Net)
    base_pool_idx, base_pool_scores = compute_global_elasticnet(adata, y_full, args)

    # Step 2: Optimization Loop for Panel Sizes
    for actual_size in sizes_to_run:
        log(f"\nProcessing Pipeline Setup for Panel Size Target: {actual_size} | Alpha: {args.alpha}")

        if actual_size > len(base_pool_idx):
            log(f"Skipping size {actual_size}: requested panel size exceeds candidate pool size ({len(base_pool_idx)}).")
            continue
            
        candidate_idx, scaler, ref_model, candidate_priors = localized_subset_refit(
            adata, y_full, base_pool_idx, base_pool_scores, args, actual_size
        )
            
        # Build Support-Calibrated Constraint Setup using the localized L2 layout
        X_scaled = scaler.transform(as_dense_array(adata[:, candidate_idx].X))
        ref_scores = ref_model.decision_function(X_scaled)

        log("Mining non-redundant margin constraints for Gurobi...")
        pair_df = select_mrmr_constraints(y_true=y_full, decision_scores=ref_scores, n_pairs_target=3000)
        
        class_counts = pd.Series(y_full).value_counts().to_dict()
        max_support = max(class_counts.values())
        penalty_map = {cls: max_support / count for cls, count in class_counts.items()}
        pair_df['support_penalty'] = pair_df['true_label'].map(penalty_map) * 5.0
        
        A, bias = build_milp_matrix(adata, candidate_idx, scaler, ref_model, pair_df)
        weights = (pair_df["macro_weight"].values.astype(float) * pair_df["support_penalty"].values.astype(float))

        # Call Gurobi with the new Dynamic Objective arguments
        log(f"Executing Dynamic Prior-Weighted MILP Allocation for Panel Size: {actual_size} (Alpha: {args.alpha})")
        try:
            # Pass args.alpha directly to solve_milp
            selected_local = solve_milp(A, bias, weights, args, actual_size, candidate_priors)
            selected_global_idx = candidate_idx[selected_local]
            
            generate_outputs(
                out_dir, str(actual_size), adata, y_full, 
                selected_global_idx, unique_types, args
            )
        except Exception as e:
            log(f"Solver routine failed for size {actual_size}: {e}")


def cli():
    try:
        main()
    except KeyboardInterrupt:
        print("\n[!] Ctrl+C detected. Assassinating orphaned background workers...")
        
        try:
            from loky import get_reusable_executor
            get_reusable_executor().shutdown(wait=True, kill_workers=True)
        except Exception:
            pass 
            
        print("[!] Cleanup complete. Exiting.")
        sys.exit(1)


if __name__ == "__main__":
    cli()
