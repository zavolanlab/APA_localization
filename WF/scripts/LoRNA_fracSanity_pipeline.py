"""
fracSanity on the U2OS density-gradient fractionation (LoRNA; Villanueva et al., Nat Methods 2024, ENA PRJEB49479):
8 iodixanol-gradient fractions and the whole-cell lysate (dTotal), 3 replicates each, DMSO (control) only; CORALL
total RNA-seq after rRNA depletion, with Lexogen SIRV spike-ins. Processed with the bulk RNA-seq workflow (run
U2OS_subcellular_fractions, NFYA_project.ipynb) and WF/scripts/count_spike_ins.sbatch (SIRV read pairs per library).
Heavy steps for LoRNA_U2OS.fracSanity.ipynb, which reads the outputs. Same design as WF/scripts/APEX_fracSanity_pipeline.py.

The spike-ins give the RNA content of every fraction only together with two per-library numbers that are not in the reads.
Per the published methods, each library was made from 1 ug of a fraction's RNA, and SIRVs were added after rRNA depletion at
"1%" of the depleted RNA (the authors' table: 0.02-0.12 ng, i.e. not the same amount in every library); the authors combine
the SIRV share of the reads with the ng SIRVs added and with the share of the fraction's RNA used for the library
(metadata/LoRNA_UPR_Villanueva2024: sirv_correction.rds of their GitHub repository). Step authors_content turns that table
into F_auth (the reference F of the notebook). The steps below first assumed equal SIRV amounts per library:
K_kr = sum_g (c_gkr / L_g) / SIRV_kr (molecules of human RNA per SIRV molecule, up to a constant) and
F_spike,kr = K_kr / sum_k K_kr, the share of all RNA molecules (of the analysed genes) of replicate r in fraction k.
That is the F of fracSanity, which is otherwise not identified without external information; the notebook tests the
fracSanity estimates of F and alpha_g against it.
That assumption does not hold (K does not differ between fractions; F_spike fits the gene-wise mass balance no better than
F = 1/8); F_spike is kept under its name for the steps that used it. F_auth agrees with the mass balance on the pooled
compartments (step mass_balance_pooled), which needs no spike-ins.

Steps (each writes its own files into --outdir and is skipped when they exist, so nothing is overwritten):
  prepare      sample metadata, gene counts (featureCounts on gencode.v42.basic, unique alignments, read pairs,
               reverse-stranded), gene lengths (union of exons), read QC
  spike_ins    SIRV read pairs per library (STAR on the spike-in sequences), RNA content K and F_spike per fraction
  sanity       C++ Sanity (v_method EAP) on genes, log2 TPM (gene lengths), Total and fractions separately
  fit_cohorts  fracSanity per fraction (average of the 3 replicates) vs the Total average, uniform prior on F
  fit_cohorts_spikeprior  the same with a tight Beta prior on F around F_spike (central 95% within F_spike / 1.25 ..
               F_spike * 1.25): the gene-level part of fracSanity when F is (nearly) known
  fit_replicates  fracSanity per fraction library vs the Total library of the same replicate, uniform prior on F
  mass_balance F from the expression data alone: if the 8 fractions together contain the lysate, then for every gene
               rho_g,Total = sum_k F_k rho_g,k (rho: TPM); non-negative least squares over the well-expressed genes
               (Total >= 10 TPM) of 1 = sum_k F_k rho_g,k / rho_g,Total, on the cohort averages (bootstrap over genes,
               95% interval) and per replicate (fractions vs the Total of the same replicate)
  fit_cohorts_mbprior  fracSanity per fraction with a broad Beta prior on F from the mass balance (centre max(F_mb, 0.005),
               95% interval a factor 4 either side): a way to constrain F without spike-ins
  fit_cohorts_lowpi, fit_cohorts_spikeprior_lowpi  as fit_cohorts and fit_cohorts_spikeprior with the outlier fraction of the
               error model bounded at LOW_PI_MAX instead of 0.5: genes flagged as outliers keep almost no information on
               alpha_g (error variance V_g + tau^2), which in the APEX-seq fits discounted the most strongly localised genes
  pooled       the fractions of one compartment (POOLS) merged per replicate by binomial thinning of the reads in proportion
               to F_spike (each pool then holds all molecules of the compartment, so a gene confined to it has alpha near 1);
               Sanity on the pooled libraries, fracSanity against Total (uniform prior; default and low pi)
  mass_balance_pooled  the mass balance on the 3 pooled compartments: rho_Total = sum_pool F_pool rho_pool; better determined
               than the 8 collinear fractions (whose F_mb is unstable between replicates)
  simulate     fracSanity on simulated D_g with F = F_spike and the real errors, for the 8 fractions and the 3 pools:
               alpha_g = the spike-in alpha^T (real shape) or drawn from the Beta fitted to alpha^T (the model holds)
  authors_content  F_auth per fraction library from the authors' table (inputs/F_authors.tsv)
  fit_cohorts_authprior, fit_cohorts_authprior_lowpi  as fit_cohorts_spikeprior(_lowpi), prior around the mean F_auth
  pooled_auth, mass_balance_pooled_auth  as pooled, mass_balance_pooled, the fractions merged in proportion to F_auth

Usage (sanity_utils environment):
    python WF/scripts/LoRNA_fracSanity_pipeline.py --outdir <dir> --threads 24 [--steps prepare,...]
"""

import argparse
import os
import re
import sys
import time

import numpy as np
import pandas as pd
from dotenv import load_dotenv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ONT_fracSanity_pipeline import dump, load, log, parse_gtf, run_jobs, run_one_sanity  # noqa: E402

PROJECT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RUN = "U2OS_subcellular_fractions"
# fraction identities as in NFYA_project.ipynb (protein markers of the paper)
FRACTION_LABELS = {1: "ER+mito (1)", 2: "ER+mito (2)", 3: "unclear (3)", 4: "nucleolus (4)", 5: "unclear (5)",
                   6: "cytosol (6)", 7: "cytosol (7)", 8: "cytosol (8)"}
SPIKE_PRIOR_FACTOR = 1.25
MB_MIN_TPM = 10              # genes used for the mass balance: Total average >= 10 TPM
MB_BOOTSTRAP = 200
MB_FLOOR, MB_PRIOR_FACTOR = 0.005, 4.0   # prior on F from the mass balance (step fit_cohorts_mbprior)
LOW_PI_MAX = 0.02            # bound on the outlier fraction in the *_lowpi fits
# compartments from the spike-in alpha of marker RNAs (mito-encoded and ER mRNAs in 1-3, snoRNAs/MALAT1/RMRP in 4-5)
POOLS = {"Membranes": ["F1", "F2", "F3"], "Nuclear": ["F4", "F5"], "Cytosol": ["F6", "F7", "F8"]}


def paths():
    load_dotenv(os.path.join(PROJECT, "APA_localization.scicore.env"), override=True)
    runs = os.getenv("WF_RUNS_DIR")
    return {"wf_out": os.path.join(runs, RUN, "output"), "spike_ins": os.path.join(runs, f"{RUN}_spike_ins"),
            "samples": os.path.join(os.getenv("METADATA_DIR"), f"start_samples.{RUN}.tsv"),
            "authors": os.path.join(os.getenv("METADATA_DIR"), "LoRNA_UPR_Villanueva2024", "sirv_correction.tsv"),
            "gtf": os.getenv("HUMAN_BASIC_ANNOTATION_FILE")}


def safe(name):
    return re.sub(r"[^A-Za-z0-9]+", "_", name)


def sample_table(cfg):
    rows = []
    for s in pd.read_csv(cfg["samples"], sep="\t")["sample"]:
        parts = s.split("-")
        if parts[0] == "dTotal" and parts[1] == "DMSO":
            rows.append({"library": s, "fraction": "Total", "fraction_number": 0, "replicate": int(parts[2]), "cohort": "Total"})
        elif parts[0] == "dLoRNA" and parts[1] == "DMSO":
            k = int(parts[3])
            rows.append({"library": s, "fraction": f"F{k}", "fraction_number": k, "replicate": int(parts[2]), "cohort": f"F{k}"})
    meta = pd.DataFrame(rows).sort_values(["fraction_number", "replicate"]).reset_index(drop=True)
    meta["label"] = meta["fraction_number"].map(FRACTION_LABELS).fillna("Total")
    meta["sample"] = meta["cohort"] + " R" + meta["replicate"].astype(str)
    assert len(meta) == 27 and meta["sample"].is_unique, meta
    return meta


def read_featurecounts(path, lengths=False):
    fc = pd.read_csv(path, sep="\t", comment="#", index_col=0)
    return fc["Length"] if lengths else fc.iloc[:, -1]


def step_prepare(out, cfg):
    d = os.path.join(out, "inputs")
    if os.path.exists(os.path.join(d, "done")):
        log("prepare: done before, skipped")
        return
    os.makedirs(d, exist_ok=True)
    meta = sample_table(cfg)
    genes, _ = parse_gtf(cfg["gtf"])
    genes.to_csv(os.path.join(d, "genes.tsv"), sep="\t", index_label="gene_id")
    fc_dir = os.path.join(cfg["wf_out"], "gene_expression_quantification", "human", "FeatureCounts_standard")
    counts = pd.DataFrame({s: read_featurecounts(os.path.join(fc_dir, f"{l}.standard.txt"))
                           for s, l in zip(meta["sample"], meta["library"])}).reindex(genes.index).fillna(0).round().astype(np.int64)
    counts.to_csv(os.path.join(d, "gene_counts.tsv.gz"), sep="\t", index_label="gene_id")
    read_featurecounts(os.path.join(fc_dir, f"{meta['library'].iloc[0]}.standard.txt"), lengths=True).rename("length").to_csv(
        os.path.join(d, "gene_lengths.tsv"), sep="\t", index_label="gene_id")

    qc = []
    for s, l in zip(meta["sample"], meta["library"]):
        with open(os.path.join(cfg["wf_out"], "mapping_stats", f"{l}.mapping_stats.txt")) as fh:
            lines = [x.strip() for x in fh]
        ms = {lines[i].lstrip(">").split(";")[0].strip(): float(lines[i + 1]) for i in range(0, len(lines) - 1, 2) if lines[i].startswith(">")}
        summ = pd.read_csv(os.path.join(fc_dir, f"{l}.standard.txt.summary"), sep="\t", index_col=0).iloc[:, 0]
        qc.append({"sample": s, "input_reads": ms.get("total number of input reads"), "mapped_reads": ms.get("total number of mapped reads"),
                   "uniquely_mapped_reads": ms.get("number of UM reads"), "assigned_read_pairs": summ.get("Assigned"),
                   "no_feature_read_pairs": summ.get("Unassigned_NoFeatures"), "ambiguous_read_pairs": summ.get("Unassigned_Ambiguity")})
    meta = meta.merge(pd.DataFrame(qc), on="sample")
    meta["use"] = True
    meta.to_csv(os.path.join(d, "metadata.tsv"), sep="\t", index=False)
    open(os.path.join(d, "done"), "w").close()
    log(f"prepare: {len(meta)} libraries, {len(counts):,} genes, {int(counts.sum().sum()):,} assigned read pairs")


def step_spike_ins(out, cfg):
    f = os.path.join(out, "inputs", "spike_ins.tsv")
    if os.path.exists(f):
        log("spike_ins: done before, skipped")
        return
    meta = pd.read_csv(os.path.join(out, "inputs", "metadata.tsv"), sep="\t")
    counts = pd.read_csv(os.path.join(out, "inputs", "gene_counts.tsv.gz"), sep="\t", index_col=0)
    lengths = pd.read_csv(os.path.join(out, "inputs", "gene_lengths.tsv"), sep="\t", index_col=0)["length"].reindex(counts.index)
    rows, per_gene = [], {}
    for s, l in zip(meta["sample"], meta["library"]):
        sd = os.path.join(cfg["spike_ins"], l)
        log_final = os.path.join(sd, "Log.final.out")
        if not os.path.exists(log_final):
            raise FileNotFoundError(f"{log_final}: run WF/scripts/count_spike_ins.sbatch first")
        lf = {}
        for line in open(log_final):
            if "|" in line:
                k, v = line.split("|")
                lf[k.strip()] = v.strip()
        rpg = pd.read_csv(os.path.join(sd, "ReadsPerGene.out.tab"), sep="\t", header=None, index_col=0)
        per_gene[s] = rpg.iloc[4:, 0]           # unstranded column: SIRV genes have antisense isoforms
        rows.append({"sample": s, "spike_input_read_pairs": int(lf["Number of input reads"]),
                     "spike_unique_pairs": int(lf["Uniquely mapped reads number"]),
                     "spike_multi_pairs": int(lf["Number of reads mapped to multiple loci"])})
    sp = pd.DataFrame(rows).set_index("sample")
    pg = pd.DataFrame(per_gene).fillna(0).astype(np.int64)
    sp["SIRV_pairs"] = sp["spike_unique_pairs"] + sp["spike_multi_pairs"]
    sp["SIRV_pairs_on_SIRV_genes"] = pg.loc[pg.index.str.startswith("SIRV")].sum()
    sp["ERCC_pairs"] = pg.loc[pg.index.str.startswith("ERCC")].sum()
    sp["SIRV_share"] = sp["SIRV_pairs"] / sp["spike_input_read_pairs"]
    # human RNA molecules (up to a constant): length-normalised assigned read pairs of all genes with a length
    sp["human_molecule_index"] = (counts / lengths.to_numpy()[:, None]).sum().reindex(sp.index)
    sp["K"] = sp["human_molecule_index"] / sp["SIRV_pairs"]
    sp = meta.set_index("sample")[["library", "cohort", "fraction_number", "replicate"]].join(sp)
    fr = sp["fraction_number"] > 0
    sp["F_spike"] = np.nan
    sp.loc[fr, "F_spike"] = sp.loc[fr, "K"] / sp.loc[fr].groupby("replicate")["K"].transform("sum")
    sp.to_csv(f, sep="\t", index_label="sample")
    pg.to_csv(os.path.join(out, "inputs", "spike_ins.per_gene.tsv"), sep="\t", index_label="spike_in")
    log("spike_ins: SIRV share " + ", ".join(f"{s} {v:.4f}" for s, v in sp["SIRV_share"].items()))
    log("spike_ins: mean F_spike per fraction " + ", ".join(f"{c} {v:.3f}" for c, v in sp.loc[fr].groupby("cohort")["F_spike"].mean().items()))


def step_sanity(out, threads):
    d = os.path.join(out, "inputs")
    meta = pd.read_csv(os.path.join(d, "metadata.tsv"), sep="\t")
    meta = meta.loc[meta["use"]].reset_index(drop=True)
    counts = pd.read_csv(os.path.join(d, "gene_counts.tsv.gz"), sep="\t", index_col=0)
    lengths = pd.read_csv(os.path.join(d, "gene_lengths.tsv"), sep="\t", index_col=0)["length"]
    for part, m in (("Total", meta.loc[meta["cohort"] == "Total"]), ("fractions", meta.loc[meta["cohort"] != "Total"])):
        od = os.path.join(out, "sanity", "genes", part)
        os.makedirs(od, exist_ok=True)
        run_one_sanity(counts, m, od, threads, gene_lengths=lengths)


def _cohort_jobs(out, sub, prior=None, pi_max=None):
    import sanity_utils as su
    tot = load(os.path.join(out, "sanity", "genes", "Total", "sanity.pkl"))
    frac = load(os.path.join(out, "sanity", "genes", "fractions", "sanity.pkl"))
    feats = tot["averages"].means.index.intersection(frac["averages"].means.index)
    od = os.path.join(out, "fracSanity", "genes", sub)
    os.makedirs(od, exist_ok=True)
    jobs = []
    for cohort in frac["averages"].means.columns:
        inp = su.prepare_frac_sanity_inputs(tot["averages"], frac["averages"], "Total", cohort, genes=feats)
        info = {"cohort": cohort, "level": "genes"}
        if prior is not None:
            info["prior_F"] = prior[cohort]
        if pi_max is not None:
            info["pi_max"] = pi_max
        jobs.append((os.path.join(od, safe(cohort) + ".pkl"), inp["ltq_total"], inp["var_total"], inp["ltq_pd"], inp["var_pd"],
                     inp["gene_ids"], info))
    return jobs


def step_fit_cohorts(out, threads):
    run_jobs(_cohort_jobs(out, "cohorts"), threads)


def _spike_prior(out):
    import sanity_utils as su
    sp = pd.read_csv(os.path.join(out, "inputs", "spike_ins.tsv"), sep="\t", index_col=0)
    f_mean = sp.loc[sp["fraction_number"] > 0].groupby("cohort")["F_spike"].mean()
    prior = {c: su.beta_prior_from_interval(f / SPIKE_PRIOR_FACTOR, min(f * SPIKE_PRIOR_FACTOR, 0.99)) for c, f in f_mean.items()}
    log("spike prior: " + ", ".join(f"{c} F {f_mean[c]:.3f} Beta({a:.1f}, {b:.1f})" for c, (a, b) in prior.items()))
    return prior


def step_fit_cohorts_spikeprior(out, threads):
    run_jobs(_cohort_jobs(out, "cohorts_spikeprior", _spike_prior(out)), threads)


def step_fit_cohorts_lowpi(out, threads):
    run_jobs(_cohort_jobs(out, "cohorts_lowpi", pi_max=LOW_PI_MAX), threads)


def step_fit_cohorts_spikeprior_lowpi(out, threads):
    run_jobs(_cohort_jobs(out, "cohorts_spikeprior_lowpi", _spike_prior(out), pi_max=LOW_PI_MAX), threads)


def step_authors_content(out, cfg):
    # RNA content per fraction from the authors' processed data (metadata/LoRNA_UPR_Villanueva2024/README.md): SIRV share x
    # ng SIRVs added to each library x total RNA of the fraction; F_auth = ng_non_rRNA / sum over the 8 fractions
    f = os.path.join(out, "inputs", "F_authors.tsv")
    if os.path.exists(f):
        log("authors_content: done before, skipped")
        return
    a = pd.read_csv(cfg["authors"], sep="\t")
    a = a.loc[(a["LoRNA"] == "dLoRNA") & (a["condition"] == "DMSO")].copy()
    a["cohort"] = "F" + a["fraction"].astype(str)
    a["sample"] = a["cohort"] + " R" + a["replicate"].astype(str)
    a["F_auth"] = a["ng_non_rRNA"] / a.groupby("replicate")["ng_non_rRNA"].transform("sum")
    assert len(a) == 24, a
    a.to_csv(f, sep="\t", index=False)
    log("authors_content: mean F_auth per fraction " + ", ".join(f"{c} {v:.3f}" for c, v in a.groupby("cohort")["F_auth"].mean().items()))


def _auth_prior(out):
    import sanity_utils as su
    a = pd.read_csv(os.path.join(out, "inputs", "F_authors.tsv"), sep="\t")
    f_mean = a.groupby("cohort")["F_auth"].mean()
    prior = {c: su.beta_prior_from_interval(f / SPIKE_PRIOR_FACTOR, min(f * SPIKE_PRIOR_FACTOR, 0.99)) for c, f in f_mean.items()}
    log("authors' prior: " + ", ".join(f"{c} F {f_mean[c]:.3f} Beta({a:.1f}, {b:.1f})" for c, (a, b) in prior.items()))
    return prior


def step_fit_cohorts_authprior(out, threads):
    run_jobs(_cohort_jobs(out, "cohorts_authprior", _auth_prior(out)), threads)


def step_fit_cohorts_authprior_lowpi(out, threads):
    run_jobs(_cohort_jobs(out, "cohorts_authprior_lowpi", _auth_prior(out), pi_max=LOW_PI_MAX), threads)


def step_fit_replicates(out, threads):
    # each fraction library against the Total library of the same replicate (same lysate)
    tot = load(os.path.join(out, "sanity", "genes", "Total", "sanity.pkl"))
    frac = load(os.path.join(out, "sanity", "genes", "fractions", "sanity.pkl"))
    genes = tot["result"].log2_expression.index.intersection(frac["result"].log2_expression.index)
    tr, rp = tot["result"], frac["result"]
    tot_of = tot["metadata"].set_index("replicate")["sample"]
    od = os.path.join(out, "fracSanity", "genes", "libraries")
    os.makedirs(od, exist_ok=True)
    jobs = []
    for _, r in frac["metadata"].iterrows():
        t = tot_of[r["replicate"]]
        jobs.append((os.path.join(od, safe(r["sample"]) + ".pkl"),
                     tr.log2_expression.loc[genes, t].to_numpy(), tr.log2_sd.loc[genes, t].to_numpy() ** 2,
                     rp.log2_expression.loc[genes, r["sample"]].to_numpy(), rp.log2_sd.loc[genes, r["sample"]].to_numpy() ** 2,
                     genes, {"cohort": r["cohort"], "sample": r["sample"], "replicate": int(r["replicate"]), "total": t, "level": "genes"}))
    run_jobs(jobs, threads)


def _nnls_F(ratio):
    from scipy.optimize import nnls
    ok = np.isfinite(ratio).all(axis=1)
    return nnls(ratio[ok], np.ones(ok.sum()))[0]


def step_mass_balance(out):
    f = os.path.join(out, "inputs", "mass_balance.tsv")
    if os.path.exists(f):
        log("mass_balance: done before, skipped")
        return
    tot = load(os.path.join(out, "sanity", "genes", "Total", "sanity.pkl"))
    frac = load(os.path.join(out, "sanity", "genes", "fractions", "sanity.pkl"))
    cohorts = sorted(frac["averages"].means.columns, key=lambda c: int(c[1:]))
    mt, mf = tot["averages"].means["Total"], frac["averages"].means[cohorts]
    genes = mt.index.intersection(mf.index)
    genes = genes[(mt.loc[genes] >= np.log2(MB_MIN_TPM)).to_numpy()]
    ratio = 2.0 ** (mf.loc[genes].to_numpy() - mt.loc[genes].to_numpy()[:, None])
    F = _nnls_F(ratio)
    rng = np.random.default_rng(0)
    boot = np.array([_nnls_F(ratio[rng.integers(0, len(genes), len(genes))]) for _ in range(MB_BOOTSTRAP)])
    rows = [{"cohort": c, "replicate": "average", "F_mb": F[i], "F_mb_CI_lower": np.quantile(boot[:, i], 0.025),
             "F_mb_CI_upper": np.quantile(boot[:, i], 0.975), "n_genes": len(genes)} for i, c in enumerate(cohorts)]
    # per replicate: fraction libraries vs the Total library of the same replicate
    lt, lf = tot["result"].log2_expression, frac["result"].log2_expression
    fm = frac["metadata"]
    for _, t in tot["metadata"].iterrows():
        libs = [fm.loc[(fm["cohort"] == c) & (fm["replicate"] == t["replicate"]), "sample"].iloc[0] for c in cohorts]
        g = lt.index.intersection(lf.index)
        g = g[(lt.loc[g, t["sample"]] >= np.log2(MB_MIN_TPM)).to_numpy()]
        Fr = _nnls_F(2.0 ** (lf.loc[g, libs].to_numpy() - lt.loc[g, t["sample"]].to_numpy()[:, None]))
        rows += [{"cohort": c, "replicate": int(t["replicate"]), "F_mb": Fr[i], "n_genes": len(g)} for i, c in enumerate(cohorts)]
    mb = pd.DataFrame(rows)
    mb.to_csv(f, sep="\t", index=False)
    a = mb.loc[mb["replicate"] == "average"]
    log(f"mass_balance ({len(genes):,} genes): " + ", ".join(f"{c} {v:.3f}" for c, v in zip(a["cohort"], a["F_mb"])) + f"; sum {a['F_mb'].sum():.3f}")


def step_fit_cohorts_mbprior(out, threads):
    import sanity_utils as su
    mb = pd.read_csv(os.path.join(out, "inputs", "mass_balance.tsv"), sep="\t")
    mb = mb.loc[mb["replicate"] == "average"].set_index("cohort")
    # broad: the bootstrap interval misses the collinearity of similar fractions (an estimate at the bound 0 has a
    # bootstrap interval of width 0); centre max(F_mb, MB_FLOOR), 95% within a factor MB_PRIOR_FACTOR either side
    prior = {}
    for c, r in mb.iterrows():
        centre = max(r["F_mb"], MB_FLOOR)
        prior[c] = su.beta_prior_from_interval(centre / MB_PRIOR_FACTOR, min(centre * MB_PRIOR_FACTOR, 0.99))
    log("mass-balance prior: " + ", ".join(f"{c} Beta({a:.1f}, {b:.1f})" for c, (a, b) in prior.items()))
    run_jobs(_cohort_jobs(out, "cohorts_mbprior", prior), threads)


def _pool_name(c):
    return next(p for p, members in POOLS.items() if c in members)


def step_pooled(out, threads, tag="", weight="F_spike"):
    # per replicate, the fraction libraries of one compartment merged in proportion to their RNA content: reads of library k
    # kept with probability proportional to F_k / (human molecule index)_k (binomial thinning), so that the pooled library
    # samples the molecules of the compartment as the cell holds them. F_k: weight "F_spike" (spike_ins.tsv, equal-SIRV
    # assumption; tag "") or "F_auth" (F_authors.tsv, the authors' RNA content; tag "_auth")
    od = os.path.join(out, "sanity", "genes", f"pooled{tag}")
    counts_f = os.path.join(out, "inputs", f"pooled{tag}_counts.tsv.gz")
    meta_f = os.path.join(out, "inputs", f"pooled{tag}_metadata.tsv")
    if not os.path.exists(counts_f):
        counts = pd.read_csv(os.path.join(out, "inputs", "gene_counts.tsv.gz"), sep="\t", index_col=0)
        sp = pd.read_csv(os.path.join(out, "inputs", "spike_ins.tsv"), sep="\t", index_col=0)
        sp = sp.loc[sp["fraction_number"] > 0]
        if weight == "F_auth":
            fa = pd.read_csv(os.path.join(out, "inputs", "F_authors.tsv"), sep="\t").set_index("sample")["F_auth"]
            sp = sp.assign(F_auth=fa.reindex(sp.index))
        rng = np.random.default_rng(1)
        pooled, rows = {}, []
        for pool, members in POOLS.items():
            for rep in (1, 2, 3):
                s = sp.loc[sp["cohort"].isin(members) & (sp["replicate"] == rep)]
                keep = s[weight] / s["human_molecule_index"]
                keep = keep / keep.max()
                name = f"{pool} R{rep}"
                pooled[name] = sum(rng.binomial(counts[lib].to_numpy(), p) for lib, p in keep.items())
                rows.append({"sample": name, "cohort": pool, "replicate": rep, weight: s[weight].sum(),
                             "libraries": ",".join(keep.index), "keep_probability": ",".join(f"{p:.3f}" for p in keep)})
        pd.DataFrame(pooled, index=counts.index).to_csv(counts_f, sep="\t")
        pd.DataFrame(rows).to_csv(meta_f, sep="\t", index=False)
        log(f"pooled{tag} counts: " + ", ".join(f"{r['sample']} keep {r['keep_probability']}" for r in rows))
    counts = pd.read_csv(counts_f, sep="\t", index_col=0)
    meta = pd.read_csv(meta_f, sep="\t")
    lengths = pd.read_csv(os.path.join(out, "inputs", "gene_lengths.tsv"), sep="\t", index_col=0)["length"]
    os.makedirs(od, exist_ok=True)
    run_one_sanity(counts, meta, od, threads, gene_lengths=lengths)
    tot = load(os.path.join(out, "sanity", "genes", "Total", "sanity.pkl"))
    pl = load(os.path.join(od, "sanity.pkl"))
    import sanity_utils as su
    feats = tot["averages"].means.index.intersection(pl["averages"].means.index)
    jobs = []
    for sub, pi_max in ((f"pooled{tag}", None), (f"pooled{tag}_lowpi", LOW_PI_MAX)):
        fd = os.path.join(out, "fracSanity", "genes", sub)
        os.makedirs(fd, exist_ok=True)
        for pool in POOLS:
            inp = su.prepare_frac_sanity_inputs(tot["averages"], pl["averages"], "Total", pool, genes=feats)
            info = {"cohort": pool, "level": "genes"}
            if pi_max is not None:
                info["pi_max"] = pi_max
            jobs.append((os.path.join(fd, safe(pool) + ".pkl"), inp["ltq_total"], inp["var_total"], inp["ltq_pd"], inp["var_pd"],
                         inp["gene_ids"], info))
    run_jobs(jobs, threads)


def step_mass_balance_pooled(out, tag=""):
    # mass balance on the 3 pooled compartments (3 nearly independent profiles instead of 8 collinear ones), as step_mass_balance
    f = os.path.join(out, "inputs", f"mass_balance_pooled{tag}.tsv")
    if os.path.exists(f):
        log(f"mass_balance_pooled{tag}: done before, skipped")
        return
    tot = load(os.path.join(out, "sanity", "genes", "Total", "sanity.pkl"))
    pl = load(os.path.join(out, "sanity", "genes", f"pooled{tag}", "sanity.pkl"))
    pools = list(POOLS)
    mt, mp = tot["averages"].means["Total"], pl["averages"].means[pools]
    genes = mt.index.intersection(mp.index)
    genes = genes[(mt.loc[genes] >= np.log2(MB_MIN_TPM)).to_numpy()]
    ratio = 2.0 ** (mp.loc[genes].to_numpy() - mt.loc[genes].to_numpy()[:, None])
    F = _nnls_F(ratio)
    rng = np.random.default_rng(0)
    boot = np.array([_nnls_F(ratio[rng.integers(0, len(genes), len(genes))]) for _ in range(MB_BOOTSTRAP)])
    rows = [{"pool": p, "replicate": "average", "F_mb": F[i], "F_mb_CI_lower": np.quantile(boot[:, i], 0.025),
             "F_mb_CI_upper": np.quantile(boot[:, i], 0.975), "n_genes": len(genes)} for i, p in enumerate(pools)]
    lt, lp = tot["result"].log2_expression, pl["result"].log2_expression
    pm = pl["metadata"]
    for _, t in tot["metadata"].iterrows():
        libs = [pm.loc[(pm["cohort"] == p) & (pm["replicate"] == t["replicate"]), "sample"].iloc[0] for p in pools]
        g = lt.index.intersection(lp.index)
        g = g[(lt.loc[g, t["sample"]] >= np.log2(MB_MIN_TPM)).to_numpy()]
        r = 2.0 ** (lp.loc[g, libs].to_numpy() - lt.loc[g, t["sample"]].to_numpy()[:, None])
        Fr = _nnls_F(r)
        bt = np.array([_nnls_F(r[rng.integers(0, len(g), len(g))]) for _ in range(MB_BOOTSTRAP)])
        rows += [{"pool": p, "replicate": int(t["replicate"]), "F_mb": Fr[i], "F_mb_CI_lower": np.quantile(bt[:, i], 0.025),
                  "F_mb_CI_upper": np.quantile(bt[:, i], 0.975), "n_genes": len(g)} for i, p in enumerate(pools)]
    mb = pd.DataFrame(rows)
    mb.to_csv(f, sep="\t", index=False)
    log(f"mass_balance_pooled{tag}: " + ", ".join(f"{p} {v:.3f}" for p, v in zip(pools, F)) + f"; sum {F.sum():.3f}")


def step_simulate(out, threads):
    # D_g = log2 alpha_g - log2 F_spike + eps with the real errors (eps ~ N(0, V_total + V_fraction)) on the genes with
    # Total >= MB_MIN_TPM, for the 8 fractions and the 3 pooled compartments; alpha_g either the spike-in alpha^T of the
    # gene (the real shape of the distribution) or drawn from the Beta fitted to the alpha^T of all genes (the model
    # holds exactly). Uniform prior on F, default error model.
    from scipy.stats import beta as beta_dist
    tot = load(os.path.join(out, "sanity", "genes", "Total", "sanity.pkl"))
    frac = load(os.path.join(out, "sanity", "genes", "fractions", "sanity.pkl"))
    sp = pd.read_csv(os.path.join(out, "inputs", "spike_ins.tsv"), sep="\t", index_col=0)
    F_spike = sp.loc[sp["fraction_number"] > 0].groupby("cohort")["F_spike"].mean()
    mT, eT = tot["averages"].means["Total"], tot["averages"].errors["Total"]
    M, E = frac["averages"].means, frac["averages"].errors
    genes = mT.index.intersection(M.index)
    genes = genes[(mT.loc[genes] >= np.log2(MB_MIN_TPM)).to_numpy()]
    alpha_T = 2.0 ** M.loc[genes].sub(mT.loc[genes], axis=0) * F_spike
    cases = {c: (alpha_T[c], F_spike[c], E.loc[genes, c] ** 2) for c in alpha_T.columns}
    for pool, members in POOLS.items():
        # error of the pooled fraction: that of the member with the smallest error (the pooled library is deeper)
        cases[pool] = (alpha_T[members].sum(axis=1), F_spike[members].sum(), (E.loc[genes, members] ** 2).min(axis=1))
    rng = np.random.default_rng(2)
    rows, jobs = [], {"sim_alpha_T": [], "sim_beta": []}
    for c, (a_true, F, var_p) in cases.items():
        a_true = a_true.clip(1e-6, 0.999)
        a_b, b_b, _, _ = beta_dist.fit(a_true.to_numpy(), floc=0, fscale=1)
        rows.append({"cohort": c, "F": F, "beta_a": a_b, "beta_b": b_b, "n_genes": len(genes)})
        var_t = (eT.loc[genes] ** 2).to_numpy()
        for sub, alpha in (("sim_alpha_T", a_true.to_numpy()), ("sim_beta", rng.beta(a_b, b_b, len(genes)))):
            ltq_p = mT.loc[genes].to_numpy() + np.log2(alpha / F) + rng.normal(0, np.sqrt(var_t + var_p.to_numpy()))
            fd = os.path.join(out, "fracSanity", "genes", sub)
            os.makedirs(fd, exist_ok=True)
            pd.DataFrame({"gene_id": genes, "alpha_sim": alpha}).to_csv(os.path.join(fd, safe(c) + ".alpha_sim.tsv.gz"), sep="\t", index=False)
            jobs[sub].append((os.path.join(fd, safe(c) + ".pkl"), mT.loc[genes].to_numpy(), var_t, ltq_p, var_p.to_numpy(), genes,
                              {"cohort": c, "level": "genes", "F_true": F, "simulation": sub, "beta_a": a_b, "beta_b": b_b}))
    pd.DataFrame(rows).to_csv(os.path.join(out, "fracSanity", "genes", "simulation_setup.tsv"), sep="\t", index=False)
    log("simulation: " + ", ".join(f"{r['cohort']} F {r['F']:.3f} Beta({r['beta_a']:.2f}, {r['beta_b']:.2f})" for r in rows))
    run_jobs(jobs["sim_alpha_T"] + jobs["sim_beta"], threads)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--steps", default="prepare,spike_ins,sanity,fit_cohorts,fit_replicates,mass_balance,fit_cohorts_mbprior,"
                                       "fit_cohorts_spikeprior,fit_cohorts_lowpi,fit_cohorts_spikeprior_lowpi,pooled,mass_balance_pooled,simulate,"
                                       "authors_content,fit_cohorts_authprior,fit_cohorts_authprior_lowpi,pooled_auth,"
                                       "mass_balance_pooled_auth")
    a = ap.parse_args()
    os.environ.setdefault("OMP_NUM_THREADS", "6")
    cfg = paths()
    os.makedirs(a.outdir, exist_ok=True)
    steps = {"prepare": lambda: step_prepare(a.outdir, cfg), "spike_ins": lambda: step_spike_ins(a.outdir, cfg),
             "sanity": lambda: step_sanity(a.outdir, a.threads), "fit_cohorts": lambda: step_fit_cohorts(a.outdir, a.threads),
             "fit_cohorts_spikeprior": lambda: step_fit_cohorts_spikeprior(a.outdir, a.threads),
             "fit_replicates": lambda: step_fit_replicates(a.outdir, a.threads), "mass_balance": lambda: step_mass_balance(a.outdir),
             "fit_cohorts_mbprior": lambda: step_fit_cohorts_mbprior(a.outdir, a.threads),
             "fit_cohorts_lowpi": lambda: step_fit_cohorts_lowpi(a.outdir, a.threads),
             "fit_cohorts_spikeprior_lowpi": lambda: step_fit_cohorts_spikeprior_lowpi(a.outdir, a.threads),
             "pooled": lambda: step_pooled(a.outdir, a.threads), "mass_balance_pooled": lambda: step_mass_balance_pooled(a.outdir),
             "simulate": lambda: step_simulate(a.outdir, a.threads),
             "authors_content": lambda: step_authors_content(a.outdir, cfg),
             "fit_cohorts_authprior": lambda: step_fit_cohorts_authprior(a.outdir, a.threads),
             "fit_cohorts_authprior_lowpi": lambda: step_fit_cohorts_authprior_lowpi(a.outdir, a.threads),
             "pooled_auth": lambda: step_pooled(a.outdir, a.threads, "_auth", "F_auth"),
             "mass_balance_pooled_auth": lambda: step_mass_balance_pooled(a.outdir, "_auth")}
    for step in a.steps.split(","):
        t0 = time.time()
        log(f"=== {step}")
        steps[step]()
        log(f"=== {step} finished in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
