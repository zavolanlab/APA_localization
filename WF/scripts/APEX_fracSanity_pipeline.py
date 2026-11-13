"""
fracSanity, tandem poly(A) site (PAS) usage and APA of the APEX-seq data of Fazal et al. 2019 (GSE116008, HEK293T,
TruSeq stranded poly(A)+, 2 x 75 nt), processed with the bulk RNA-seq workflow (WF/Snakefile-prepare-faster,
-quantification-faster, -PAQR-quantify; run APEX_seq_Fazal2019). Heavy steps for APEX_seq.fracSanity.ipynb and
APEX_seq_vs_Control.ipynb, which read the outputs. Same design as WF/scripts/ONT_fracSanity_pipeline.py.

Reference ("Total"): the 18 unlabelled controls (no H2O2, otherwise processed like the targets, incl. the streptavidin
enrichment), as in the 18-control DESeq2 analysis of the paper. Fractions: the labelled targets of each APEX construct
(2 biological replicates; the technical replicate of MITO T2 is added to MITO T2), and the nuclear and cytosolic
fractionation libraries (2 each; not streptavidin-enriched, compared with the same reference).

Steps (each writes its own files into --outdir and is skipped when they exist, so nothing is overwritten):
  prepare      sample metadata, gene counts (featureCounts on gencode.v42.basic, unique alignments, read pairs) and gene
               lengths (featureCounts: length of the union of the gene's exons),
               tandem PAS "counts" (PAQR mean coverage per PAS, rounded, as in NFYA_project.ipynb; terminal exons with
               a missing value (-1) in any library are left out), PAS annotation (3' UTR length)
  complexity   read pairs and distinct fragment coordinates per library (uniquely mapped BAM, first mates)
  sanity       C++ Sanity (v_method EAP): genes and PAS, Total and fractions separately and jointly (PAS). Genes are
               length-corrected (log2 TPM): reads come from along the whole transcript, so TPM here corresponds to CPM of
               the ONT data (one read per molecule). PAS are not: PAQR mean coverage is already per nucleotide.
  fit_cohorts  fracSanity per cohort (fraction average vs Total average), genes and PAS
  fit_replicates  fracSanity per fraction library vs the Total average (genes): input of the gene-set analysis
  usage        PAS relative usage (each PAS vs the other PAS of its terminal exon), 3' UTR length index, most proximal
               vs distal PAS, per cohort vs Total

Usage (sanity_utils environment):
    python WF/scripts/APEX_fracSanity_pipeline.py --outdir <dir> --threads 24 [--steps prepare,...]
"""

import argparse
import glob
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from dotenv import load_dotenv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ONT_fracSanity_pipeline import PI_MAX, dump, load, log, parse_gtf, run_jobs, run_one_sanity  # noqa: E402

PROJECT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RUN = "APEX_seq_Fazal2019"
CONSTRUCTS = {"NLS": "nucleus", "NIK": "nucleolus", "LMA": "nuclear lamina", "NucPore": "nuclear pore", "NES": "cytosol",
              "ERM": "ER membrane (cytosolic face)", "KDEL": "ER lumen", "OMM": "outer mitochondrial membrane",
              "MITO": "mitochondrial matrix"}
MIN_PAS_COUNT = 5            # PAS "reads" (rounded mean coverage) summed over the analysed libraries
SAMTOOLS = os.path.expanduser("~/miniforge3_Dec2024/envs/samtools/bin/samtools")
LEVEL_PAS = "pas"


def paths():
    load_dotenv(os.path.join(PROJECT, "APA_localization.scicore.env"), override=True)
    out = os.path.join(os.getenv("WF_RUNS_DIR"), RUN, "output")
    return {"wf_out": out, "samples": os.path.join(os.getenv("METADATA_DIR"), f"{RUN}.sample_annotation.tsv"),
            "gtf": os.getenv("HUMAN_BASIC_ANNOTATION_FILE")}


def sample_table(cfg):
    ann = pd.read_csv(cfg["samples"], sep="\t")
    rows = []
    for _, r in ann.iterrows():
        s = r["sample"]
        if s.startswith("APEX-"):
            _, construct, kind = s.split("-")
            rep = kind[:2]
            if kind.startswith("C"):
                rows.append({"library": s, "sample": f"Control {construct} {rep}", "batch": "Controls", "fraction": "Total",
                             "construct": construct, "replicate": f"{construct} {rep}", "kind": "control"})
            else:
                rows.append({"library": s, "sample": f"APEX: {construct} {rep}", "batch": "APEX", "fraction": construct,
                             "construct": construct, "replicate": rep, "kind": "target"})
        else:
            _, frac, rep = s.split("-")
            rows.append({"library": s, "sample": f"Fractionation: {frac} {rep}", "batch": "Fractionation", "fraction": frac,
                         "construct": "", "replicate": rep, "kind": "fractionation"})
    lib = pd.DataFrame(rows).merge(ann[["sample", "gsm", "srr", "title"]].rename(columns={"sample": "library"}), on="library")
    lib["cohort"] = np.where(lib["fraction"] == "Total", "Total", lib["batch"] + ": " + lib["fraction"])
    return lib


def read_featurecounts(path, lengths=False):
    fc = pd.read_csv(path, sep="\t", comment="#", index_col=0)
    return fc["Length"] if lengths else fc.iloc[:, -1]


def read_paqr(wf_out, libraries):
    # tandem PAS expression (mean coverage) of all chromosomes; columns <library>.<chr>
    parts = []
    for f in sorted(glob.glob(os.path.join(wf_out, "PAQR", "human", "quantifications", "tandem_pas_expression-mean.chr*.tsv"))):
        d = pd.read_csv(f, sep="\t")
        d = d.rename(columns={c: c.rsplit(".chr", 1)[0] for c in d.columns[10:]})
        parts.append(d[list(d.columns[:10]) + libraries])
    return pd.concat(parts, ignore_index=True)


def step_prepare(out, cfg):
    d = os.path.join(out, "inputs")
    if os.path.exists(os.path.join(d, "done")):
        log("prepare: done before, skipped")
        return
    os.makedirs(d, exist_ok=True)
    lib = sample_table(cfg)
    genes, _ = parse_gtf(cfg["gtf"])
    genes.to_csv(os.path.join(d, "genes.tsv"), sep="\t", index_label="gene_id")

    # gene counts per library (read pairs); MITO T2 + its technical replicate
    fc_dir = os.path.join(cfg["wf_out"], "gene_expression_quantification", "human", "FeatureCounts_standard")
    lc = pd.DataFrame({l: read_featurecounts(os.path.join(fc_dir, f"{l}.standard.txt")) for l in lib["library"]}).reindex(genes.index).fillna(0)
    lc.round().astype(np.int64).to_csv(os.path.join(d, "gene_counts.per_library.tsv.gz"), sep="\t", index_label="gene_id")
    read_featurecounts(os.path.join(fc_dir, f"{lib['library'].iloc[0]}.standard.txt"), lengths=True).rename("length").to_csv(
        os.path.join(d, "gene_lengths.tsv"), sep="\t", index_label="gene_id")
    merged_of = lib.set_index("library")["sample"]
    counts = lc.T.groupby(merged_of.reindex(lc.columns).to_numpy()).sum().T.round().astype(np.int64)
    meta = lib.groupby("sample", sort=False).agg({c: "first" for c in ["batch", "fraction", "construct", "replicate", "kind", "cohort"]}
                                                  ).join(lib.groupby("sample")["library"].agg(",".join)).reset_index()
    meta["assigned_read_pairs"] = counts.sum().reindex(meta["sample"]).to_numpy()
    meta["use"] = True
    meta.to_csv(os.path.join(d, "metadata.tsv"), sep="\t", index=False)
    lib.to_csv(os.path.join(d, "libraries.tsv"), sep="\t", index=False)
    counts[meta["sample"]].to_csv(os.path.join(d, "gene_counts.tsv.gz"), sep="\t", index_label="gene_id")

    # QC: mapping statistics and featureCounts summary per library
    qc = []
    for l in lib["library"]:
        with open(os.path.join(cfg["wf_out"], "mapping_stats", f"{l}.mapping_stats.txt")) as fh:
            v = [int(x) for x in fh.read().split("\n") if x.strip() and not x.startswith(">")]
        sm = pd.read_csv(os.path.join(fc_dir, f"{l}.standard.txt.summary"), sep="\t", index_col=0).iloc[:, 0]
        qc.append({"library": l, "uniquely_mapped_reads": v[0], "mapped_reads": v[1], "input_reads": v[2],
                   **{f"featureCounts_{k}": int(x) for k, x in sm.items() if x > 0}})
    pd.DataFrame(qc).to_csv(os.path.join(d, "read_qc.per_library.tsv"), sep="\t", index=False)

    # tandem PAS of terminal exons (PAQR)
    paqr = read_paqr(cfg["wf_out"], lib["library"].tolist())
    vals = paqr[lib["library"]]
    bad = paqr.loc[(vals < 0).any(axis=1), "exon"].unique()
    paqr = paqr.loc[~paqr["exon"].isin(bad)].reset_index(drop=True)
    pas_lib = paqr.set_index("pas")[lib["library"]].round().astype(np.int64)
    pas = pas_lib.T.groupby(merged_of.reindex(pas_lib.columns).to_numpy()).sum().T[meta["sample"]]
    tab = paqr.set_index("pas")[["chrom", "start", "end", "strand", "polyAsite_exon_idx", "nr_polyAsites_on_exon", "exon", "gene"]].rename(
        columns={"gene": "gene_id"})
    tab["pas_pos"] = tab.index.str.split(":").str[1].astype(np.int64)
    ex = tab["exon"].str.split(":", expand=True)
    tab["te_start"], tab["te_end"] = ex[3].astype(np.int64), ex[4].astype(np.int64)
    # 3' UTR length: from the stop codon of the Ensembl canonical transcript when it lies in the terminal exon, else from
    # the 5' end of the terminal exon (as for the ONT data)
    tab = tab.join(genes[["gene_name", "gene_type", "stop_pos"]], on="gene_id")
    plus, p, stop = (tab["strand"] == "+").to_numpy(), tab["pas_pos"].to_numpy(), tab["stop_pos"].to_numpy(float)
    ts, te = tab["te_start"].to_numpy(), tab["te_end"].to_numpy()
    stop_in = ~np.isnan(stop) & (stop >= ts) & (stop <= te)
    origin = np.where(plus, np.where(stop_in, stop, ts - 1), np.where(stop_in, stop, te + 1))
    tab["utr_length"] = np.where(plus, p - origin, origin - p)
    tab["utr_origin"] = np.where(stop_in, "stop_codon", "exon_start")
    tab["reads"] = pas.sum(axis=1).reindex(tab.index).to_numpy()
    tab["passes_min_reads"] = (tab["reads"] >= MIN_PAS_COUNT) & (tab["utr_length"] > 0)
    tab.to_csv(os.path.join(d, "pas_table.tsv.gz"), sep="\t", index_label="pas")
    pas.loc[tab.index[tab["passes_min_reads"]]].to_csv(os.path.join(d, "pas_counts.tsv.gz"), sep="\t", index_label="pas")
    log(f"prepare: {len(meta)} samples ({len(lib)} libraries); {len(counts)} genes; {len(tab):,} tandem PAS in {tab['exon'].nunique():,} "
        f"terminal exons ({len(bad):,} exons with missing values left out), {int(tab['passes_min_reads'].sum()):,} with >= {MIN_PAS_COUNT}")
    open(os.path.join(d, "done"), "w").close()


# ---------------------------------------------------------------------------------------------------------------
def _fragments(job):
    # read pairs (first mates of proper pairs, uniquely mapped) and their distinct (chrom, position, strand, insert size)
    lib, bam, tmp = job
    cmd = (f"{SAMTOOLS} view -f 66 -F 2316 {bam} | awk -v OFS='\\t' '{{print $3,$4,and($2,16),$9}}' | "
           f"LC_ALL=C sort -S 6G --parallel=4 -T {tmp} | uniq -c | awk '{{n+=$1; u+=1; if ($1==1) s+=1}} END {{print n, u, s}}'")
    n, u, s = subprocess.run(["bash", "-c", "set -o pipefail; " + cmd], check=True, capture_output=True, text=True).stdout.split()
    return {"library": lib, "read_pairs": int(n), "distinct_fragments": int(u), "singleton_fragments": int(s)}


def step_complexity(out, cfg, threads):
    f = os.path.join(out, "inputs", "fragment_duplicates.per_library.tsv")
    if os.path.exists(f):
        log("complexity: done before, skipped")
        return
    lib = pd.read_csv(os.path.join(out, "inputs", "libraries.tsv"), sep="\t")
    tmp = os.environ.get("TMPDIR", "/tmp")
    jobs = [(l, os.path.join(cfg["wf_out"], "samples", l, "map_genome", f"{l}.Aligned.sortedByCoord.uniquely_mapped.bam"), tmp)
            for l in lib["library"]]
    with ProcessPoolExecutor(max_workers=max(1, threads // 6)) as ex:
        rows = list(ex.map(_fragments, jobs))
    pd.DataFrame(rows).to_csv(f, sep="\t", index=False)
    log(f"complexity: {len(rows)} libraries")


# ---------------------------------------------------------------------------------------------------------------
def step_sanity(out, threads):
    d = os.path.join(out, "inputs")
    meta = pd.read_csv(os.path.join(d, "metadata.tsv"), sep="\t")
    meta = meta.loc[meta["use"]].reset_index(drop=True)
    genes = pd.read_csv(os.path.join(d, "gene_counts.tsv.gz"), sep="\t", index_col=0)
    pas = pd.read_csv(os.path.join(d, "pas_counts.tsv.gz"), sep="\t", index_col=0)
    lengths = pd.read_csv(os.path.join(d, "gene_lengths.tsv"), sep="\t", index_col=0)["length"]
    tot, frac = meta.loc[meta["fraction"] == "Total"], meta.loc[meta["fraction"] != "Total"]
    for level, counts, lens in (("genes", genes, lengths), (LEVEL_PAS, pas, None)):
        for part, m in (("Total", tot), ("fractions", frac)) + ((("joint", meta),) if level == LEVEL_PAS else ()):
            od = os.path.join(out, "sanity", level, part)
            os.makedirs(od, exist_ok=True)
            run_one_sanity(counts, m, od, threads, gene_lengths=lens)


def step_fit_cohorts(out, threads):
    import sanity_utils as su
    jobs = []
    for level in ("genes", LEVEL_PAS):
        tot = load(os.path.join(out, "sanity", level, "Total", "sanity.pkl"))
        frac = load(os.path.join(out, "sanity", level, "fractions", "sanity.pkl"))
        feats = tot["averages"].means.index.intersection(frac["averages"].means.index)
        od = os.path.join(out, "fracSanity", level, "cohorts")
        os.makedirs(od, exist_ok=True)
        for cohort in frac["averages"].means.columns:
            inp = su.prepare_frac_sanity_inputs(tot["averages"], frac["averages"], "Total", cohort, genes=feats)
            jobs.append((os.path.join(od, re.sub(r"[^A-Za-z0-9]+", "_", cohort) + ".pkl"), inp["ltq_total"], inp["var_total"],
                         inp["ltq_pd"], inp["var_pd"], inp["gene_ids"], {"cohort": cohort, "level": level}))
    run_jobs(jobs, threads)


def step_fit_replicates(out, threads):
    jobs = []
    tot = load(os.path.join(out, "sanity", "genes", "Total", "sanity.pkl"))
    frac = load(os.path.join(out, "sanity", "genes", "fractions", "sanity.pkl"))
    genes = tot["averages"].means.index.intersection(frac["averages"].means.index)
    ta, rp = tot["averages"], frac["result"]
    od = os.path.join(out, "fracSanity", "genes", "libraries")
    os.makedirs(od, exist_ok=True)
    for _, r in frac["metadata"].iterrows():
        jobs.append((os.path.join(od, re.sub(r"[^A-Za-z0-9]+", "_", r["sample"]) + ".pkl"),
                     ta.means.loc[genes, "Total"].to_numpy(), ta.errors.loc[genes, "Total"].to_numpy() ** 2,
                     rp.log2_expression.loc[genes, r["sample"]].to_numpy(), rp.log2_sd.loc[genes, r["sample"]].to_numpy() ** 2,
                     genes, {"cohort": r["cohort"], "sample": r["sample"], "replicate": r["replicate"], "level": "genes"}))
    run_jobs(jobs, threads)


def step_usage(out, threads):
    import sanity_utils as su
    od = os.path.join(out, f"usage_{LEVEL_PAS}")
    os.makedirs(od, exist_ok=True)
    joint = load(os.path.join(out, "sanity", LEVEL_PAS, "joint", "sanity.pkl"))
    res, meta = joint["result"], joint["metadata"]
    tab = pd.read_csv(os.path.join(out, "inputs", "pas_table.tsv.gz"), sep="\t", index_col=0).loc[res.log2_expression.index]
    cohorts = [c for c in meta["cohort"].unique() if c != "Total"]
    cond_groups = {c: meta.loc[meta["cohort"] == c, "sample"].tolist() for c in meta["cohort"].unique()}

    # relative usage of each PAS against the other PAS of its terminal exon
    p = os.path.join(od, "pas_usage.per_sample.pkl")
    groups = su.isoform_groups(tab["exon"].to_dict())
    if not os.path.exists(p):
        y, s = su.grouped_log2_ratios(res, groups)
        avg = su.average_replicates(y, s, cond_groups, method="reml")
        dump({"groups": groups, "y": y, "sd": s, "averages": avg}, p)
        log(f"usage: {len(groups['group'].unique()):,} PAS tests in {tab.loc[groups['feature'].unique(), 'exon'].nunique():,} terminal exons")
    for c in cohorts:
        f = os.path.join(od, "pas_usage." + re.sub(r"[^A-Za-z0-9]+", "_", c) + "_vs_Total.tsv.gz")
        if os.path.exists(f):
            continue
        du = su.differential_relative_usage(res, groups, meta, c, "Total", sample_col="sample", cond_col="cohort", method="reml",
                                            n_draws=500)
        du.to_csv(f, sep="\t", index_label="group")
        log(f"usage {c}: {len(du):,} tests, padj < 0.05: {int((du['padj'] < 0.05).sum()):,}")

    # 3' UTR length index (one terminal exon per gene: the one with the most reads)
    best = tab.groupby(["gene_id", "exon"])["reads"].sum().reset_index().sort_values("reads").drop_duplicates("gene_id", keep="last")["exon"]
    sites = tab.loc[tab["exon"].isin(set(best))].reset_index().rename(columns={"pas": "feature", "gene_id": "gene"})[["gene", "feature", "utr_length"]]
    p = os.path.join(od, "utr_length_index.per_sample.pkl")
    if not os.path.exists(p):
        y, s = su.utr_length_index(res, sites)
        avg = su.average_replicates(y, s, cond_groups, method="reml")
        dump({"sites": sites, "y": y, "sd": s, "averages": avg}, p)
        log(f"APA: log2 mean 3' UTR length for {len(y):,} genes")
    for c in cohorts:
        f = os.path.join(od, "utr_length." + re.sub(r"[^A-Za-z0-9]+", "_", c) + "_vs_Total.tsv.gz")
        if os.path.exists(f):
            continue
        du = su.differential_utr_length(res, sites, meta, c, "Total", sample_col="sample", cond_col="cohort", method="reml", n_draws=500)
        du.to_csv(f, sep="\t", index_label="gene")
        log(f"APA {c}: {len(du):,} genes, {du['call'].value_counts().to_dict()}")
        pdd = su.differential_proximal_distal(res, sites, meta, c, "Total", splits="first", sample_col="sample", cond_col="cohort",
                                              method="reml", n_draws=500)
        pdd.to_csv(f.replace("utr_length.", "proximal_vs_distal."), sep="\t", index_label="group")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--steps", default="prepare,complexity,sanity,fit_cohorts,fit_replicates,usage")
    a = ap.parse_args()
    os.environ.setdefault("OMP_NUM_THREADS", "6")
    cfg = paths()
    os.makedirs(a.outdir, exist_ok=True)
    for step in a.steps.split(","):
        t0 = time.time()
        log(f"=== {step}")
        if step == "prepare":
            step_prepare(a.outdir, cfg)
        elif step == "complexity":
            step_complexity(a.outdir, cfg, a.threads)
        elif step == "sanity":
            step_sanity(a.outdir, a.threads)
        elif step == "fit_cohorts":
            step_fit_cohorts(a.outdir, a.threads)
        elif step == "fit_replicates":
            step_fit_replicates(a.outdir, a.threads)
        elif step == "usage":
            step_usage(a.outdir, a.threads)
        else:
            raise ValueError(step)
        log(f"=== {step} finished in {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
