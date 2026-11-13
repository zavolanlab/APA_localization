"""
fracSanity, cleavage-site (CS) usage and APA analysis of the ONT subcellular fractions (siCtrl only) against the
siCtrl HCT116 Total, with Sanity_utils. Heavy steps for ONT_organelles.fracSanity.ipynb, which reads the outputs.

Steps (each writes its own files into --outdir and is skipped when they exist, so nothing is overwritten):
  prepare      sample metadata, gene counts (featureCounts, all genes of the annotation), gene annotation
  prepare_cs   CS counts from the read tags (XT = gene, XF = cleavage site): every position is its own CS (no
               clustering of nearby positions); CS annotation (site class, PAQR3 terminal exon and 3' UTR length in
               it, PolyASite atlas overlap)
  sanity       C++ Sanity (v_method EAP): genes and CS, Total and fractions separately (fracSanity inputs) and
               jointly (CS usage and APA, fraction vs Total); ISMARA replicate averages per cohort
  fit_cohorts  fracSanity per cohort (fraction average vs Total average), genes and CS
  fit_replicates  fracSanity per fraction library vs the Total average (genes): input of the gene-set analysis
  usage        CS relative usage (each CS vs the other CS of its gene) and 3' UTR length index (tandem APA: CS in
               the PAQR3 terminal exon of the gene), per cohort vs Total

The CS level is "cs_pos" (inputs/cs_pos_*, sanity/cs_pos, fracSanity/cs_pos, usage_cs_pos); the outputs of the
earlier version with CS clustered within 12 nt ("cs", usage) are left in place but no longer used.

Usage (sanity_utils environment):
    python WF/scripts/ONT_fracSanity_pipeline.py --outdir <dir> --threads 48 [--steps prepare,sanity,...]
"""

import argparse
import glob
import gzip
import os
import pickle
import re
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from dotenv import load_dotenv

PROJECT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# ---------------------------------------------------------------------------------------------------------------
# samples (as in ONT_organelles_vs_Total.ipynb); siCtrl only
TOTAL_SICTRL = {"barcode02": "R1", "barcode05": "R2", "barcode08": "R3"}
POOLING_SHEETS = {
    "p40713_o43608_1": {1: "14.07.25 siCTRL Free", 2: "14.07.25 siCTRL Org2", 3: "14.07.25 siCTRL CANX", 4: "14.07.25 siCTRL TIS",
                        5: "14.07.25 si5 Free", 6: "14.07.25 si5 Org2", 7: "14.07.25 si5 CANX", 8: "14.07.25 si5 TIS",
                        9: "18.07.25 siCTRL Free", 10: "18.07.25 siCTRL Org2", 11: "18.07.25 siCTRL CANX", 12: "18.07.25 siCTRL TIS",
                        13: "18.07.25 si5 Free", 14: "18.07.25 si5 Org2", 15: "18.07.25 si5 CANX", 16: "18.07.25 si5 TIS",
                        22: "23.01.26 si5 Org2"},
    "p40713_o43608_2": {1: "29.09 siCTRL Free", 2: "29.09 siCTRL Org2", 3: "29.09 siCTRL Golgi", 4: "29.09 siCTRL Perox",
                        5: "29.09 si5 Free", 6: "29.09 si5 Org2", 7: "29.09 si5 Golgi", 8: "29.09 si5 Perox",
                        9: "2.10 siCTRL Free", 10: "2.10 siCTRL Org2", 11: "2.10 siCTRL Golga2", 12: "2.10 siCTRL Pex14",
                        13: "2.10 si5 Free", 14: "2.10 si5 Org2", 15: "2.10 si5 Golga2", 16: "2.10 si5 Pex14"}}
BATCH_OF_FLOWCELL = {"p40713_o43608_1": "ER/TIS", "p40713_o43608_2": "Golgi/PEX"}
END_TJ_LIBRARIES = {f"barcode{n + 1:02d}": name for n, name in enumerate(
    f"{d} {c} {f}" for d in ("11.09", "22.09", "25.09") for c in ("siCTRL", "si5") for f in ("Free", "Org2", "Eea1", "Occl"))}
FRACTION_NAMES = {"Free": "Free", "Org2": "Org2", "CANX": "ER", "TIS": "TIS", "Golgi": "Golgi", "Golga2": "Golgi", "Perox": "PEX",
                  "Pex14": "PEX", "Eea1": "Endo", "Occl": "TJ"}
MIN_READS = 100_000          # reads on nuclear protein-coding genes with >= 10 CPM in the pooled Total (as in ONT_organelles_vs_Total)
MIN_CS_READS = 5             # reads of a CS (one position) summed over the analysed libraries
TERMINAL_EXTENSION = 5000    # CS up to this distance downstream of the annotated end of the last exon count as terminal
ATLAS_WINDOW = 10            # CS within this distance of a PolyASite cluster count as in the atlas
TE_EXTENSION = 10            # CS up to this distance downstream of the PAQR3 terminal exon end count as in it
CS_LEVEL = "cs_pos"
PI_MAX = 0.5                 # bound on the outlier fraction of fracSanity (as in the P-body analysis)


def dump(obj, path):
    tmp = path + ".tmp"
    with open(tmp, "wb") as fh:
        pickle.dump(obj, fh, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp, path)


def load(path):
    with open(path, "rb") as fh:
        return pickle.load(fh)


def log(msg):
    print(time.strftime("%H:%M:%S"), msg, flush=True)


def paths():
    load_dotenv(os.path.join(PROJECT, "APA_localization.scicore.env"), override=True)
    runs_dir, group = os.getenv("WF_RUNS_DIR"), os.getenv("LAB_GROUP_DIR")
    return {"runs": {"Total": os.path.join(group, "APA_AI_data", "wf_runs", "v2p0_ONT_HCT116_v1", "results"),
                     "ER_TIS_Golgi_PEX": os.path.join(runs_dir, "v2p0_ONT_FACSorganelles_ER_TIS_PEX_Golgi_Aleksei_v1",
                                                      "results_alignment_with_SpikeIns"),
                     "End_TJ": os.path.join(runs_dir, "v2p0_ONT_FACSorganelles_End_TJ_Aleksei_v1", "results_alignment_with_SpikeIns")},
            "gtf": os.getenv("HUMAN_BASIC_ANNOTATION_FILE"), "atlas": os.getenv("HUMAN_POLYASITE_ATLAS"),
            "paqr3": os.path.join(group, "APA_AI_data", "paqr3_segments_results_atlasv3_min77", "paqr3_result_atlasv3_min77")}


# ---------------------------------------------------------------------------------------------------------------
# prepare
def sample_table(runs):
    rows = []
    for bc, rep in TOTAL_SICTRL.items():
        rows.append({"source": "Total", "sample_id": bc, "batch": "HCT116 Total", "fraction": "Total", "replicate": rep})
    for flowcell, sheet in POOLING_SHEETS.items():
        for n, name in sheet.items():
            date, cond, frac = name.split()
            if cond == "siCTRL":
                rows.append({"source": "ER_TIS_Golgi_PEX", "sample_id": f"barcode{n:02d}_{flowcell}", "batch": BATCH_OF_FLOWCELL[flowcell],
                             "fraction": FRACTION_NAMES[frac], "replicate": date})
    for bc, name in END_TJ_LIBRARIES.items():
        date, cond, frac = name.split()
        if cond == "siCTRL":
            rows.append({"source": "End_TJ", "sample_id": bc, "batch": "End/TJ", "fraction": FRACTION_NAMES[frac], "replicate": date})
    meta = pd.DataFrame(rows)
    meta["cohort"] = np.where(meta["fraction"] == "Total", "Total", meta["batch"] + ": " + meta["fraction"])
    meta["sample"] = np.where(meta["fraction"] == "Total", "Total " + meta["replicate"], meta["cohort"] + " " + meta["replicate"])
    meta["featurecounts"] = [os.path.join(runs[s], "gene_assignments", f"{i}.featureCounts.txt") for s, i in zip(meta["source"], meta["sample_id"])]
    meta["read_tags"] = [os.path.join(runs[s], "read_tags", f"{i}.tags.tsv.gz") for s, i in zip(meta["source"], meta["sample_id"])]
    assert meta["sample"].is_unique
    for f in list(meta["featurecounts"]) + list(meta["read_tags"]):
        assert os.path.exists(f), f
    return meta


def parse_gtf(gtf):
    # genes, exons of all transcripts, and stop codon + last exon of the Ensembl canonical transcript of each gene
    rows = []
    for chunk in pd.read_csv(gtf, sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4, 6, 8],
                             names=["chrom", "type", "start", "end", "strand", "attr"], chunksize=500_000):
        rows.append(chunk.loc[chunk["type"].isin(["gene", "transcript", "exon", "stop_codon"])])
    g = pd.concat(rows, ignore_index=True)
    g["gene_id"] = g["attr"].str.extract(r'gene_id "([^"]+)"')[0]
    genes = g.loc[g["type"] == "gene"].copy()
    genes["gene_name"] = genes["attr"].str.extract(r'gene_name "([^"]+)"')[0]
    genes["gene_type"] = genes["attr"].str.extract(r'gene_type "([^"]+)"')[0]
    genes = genes.set_index("gene_id")[["gene_name", "gene_type", "chrom", "start", "end", "strand"]]
    tx = g.loc[g["type"] != "gene"].copy()
    tx["transcript_id"] = tx["attr"].str.extract(r'transcript_id "([^"]+)"')[0]
    canonical = set(tx.loc[(tx["type"] == "transcript") & tx["attr"].str.contains('tag "Ensembl_canonical"'), "transcript_id"])
    exons = tx.loc[tx["type"] == "exon", ["gene_id", "transcript_id", "chrom", "start", "end", "strand"]].reset_index(drop=True)
    stops = tx.loc[(tx["type"] == "stop_codon") & tx["transcript_id"].isin(canonical)]
    # last exon of the canonical transcript: the most downstream exon on its strand
    can_ex = exons.loc[exons["transcript_id"].isin(canonical)]
    plus = can_ex["strand"] == "+"
    last = can_ex.loc[np.concatenate([can_ex.loc[plus].groupby("transcript_id")["end"].idxmax().to_numpy(),
                                      can_ex.loc[~plus].groupby("transcript_id")["start"].idxmin().to_numpy()])]
    term = last.set_index("transcript_id")[["gene_id", "start", "end"]].rename(
        columns={"start": "last_exon_start", "end": "last_exon_end"})
    st = stops.groupby("transcript_id").agg(strand=("strand", "first"), stop_start=("start", "min"), stop_end=("end", "max"))
    # stop codon position: last nucleotide of the stop codon on the transcript (+: end, -: start)
    st["stop_pos"] = np.where(st["strand"] == "+", st["stop_end"], st["stop_start"])
    term = term.join(st[["stop_pos"]], how="left")
    # exonic length of the 3' UTR upstream of the last exon (spliced 3' UTRs), so that the 3' UTR length of a CS in the
    # last exon is counted in transcript coordinates
    offset = {}
    for tid, e in can_ex.loc[can_ex["transcript_id"].isin(term.index[term["stop_pos"].notna()])].groupby("transcript_id"):
        stop, le = term.loc[tid, "stop_pos"], term.loc[tid]
        rest = e.loc[e["start"] != le["last_exon_start"]]
        if e["strand"].iloc[0] == "+":
            offset[tid] = int(np.clip(rest["end"] - np.maximum(rest["start"] - 1, stop), 0, None).sum())
        else:
            offset[tid] = int(np.clip(np.minimum(rest["end"] + 1, stop) - rest["start"], 0, None).sum())
    term["utr_offset"] = pd.Series(offset)
    term = term.reset_index().drop_duplicates("gene_id").set_index("gene_id")
    genes = genes.join(term[["transcript_id", "last_exon_start", "last_exon_end", "stop_pos", "utr_offset"]].rename(
        columns={"transcript_id": "canonical_transcript"}), how="left")
    return genes, exons


def merged_exons(exons):
    # gene -> (starts, ends) of the union of its exons over all transcripts
    out = {}
    for gid, e in exons.groupby("gene_id"):
        iv = e[["start", "end"]].sort_values("start").to_numpy()
        merged = [list(iv[0])]
        for s, t in iv[1:]:
            if s <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], t)
            else:
                merged.append([s, t])
        m = np.array(merged)
        out[gid] = (m[:, 0], m[:, 1])
    return out


def read_cs_positions(path):
    # reads per (chrom, strand, gene, cleavage site) of one library: reads with an estimated poly(A) tail (pt > 0) and
    # assigned to one gene (XT without ','), as in CD47_analysis.ipynb
    agg = []
    for chunk in pd.read_csv(path, sep="\t", usecols=["chrom", "strand", "pt", "XT", "XF"], chunksize=2_000_000,
                             dtype={"chrom": str, "strand": str, "XT": str}):
        chunk = chunk.loc[(chunk["pt"] > 0) & chunk["XT"].notna() & ~chunk["XT"].str.contains(",", regex=False) & chunk["XF"].notna()]
        agg.append(chunk.astype({"XF": np.int64}).groupby(["chrom", "strand", "XT", "XF"]).size())
    s = pd.concat(agg)
    return s.groupby(level=[0, 1, 2, 3]).sum()


def paqr3_terminal_exons(prefix):
    # PAQR3 segments (<gene index>.<segment index>, segment index increasing from the first to the last exon on either
    # strand): the terminal exon of a gene is its segment with the largest index; gene ids without version
    genes = pd.read_csv(prefix + "_genes.bed", sep="\t", header=None, usecols=[3, 4], names=["gene_idx", "gene_id"], dtype=str)
    seg = pd.read_csv(prefix + "_segments.bed", sep="\t", header=None, usecols=[0, 1, 2, 3, 5],
                      names=["chrom", "start", "end", "segment", "strand"], dtype={"segment": str})
    seg["gene_idx"] = seg["segment"].str.split(".").str[0]
    seg["seg_idx"] = seg["segment"].str.split(".").str[1].astype(int)
    te = seg.loc[seg.groupby("gene_idx")["seg_idx"].idxmax()].merge(genes, on="gene_idx")
    te["gene"] = te["gene_id"].str.split(".").str[0]
    te["te_start"] = te["start"] + 1                                    # BED -> 1-based
    return te.drop_duplicates("gene").set_index("gene")[["segment", "chrom", "strand", "te_start", "end"]].rename(
        columns={"end": "te_end", "segment": "te_segment"})


def annotate_sites(cl, genes, exons, atlas_path, te):
    cl = cl.join(genes[["gene_name", "gene_type", "start", "end", "strand", "canonical_transcript", "last_exon_start",
                        "last_exon_end", "stop_pos"]].rename(columns={"start": "gene_start", "end": "gene_end",
                                                                       "strand": "gene_strand"}), on="gene_id")
    p, plus = cl["cs_pos"].to_numpy(), (cl["strand"] == "+").to_numpy()
    les, lee, stop = cl["last_exon_start"].to_numpy(float), cl["last_exon_end"].to_numpy(float), cl["stop_pos"].to_numpy(float)
    has = ~np.isnan(stop)
    # site class from the GTF (overview only): last exon of the Ensembl canonical transcript downstream of the stop codon
    terminal = has & np.where(plus, (p >= les) & (p > stop) & (p <= lee + TERMINAL_EXTENSION),
                              (p <= lee) & (p < stop) & (p >= les - TERMINAL_EXTENSION))
    mex = merged_exons(exons.loc[exons["gene_id"].isin(set(cl["gene_id"]))])
    exonic = np.zeros(len(cl), bool)
    for i, (gid, pos) in enumerate(zip(cl["gene_id"], p)):
        if gid in mex:
            s, e = mex[gid]
            k = np.searchsorted(s, pos, side="right") - 1
            exonic[i] = k >= 0 and pos <= e[k]
    in_gene = (p >= cl["gene_start"].to_numpy()) & (p <= cl["gene_end"].to_numpy())
    cl["site_class"] = np.select([terminal, exonic, in_gene], ["terminal", "exonic", "intronic"], "outside")

    # tandem APA: CS in the PAQR3 terminal exon of its gene (same strand, up to TE_EXTENSION nt past its 3' end).
    # 3' UTR length = distance from the stop codon of the canonical transcript (+1) when the stop codon lies in the
    # terminal exon, otherwise from the 5' end of the terminal exon (+1); only the within-gene differences matter
    t = te.reindex(cl["gene_id"].str.split(".").str[0].to_numpy())
    ts, tend = t["te_start"].to_numpy(float), t["te_end"].to_numpy(float)
    ok = (t["chrom"].to_numpy() == cl["chrom"].to_numpy()) & (t["strand"].to_numpy() == cl["strand"].to_numpy())
    in_te = ok & np.where(plus, (p >= ts) & (p <= tend + TE_EXTENSION), (p <= tend) & (p >= ts - TE_EXTENSION))
    stop_in_te = has & (stop >= ts) & (stop <= tend)
    origin = np.where(plus, np.where(stop_in_te, stop, ts - 1), np.where(stop_in_te, stop, tend + 1))
    utr = np.where(plus, p - origin, origin - p)
    cl["te_segment"] = np.where(in_te, t["te_segment"].to_numpy(), None)
    cl["te_utr_origin"] = np.where(in_te, np.where(stop_in_te, "stop_codon", "exon_start"), None)
    cl["utr_length"] = np.where(in_te & (utr > 0), utr, np.nan)
    cl["in_terminal_exon"] = in_te & (utr > 0)

    # PolyASite atlas (the env file names a .bed.gz; the uncompressed .bed is used if only that exists)
    if not os.path.exists(atlas_path) and os.path.exists(re.sub(r"\.gz$", "", atlas_path)):
        atlas_path = re.sub(r"\.gz$", "", atlas_path)
    atlas = pd.read_csv(atlas_path, sep="\t", header=None, usecols=[0, 1, 2, 5], names=["chrom", "start", "end", "strand"],
                        dtype={"chrom": str})
    atlas["chrom"] = np.where(atlas["chrom"].str.startswith("chr"), atlas["chrom"], "chr" + atlas["chrom"])
    in_atlas = np.zeros(len(cl), bool)
    for (c, s), idx in cl.groupby(["chrom", "strand"]).groups.items():
        a = atlas.loc[(atlas["chrom"] == c) & (atlas["strand"] == s)].sort_values("start")
        if a.empty:
            continue
        st, en = a["start"].to_numpy() + 1, a["end"].to_numpy()       # BED -> 1-based
        pos = cl.loc[idx, "cs_pos"].to_numpy()
        k = np.searchsorted(st, pos + ATLAS_WINDOW, side="right") - 1
        ok = (k >= 0) & (pos - ATLAS_WINDOW <= en[np.clip(k, 0, None)])
        in_atlas[cl.index.get_indexer(idx)] = ok
    cl["in_PolyASite_atlas"] = in_atlas
    return cl


def step_prepare(out, cfg):
    d = os.path.join(out, "inputs")
    if os.path.exists(os.path.join(d, "done")):
        log("prepare: done before, skipped")
        return
    os.makedirs(d, exist_ok=True)
    meta = sample_table(cfg["runs"])
    log(f"prepare: {len(meta)} siCtrl libraries; parsing {cfg['gtf']}")
    genes, exons = parse_gtf(cfg["gtf"])
    genes.to_csv(os.path.join(d, "genes.tsv"), sep="\t", index_label="gene_id")

    counts = {}
    for _, r in meta.iterrows():
        fc = pd.read_csv(r["featurecounts"], sep="\t", comment="#", index_col=0).iloc[:, -1]
        counts[r["sample"]] = fc.reindex(genes.index)
    counts = pd.DataFrame(counts).fillna(0).astype(np.int64)
    # depth filter on the libraries, as in ONT_organelles_vs_Total.ipynb
    nuc_pc = genes.index[(genes["gene_type"] == "protein_coding") & (genes["chrom"] != "chrM")]
    pooled = counts[meta.loc[meta["fraction"] == "Total", "sample"]].sum(axis=1)[nuc_pc]
    universe = pooled.index[pooled / pooled.sum() * 1e6 >= 10]
    meta["reads_on_genes"] = counts.loc[universe].sum().reindex(meta["sample"]).to_numpy()
    meta["assigned_reads"] = counts.sum().reindex(meta["sample"]).to_numpy()
    meta["use"] = (meta["fraction"] == "Total") | (meta["reads_on_genes"] >= MIN_READS)
    meta.to_csv(os.path.join(d, "metadata.tsv"), sep="\t", index=False)
    counts.to_csv(os.path.join(d, "gene_counts.tsv.gz"), sep="\t", index_label="gene_id")
    log(f"prepare: {int(meta['use'].sum())} libraries used; {len(counts)} genes")

    open(os.path.join(d, "done"), "w").close()


def step_prepare_cs(out, cfg):
    # every cleavage position (XF) of a gene (XT) is one CS, named <gene>_<position>
    d = os.path.join(out, "inputs")
    if os.path.exists(os.path.join(d, f"{CS_LEVEL}_done")):
        log("prepare_cs: done before, skipped")
        return
    meta = pd.read_csv(os.path.join(d, "metadata.tsv"), sep="\t")
    genes = pd.read_csv(os.path.join(d, "genes.tsv"), sep="\t", index_col=0)
    _, exons = parse_gtf(cfg["gtf"])
    pos = {}
    for _, r in meta.loc[meta["use"]].iterrows():
        pos[r["sample"]] = read_cs_positions(r["read_tags"])
        log(f"  read tags {r['sample']}: {pos[r['sample']].sum():,} reads at {len(pos[r['sample']]):,} positions")
    pos = pd.DataFrame(pos).fillna(0).astype(np.int64)
    pos.index.names = ["chrom", "strand", "gene_id", "cs_pos"]
    # 3' ends on the chromosome and strand of the assigned gene only (a gene and position on two strands would give
    # two CS of one name)
    idx = pos.index.to_frame(index=False)
    g = genes.reindex(idx["gene_id"])
    same = (g["chrom"].to_numpy() == idx["chrom"].to_numpy()) & (g["strand"].to_numpy() == idx["strand"].to_numpy())
    log(f"prepare_cs: {pos.to_numpy()[~same].sum():,} of {pos.to_numpy().sum():,} reads at positions off the gene's chromosome or strand, left out")
    pos = pos.loc[same]
    cl = pos.sum(axis=1).rename("reads").reset_index()
    cl["cs_id"] = cl["gene_id"] + "_" + cl["cs_pos"].astype(str)
    assert cl["cs_id"].is_unique
    cs_counts = pos.set_axis(cl["cs_id"].to_numpy())
    cl = cl.set_index("cs_id")
    cl = cl.loc[cl["gene_id"].isin(genes.index)]
    cl = annotate_sites(cl, genes, exons, cfg["atlas"], paqr3_terminal_exons(cfg["paqr3"]))
    cl["passes_min_reads"] = cl["reads"] >= MIN_CS_READS
    cs_counts = cs_counts.loc[cl.index]
    cl.to_csv(os.path.join(d, f"{CS_LEVEL}_table.tsv.gz"), sep="\t", index_label="cs_id")
    cs_counts.loc[cl.index[cl["passes_min_reads"]]].to_csv(os.path.join(d, f"{CS_LEVEL}_counts.tsv.gz"), sep="\t", index_label="cs_id")
    k = cl["passes_min_reads"]
    log(f"prepare_cs: {len(cl):,} CS ({int(k.sum()):,} with >= {MIN_CS_READS} reads, {cl.loc[k, 'reads'].sum() / cl['reads'].sum():.1%} "
        f"of the reads) in {cl['gene_id'].nunique():,} genes; classes {cl.loc[k, 'site_class'].value_counts().to_dict()}; "
        f"in PAQR3 terminal exons: {int(cl.loc[k, 'in_terminal_exon'].sum()):,}")
    open(os.path.join(d, f"{CS_LEVEL}_done"), "w").close()


# ---------------------------------------------------------------------------------------------------------------
# sanity
def run_one_sanity(counts, meta, outdir, threads, cond_col="cohort", gene_lengths=None):
    # gene_lengths: feature -> length for short-read data (log2 TPM instead of log2 CPM); None for ONT (one read per molecule)
    import sanity_utils as su
    if os.path.exists(os.path.join(outdir, "sanity.pkl")):
        log(f"sanity {outdir}: done before, skipped")
        return
    sub = counts[meta["sample"].tolist()]
    sub = sub.loc[sub.sum(axis=1) > 0]
    t0 = time.time()
    if gene_lengths is not None:
        sub = sub.loc[sub.index.isin(gene_lengths.index)]
    res = su.sanity_normalization(sub, v_method="EAP", n_threads=threads, output_dir=os.path.join(outdir, "sanity_run"),
                                  gene_lengths=None if gene_lengths is None else gene_lengths.reindex(sub.index))
    avg = su.average_conditions(res, meta, sample_col="sample", cond_col=cond_col, method="reml")
    dump({"result": res, "averages": avg, "metadata": meta}, os.path.join(outdir, "sanity.pkl"))
    log(f"sanity {outdir}: {res.log2_expression.shape}, {res.params.get('units')}, {time.time() - t0:.0f} s")


def step_sanity(out, threads):
    d = os.path.join(out, "inputs")
    meta = pd.read_csv(os.path.join(d, "metadata.tsv"), sep="\t")
    meta = meta.loc[meta["use"]].reset_index(drop=True)
    genes = pd.read_csv(os.path.join(d, "gene_counts.tsv.gz"), sep="\t", index_col=0)
    cs = pd.read_csv(os.path.join(d, f"{CS_LEVEL}_counts.tsv.gz"), sep="\t", index_col=0)
    tot, frac = meta.loc[meta["fraction"] == "Total"], meta.loc[meta["fraction"] != "Total"]
    for level, counts in (("genes", genes), (CS_LEVEL, cs)):
        for part, m in (("Total", tot), ("fractions", frac)) + ((("joint", meta),) if level == CS_LEVEL else ()):
            od = os.path.join(out, "sanity", level, part)
            os.makedirs(od, exist_ok=True)
            run_one_sanity(counts, m, od, threads)


# ---------------------------------------------------------------------------------------------------------------
# fracSanity
def _fit(job):
    import sanity_utils as su
    out_path, ltq_t, var_t, ltq_p, var_p, ids, info, threads = job
    if os.path.exists(out_path):
        return f"{os.path.basename(out_path)}: done before, skipped"
    t0 = time.time()
    fit = su.fit_frac_sanity(ltq_t, var_t, ltq_p, var_p, prior_F=info.get("prior_F"), pi_max=info.get("pi_max", PI_MAX), n_threads=threads,
                             keep_gene_posterior=False, verbose=False)
    fit.update(info)
    fit["gene_ids"] = list(ids)
    dump(fit, out_path)
    return (f"{os.path.basename(out_path)}: F = {fit['F']:.4g} [{fit['F_CI_lower']:.4g}, {fit['F_CI_upper']:.4g}], "
            f"identifiable = {fit['F_identifiable']}, pi = {fit['pi_outlier']:.3f}, n = {len(ids)}, {time.time() - t0:.0f} s")


def run_jobs(jobs, threads):
    workers = max(1, min(len(jobs), threads // 6))
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for msg in ex.map(_fit, [j + (max(1, threads // workers),) for j in jobs]):
            log(msg)


def step_fit_cohorts(out, threads):
    import sanity_utils as su
    jobs = []
    for level in ("genes", CS_LEVEL):
        tot = load(os.path.join(out, "sanity", level, "Total", "sanity.pkl"))
        frac = load(os.path.join(out, "sanity", level, "fractions", "sanity.pkl"))
        genes = tot["averages"].means.index.intersection(frac["averages"].means.index)
        od = os.path.join(out, "fracSanity", level, "cohorts")
        os.makedirs(od, exist_ok=True)
        for cohort in frac["averages"].means.columns:
            inp = su.prepare_frac_sanity_inputs(tot["averages"], frac["averages"], "Total", cohort, genes=genes)
            jobs.append((os.path.join(od, re.sub(r"[^A-Za-z0-9]+", "_", cohort) + ".pkl"), inp["ltq_total"], inp["var_total"],
                         inp["ltq_pd"], inp["var_pd"], inp["gene_ids"], {"cohort": cohort, "level": level}))
    run_jobs(jobs, threads)


def step_fit_replicates(out, threads):
    # each fraction library against the Total average (replicates of fractions and Total are not paired)
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


# ---------------------------------------------------------------------------------------------------------------
# CS usage and tandem APA (joint Sanity run of Total and fractions)
def step_usage(out, threads):
    import sanity_utils as su
    od = os.path.join(out, f"usage_{CS_LEVEL}")
    os.makedirs(od, exist_ok=True)
    joint = load(os.path.join(out, "sanity", CS_LEVEL, "joint", "sanity.pkl"))
    res, meta = joint["result"], joint["metadata"]
    tab = pd.read_csv(os.path.join(out, "inputs", f"{CS_LEVEL}_table.tsv.gz"), sep="\t", index_col=0).loc[res.log2_expression.index]
    cohorts = [c for c in meta["cohort"].unique() if c != "Total"]
    cond_groups = {c: meta.loc[meta["cohort"] == c, "sample"].tolist() for c in meta["cohort"].unique()}

    # relative usage of each CS against the other CS of its gene (all site classes)
    p = os.path.join(od, "cs_usage.per_sample.pkl")
    groups = su.isoform_groups(tab["gene_id"].to_dict())
    if not os.path.exists(p):
        y, s = su.grouped_log2_ratios(res, groups)
        avg = su.average_replicates(y, s, cond_groups, method="reml")
        dump({"groups": groups, "y": y, "sd": s, "averages": avg}, p)
        log(f"usage: {len(groups['group'].unique()):,} CS tests in {tab.loc[groups['feature'].unique(), 'gene_id'].nunique():,} genes")
    for c in cohorts:
        f = os.path.join(od, "cs_usage." + re.sub(r"[^A-Za-z0-9]+", "_", c) + "_vs_Total.tsv.gz")
        if os.path.exists(f):
            continue
        du = su.differential_relative_usage(res, groups, meta, c, "Total", sample_col="sample", cond_col="cohort", method="reml",
                                            n_draws=500)
        du.to_csv(f, sep="\t", index_label="group")
        log(f"usage {c}: {len(du):,} tests, padj < 0.05: {int((du['padj'] < 0.05).sum()):,}")

    # 3' UTR length index of the tandem CS in the PAQR3 terminal exon
    sites = tab.loc[tab["in_terminal_exon"]].reset_index().rename(columns={"cs_id": "feature", "gene_id": "gene"})
    sites = sites[["gene", "feature", "utr_length"]]
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
        du = su.differential_utr_length(res, sites, meta, c, "Total", sample_col="sample", cond_col="cohort", method="reml",
                                        n_draws=500)
        du.to_csv(f, sep="\t", index_label="gene")
        log(f"APA {c}: {len(du):,} genes, {du['call'].value_counts().to_dict()}")
        f2 = f.replace("utr_length.", "proximal_vs_distal.")
        pdd = su.differential_proximal_distal(res, sites, meta, c, "Total", splits="first", sample_col="sample", cond_col="cohort",
                                              method="reml", n_draws=500)
        pdd.to_csv(f2, sep="\t", index_label="group")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--steps", default="prepare,prepare_cs,sanity,fit_cohorts,fit_replicates,usage")
    a = ap.parse_args()
    os.environ.setdefault("OMP_NUM_THREADS", "6")
    cfg = paths()
    os.makedirs(a.outdir, exist_ok=True)
    for step in a.steps.split(","):
        t0 = time.time()
        log(f"=== {step}")
        if step == "prepare":
            step_prepare(a.outdir, cfg)
        elif step == "prepare_cs":
            step_prepare_cs(a.outdir, cfg)
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
