# Status of the current analysis work

Working notes for the ongoing analysis (plan, state of the running jobs, how to resume). Last update: 2026-11-12 01:30.

## Plan

### A. ONT fractions (`ONT_organelles.fracSanity.ipynb`): revised cleavage-site analysis
Requested changes to the first version (2026-11-10):
1. **No clustering of cleavage sites.** Every read 3' end position (`XF`) of a gene (`XT`) is its own CS (`<gene>_<position>`). Nearby positions can be biological (poly(A) tail length, splice isoforms), and a locus usually has one dominant position that should be quantified on its own; minor positions get expression close to the prior after Sanity.
2. **3' UTR length (APA) only on terminal exons (tandem APA).** Terminal exon = PAQR3 segment with the largest segment index of the gene (`APA_AI_data/paqr3_segments_results_atlasv3_min77/paqr3_result_atlasv3_min77_segments.bed`, column 4 = `<gene index>.<segment index>`, segment index increases from the first to the last exon on either strand; gene index → `gene_id` in `..._genes.bed`, column 5). 3' UTR length = distance from the stop codon of the canonical transcript when it is in the terminal exon, otherwise from the 5' end of the terminal exon.
3. **No intronic polyadenylation** for now (the assignment of reads to transcript isoforms by `nanoflowz` is to be checked first).
4. Sanity `v_method='EAP'` stays.

### B. APEX-seq (Fazal et al. 2019, GSE116008, HEK293T) as an external reference
1. Download the FASTQ files with the SRA toolkit on the vscode node (compute nodes have no internet).
2. Process them with the bulk RNA-seq Snakemake workflows (`WF/Snakefile-prepare-faster`, `-quantification-faster`, `-PAQR-quantify`; config as in `NFYA_project.ipynb`).
3. Two notebooks for APEX-seq, as for the ONT data:
   - `APEX_seq_vs_Control.ipynb` (like `ONT_organelles_vs_Total.ipynb`): read QC, nuclear RNA, log2 CPM slopes, depth, library complexity, APEX vs biochemical fractionation;
   - `APEX_seq.fracSanity.ipynb` (like `ONT_organelles.fracSanity.ipynb`): fracSanity, ISMARA gene sets, tandem PAS usage (PAQR), APA, PAS-level fracSanity.
4. Read mapping and gene-biotype content of all ONT libraries (incl. Totals) and APEX-seq libraries, as `pbody_RNA_Spang.ReadMapping.ipynb` (user, 2026-11-11): `ONT_organelles_and_APEX.ReadMapping.ipynb` (builder `build_readmapping_nb.py`). The ONT part has been executed and interpreted; the APEX part runs in `APEX_fracSanity.sbatch`.
5. Cross-comparison of the fracSanity patterns of the ONT fractions and the APEX-seq constructs for matching fractions (e.g. ER ↔ ERM / KDEL, Free ↔ NES): `APEX_vs_ONT_organelles.fracSanity.ipynb`.

### C. U2OS density-gradient fractionation with spike-ins (LoRNA) as a test of fracSanity (user, 2026-11-11)
Data: Villanueva et al., Nat Methods 2024 (ENA PRJEB49479), already processed by the bulk workflow as `U2OS_subcellular_fractions` (`NFYA_project.ipynb`): 8 iodixanol-gradient fractions + whole-cell lysate (`dTotal`), 3 replicates, DMSO only; CORALL total RNA-seq with Lexogen SIRV spike-ins. The paper itself was not readable here (attachment did not come through, Nature page behind login); methods from the bioRxiv preprint: SIRVs "1%" at library preparation, RNA content per fraction from SIRV vs human reads and the share of the fraction used.
1. Spike-in counts: STAR on SIRV + ERCC92 sequences only (`WF/scripts/count_spike_ins.sbatch`, reference `$EXTERNAL_DATA_DIR/spike_ins`; SIRV sequences + annotation from a public CRG mirror of the Lexogen files, ERCC92 from Thermo Fisher). Test on 2 M pairs: 0.23% SIRV, no ERCC; Totals also contain SIRVs.
2. `F_spike` = share of the gradient's RNA molecules per fraction: K = sum_g(c_g / L_g) / SIRV pairs, normalised over the 8 fractions per replicate. Assumes equal SIRV amounts per fraction library; checked against the **mass balance** (NNLS: rho_Total = sum_k F_k rho_k over genes ≥ 10 TPM), which needs no spike-ins.
3. fracSanity per fraction vs Total (uniform prior; prior from the mass balance; tight prior around F_spike) and per library vs the Total of the same replicate; compare F, alpha_g (posterior, CI coverage, r_g), the Beta assumption, sum of alpha over fractions.
4. Compare r_g of the gradient fractions with the ONT fractions and APEX-seq.
Notebook `LoRNA_U2OS.fracSanity.ipynb` (builder `$SCRATCH/build_lorna_nb.py`), pipeline `WF/scripts/LoRNA_fracSanity_pipeline.py`, job `LoRNA_fracSanity.sbatch` (outputs `wf_runs/Sanity_utils_LoRNA_U2OS/`).

**Design choices for APEX-seq.**
- **Samples:** the unperturbed libraries only: 9 constructs × 2 targets and 2 unlabelled controls, the MITO T2 technical replicate (added to T2), and the nuclear and cytosolic fractionation libraries (2 each). This makes 41 runs. Not used: the drug and nocodazole perturbations, and the total-RNA MITO libraries.
- **Reference ("Total"):** the 18 unlabelled controls together, as in the 18-control DESeq2 analysis of the paper. They went through the same streptavidin enrichment.
- **Gene counts:** featureCounts on GENCODE v42 basic, all genes, read pairs, reverse-stranded.
- **PAS quantification:** PAQR `tandem_pas_expression-mean`, rounded as in `NFYA_project.ipynb`. Terminal exons with −1 in any library are left out. `PAQ_coverage_unstranded: "no"`, because TruSeq stranded is read 1 reverse.
- **Library complexity:** without UMIs, estimated from distinct fragment coordinates of first mates (Lander–Waterman).
- **Expression analysis with Sanity_utils, as for ONT, not DESeq2** (user, 2026-11-11). Genes are length-corrected: `sanity_normalization(..., gene_lengths=...)` uses the exonic union length from featureCounts, which gives log2 TPM. ONT is unfragmented (one read per molecule), so ONT CPM corresponds to APEX TPM. The QC notebook also uses log2 TPM. PAS (PAQR mean coverage) are not length-corrected, because coverage is already per nucleotide.
- **ISMARA gene sets (ONT and APEX):** every database fit includes the gene biotypes as covariates, as in `pbody_RNA_Spang.TotalDE.ipynb`: `BIOTYPE:` sets relative to protein-coding genes, plus `mitochondrially_encoded`. The chrM genes are taken out of their own biotypes, so the sets form a partition. The fits also include the deciles of reference expression (`EXPR:`).

## Current status (2026-11-12 14:30)

| step | state |
|:---|:---|
| A1–A3 pipeline code (`WF/scripts/ONT_fracSanity_pipeline.py`, steps `prepare_cs`, `usage` rewritten) | done |
| A: `prepare_cs` | done: 596,094 CS positions, 163,256 with ≥ 5 reads (98.8% of reads) in 23,680 genes; 87,340 in PAQR3 terminal exons; 19,436 of 57.3 M reads off the gene's chromosome/strand left out |
| A: Sanity on `cs_pos` (Total, fractions, joint) | done |
| A: `fit_cohorts` (fracSanity on `cs_pos`), `usage` | done, job 31201273 (4 h 55 min, 21 GB): 160,721 CS usage tests in 13,065 genes; APA on 7,512 genes |
| A: notebook edits | done; copied to `ONT_organelles.fracSanity.ipynb` at 16:00 (backup of the previous notebook: `$SCRATCH/ONT_organelles.fracSanity.before_cs_pos.ipynb`). Old CS/APA figures and tables moved to `_archive_clustered_CS_v1/` |
| A: execute the notebook, write interpretations | done (job 31373322; contamination plot now with the mean $D_g$) |
| B1 downloads | done, 41 runs (pass 1 died at 13:12 when the script was edited while running; the 4 missing runs were fetched by `download_then_workflow.v2.sh`, 16:11) |
| B2 workflow | done 17:11 (prepare 411, quantification 1,992, PAQR 24 steps; `workflow_head.log`). PAQR: 11,219 tandem PAS in 3,716 terminal exons, as in the other runs (10–13k); 1,835 exons have −1 (not quantified) in at least one library (most in APEX-NIK-T1, OMM-C2, MITO-T2tech) and are left out → 6,034 PAS in 1,881 exons |
| B3 pipeline + notebooks | done: job 31372648 (3 h 24 min, 14 GB; all 4 notebooks status 0). Bug fixed afterwards: the QC and ReadMapping notebooks halved the STAR read numbers (STAR / `mapping_stats` count pairs, not mates); builders and notebooks patched, both re-executed by job 31412850 |
| B interpretations | done for all four APEX notebooks (`$SCRATCH/interpretations_apex.py`, inserted into the executed notebooks); QC and ReadMapping re-executed after the read-pair fix by job 31412850 (31409277 failed: `--wrap` runs `sh`) |
| C1 spike-in counts | done (all 27 libraries; the 5 corrupt raw FASTQs counted from `fastq_sanitized`, array 31408749) |
| C2 pipeline `prepare`, `sanity`, `mass_balance` | done (Sanity: 41,976 genes in Total, 48,998 in fractions, log2 TPM). Mass balance (8,494 genes ≥ 10 TPM): F3 0.20, F5 0.12, F6 0.36, F7 0.37, F1 / F2 / F8 at 0 (collinear fractions; per replicate unstable); sum 1.05 |
| C3 fits | done: 31375302 (fit_cohorts, fit_replicates, fit_cohorts_mbprior), 31410284 (fit_cohorts_lowpi), 31408750 (spike_ins, spike-prior fits) |
| C4 troubleshooting steps (`pooled`, `mass_balance_pooled`, `simulate`) | done: job 31446656 (1 h 11 min) and 31471167. Pooled libraries = per-replicate binomial thinning of the fraction libraries in proportion to F_spike (Membranes 1–3, Nuclear 4–5, Cytosol 6–8), Sanity + fracSanity (default and low pi) vs Total; pooled mass balance; simulation with known F (real α^T shape vs fitted Beta) for 8 fractions + 3 pools |
| C5 authors' RNA content (`authors_content`, `fit_cohorts_authprior(_lowpi)`, `pooled_auth`, `mass_balance_pooled_auth`) | done: job 31507929 (2 h 11 min, 5 GB). Notebook rebuilt with F_auth as the reference, the section on the authors' method vs equal spike-ins before RNA isolation, the mass→molecule conversion (`spike_ins/F_spike_mass_vs_molecules.tsv`), executed by job 31535707, interpretations and summary rewritten (`scratchpad/patch_lorna_interp.py` → builder; backups `$SCRATCH/build_lorna_nb.before_interp_auth.py`, `$SCRATCH/LoRNA_U2OS.fracSanity.executed_31535707_no_interp.ipynb`) (previous outputs archived in `tables|figures/gene_expression_analysis/LoRNA_U2OS_fracSanity_archive_equal_SIRV_2026-11-12/`, executed notebook and builder in `$SCRATCH/*equal_SIRV*`) |
| D pairwise slopes (user, 2026-11-12: does averaging the reference bias the buffering slopes?) | done: section "Pairwise comparisons" added to `ONT_organelles_vs_Total.ipynb` and `APEX_seq_vs_Control.ipynb` (`$SCRATCH/add_pairwise_buffering.py`; APEX also in `build_apex_qc_nb.py` + `interpretations_apex.py`), executed by job 31510170 (`QC_notebooks.sbatch`). ONT: vs single Total all slopes −0.022 to −0.030 lower, ranking and MALAT1 association unchanged. APEX: single controls inconsistent (control vs control 0.62–1.62); matched-control slopes −0.11 lower, complexity association −0.46 (low complexity → wider). The pbody `SanityVmethod` complexity analysis was already per replicate pair. Sanity cohort averages (log means) sum 0.01–0.11 log2 below 10^6; offset fraction vs Total ≤ 0.07 log2 (affects F scale ≤ 5%, not slopes) |
| C notebook | done: 37 cells (F_auth version, see C5); the equal-SIRV version was executed by 31473189. Previous execution kept as `$SCRATCH/LoRNA_U2OS.fracSanity.executed_31471167.ipynb` |
| README, `make_html_reports.py` (NOTEBOOKS list) | updated (incl. LoRNA) |

`$SCRATCH` = `/scicore/home/zavolan/GROUP/APA_localization/temp_dir/notebook_builders_2026-11` (persistent copy of the notebook builders, the edited ONT notebook and its backup; working copies in the session scratchpad `/scratch/claude-100158/.../scratchpad`). The builders write to the path given as argument; `make_apex_fracsanity_builder.py` reads `build_fracsanity_nb.py` from its own folder.

## Findings that change the plan

- **LoRNA spike-ins: corrected (2026-11-12, after the user's check of the paper).** My first claim ("K constant by design, yields not public") was wrong. The published methods: 1 µg RNA per library, SIRVs added after rRNA depletion at ~1% of the measured depleted RNA (0.02–0.12 ng, not the same amount per library). The authors convert to RNA content with ng SIRVs added × SIRV read share × total RNA of the fraction; their processed table (GitHub LoRNA_UPR `sirv_correction.rds`, commit 054f938) is in `metadata/LoRNA_UPR_Villanueva2024/` (parsed without running R). F_auth per compartment: Membranes 0.18, Nuclear 0.12, Cytosol 0.69, matching the pooled mass balance normalised to 1 (0.18 / 0.13 / 0.69). The read-only F^K (equal SIRV amounts assumed) was wrong; it is kept as a comparison. F_auth is now the reference (pipeline step `authors_content`; new fits `cohorts_authprior(_lowpi)`, pools `pooled_auth`).
- **fracSanity fails on F because of the Beta shape of the α prior (2026-11-12).** Fractions with an enriched minority over a background (gradient 1–5, Membranes, Nuclear) have a right-skewed log α, which no Beta can produce → Gamma limit, F → 0, enriched genes become outliers. Narrow cytosolic fractions (6–8) put the edge of D_g at α = 1 → F 1.8–5× too high vs F_auth (0.64 / 0.55 / 0.69 vs 0.12 / 0.31 / 0.26; sum 1.88 vs 0.69), and the data reject F_auth by 25–437 nats, so a tight prior does not hold. Simulation with the same data and known F: Beta-shaped α → F recovered (10/11 in the 95% CI, coverage 0.91–1.00); real α shape → the real failure exactly. Only the complete cytosol pool (weighted by F_auth) gives the right F: 0.684 / low pi 0.637 vs F_auth 0.692, mass balance 0.729 (the F^K-weighted pool gave 0.82). With a tight F_auth prior in fractions 1–5 (likelihood shallow there, 8–15 nats) the absolute α is nearly right (bias −0.07 to −0.25 log2, coverage 0.79–0.95), but r_g still shrunk (slope 0.20–0.50); with low pi the likelihood pulls F to 0.014–0.036. Gene rankings are right everywhere (Spearman 0.97–0.99). Low pi removes the shrinkage (slope 0.2–0.5 → 0.61–0.67) but not the F error. Same caveat for the ONT and APEX F values. Proposed Sanity_utils changes: Beta-mixture or free-α_max prior; F from a mass balance; low default pi_max / posterior from the non-outlier component.
- **Gradient ↔ ONT ↔ APEX:** the gradient membrane set (top 5% α(1–3)/α(6–8), 343 genes; secreted / ER proteins) is enriched in APEX OMM 0.91 / ERM 0.83 and ONT Golgi / Endo / TJ / PEX / ER 0.65–0.72 (AUC), depleted in ONT Free 0.15–0.23. Fractions 7–8 anticorrelate with APEX ERM (Spearman −0.57 / −0.67), partly via transcript length.

- **fracSanity discards the most strongly localised genes (APEX-seq).** With `pi_max` 0.5, the contaminated-Gaussian error model flags the tails of $D_g$ as outliers (error variance $V_g + \tau^2$, $\tau$ = 2–4), and their posterior falls back to the Beta: mitochondrially encoded genes in APEX MITO (+5 log2 above the median gene) get $r_g$ = 0.007, snoRNAs in the nuclear fraction (+4) 0.011, MALAT1 in NLS −0.02. 40% of genes are "outliers" in the nuclear fraction. Also, with $F \to 0$ the Beta becomes a scaled Gamma with light log-tails, which compresses the rest. Tested on LoRNA with the low-pi fits (`pi_max` 0.02); possible fixes in Sanity_utils: per-gene posterior from the non-outlier component, a heavier-tailed $\alpha$ prior, or $F$ from spike-ins.
- **APEX APA follows library complexity** (median $D_g$ vs log2 complexity: Spearman 0.80 over 9 constructs), not nuclear RNA (0.18): not interpretable as compartment APA without de-duplication.
- **ONT ↔ APEX:** ONT Golgi / PEX / Endo / TJ share their over- and under-represented mRNAs with APEX ERM / OMM (Spearman up to 0.44): ER-translated mRNAs co-sort with these membranes.

## Next steps (in order)

1. Decide with the user how to change fracSanity (α prior that allows an enriched minority, F from a mass balance, default pi_max), then re-fit LoRNA (test: F against the pooled mass balance, and the simulation with the real α shape) and the ONT / APEX fractions.
2. (done 2026-11-12) The authors' per-fraction RNA content is used as the reference F.

## Outputs

| what | where |
|:---|:---|
| ONT pipeline outputs | `wf_runs/Sanity_utils_ONT_organelles_siCtrl/`. New CS level: `inputs/cs_pos_*`, `sanity/cs_pos`, `fracSanity/cs_pos`, `usage_cs_pos`. The old clustered level (`cs`, `usage`) is kept, unused |
| APEX FASTQ | `input_fastq/APEX_seq_Fazal2019/<SRR>_{1,2}.fastq.gz` (temporary files in `/scratch/mirono0000/sra_download`) |
| APEX metadata | `metadata/APEX_seq_Fazal2019.{GSE116008_runs.tsv,selected_runs.txt,sample_annotation.tsv}`, `metadata/start_samples.APEX_seq_Fazal2019.tsv` |
| APEX workflow | `wf_runs/APEX_seq_Fazal2019/` (`config.yaml`, `output/`) |
| APEX pipeline outputs | `wf_runs/Sanity_utils_APEX_seq_Fazal2019/` |
| HTML reports (all 7 notebooks, without code) | `reports/html_2026-11-12/index.html` (before the F_auth / pairwise updates); `reports/html_2026-11-12b/index.html` (current) |
| LoRNA pipeline outputs | `wf_runs/Sanity_utils_LoRNA_U2OS/` (new: `inputs/pooled_*`, `inputs/mass_balance_pooled.tsv`, `sanity/genes/pooled`, `fracSanity/genes/{pooled,pooled_lowpi,sim_alpha_T,sim_beta}`; F_auth: `inputs/F_authors.tsv`, `inputs/pooled_auth_*`, `inputs/mass_balance_pooled_auth.tsv`, `sanity/genes/pooled_auth`, `fracSanity/genes/{cohorts_authprior,cohorts_authprior_lowpi,pooled_auth,pooled_auth_lowpi}`) |
| notebook tables and figures | `tables/` and `figures/gene_expression_analysis/{ONT_organelles_fracSanity,APEX_seq_vs_Control,APEX_seq_fracSanity,APEX_vs_ONT_organelles_fracSanity,LoRNA_U2OS_fracSanity}/` |

(`wf_runs`, `input_fastq`, `metadata`, `tables` and `figures` are under `/scicore/home/zavolan/GROUP/APA_localization/`.)

## Notes for long-running jobs

- The vscode node has 6 cores and 10 GB for the user (cgroup), and compute nodes have no internet: downloads and the Snakemake head run on the vscode node (detached with `nohup setsid`), and everything heavy runs on SLURM.
- Notebooks with gene-set fits: BLAS limited to 6–8 threads (set in the sbatch files); more threads were much slower.
- Every pipeline step skips if its outputs exist, so a failed job can be resubmitted as is. Nothing existing is overwritten; new analysis levels get new names.
- `download_sra_fastq.sh` can be resumed: it skips finished runs, reuses `.sra` files and completed `fasterq-dump` outputs, and clears stale prefetch locks.
- Never `mv` into an archive folder that already holds files of the same name (it happened once, 17:35, to `apa.median_D_vs_nuclear_contamination`; the clustered-CS version was regenerated from the archived tables).
