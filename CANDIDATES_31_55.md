# Extended-BioAgentBench: Task Candidates 31-55

25 gold candidates for tasks 31-55 (task 31-36 assigned to this session, 37-55 assigned to parallel sessions).

## Session Assignments

| Session | Tasks | Candidates |
|---------|-------|------------|
| **This session** | 35-40 | #1, #2, #6, #8, #9, #15 |
| Session B | 41-45 | #3, #4, #5, #7, #10 |
| Session C | 46-50 | #11, #12, #13, #14, #16 |
| Session D | 51-55 | #17, #18, #19, #20, #21 |
| Session E | 56-59 | #22, #23, #24, #25 |

---

## MEDIUM (5 candidates, depth 5-7, 2 convergence points)

### C1. Small RNA-seq / miRNA Discovery
- **Domain:** Non-coding RNA
- **Depth:** 7 | **Convergence:** 2 | **Tools:** 9
- **DAG:** fastp → bowtie(contam) → [bowtie(mature) || bowtie(hairpin)] → mirtop CONVERGE → [miRDeep2(novel) || counts] → edgeR CONVERGE → report
- **Conda tools:** `fastp`, `bowtie`, `samtools`, `mirtop`, `mirdeep2`, `mirtrace`, `bioconductor-edger`, `multiqc`
- **Data:** nf-core/smrnaseq test (~50MB)
- **Runtime:** ~30 min
- **Key trap:** Adapter trimming is critical — miRNAs are 18-25nt, adapters dominate read length

### C2. CUT&RUN Epigenomic Profiling
- **Domain:** Epigenomics (antibody-free chromatin)
- **Depth:** 7 | **Convergence:** 2 | **Tools:** 10
- **DAG:** fastp → [bowtie2(genome) || bowtie2(spike-in)] → dedup → [SEACR(IgG-ctrl) || MACS2(statistical)] CONVERGE → consensus → [deeptools heatmap || HOMER motif] CONVERGE → report
- **Conda tools:** `fastp`, `bowtie2`, `samtools`, `picard`, `seacr`, `macs2`, `bedtools`, `deeptools`, `multiqc`
- **Data:** nf-core/cutandrun test (~100MB)
- **Runtime:** ~45 min
- **Key trap:** Spike-in normalization (E. coli DNA) — must align to both genomes

### C3. GC-MS Metabolomics Profiling
- **Domain:** Metabolomics
- **Depth:** 8 | **Convergence:** 2 | **Tools:** 7
- **DAG:** mzML → MSnbase → XCMS(peaks) → [XCMS(align) || XCMS(group)] CONVERGE → fillPeaks → [RAMClustR(deconv) || RIAssigner(RI)] CONVERGE → matchms(library) → stats
- **Conda tools:** `bioconductor-xcms`, `bioconductor-msnbase`, `bioconductor-camera`, `r-ramclustr`, `matchms`
- **Data:** Galaxy GC-MS tutorial (Zenodo 3244991, ~150MB)
- **Runtime:** ~40 min
- **Key trap:** Retention time alignment must happen before peak grouping

### C4. RADseq Population Genetics
- **Domain:** Population Genomics / Ecology
- **Depth:** 9 | **Convergence:** 2 | **Tools:** 7
- **DAG:** process_radtags(demux) → ustacks(per-sample) → cstacks(catalog) CONVERGE → sstacks(rematch) → gstacks(genotype) → populations CONVERGE → [VCF || STRUCTURE exports]
- **Conda tools:** `stacks`, `fastqc`, `samtools`, `bwa`, `vcftools`, `plink`
- **Data:** Stacks tutorial (~200MB)
- **Runtime:** ~1h
- **Key trap:** N-to-1-to-N-to-1 pattern — per-sample loci must converge to catalog, then re-diverge

### C5. DIA Proteomics Quantification
- **Domain:** Proteomics (data-independent acquisition)
- **Depth:** 7 | **Convergence:** 2 | **Tools:** 7
- **DAG:** RAW → ThermoRawFileParser → DIA-NN/OpenSWATH → pyprophet(FDR) CONVERGE(library+search) → TRIC(alignment) → MSstats(quant) CONVERGE(design+data) → report
- **Conda tools:** `thermorawfileparser`, `openms`, `pyprophet`, `msstats`
- **Data:** PRIDE PXD014414 (~500MB subset)
- **Runtime:** ~2h
- **Key trap:** Spectral library generation is a separate branch that feeds into search

---

## MEDIUM-HARD (10 candidates, depth 6-10, 3 convergence points)

### C6. RNA Fusion Detection (Cancer)
- **Domain:** Cancer Diagnostics
- **Depth:** 8 | **Convergence:** 3 | **Tools:** 10
- **DAG:** fastp → STAR → [Arriba || STAR-Fusion || FusionCatcher] CONVERGE(fusion-report) → FusionInspector(validate) CONVERGE(BAM+fusions) → [StringTie || Salmon] CONVERGE(expression+fusions) → report
- **Conda tools:** `star`, `arriba`, `star-fusion`, `fusioncatcher`, `salmon`, `stringtie`, `fastp`, `samtools`, `multiqc`
- **Data:** nf-core/rnafusion test (~500MB)
- **Runtime:** ~2.5h
- **Key trap:** Three independent fusion callers with different output formats must be reconciled; FusionInspector needs BOTH the merged list AND the original BAM

### C7. Circular RNA Discovery + miRNA Target
- **Domain:** Non-coding RNA / Cancer
- **Depth:** 8 | **Convergence:** 3 | **Tools:** 11
- **DAG:** TrimGalore → STAR(chimeric) → [CIRCexplorer2 || CIRIquant || find_circ || DCC || circRNA_finder] CONVERGE(≥2 tools) → [psirc-quant(counts) || getfasta(seqs)] → [miRanda(targets) || CircTest(DE)] CONVERGE → report
- **Conda tools:** `star`, `circexplorer2`, `ciriquant`, `trim-galore`, `miranda`, `bioconductor-deseq2`, `samtools`
- **Data:** nf-core/circrna test (~400MB)
- **Runtime:** ~1.5h
- **Key trap:** Width-5 fan-out at BSJ detection; consensus requires ≥2 tools agreeing

### C8. Immune Receptor Repertoire (BCR/TCR)
- **Domain:** Immunology / Clinical
- **Depth:** 10 | **Convergence:** 3 | **Tools:** 10
- **DAG:** fastp → pRESTO(Filter→Mask→Pair→Cluster→Consensus→Assemble→Collapse→Split) → IgBLAST → Change-O(MakeDB→ParseDB) → [SHazaM(threshold) || TIgGER(alleles)] CONVERGE → SCOPer(clones) → [Dowser(trees) || Alakazam(diversity)] CONVERGE → EnchantR CONVERGE → report
- **Conda tools:** `fastp`, `presto`, `igblast`, `changeo`, `r-shazam`, `r-tigger`, `r-scoper`, `r-dowser`, `r-alakazam`
- **Data:** nf-core/airrflow test (~300MB); Immcantation examples
- **Runtime:** ~1h
- **Key trap:** 8-step sequential pRESTO chain; threshold + novel allele branches must converge before clonal assignment

### C9. Ribosome Profiling + Translational Efficiency
- **Domain:** Translation / Gene Regulation
- **Depth:** 9 | **Convergence:** 3 | **Tools:** 10
- **DAG:** [Ribo-seq: fastp→SortMeRNA→STAR→filter(28-32nt)] || [RNA-seq: fastp→STAR→Salmon] → [riboWaltz || Ribo-TISH] CONVERGE(P-site) → [Ribotricer || RiboCode || ORFquant] CONVERGE(ORFs) → anota2seq CONVERGE(Ribo+RNA) → report
- **Conda tools:** `fastp`, `bowtie2`, `star`, `salmon`, `sortmerna`, `samtools`, `multiqc`
- **Data:** nf-core/riboseq test (~500MB)
- **Runtime:** ~2h
- **Key trap:** Dual-input (Ribo-seq + RNA-seq); P-site offset calibration is domain-specific; rRNA removal is critical (80-90% of Ribo-seq reads are rRNA)

### C10. Quantitative Proteomics (DDA-LFQ, dual search engine)
- **Domain:** Proteomics / Clinical
- **Depth:** 10 | **Convergence:** 3 | **Tools:** 7
- **DAG:** RAW → ThermoRawFileParser → PeakPicker → [Comet || MSGF+] CONVERGE(PeptideIndexer) → PSMFeatureExtractor → Percolator(FDR) → [Luciphor(PTM) || ProteomicsLFQ(quant)] CONVERGE → MSstats CONVERGE(design) → report
- **Conda tools:** `thermorawfileparser`, `openms`, `comet-ms`, `msgf_plus`, `percolator`, `msstats`
- **Data:** nf-core/quantms test; PRIDE PXD000001 (~300MB)
- **Runtime:** ~2h
- **Key trap:** Dual search engines with different scoring — must be properly merged before FDR; PTM localization parallel to quantification

### C11. MHC Immunopeptidomics
- **Domain:** Cancer Immunology
- **Depth:** 10 | **Convergence:** 3 | **Tools:** 7
- **DAG:** RAW → ThermoRawFileParser → PeakPicker → Comet(no enzyme) → PSMFeatureExtractor → [Percolator || DeepLC(RT) || MS2PIP(spectra)] CONVERGE(rescoring) → IDFilter → MapAligner → FeatureFinder → report
- **Conda tools:** `thermorawfileparser`, `openms`, `comet-ms`, `percolator`, `deeplc`, `ms2pip`
- **Data:** nf-core/mhcquant test (~200MB)
- **Runtime:** ~1.5h
- **Key trap:** No enzymatic cleavage specificity — search space is enormous; ML rescoring integrates 3 tools

### C12. HiCAR Multi-omic Chromatin Interaction
- **Domain:** 3D Genomics / Epigenomics
- **Depth:** 10 | **Convergence:** 3 | **Tools:** 11
- **DAG:** FASTQ(R1=prox,R2=ATAC) → cutadapt → bwa → [pairtools(pairs) || samtools(R2→MACS2 peaks)] CONVERGE(pairs+peaks) → cooler(matrix) → [MAPS(interactions) || cooltools(compartments)] → [edgeR(diff) || ChIPpeakAnno] CONVERGE → viz
- **Conda tools:** `bwa`, `cutadapt`, `pairtools`, `cooler`, `macs2`, `bioconductor-edger`, `bioconductor-chippeakanno`, `samtools`
- **Data:** nf-core/hicar test (~600MB)
- **Runtime:** ~2h
- **Key trap:** Same reads serve dual purpose (proximity ligation + accessibility); must split by mate role

### C13. Multi-classifier Taxonomic Profiling
- **Domain:** Clinical Metagenomics
- **Depth:** 8 | **Convergence:** 3 | **Tools:** 10+
- **DAG:** fastp → BBDuk → bowtie2(host) → [Kraken2→Bracken || MetaPhlAn4 || mOTUs || Centrifuge || DIAMOND] CONVERGE(Taxpasta) → [Krona(viz) || diversity(stats)] CONVERGE → report
- **Conda tools:** `fastp`, `bowtie2`, `kraken2`, `bracken`, `metaphlan`, `centrifuge`, `kaiju`, `diamond`, `taxpasta`
- **Data:** nf-core/taxprofiler test (~500MB + DBs)
- **Runtime:** ~2h
- **Key trap:** Width-5 classifier fan-out; each tool has different DB format and output schema; Taxpasta standardization is the key convergence

### C14. 16S Amplicon Microbiome (full DADA2→diversity)
- **Domain:** Microbiome / Ecology
- **Depth:** 9 | **Convergence:** 3 | **Tools:** 8
- **DAG:** cutadapt(primers) → DADA2(filter→learn→denoise→merge→chimera) → [classify(taxonomy) || MAFFT→FastTree(phylogeny)] CONVERGE → [PICRUSt2(function) || UniFrac(diversity)] CONVERGE → ANCOM-BC(DA) CONVERGE(metadata) → report
- **Conda tools:** `cutadapt`, `bioconductor-dada2`, `mafft`, `fasttree`, `picrust2`, `biom-format`
- **Data:** QIIME2 Moving Pictures (~50MB); nf-core/ampliseq test
- **Runtime:** ~1h
- **Key trap:** Primer removal is make-or-break (wrong primers = 0 ASVs); UniFrac requires BOTH tree and ASV table

### C15. Functional Screening (AMR + BGC)
- **Domain:** AMR / Natural Products
- **Depth:** 8 | **Convergence:** 3 | **Tools:** 10
- **DAG:** contigs → Prodigal → [ABRicate || AMRFinderPlus || fARGene] CONVERGE(argNorm) + [antiSMASH || DeepBGC || GECCO] CONVERGE(comBGC) → hAMRonization CONVERGE(ARG+BGC) → report
- **Conda tools:** `prodigal`, `abricate`, `amrfinderplus`, `fargene`, `antismash`, `deepbgc`, `gecco`, `hmmer`
- **Data:** nf-core/funcscan test (~300MB)
- **Runtime:** ~1.5h
- **Key trap:** Two parallel nested fan-outs (3 ARG tools + 3 BGC tools), each with internal convergence, then outer convergence

---

## HARD (10 candidates, depth 8-12, 4+ convergence points)

### C16. Somatic Variant Calling (Tumor-Normal)
- **Domain:** Cancer Genomics (Clinical)
- **Depth:** 12 | **Convergence:** 4 | **Tools:** 14
- **DAG:** [T: fastp→BWA→dedup→BQSR] || [N: fastp→BWA→dedup→BQSR] CONVERGE(shared-sites) → [Mutect2 || Strelka2 || FreeBayes || Manta(SV)] CONVERGE(multi-caller) → [FilterMutect || ASCAT(CNA) || MSIsensor] CONVERGE(variant+CNA+MSI) → [VEP || snpEff] CONVERGE(dual-annot) → report
- **Conda tools:** `bwa-mem2`, `gatk4`, `samtools`, `strelka2`, `freebayes`, `manta`, `bcftools`, `ensembl-vep`, `snpsift`, `ascat`, `msisensor-pro`, `mosdepth`, `fastp`, `multiqc`
- **Data:** nf-core/sarek test (~700MB)
- **Runtime:** ~3h
- **Key trap:** T+N mandatory paired processing; Mutect2 needs 3-step internal chain (LearnROM+GetPileup+Filter); 4 convergence points with nested diamonds

### C17. Neoantigen Prediction (Triple-input)
- **Domain:** Cancer Immunotherapy
- **Depth:** 11 | **Convergence:** 4+ | **Tools:** 12
- **DAG:** [T-WES → BWA → GATK] || [N-WES → BWA → GATK] || [T-RNA → STAR → Salmon] → Mutect2 CONVERGE(T+N) → VEP → [OptiType(HLA-DNA) || HLA-HD(HLA-RNA)] CONVERGE(HLA) → pVACseq CONVERGE(variants+HLA+expression) → [CNVkit || MiXCR(TCR)] CONVERGE(TME) → scoring
- **Conda tools:** `bwa`, `gatk4`, `star`, `salmon`, `optitype`, `pvactools`, `ensembl-vep`, `cnvkit`, `samtools`
- **Data:** nextNEOpi test (~900MB)
- **Runtime:** ~3.5h
- **Key trap:** Triple-input (tumor DNA + normal DNA + tumor RNA); 5 convergence points; HLA typing from both DNA and RNA

### C18. Single-Cell RNA-seq (Full Pipeline)
- **Domain:** Single-cell Genomics
- **Depth:** 12 | **Convergence:** 4 | **Tools:** 10
- **DAG:** FASTQ → [STARsolo || Alevin-fry || Kallisto+BUS] CONVERGE(count matrix) → emptyDrops → scanpy(QC→norm→HVG→PCA) → [Harmony(batch) || Scrublet(doublets)] CONVERGE → UMAP→Leiden → [CellTypist(types) || PAGA(trajectory)] CONVERGE → pseudobulk-DE CONVERGE(metadata) → report
- **Conda tools:** `star`, `salmon`, `kallisto`, `bustools`, `scanpy`, `scvi-tools`, `scrublet`, `celltypist`, `harmonypy`
- **Data:** 10x PBMC 3k (~100MB); nf-core/scrnaseq test
- **Runtime:** ~2.5h
- **Key trap:** 3 parallel quantification engines; batch correction + doublet detection convergence; cell typing + trajectory convergence

### C19. Clinical Metaproteomics
- **Domain:** Clinical Proteomics / Microbiome
- **Depth:** 11 | **Convergence:** 4 | **Tools:** 8
- **DAG:** RAW → msconvert → PeakPicker → SearchGUI(X!Tandem+Comet+OMSSA) CONVERGE → PeptideShaker → [Unipept(taxonomy) || filtering || verification] CONVERGE → [FlashLFQ || metaQuantome-prep] CONVERGE → [metaQuantome(function) || metaQuantome(taxonomy) || metaQuantome(interaction)] CONVERGE → report
- **Conda tools:** `openms`, `searchgui`, `peptideshaker`, `r-base`, `bioconductor-msnbase`
- **Data:** Galaxy metaproteomics tutorial (Zenodo ~400MB)
- **Runtime:** ~3h
- **Key trap:** Multi-search-engine convergence; three parallel metaQuantome modes converge

### C20. Somatic + Germline Dual Analysis
- **Domain:** Hereditary Cancer / Clinical
- **Depth:** 12 | **Convergence:** 4 | **Tools:** 10
- **DAG:** [T → BWA → GATK-preprocess] || [N → BWA → GATK-preprocess] CONVERGE(shared-BQSR) → [Mutect2(somatic) || HaplotypeCaller→GVCF→GenotypeGVCFs(germline)] → [FilterMutect || VQSR(SNP)+VQSR(Indel)] CONVERGE → bcftools-isec CONVERGE(somatic∩germline) → [VEP || SnpSift] CONVERGE → clinical-report
- **Conda tools:** `bwa-mem2`, `gatk4`, `samtools`, `bcftools`, `ensembl-vep`, `snpsift`, `mosdepth`, `fastp`, `picard`, `multiqc`
- **Data:** nf-core/sarek test; GATK bundle
- **Runtime:** ~3.5h
- **Key trap:** Somatic-germline intersection is unique convergence; dual VQSR models

### C21. Spatial Transcriptomics (Visium)
- **Domain:** Spatial Biology / Oncology
- **Depth:** 11 | **Convergence:** 4 | **Tools:** 7
- **DAG:** [FASTQ + H&E + slide-layout] CONVERGE(SpaceRanger) → scanpy(QC→norm→PCA) → [Leiden || squidpy(spatial-autocorr)] CONVERGE → [cell2location(deconv) || CellPhoneDB(LR)] CONVERGE → scVI-spatial(integration) CONVERGE(replicates) → report
- **Conda tools:** `scanpy`, `squidpy`, `cell2location`, `scvi-tools`, `anndata`
- **Data:** 10x Visium demo (~500MB); nf-core/spatialvi test
- **Runtime:** ~2.5h
- **Key trap:** Multi-modal input (FASTQ + image + spatial); SpaceRanger is proprietary (agent must find alternative)

### C22. CRISPR Screen (Dual-Library)
- **Domain:** Functional Genomics / Cancer
- **Depth:** 10 | **Convergence:** 4 | **Tools:** 9
- **DAG:** [plasmid-FASTQ || treated-FASTQ || control-FASTQ] → cutadapt → MAGeCK-count CONVERGE(all-samples) → [MAGeCK-test(RRA) || MAGeCK-mle(MLE) || BAGEL2] CONVERGE(multi-method) → MAGeCKFlute(QC+CNV) → [GSEA(pathways) || VISPR(viz)] CONVERGE → hit-ranking CONVERGE(CNV-data) → report
- **Conda tools:** `mageck`, `mageck-vispr`, `cutadapt`, `fastqc`, `bowtie2`, `samtools`
- **Data:** MAGeCK test data; DepMap
- **Runtime:** ~1h
- **Key trap:** Three sample types must converge for counting; three analysis algorithms converge

### C23. LC-MS Untargeted Metabolomics
- **Domain:** Metabolomics / Clinical Chemistry
- **Depth:** 10 | **Convergence:** 4 | **Tools:** 7
- **DAG:** [study + QC + blanks] CONVERGE(MSnbase) → XCMS(peaks) → [XCMS(align) || XCMS(group)] CONVERGE → fillPeaks → [CAMERA(annotation) || blank-subtract || QC-filter] CONVERGE → [SIRIUS(structure) || matchms(library)] CONVERGE → MetaboAnalyst(pathways) → report
- **Conda tools:** `bioconductor-xcms`, `bioconductor-camera`, `bioconductor-msnbase`, `matchms`, `sirius-ms`
- **Data:** Workflow4Metabolomics (Zenodo 3757956, ~300MB)
- **Runtime:** ~2h
- **Key trap:** QC samples and blanks create parallel tracks; 4 convergence points

### C24. Variant Annotation + Clinical Interpretation (Trio)
- **Domain:** Rare Disease / Clinical
- **Depth:** 12 | **Convergence:** 4 | **Tools:** 10
- **DAG:** [proband → BWA → GATK] || [mother → BWA → GATK] || [father → BWA → GATK] → GenomicsDBImport CONVERGE(joint-genotyping) → [VQSR(SNP) || VQSR(Indel)] CONVERGE → [VEP || CADD || SpliceAI] CONVERGE(multi-score) → genmod CONVERGE(pedigree) → scout-ready VCF
- **Conda tools:** `bwa-mem2`, `gatk4`, `samtools`, `ensembl-vep`, `bcftools`, `genmod`, `picard`, `mosdepth`
- **Data:** GATK bundle; GIAB NA12878 trio
- **Runtime:** ~3.5h
- **Key trap:** Trio joint genotyping; dual VQSR; 3 annotation tools; pedigree-aware inheritance scoring

### C25. MAG Recovery (Hybrid Assembly)
- **Domain:** Environmental Metagenomics
- **Depth:** 12 | **Convergence:** 4 | **Tools:** 14
- **DAG:** [short: fastp→bowtie2(host)] || [long: Porechop→minimap2(host)] → [MEGAHIT || SPAdes-hybrid] CONVERGE → bowtie2(map-back) → [MetaBAT2 || MaxBin2 || CONCOCT || SemiBin2] CONVERGE(DAS_Tool) → [CheckM2 || GTDB-Tk || geNomad] CONVERGE → Prokka/Bakta CONVERGE → report
- **Conda tools:** `fastp`, `megahit`, `spades`, `bowtie2`, `minimap2`, `metabat2`, `maxbin2`, `concoct`, `das_tool`, `busco`, `gtdbtk`, `checkm2`, `quast`
- **Data:** nf-core/mag test; CAMI challenge
- **Runtime:** ~3h
- **Key trap:** Dual-input (short+long reads); 4 parallel binners → DAS_Tool merge; QC+taxonomy+viral screening converge; most complex nested diamonds
