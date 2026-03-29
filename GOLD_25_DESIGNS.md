
> **WARNING: COLLISION ALERT** (added by Session A)
> These task IDs already exist from the CANDIDATES phase. Sessions F-J MUST SKIP them:
> - `circrna-detection` — ALREADY BUILT (skip in Session F)
> - `dda-lfq-proteomics` — ALREADY BUILT (skip in Session I)
> - `mag-recovery` — ALREADY BUILT (skip in Session H)
> - `immune-repertoire` — BLOCKED by NCBI FTP (Session F should try, but may fail)

# 25 GOLD Bioinformatics Benchmark Task Designs

> All tasks: HARD-level DAG complexity (depth 8+, 3+ convergence points, nested diamonds), 8-16 CLI tools, real public data (<1 GB), runtime <4h on 8 CPUs, all tools conda-installable.

---

## HIGH PRIORITY (Clinical/Translational): Tasks 1-8

---

### TASK 1: Hi-C 3D Genome Conformation Analysis

**NOTE:** `hicar-chromatin` exists but covers HiCAR (multi-omic Hi-C + ATAC hybrid). This task covers STANDARD Hi-C (pure proximity ligation), which is a fundamentally different protocol with different tools (no ATAC component, full contact matrix analysis, A/B compartments, TADs, loops). No overlap.

**Name:** `hic-3d-conformation`
**Description:** Analyze Hi-C chromatin interaction data to map 3D genome organization, identify compartments, TADs, and significant chromatin loops.

**DAG Diagram (depth=10, convergence=4, tools=10):**
```
 R1.fastq.gz    R2.fastq.gz
     │               │
     └───────┬───────┘
             │
        [fastp] ─────────────────────────────── Level 1
             │
     ┌───────┴───────┐
     │               │
 [bwa-mem2        [bwa-mem2                     Level 2
  align R1]        align R2]
     │               │
 [samtools         [samtools                    Level 3
  view/filter]      view/filter]
     │               │
     └───────┬───────┘
             │
     [pairtools parse]  ◄── [CONVERGENCE 1]     Level 4
             │              (mate pairing)
     [pairtools sort]                            Level 5
             │
     ┌───────┼───────────────┐
     │       │               │
 [pairtools [cooler       [pairtools              Level 6
  dedup]     cload]        stats]
     │       │               │
     │   [cooler          ┌──┘
     │    balance]        │
     │       │            │
     ├───────┼────────────┘
     │       │
     │   ┌───┴────────────┐
     │   │                │
     │ [cooltools       [cooltools               Level 7
     │  eigs]            insulation]
     │   │                │
     │   │(A/B compart)   │(TAD boundaries)
     │   │                │
     └───┴───────┬────────┘
                 │
         [CONVERGENCE 2] ◄──────────────────     Level 8
         (compartments + TADs + pairs)
                 │
     ┌───────────┼───────────┐
     │           │           │
 [cooltools   [chromosight [fithic               Level 9
  saddle]      detect]      (sig loops)]
     │           │           │
     └───────────┴─────┬─────┘
                       │
               [CONVERGENCE 3] ◄─────────────   Level 10
               (saddle + loops + dots)
                       │
               [python report]
               [CONVERGENCE 4] ◄── stats
```

**Longest path:** fastp -> bwa R1 -> filter -> pairtools parse -> sort -> cooler cload -> balance -> cooltools insulation -> fithic -> report (depth 10)

**Tool List (10 tools):**
| Tool | Conda package | Purpose |
|------|--------------|---------|
| fastp | `fastp` | Read QC and adapter trimming |
| bwa-mem2 | `bwa-mem2` | Alignment of Hi-C reads |
| samtools | `samtools` | BAM filtering (MAPQ, chimeric) |
| pairtools | `pairtools` | Parse, sort, dedup contact pairs |
| cooler | `cooler` | Create and balance contact matrices |
| cooltools | `cooltools` | Eigenvectors (A/B), insulation (TADs) |
| chromosight | `chromosight` | Pattern detection (loops, borders) |
| FitHiC | `fithic` | Statistical significance of interactions |
| bedtools | `bedtools` | Interval operations |
| multiqc | `multiqc` | Aggregate QC |

**Data Source:** 4DN Data Portal - GM12878 Hi-C (4DNESUL2GOMU), subset to chr22
- Alternative: Rao et al. 2014 (SRR1658570), downsample to 5M read pairs
- nf-core/hic test data: https://github.com/nf-core/test-datasets/tree/hic

**Data Size:** ~400 MB (5M PE reads + hg38 chr22 reference)
**Estimated Runtime:** ~2.5h on 8 CPUs

**Key Domain-Specific Traps:**
1. **Ligation junction handling**: Hi-C reads contain chimeric sequences at ligation junctions. Must parse with `--walks-policy mask` or similar; naive alignment loses 40-60% of informative contacts
2. **MAPQ filtering**: Must filter MAPQ >= 30 for EACH mate independently BEFORE pairing, not after
3. **Matrix balancing**: ICE/KR balancing is mandatory before compartment/TAD calling; unbalanced matrices give garbage compartments
4. **Resolution selection**: TADs need 10-40 kb resolution, loops need 5-10 kb; wrong resolution = no signal
5. **Dangling ends vs self-circles**: Must identify and remove dangling-end and self-circle artifacts using restriction enzyme cut site info

**Why the DAG is genuinely complex:**
- Nested diamond 1: R1/R2 parallel alignment converges at pairtools parse
- Nested diamond 2: dedup/matrix/stats branches converge for feature calling
- Nested diamond 3: compartments/TADs/loops three-way convergence for integrated report
- Cross-branch dependency: dedup stats feed into loop significance thresholds
- The cooler balance step depends on cload output AND the dedup filtering, creating a non-trivial dependency chain within the second diamond

---

### TASK 2: Circular RNA Detection and Quantification

**Name:** `circrna-detection`
**Description:** Detect and quantify circular RNAs from RNA-seq data using multiple back-splice junction callers, filter by consensus, and predict miRNA binding sites.

**DAG Diagram (depth=9, convergence=4, tools=11):**
```
   R1.fastq.gz     R2.fastq.gz
       │                │
       └───────┬────────┘
               │
         [trim-galore] ──────────────────────── Level 1
               │
         [STAR align] ────────────────────────── Level 2
         (--chimSegmentMin 10
          --chimOutType Junctions)
               │
       ┌───────┼────────────┬─────────────┐
       │       │            │             │
  [CIRCexplorer2] [CIRIquant]  [find_circ]  [DCC]    Level 3
  (BSJ detect)     (BSJ)       (BSJ)       (BSJ)
       │       │            │             │
       └───────┴──────┬─────┴─────────────┘
                      │
              [CONVERGENCE 1] ◄──────────────── Level 4
              (consensus: require >= 2 tools)
              [python merge_bsj.py]
                      │
           ┌──────────┼──────────────┐
           │          │              │
     [bedtools     [CIRIquant     [samtools        Level 5
      getfasta]     quant]         flagstat]
      (circRNA      (BSJ reads)    (linear reads)
       sequences)       │              │
           │          └──────┬─────────┘
           │                 │
           │          [CONVERGENCE 2] ◄──────── Level 6
           │          (circular/linear ratio)
           │                 │
     ┌─────┴────┐     [CircTest DE] ◄── design  Level 7
     │          │            │
  [miranda   [RNAhybrid]    │
   targets]   targets]      │
     │          │            │
     └────┬─────┘            │
          │                  │
   [CONVERGENCE 3] ◄────────┘                   Level 8
   (targets + DE results)
          │
   [python report]
   [CONVERGENCE 4] ◄── QC stats                 Level 9
```

**Longest path:** trim-galore -> STAR -> CIRIquant -> consensus -> CIRIquant quant -> ratio -> CircTest -> report (depth 9)

**Tool List (11 tools):**
| Tool | Conda package | Purpose |
|------|--------------|---------|
| Trim Galore | `trim-galore` | Adapter + quality trimming |
| STAR | `star` | Chimeric alignment |
| CIRCexplorer2 | `circexplorer2` | BSJ detection (annotation-based) |
| CIRIquant | `ciriquant` | BSJ detection + quantification |
| find_circ | `find_circ` | BSJ detection (unmapped reads) |
| DCC | `dcc` | BSJ detection (mate pair) |
| samtools | `samtools` | BAM operations |
| bedtools | `bedtools` | Sequence extraction |
| miRanda | `miranda` | miRNA target prediction |
| DESeq2 | `bioconductor-deseq2` | Differential circRNA expression |
| Python/pandas | `pandas` | Merge, consensus, reporting |

**Data Source:** nf-core/circrna test dataset
- https://github.com/nf-core/test-datasets/tree/circrna
- SRA alternative: SRP098984 (HeLa circRNA profiling, ~300 MB subset)

**Data Size:** ~400 MB (PE RNA-seq reads + genome + annotation)
**Estimated Runtime:** ~2h on 8 CPUs

**Key Domain-Specific Traps:**
1. **STAR chimeric parameters**: Must set `--chimSegmentMin 10 --chimJunctionOverhangMin 10 --chimOutType Junctions SeparateSAMold`; default STAR produces no chimeric reads
2. **Width-5 fan-out consensus**: Each tool produces different BSJ coordinate formats (0-based vs 1-based, strand conventions); must normalize before intersection
3. **Circular-to-linear ratio (CLR)**: Critical metric requires BOTH BSJ counts AND linear junction counts from the same locus; agents often forget the linear read quantification branch
4. **Gene annotation dependency**: CIRCexplorer2 and DCC require specific annotation formats (refFlat vs GTF); wrong format = silent 0 circRNAs

**Why the DAG is genuinely complex:**
- Width-4 fan-out (4 BSJ callers) into consensus convergence
- Second diamond: sequence extraction vs quantification, both needed for downstream
- Third convergence: miRNA targets require circRNA sequences, DE requires counts -- independent branches merging
- Cross-branch: linear read counts from flagstat feed into CLR calculation alongside BSJ counts from a different branch

---

### TASK 3: CNV Detection from WGS/WES

**Name:** `cnv-detection-wes`
**Description:** Detect copy number variants from whole-exome sequencing data using multiple CNV callers, merge calls, and annotate with gene impact.

**DAG Diagram (depth=10, convergence=4, tools=10):**
```
 tumor.fastq.gz         normal.fastq.gz
     │                       │
 [fastp QC] ─────────── [fastp QC]              Level 1
     │                       │
 [bwa-mem2 align] ────── [bwa-mem2 align]        Level 2
     │                       │
 [picard MarkDup] ────── [picard MarkDup]        Level 3
     │                       │
 [mosdepth coverage] ─── [mosdepth coverage]     Level 4
     │                       │
     └───────────┬───────────┘
                 │
         [CONVERGENCE 1] ◄──────────────────     Level 5
         (T+N BAMs + coverage)
                 │
     ┌───────────┼───────────┬──────────┐
     │           │           │          │
 [CNVkit       [Control-  [gatk4      [cnvnator  Level 6
  batch]        FREEC]     CollectRC→   (bin+
                           DenoiseRC→   partition+
                           ModelSeg→    call)]
                           CallSeg]
     │           │           │          │
     └───────────┴─────┬─────┴──────────┘
                       │
               [CONVERGENCE 2] ◄─────────────   Level 7
               (4-caller merge)
               [SURVIVOR merge OR
                bedtools intersect]
                       │
               ┌───────┼───────┐
               │       │       │
          [bcftools  [bedtools [python          Level 8
           annotate]  intersect gene-level
           (dbVar)]   w/ genes] stats]
               │       │       │
               └───────┼───────┘
                       │
               [CONVERGENCE 3] ◄─────────────   Level 9
               (annotated + gene-impact + stats)
                       │
               [python report]
               [CONVERGENCE 4] ◄── QC + cov     Level 10
```

**Longest path:** fastp -> bwa -> markdup -> mosdepth -> convergence1 -> GATK CNV (4 sub-steps) -> merge -> annotate -> report (depth 10)

**Tool List (10 tools):**
| Tool | Conda package | Purpose |
|------|--------------|---------|
| fastp | `fastp` | Read QC |
| bwa-mem2 | `bwa-mem2` | Alignment |
| picard | `picard` | Duplicate marking |
| mosdepth | `mosdepth` | Coverage calculation |
| CNVkit | `cnvkit` | Read-depth CNV caller |
| Control-FREEC | `control-freec` | BAF+depth CNV caller |
| GATK4 | `gatk4` | Model-based CNV caller |
| CNVnator | `cnvnator` | Histogram-based CNV caller |
| bcftools | `bcftools` | VCF operations |
| bedtools | `bedtools` | Interval operations and gene overlap |

**Data Source:** GIAB HG002 (NA24385) WES data, chr22 subset
- Tumor-like: HG002 with in silico CNV spikes (or use TCGA test data)
- nf-core/sarek test data includes T/N WES pairs
- Alternative: SEQC2 somatic reference samples (SRR7890824/SRR7890825)

**Data Size:** ~500 MB (T+N WES BAMs for chr22 + reference + targets BED)
**Estimated Runtime:** ~3h on 8 CPUs

**Key Domain-Specific Traps:**
1. **Target/bait BED file**: CNVkit and GATK4 CNV require exome capture intervals; Control-FREEC needs different format (chr\tstart\tend\tgc); forgetting this = whole-genome mode on WES data = garbage
2. **GC correction**: All callers need GC content correction; Control-FREEC computes internally but CNVkit needs explicit `--access` BED
3. **Normal panel**: GATK4 CNV best practice uses a panel of normals (PoN); with single normal, must use matched mode which has different command syntax
4. **Merge strategy**: SURVIVOR merge requires VCF format, but CNVkit outputs `.cns` segments -- format conversion is a hidden step
5. **Ploidy assumption**: Control-FREEC needs ploidy parameter (default=2); tumor samples with genome doubling will fail silently

**Why the DAG is genuinely complex:**
- Parallel T/N preprocessing (2 branches, depth 4 each) converge at caller stage
- Width-4 caller fan-out (CNVkit, Control-FREEC, GATK4 CNV, CNVnator) with different input formats
- GATK4 CNV is internally a 4-step chain (CollectReadCounts -> DenoiseReadCounts -> ModelSegments -> CallCopyRatioSegments)
- Annotation diamond: dbVar annotation vs gene overlap vs summary stats
- Cross-dependency: mosdepth coverage from preprocessing feeds both callers AND final QC

---

### TASK 4: Immune Repertoire Analysis (TCR/BCR-seq)

**Name:** `immune-repertoire`
**Description:** Analyze B cell receptor sequencing data through UMI-based consensus assembly, V(D)J annotation, clonal grouping, and repertoire diversity analysis.

**DAG Diagram (depth=10, convergence=4, tools=10):**
```
   R1.fastq.gz (V-region)    R2.fastq.gz (C-region + UMI)
       │                          │
   [fastp QC] ──────────── [fastp QC]              Level 1
       │                          │
       └──────────┬───────────────┘
                  │
   [pRESTO FilterSeq] ◄── quality filter           Level 2
                  │
   [pRESTO MaskPrimers] ◄── primer masking          Level 3
                  │
   [pRESTO PairSeq] ◄── coordinate pairing          Level 4
                  │
   [pRESTO ClusterSets] ◄── UMI clustering           Level 5
                  │
   [pRESTO BuildConsensus] ◄── consensus per UMI     Level 6
                  │
   [pRESTO AssemblePairs]                            Level 7
                  │
         [IgBLAST] ◄── V(D)J germline DB            Level 8
                  │
         [Change-O MakeDB]                           Level 8 (cont.)
                  │
         ┌────────┼────────┐
         │        │        │
   [SHazaM     [TIgGER   [Change-O                  Level 9
    threshold]  findNovel] ParseDB]
    (clonal      (novel     (filtering)
     distance)   alleles)
         │        │        │
         └────┬───┴────────┘
              │
      [CONVERGENCE 1] ◄─────────────────
      (threshold + alleles + filtered DB)
              │
      [SCOPer defineClonesScoper]                    Level 10
      (clonal assignment)
              │
      ┌───────┼───────────┐
      │       │           │
 [Dowser    [Alakazam    [Alakazam                  Level 11
  trees]     diversity]   geneUsage]
  (lineage   (Hill/       (V/J usage)
   trees)    Chao/Shan)
      │       │           │
      └───────┴─────┬─────┘
                    │
            [CONVERGENCE 2] ◄────────────
            (trees + diversity + usage)
                    │
       ┌────────────┼────────────┐
       │            │            │
  [Alakazam      [python       [python              Level 12
   mutationAnal]  CDR3-stats]   SHM-stats]
       │            │            │
       └────────────┼────────────┘
                    │
            [CONVERGENCE 3] ◄────────────
            (mutation + CDR3 + SHM)
                    │
            [python report]
            [CONVERGENCE 4] ◄── QC stats
```

**Longest path:** fastp -> FilterSeq -> MaskPrimers -> PairSeq -> ClusterSets -> BuildConsensus -> AssemblePairs -> IgBLAST -> SHazaM -> SCOPer -> Dowser -> report (depth 12)

**Tool List (10 tools):**
| Tool | Conda package | Purpose |
|------|--------------|---------|
| fastp | `fastp` | Read QC |
| pRESTO | `presto` | UMI processing pipeline (7 sub-steps) |
| IgBLAST | `igblast` | V(D)J gene assignment |
| Change-O | `changeo` | Database creation and parsing |
| SHazaM | `r-shazam` | Clonal distance threshold |
| TIgGER | `r-tigger` | Novel allele detection |
| SCOPer | `r-scoper` | Clonal grouping |
| Dowser | `r-dowser` | Lineage tree construction |
| Alakazam | `r-alakazam` | Diversity and gene usage |
| samtools | `samtools` | Sequence manipulation |

**Data Source:** Immcantation tutorial dataset (10x Genomics PBMC BCR)
- nf-core/airrflow test data: https://github.com/nf-core/test-datasets/tree/airrflow
- Alternative: Gupta et al. 2017 (SRP109035), ~200 MB subset

**Data Size:** ~300 MB (BCR-seq FASTQ + IMGT germline database)
**Estimated Runtime:** ~1.5h on 8 CPUs

**Key Domain-Specific Traps:**
1. **8-step sequential pRESTO chain**: Each step has specific parameter requirements and produces a specific output format; skipping MaskPrimers or ClusterSets corrupts downstream
2. **IMGT germline database format**: IgBLAST requires IMGT-gapped germline sequences in specific directory structure; wrong format = 0 annotations
3. **Clonal threshold**: SHazaM distance-to-nearest must be computed BEFORE clonal assignment; using arbitrary threshold gives wrong clonotypes
4. **Novel allele detection**: TIgGER must run BEFORE clonal grouping because novel alleles change V gene assignments
5. **R package interop**: SHazaM, TIgGER, SCOPer, Dowser, Alakazam share an `airr` data format but each expects specific columns; column name mismatches cause silent failures

**Why the DAG is genuinely complex:**
- Deep sequential chain (pRESTO 7-step) before any branching
- Three-way branch after Change-O (threshold/alleles/filtering) must ALL converge before clonal assignment
- Second three-way branch (trees/diversity/usage) after clonal assignment
- Third branch (mutation/CDR3/SHM analysis) creates triple-nested diamond
- N-to-1-to-N pattern: many sequences -> per-clone grouping -> per-clone analysis branches

---

### TASK 5: GWAS / Population Association Testing

**Name:** `gwas-association`
**Description:** Perform genome-wide association testing on genotype array data, including QC, population stratification correction, association testing with multiple methods, and annotation.

**DAG Diagram (depth=10, convergence=4, tools=10):**
```
   genotypes.bed/bim/fam    phenotype.txt    covariates.txt
           │                     │                │
   [plink2 --geno --mind         │                │
    --maf --hwe QC] ────────────────────────────── Level 1
           │
   ┌───────┼───────────┐
   │       │           │
 [plink2 [plink2     [plink2                       Level 2
  --ld     --pca]      --het
  prune]               (inbreeding)]
   │       │           │
   │   [CONVERGENCE 1] │
   │   (pruned + PCA   │
   │    eigenvecs)     │
   │       │           │
   │       │    ┌──────┘
   │       │    │
   │  [plink2 outlier removal]                     Level 3
   │       │
   └───────┴───────────┐
                       │
           [CONVERGENCE 2] ◄─────────────────      Level 4
           (QC'd genotypes + PCs + cleaned samples)
                       │
           ┌───────────┼───────────┐
           │           │           │
     [plink2        [regenie      [plink2           Level 5
      --glm          step1→step2]  --assoc
      (linear/                     (Fisher)]
      logistic)]
           │           │           │
           └───────────┼───────────┘
                       │
               [CONVERGENCE 3] ◄─────────────      Level 6
               (3-method results merge)
                       │
           ┌───────────┼───────────┐
           │           │           │
     [python        [bcftools    [plink2             Level 7
      manhattan+     annotate     --clump
      QQ plot]       (rsIDs)]     (LD clumping)]
           │           │           │
           └───────────┼───────────┘
                       │
               [CONVERGENCE 4] ◄─────────────      Level 8
               (plots + annotations + clumps)
                       │
               ┌───────┼───────┐
               │       │       │
          [SnpSift  [bedtools [python                Level 9
           annotate  closest   gene
           ClinVar]  (genes)]  enrichment]
               │       │       │
               └───────┼───────┘
                       │
               [python report]                       Level 10
```

**Longest path:** QC -> LD prune -> PCA -> outlier removal -> convergence2 -> regenie step1+step2 -> convergence3 -> clump -> convergence4 -> annotate -> report (depth 10)

**Tool List (10 tools):**
| Tool | Conda package | Purpose |
|------|--------------|---------|
| PLINK2 | `plink2` | Genotype QC, PCA, association, clumping |
| PLINK | `plink` | LD pruning, format conversion |
| REGENIE | `regenie` | Whole-genome regression association |
| bcftools | `bcftools` | VCF annotation |
| SnpSift | `snpsift` | ClinVar annotation |
| bedtools | `bedtools` | Gene proximity |
| R/qqman | `r-qqman` | Manhattan and QQ plots |
| samtools | `samtools` | Reference indexing |
| Python/scipy | `scipy` | Statistical tests, enrichment |
| multiqc | `multiqc` | QC aggregation |

**Data Source:** 1000 Genomes Project Phase 3, chr22 subset
- PLINK tutorial data: https://zzz.bwh.harvard.edu/plink/tutorial.shtml
- OpenSNP subset: https://opensnp.org/
- Alternative: GCTA tutorial data (simulated, ~100 MB)
- nf-core/gwas does not exist, but regenie tutorial data is available

**Data Size:** ~200 MB (genotype BED/BIM/FAM for chr22 + phenotype + ClinVar VCF)
**Estimated Runtime:** ~1.5h on 8 CPUs

**Key Domain-Specific Traps:**
1. **Population stratification**: Must compute PCA on LD-pruned SNPs and include top PCs as covariates; omitting this inflates false positives catastrophically (genomic inflation lambda >> 1)
2. **LD pruning parameters**: `--indep-pairwise 50 5 0.2` is standard; using wrong window or r2 threshold leads to too many or too few tag SNPs
3. **HWE filtering**: Must apply `--hwe 1e-6` for controls only in case-control studies; applying to cases removes true disease-associated variants
4. **REGENIE two-step**: Step 1 produces ridge regression predictions, step 2 uses those as offset. Agents often run step 2 without step 1 output = crash
5. **Multiple testing**: Must apply Bonferroni (5e-8) or FDR correction; reporting raw p-values without correction is a fundamental error

**Why the DAG is genuinely complex:**
- QC phase has 3-way branch (LD prune, PCA, heterozygosity) with cross-dependencies
- PCA requires LD-pruned set; outlier removal requires both PCA and heterozygosity
- Association phase: 3 parallel methods with different input requirements
- Post-association: clumping, annotation, and visualization are interdependent (clumping needs LD, annotation needs positions, visualization needs both)
- REGENIE has internal 2-step dependency that adds hidden depth

---

### TASK 6: Germline WES Variant Calling (GATK Best Practices)

**Name:** `germline-wes-gatk`
**Description:** Call germline variants from whole-exome sequencing following GATK best practices, including BQSR, HaplotypeCaller, variant filtration, and clinical annotation.

**DAG Diagram (depth=12, convergence=4, tools=12):**
```
 sample.R1.fq.gz    sample.R2.fq.gz
       │                  │
   [fastp QC] ──────── [fastp QC]                  Level 1
       │                  │
       └────────┬─────────┘
                │
        [bwa-mem2 align]                            Level 2
                │
        [samtools sort]                             Level 3
                │
        [picard MarkDuplicates]                     Level 4
                │
        ┌───────┼───────┐
        │       │       │
   [gatk BQSR [mosdepth [picard                     Level 5
    BaseRecal]  coverage] CollectHsMetrics]
        │       │       │
   [gatk BQSR  │       │
    ApplyBQSR] │       │
        │       │       │
        └───────┼───────┘
                │
        [CONVERGENCE 1] ◄── (recal BAM + QC)       Level 6
                │
        [gatk HaplotypeCaller]                       Level 7
                │
        [gatk GenotypeGVCFs]                         Level 8
                │
        ┌───────┼───────────┐
        │       │           │
   [gatk Select [gatk Select                        Level 9
    SNPs]        Indels]
        │       │
   [gatk Variant [gatk Variant
    Filtration    Filtration
    (SNP filters)] (Indel filters)]
        │       │
        └───────┼───────────┘
                │
        [CONVERGENCE 2] ◄── (filtered SNPs+Indels)  Level 10
        [gatk MergeVcfs]
                │
        ┌───────┼───────────┐
        │       │           │
   [SnpSift   [bcftools   [python                   Level 11
    annotate   stats]       Ti/Tv
    ClinVar]               het/hom]
        │       │           │
        └───────┼───────────┘
                │
        [CONVERGENCE 3] ◄── (annotations + stats)
                │
        [bcftools filter PASS]
                │
        [python clinical report]
        [CONVERGENCE 4] ◄── (coverage + Hs metrics) Level 12
```

**Longest path:** fastp -> bwa -> sort -> markdup -> BQSR(BaseRecal+Apply) -> HaplotypeCaller -> GenotypeGVCFs -> SelectSNPs -> FilterSNPs -> MergeVcfs -> SnpSift -> report (depth 12)

**Tool List (12 tools):**
| Tool | Conda package | Purpose |
|------|--------------|---------|
| fastp | `fastp` | Read QC |
| bwa-mem2 | `bwa-mem2` | Alignment |
| samtools | `samtools` | Sort, index |
| picard | `picard` | MarkDuplicates, CollectHsMetrics |
| GATK4 | `gatk4` | BQSR, HaplotypeCaller, genotyping, filtration |
| mosdepth | `mosdepth` | Coverage depth |
| bcftools | `bcftools` | VCF stats, filtering |
| SnpSift | `snpsift` | ClinVar annotation |
| bedtools | `bedtools` | Target interval operations |
| multiqc | `multiqc` | QC aggregation |
| Python/pandas | `pandas` | Report generation |
| R/ggplot2 | `r-ggplot2` | Coverage plots |

**Data Source:** GIAB HG001 (NA12878) WES data, chr22 subset
- Garvan Institute WES: ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/
- GATK resource bundle: gs://gatk-test-data/
- Alternative: SRA SRR098401 (NA12878 WES, subset to chr22)

**Data Size:** ~600 MB (WES FASTQ chr22 + hg38 chr22 + dbSNP + ClinVar + known indels + targets BED)
**Estimated Runtime:** ~3h on 8 CPUs

**Key Domain-Specific Traps:**
1. **Read group (@RG)**: bwa-mem2 MUST include `-R @RG\tID:...\tSM:...\tPL:ILLUMINA`; missing RG causes GATK to crash at MarkDuplicates
2. **BQSR two-step**: BaseRecalibrator creates table, ApplyBQSR applies it. Many agents try to run ApplyBQSR without the recal table
3. **Known sites for BQSR**: Must provide BOTH dbSNP AND known indels; using only dbSNP under-corrects indel qualities
4. **SNP vs Indel filtration**: Different hard filter thresholds (SNP: QD<2, FS>60, MQ<40; Indel: QD<2, FS>200); applying SNP filters to indels removes good calls
5. **CollectHsMetrics**: Requires BOTH bait and target interval files; WES without this step gives no capture efficiency metrics
6. **Interval padding**: HaplotypeCaller should pad target intervals by 100-150bp to capture variants near exon boundaries

**Why the DAG is genuinely complex:**
- Long preprocessing chain (7 steps) before any branching
- BQSR is internally a 2-step dependency (BaseRecal -> Apply)
- SNP/Indel split-filter-merge creates a classic nested diamond
- Three parallel annotation/stats branches converge
- Coverage metrics from early preprocessing feed into final report (cross-branch dependency spanning 8 levels)
- GATK4 internal state management (`.dict` files, `.fai` indexes) creates implicit dependencies

---

### TASK 7: Human Structural Variant Detection (Multi-caller)

**NOTE:** `sv-detection` exists but covers BACTERIAL SVs (MRSA, depth 6, Delly only). This task covers HUMAN WGS SVs with 4 callers (Manta+Delly+smoove+GRIDSS), SURVIVOR merge, gnomAD population frequency annotation, and ClinVar clinical interpretation. Completely different scale, tools, and complexity. No overlap.

**Name:** `structural-variant-multi`
**Description:** Detect structural variants from human WGS data using four orthogonal SV callers, merge with SURVIVOR, genotype, and annotate.

**DAG Diagram (depth=10, convergence=4, tools=12):**
```
 sample.R1.fq.gz    sample.R2.fq.gz
       │                  │
   [fastp QC] ──────── [fastp QC]                  Level 1
       │                  │
       └────────┬─────────┘
                │
        [bwa-mem2 align]                            Level 2
                │
        [samtools sort + index]                     Level 3
                │
        [sambamba markdup]                          Level 4
                │
     ┌──────────┼──────────────┬──────────┐
     │          │              │          │
 [Manta      [Delly        [smoove     [GRIDSS     Level 5
  SV call]    SV call]      SV call]    SV call]
  (PE+SR)     (PE signal)   (lumpy      (assembly
                             wrapper)    -based)
     │          │              │          │
     └──────────┴──────┬───────┴──────────┘
                       │
               [CONVERGENCE 1] ◄─────────────      Level 6
               [SURVIVOR merge]
               (require >= 2 callers)
                       │
               ┌───────┼───────────┐
               │       │           │
         [SURVIVOR  [bcftools    [python             Level 7
          stats]     annotate     SV size
                     (gnomAD-SV)] distribution]
               │       │           │
               └───────┼───────────┘
                       │
               [CONVERGENCE 2] ◄─────────────      Level 8
               (stats + freq + size dist)
                       │
               ┌───────┼───────┐
               │       │       │
          [bedtools [SnpSift  [bcftools              Level 9
           intersect annotate  query
           (genes)]  (ClinVar)] (PASS filter)]
               │       │       │
               └───────┼───────┘
                       │
               [CONVERGENCE 3] ◄─────────────
               (gene overlap + clinical + filtered)
                       │
               [python clinical report]
               [CONVERGENCE 4] ◄── QC + caller stats Level 10
```

**Longest path:** fastp -> bwa -> sort -> markdup -> GRIDSS (assembly, slowest) -> SURVIVOR -> annotate gnomAD -> convergence2 -> ClinVar -> report (depth 10)

**Tool List (12 tools):**
| Tool | Conda package | Purpose |
|------|--------------|---------|
| fastp | `fastp` | Read QC |
| bwa-mem2 | `bwa-mem2` | Alignment |
| samtools | `samtools` | BAM operations |
| sambamba | `sambamba` | Fast duplicate marking |
| Manta | `manta` | PE+SR SV caller |
| Delly | `delly` | PE-based SV caller |
| smoove | `smoove` | Lumpy wrapper (SR+PE) |
| GRIDSS | `gridss` | Assembly-based SV caller |
| SURVIVOR | `survivor` | SV merging and stats |
| bcftools | `bcftools` | VCF operations |
| SnpSift | `snpsift` | ClinVar annotation |
| bedtools | `bedtools` | Gene overlap |

**Data Source:** GIAB HG002 (NA24385) WGS, chr22 subset (or chr20)
- Tier 1 SV benchmark: ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002/
- nf-core/sarek test data (includes SV calling)
- Alternative: 1000G high-coverage WGS subset

**Data Size:** ~700 MB (WGS FASTQ chr22 + hg38 chr22 + SV truth set)
**Estimated Runtime:** ~3.5h on 8 CPUs (GRIDSS is compute-intensive)

**Key Domain-Specific Traps:**
1. **Caller-specific VCF formats**: Manta uses `<INV:3>/<INV:5>`, Delly uses `<INV>`, smoove uses different INFO fields. SURVIVOR merge requires normalized VCF headers
2. **GRIDSS setup**: Requires BWA index AND a blacklist BED AND exact JVM memory settings (`-Xmx4g`); default settings OOM on human data
3. **SURVIVOR merge parameters**: `SURVIVOR merge sample_files 1000 2 1 1 0 50` -- the `1000` (max distance) and `2` (min callers) are critical; wrong values = too many or too few merged calls
4. **SV type concordance**: Must merge by SV type (DEL with DEL, not DEL with INV). SURVIVOR handles this but only if VCF SVTYPE fields are consistent
5. **Breakpoint precision**: Different callers report breakpoints at different resolutions (Manta: bp-level, Delly: ~100bp). Merge distance must accommodate this

**Why the DAG is genuinely complex:**
- Width-4 SV caller fan-out, each with different algorithmic approaches (PE, SR, assembly)
- SURVIVOR merge is a genuine convergence requiring ALL 4 caller outputs
- Post-merge splits into 3 parallel annotation branches
- Gene overlap + ClinVar + population frequency all converge for clinical interpretation
- GRIDSS has internal assembly step that makes it 3x slower than others, creating an asymmetric diamond where one branch dominates runtime
- Cross-branch: per-caller stats feed into final report alongside merged results

---

### TASK 8: Clinical WGS Interpretation Pipeline

**Name:** `clinical-wgs-interpretation`
**Description:** Full clinical WGS pipeline from FASTQ to clinically actionable variant report, including SNV/indel calling, structural variant detection, coverage analysis, and multi-source annotation.

**DAG Diagram (depth=12, convergence=5, tools=14):**
```
 sample.R1.fq.gz    sample.R2.fq.gz
       │                  │
   [fastp QC] ──────── [fastp QC]                  Level 1
       │                  │
       └────────┬─────────┘
                │
        [bwa-mem2 align]                            Level 2
                │
        [samtools sort]                             Level 3
                │
        [picard MarkDuplicates]                     Level 4
                │
        [gatk BQSR (BaseRecal + Apply)]             Level 5
                │
        ┌───────┼───────────────────┐
        │       │                   │
   [gatk HC]  [Manta            [mosdepth            Level 6
    (SNV/      (SV call)]        (coverage +
     Indel)]                      per-gene)]
        │       │                   │
   [gatk       [Delly              │
    GenotypeGVCFs] (SV call)]      │
        │       │                   │
   ┌────┴────┐  │                   │
   │         │  │                   │
 [Select  [Select                   │               Level 7
  SNP]     Indel]                   │
   │         │                      │
 [Filter  [Filter                   │               Level 8
  SNP]     Indel]                   │
   │         │                      │
   └────┬────┘                      │
        │                           │
   [MergeVcfs]                      │               Level 9
        │       ┌───────────────────┘
        │       │
   [CONVERGENCE 1] ◄── (SNV/Indel + coverage)
        │
        │       [SURVIVOR merge SVs]
        │               │
        │       [CONVERGENCE 2] ◄── (merged SVs)    Level 10
        │               │
        └───────┬───────┘
                │
        [CONVERGENCE 3] ◄── (all variants merged)
                │
        ┌───────┼───────────────┐
        │       │               │
   [SnpSift   [bcftools       [genmod                Level 11
    annotate   annotate        score]
    ClinVar]   (gnomAD freq)]  (inheritance)
        │       │               │
        └───────┴───────┬───────┘
                        │
                [CONVERGENCE 4] ◄── (all annotations)
                        │
                ┌───────┼───────┐
                │       │       │
           [bcftools  [python  [python               Level 12
            filter     gene    clinical
            PASS]      panel]  report]
                │       │       │
                └───────┼───────┘
                        │
                [CONVERGENCE 5] ◄── coverage + QC
                [final clinical report]
```

**Longest path:** fastp -> bwa -> sort -> markdup -> BQSR -> HC -> GenotypeGVCFs -> SelectSNP -> FilterSNP -> MergeVcfs -> conv1 -> conv3 -> SnpSift -> conv4 -> report (depth 12)

**Tool List (14 tools):**
| Tool | Conda package | Purpose |
|------|--------------|---------|
| fastp | `fastp` | Read QC |
| bwa-mem2 | `bwa-mem2` | Alignment |
| samtools | `samtools` | BAM operations |
| picard | `picard` | Duplicate marking |
| GATK4 | `gatk4` | BQSR, HaplotypeCaller, filtration |
| Manta | `manta` | SV detection (PE+SR) |
| Delly | `delly` | SV detection (PE) |
| SURVIVOR | `survivor` | SV merging |
| mosdepth | `mosdepth` | Coverage analysis |
| bcftools | `bcftools` | VCF operations, frequency annotation |
| SnpSift | `snpsift` | ClinVar annotation |
| genmod | `genmod` | Inheritance scoring |
| bedtools | `bedtools` | Gene panel filtering |
| multiqc | `multiqc` | QC aggregation |

**Data Source:** GIAB HG001 (NA12878) WGS, chr22 subset
- GIAB truth VCF: ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/
- ClinVar VCF: ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/
- gnomAD sites VCF (chr22): https://gnomad.broadinstitute.org/

**Data Size:** ~800 MB (WGS FASTQ chr22 + ref + dbSNP + ClinVar + gnomAD + gene panel BED)
**Estimated Runtime:** ~3.5h on 8 CPUs

**Key Domain-Specific Traps:**
1. **Three-track variant calling**: SNV/Indel (GATK), SV (Manta+Delly), and coverage must all be produced from same BAM but processed independently
2. **genmod requires PED file**: Even for singleton analysis, must provide a minimal PED file with `#family_id individual_id paternal_id maternal_id sex phenotype`
3. **Gene panel filtering**: Clinical WGS reports filter to clinically relevant genes; must intersect variants with gene panel BED, not just annotate
4. **ClinVar version matters**: ClinVar annotations must match genome build (GRCh37 vs GRCh38); wrong build = 0 annotations but no error
5. **Coverage thresholds**: Clinical WGS requires >=30x median coverage; must flag regions below 20x as "uncallable" in report
6. **Variant normalization**: Must left-align and normalize indels (bcftools norm) before annotation; unnormalized indels miss ClinVar matches

**Why the DAG is genuinely complex:**
- 5 convergence points (the maximum in any task)
- Three independent tracks (SNV, SV, coverage) from single BAM
- SNV track contains nested SNP/Indel diamond
- SV track contains dual-caller merge diamond
- All three tracks converge at annotation stage
- Triple annotation branch (ClinVar + gnomAD + inheritance)
- Final report requires ALL of: filtered variants, coverage stats, QC metrics, gene panel intersection

---

## MEDIUM PRIORITY (Research Genomics): Tasks 9-15

---

### TASK 9: RNA Editing Detection (A-to-I)

**Name:** `rna-editing-detection`
**Description:** Detect A-to-I RNA editing events from matched RNA-seq and DNA-seq (or RNA-only with stringent filtering), using multiple detection tools and database validation.

**DAG Diagram (depth=10, convergence=4, tools=10):**
```
 RNA_R1.fq.gz  RNA_R2.fq.gz    DNA_R1.fq.gz  DNA_R2.fq.gz
     │              │               │              │
 [fastp] ────── [fastp]         [fastp] ────── [fastp]        Level 1
     │              │               │              │
     └──────┬───────┘               └──────┬───────┘
            │                              │
     [STAR 2-pass align]            [bwa-mem2 align]           Level 2
            │                              │
     [picard MarkDup]               [picard MarkDup]           Level 3
            │                              │
     [gatk SplitNCigarReads]        [gatk BQSR]                Level 4
            │                              │
     [gatk BQSR]                           │                   Level 5
            │                              │
            └──────────────┬───────────────┘
                           │
                   [CONVERGENCE 1] ◄──────────────             Level 6
                   (RNA BAM + DNA BAM)
                           │
            ┌──────────────┼──────────────┐
            │              │              │
      [JACUSA2        [bcftools      [gatk HC                  Level 7
       call-2          mpileup→call   RNA-mode
       (RNA vs DNA)]   (paired)]      (RNA-only)]
            │              │              │
            └──────────────┼──────────────┘
                           │
                   [CONVERGENCE 2] ◄──────────────             Level 8
                   (3-caller intersection)
                   [bcftools isec: >= 2 agree]
                           │
            ┌──────────────┼──────────────┐
            │              │              │
      [filter A>G/T>C  [bcftools     [bedtools                 Level 9
       on strand]       annotate      intersect
                        (REDIportal    (Alu regions,
                         database)]    repeat masker)]
            │              │              │
            └──────────────┼──────────────┘
                           │
                   [CONVERGENCE 3] ◄──────────────
                   (strand-filtered + DB + repeats)
                           │
                   [python editing report]
                   [CONVERGENCE 4] ◄── RNA/DNA QC              Level 10
```

**Longest path:** RNA fastp -> STAR 2-pass -> MarkDup -> SplitNCigar -> BQSR -> convergence1 -> JACUSA2 -> convergence2 -> REDIportal -> convergence3 -> report (depth 10)

**Tool List (10 tools):**
| Tool | Conda package | Purpose |
|------|--------------|---------|
| fastp | `fastp` | Read QC |
| STAR | `star` | RNA-seq alignment (2-pass) |
| bwa-mem2 | `bwa-mem2` | DNA alignment |
| picard | `picard` | Duplicate marking |
| GATK4 | `gatk4` | SplitNCigar, BQSR, HaplotypeCaller |
| JACUSA2 | `jacusa2` | RNA vs DNA variant calling |
| bcftools | `bcftools` | Variant calling, intersection, annotation |
| samtools | `samtools` | BAM operations |
| bedtools | `bedtools` | Repeat region intersection |
| Python/pandas | `pandas` | Filtering and reporting |

**Data Source:** ENCODE K562 RNA-seq + matched WGS (chr22 subset)
- RNA-seq: ENCSR000AEL (SRR521448)
- WGS: ENCSR765JPC
- REDIportal database: http://srv00.recas.ba.infn.it/atlas/download.html
- Alternative: GTEx sample with matched RNA+WGS

**Data Size:** ~500 MB (RNA FASTQ + DNA FASTQ + reference + REDIportal + RepeatMasker BED)
**Estimated Runtime:** ~2.5h on 8 CPUs

**Key Domain-Specific Traps:**
1. **STAR 2-pass alignment**: RNA editing detection requires 2-pass mode for novel splice junction discovery; 1-pass misaligns reads at editing-enriched junctions
2. **SplitNCigarReads**: Mandatory for RNA BAM before variant calling; without it, GATK HC misidentifies splice junction N-cigar as deletions
3. **Strand-specific filtering**: A-to-I editing appears as A>G on sense strand and T>C on antisense; must filter by strand orientation, not just base change
4. **SNP contamination**: Must remove known SNPs (dbSNP) from editing candidates; DNA-seq matching is the gold standard but requires careful alignment
5. **Alu enrichment**: ~90% of A-to-I editing occurs in Alu repeats; zero overlap with Alu regions signals pipeline failure

**Why the DAG is genuinely complex:**
- Dual-input (RNA + DNA) with independent preprocessing chains of different depth (RNA needs SplitNCigar, DNA does not)
- Three parallel callers (JACUSA2 for paired, bcftools for generic, GATK HC for RNA-only) converge via intersection
- Triple annotation branch (strand filter, database, repeat region) creates third diamond
- Cross-input dependency: DNA BAM feeds into RNA variant filtering
- RNA preprocessing is 5 steps deep before first convergence; DNA is 3 steps -- asymmetric branches

---

### TASK 10: Nascent Transcription (GRO-seq/PRO-seq)

**Name:** `nascent-transcription`
**Description:** Analyze PRO-seq/GRO-seq data to map RNA polymerase positions, identify transcription start sites, active enhancers, and nascent transcript quantification.

**DAG Diagram (depth=9, convergence=4, tools=10):**
```
 R1.fastq.gz    R2.fastq.gz
     │               │
     └───────┬───────┘
             │
     [fastp + UMI extract] ──────────────────── Level 1
             │
     [bowtie2 align to rRNA] ◄── remove rRNA    Level 2
     (keep unmapped)
             │
     [STAR align to genome]                      Level 3
             │
     [samtools dedup + filter]                   Level 4
             │
     ┌───────┼───────────────────┐
     │       │                   │
 [bedtools [bedtools          [samtools            Level 5
  genomecov  genomecov          flagstat]
  (+ strand)] (- strand)]
     │       │                   │
     │  [CONVERGENCE 1]          │
     │  (strand-specific         │
     │   bigWig via              │
     │   bedGraphToBigWig)       │
     │       │                   │
     ├───────┤                   │
     │       │                   │
 ┌───┴───┐   │                   │
 │       │   │                   │
[HOMER  [HOMER  [deeptools                          Level 6
 findPeaks findPeaks  computeMatrix
 -style   -style     (TSS profile)]
 tss]     groseq]
 │       │           │
 │(TSSs) │(TREs)     │(heatmap)
 │       │           │
 └───┬───┘           │
     │               │
 [CONVERGENCE 2]     │                              Level 7
 (TSS + enhancer)    │
     │               │
 ┌───┼───────┐       │
 │   │       │       │
[HOMER [bedtools [deeptools                         Level 8
 annotate intersect plotHeatmap]
 Peaks]  enhancers
         w/ H3K27ac]
 │   │       │       │
 └───┴───────┼───────┘
             │
     [CONVERGENCE 3] ◄──────────────────
     (annotated TSS + enhancers + heatmaps)
             │
     [CONVERGENCE 4] ◄── QC + flagstat              Level 9
     [python report]
```

**Longest path:** fastp -> bowtie2(rRNA) -> STAR -> dedup -> genomecov -> HOMER findPeaks(tss) -> convergence2 -> HOMER annotate -> convergence3 -> report (depth 9)

**Tool List (10 tools):**
| Tool | Conda package | Purpose |
|------|--------------|---------|
| fastp | `fastp` | Read QC + UMI extraction |
| bowtie2 | `bowtie2` | rRNA depletion alignment |
| STAR | `star` | Genome alignment |
| samtools | `samtools` | Dedup, filter, flagstat |
| bedtools | `bedtools` | Strand-specific coverage |
| HOMER | `homer` | TSS and enhancer detection |
| deeptools | `deeptools` | TSS profiles and heatmaps |
| UCSC tools | `ucsc-bedgraphtobigwig` | BigWig conversion |
| Python/pandas | `pandas` | Reporting |
| multiqc | `multiqc` | QC aggregation |

**Data Source:** ENCODE PRO-seq data (K562)
- ENCSR916GXK (PRO-seq K562), subset to chr22
- GRO-seq alternative: Core et al. 2014 (SRR1552484)
- nf-core/nascent test data: https://github.com/nf-core/test-datasets/tree/nascent

**Data Size:** ~300 MB (PE FASTQ + genome + rRNA index + gene annotation)
**Estimated Runtime:** ~1.5h on 8 CPUs

**Key Domain-Specific Traps:**
1. **3' end mapping**: PRO-seq maps the 3' end of nascent RNA (active site of Pol II); must use R1 (or R2 depending on library) for strand assignment -- getting this backwards inverts all results
2. **rRNA depletion**: 30-50% of reads map to rRNA; must remove BEFORE genome alignment or all peak callers are overwhelmed
3. **Strand-specific coverage**: Must generate SEPARATE plus and minus strand BigWig files; combined-strand coverage obscures bidirectional transcription at enhancers
4. **HOMER mode selection**: `-style tss` for promoters vs `-style groseq` for gene body -- using wrong style gives incorrect features
5. **UMI deduplication**: PRO-seq uses UMI-based dedup, not position-based; wrong dedup method removes 80% of real signal at high-expression genes

**Why the DAG is genuinely complex:**
- Sequential rRNA depletion -> genome alignment creates a dependent chain (not parallelizable)
- Strand-specific split creates first diamond (plus/minus strand coverage)
- TSS vs enhancer detection from same coverage creates second diamond
- Annotation + enhancer validation + heatmap creates third diamond
- Cross-branch: deeptools heatmap needs the BigWig files from strand-specific branch AND TSS positions from HOMER branch

---

### TASK 11: Long-read RNA Isoform (Nanopore Direct RNA)

**Name:** `longread-rna-isoform`
**Description:** Analyze Nanopore direct RNA-seq data for full-length transcript isoform discovery, novel isoform detection, differential isoform usage, and poly(A) tail estimation.

**DAG Diagram (depth=10, convergence=4, tools=11):**
```
 reads.fastq.gz (Nanopore direct RNA)
          │
    [NanoPlot QC] ─────────────────────────── Level 1
          │
    [minimap2 -ax splice] ◄── splice-aware    Level 2
          │
    [samtools sort + index]                    Level 3
          │
    ┌─────┼──────────────────┬──────────┐
    │     │                  │          │
 [StringTie [IsoQuant      [FLAIR     [samtools    Level 4
  -L        (known+novel)]  correct]   flagstat]
  long-read]                    │
    │     │               [FLAIR
    │     │                align]
    │     │                    │
    │     │               [FLAIR
    │     │                collapse]
    │     │                    │
    │     │                    │
    └─────┴──────────┬─────────┘
                     │
             [CONVERGENCE 1] ◄───────────────  Level 5
             (3 transcript assemblers merged)
             [gffcompare consensus]
                     │
             ┌───────┼───────────┐
             │       │           │
       [gffcompare [SQANTI3-   [python            Level 6
        vs ref      like         novel isoform
        annotation] classify]    stats]
             │       │           │
             └───────┼───────────┘
                     │
             [CONVERGENCE 2] ◄───────────────  Level 7
             (ref comparison + classification + stats)
                     │
             ┌───────┼───────────┐
             │       │           │
       [salmon      [python    [nanopolish/         Level 8
        quant        DTE/DTU    tailfindr
        --ont]       analysis]  (poly-A tail)]
             │       │           │
             └───────┼───────────┘
                     │
             [CONVERGENCE 3] ◄───────────────  Level 9
             (quant + DTE + poly-A)
                     │
             [python report]
             [CONVERGENCE 4] ◄── NanoPlot QC    Level 10
```

**Longest path:** NanoPlot -> minimap2 -> sort -> FLAIR(correct+align+collapse) -> gffcompare consensus -> classification -> convergence2 -> DTE/DTU -> convergence3 -> report (depth 10)

**Tool List (11 tools):**
| Tool | Conda package | Purpose |
|------|--------------|---------|
| NanoPlot | `nanoplot` | Read QC |
| minimap2 | `minimap2` | Splice-aware long-read alignment |
| samtools | `samtools` | BAM operations |
| StringTie | `stringtie` | Transcript assembly (long-read mode) |
| IsoQuant | `isoquant` | Isoform detection (known+novel) |
| FLAIR | `flair` | Full-length alternative isoform analysis |
| gffcompare | `gffcompare` | Transcript comparison/consensus |
| Salmon | `salmon` | Transcript quantification |
| Python/pandas | `pandas` | DTE/DTU analysis, reporting |
| bedtools | `bedtools` | Genomic interval operations |
| multiqc | `multiqc` | QC aggregation |

**Data Source:** SG-NEx Nanopore direct RNA (SGNex_A549)
- https://github.com/GoekeLab/sg-nex-data
- Alternative: ENCODE direct RNA-seq (ENCSR706NOW)
- Subset to chr22 (~200 MB)

**Data Size:** ~400 MB (Nanopore FASTQ + reference + annotation + barcode files)
**Estimated Runtime:** ~2h on 8 CPUs

**Key Domain-Specific Traps:**
1. **minimap2 splice preset**: Must use `-ax splice` (not `-ax map-ont`); genomic alignment mode does not handle introns and breaks all downstream tools
2. **StringTie long-read mode**: Must use `-L` flag; without it, StringTie uses short-read model and fragments long transcripts
3. **FLAIR three-step**: correct -> align -> collapse is a mandatory sequential chain; skipping correct produces 60% fewer valid isoforms
4. **Nanopore-specific error profile**: Direct RNA has ~8-12% error rate, concentrated in homopolymers; tools need error-tolerant parameters
5. **Novel isoform classification**: Must distinguish FSM (full splice match), ISM (incomplete splice match), NIC (novel in catalog), NNC (novel not in catalog); mixing these categories invalidates isoform discovery claims

**Why the DAG is genuinely complex:**
- Three parallel transcript assemblers (StringTie, IsoQuant, FLAIR) each with different internal depth
- FLAIR has 3 internal sequential steps, making it the longest branch in the first diamond
- gffcompare consensus is a genuine convergence requiring ALL assembler outputs
- Second diamond: reference comparison vs classification vs novel stats
- Third diamond: quantification vs differential vs poly(A) -- all require different inputs from previous convergence
- Cross-branch: NanoPlot QC metrics feed into final report alongside isoform analysis results

---

### TASK 12: MAG Recovery from Metagenomes

**Name:** `mag-recovery`
**Description:** Recover metagenome-assembled genomes from short-read metagenomics data using multiple assemblers, multiple binners, bin refinement, and taxonomic/quality assessment.

**DAG Diagram (depth=12, convergence=5, tools=14):**
```
 R1.fastq.gz     R2.fastq.gz
     │                │
 [fastp QC] ────── [fastp QC]                     Level 1
     │                │
     └───────┬────────┘
             │
 [bowtie2 host removal] ◄── human ref             Level 2
             │
     ┌───────┼───────┐
     │       │       │
 [MEGAHIT] [SPAdes  [MEGAHIT                       Level 3
             meta]   --k-list]
     │       │       │
     └───────┴───┬───┘
                 │
         [CONVERGENCE 1] ◄──────────────           Level 4
         (assembly merge: longest contigs)
         [quast assessment]
                 │
         [bowtie2 map-back reads]                   Level 5
                 │
         [samtools sort + depth]                    Level 6
                 │
     ┌───────────┼───────────────┬──────────┐
     │           │               │          │
 [MetaBAT2]  [MaxBin2]       [SemiBin]  [CONCOCT   Level 7
                                          (coverage
                                           + comp)]
     │           │               │          │
     └───────────┴───────┬───────┴──────────┘
                         │
                 [CONVERGENCE 2] ◄──────────       Level 8
                 [DAS Tool] (bin refinement)
                         │
                 ┌───────┼───────────┐
                 │       │           │
           [CheckM2   [GTDB-Tk   [QUAST              Level 9
            (quality)]  (taxonomy)] (assembly stats)]
                 │       │           │
                 └───────┴─────┬─────┘
                               │
                       [CONVERGENCE 3] ◄─────       Level 10
                       (quality + taxonomy + stats)
                               │
                       ┌───────┼───────┐
                       │       │       │
                  [filter   [Prokka  [python          Level 11
                   HQ bins]  annotate  diversity]
                       │       │       │
                       └───────┼───────┘
                               │
                       [CONVERGENCE 4] ◄──────       Level 11
                       (HQ MAGs + annotation + diversity)
                               │
                       [python report]
                       [CONVERGENCE 5] ◄── QC        Level 12
```

**Longest path:** fastp -> host removal -> MEGAHIT -> convergence1 -> map-back -> sort -> MetaBAT2 -> DAS Tool -> CheckM2 -> convergence3 -> Prokka -> convergence4 -> report (depth 12)

**Tool List (14 tools):**
| Tool | Conda package | Purpose |
|------|--------------|---------|
| fastp | `fastp` | Read QC |
| bowtie2 | `bowtie2` | Host removal, read mapping |
| MEGAHIT | `megahit` | Metagenomic assembly |
| SPAdes | `spades` | Metagenomic assembly (meta mode) |
| samtools | `samtools` | BAM operations |
| MetaBAT2 | `metabat2` | Composition+coverage binning |
| MaxBin2 | `maxbin2` | EM-based binning |
| SemiBin | `semibin` | Deep learning binning |
| CONCOCT | `concoct` | Variational inference binning |
| DAS Tool | `das_tool` | Bin refinement |
| CheckM2 | `checkm2` | Bin quality (completeness/contamination) |
| GTDB-Tk | `gtdbtk` | Taxonomic classification |
| QUAST | `quast` | Assembly quality |
| multiqc | `multiqc` | QC aggregation |

**Data Source:** CAMI2 challenge marine dataset (subset)
- https://data.cami-challenge.org/
- nf-core/mag test data: https://github.com/nf-core/test-datasets/tree/mag
- Alternative: HMP2 metagenome subset

**Data Size:** ~600 MB (PE metagenome reads + human ref for decontamination)
**Estimated Runtime:** ~3h on 8 CPUs

**Key Domain-Specific Traps:**
1. **MEGAHIT output directory**: MEGAHIT refuses to run if output directory exists (even empty); must delete before re-run
2. **Coverage depth file format**: MetaBAT2 needs `jgi_summarize_bam_contig_depths` output; MaxBin2 needs per-contig coverage; CONCOCT needs its own coverage table -- three different formats from same BAM
3. **DAS Tool input format**: Requires `contig_id\tbin_id` TSV per binner; must convert each binner's native format
4. **GTDB-Tk database**: Requires ~70 GB database; for benchmark, pre-download and provide path. Agent must use `--mash_db` for fast ANI
5. **SemiBin requires**: Must specify `--environment` (human_gut, dog_gut, ocean, etc.) or use `--self-supervised` mode
6. **CheckM2 vs CheckM1**: CheckM2 uses ML, no marker gene database needed; CheckM1 needs 1.4 GB database. Different APIs

**Why the DAG is genuinely complex:**
- Dual assembler convergence (assembly merge/selection)
- Width-4 binner fan-out with DAS Tool genuine convergence (requires ALL binner outputs)
- Triple QC convergence (CheckM2 + GTDB-Tk + QUAST)
- Per-bin annotation creates fan-out after quality filtering
- N-to-M-to-N pattern: reads -> multiple assemblies -> multiple binners -> refined bins -> per-bin analysis
- Cross-branch: read depth from mapping feeds both binners and final coverage stats
- 5 convergence points is the maximum complexity tier

---

### TASK 13: DIA Proteomics Quantification

**Name:** `dia-proteomics`
**Description:** Quantify proteins from data-independent acquisition mass spectrometry data using library-based and library-free approaches with statistical differential analysis.

**DAG Diagram (depth=10, convergence=4, tools=9):**
```
 sample1.raw  sample2.raw  ...  library.raw (DDA)
     │            │                 │
 [ThermoRaw   [ThermoRaw       [ThermoRaw              Level 1
  FileParser]   FileParser]      FileParser]
     │            │                 │
     │            │          [OpenMS FeatureFinder       Level 2
     │            │           + ID (library build)]
     │            │                 │
     │            │          [OpenMS spectral            Level 3
     │            │           library generation]
     │            │                 │
     └────────────┼─────────────────┘
                  │
          [CONVERGENCE 1] ◄──────────────────           Level 4
          (DIA mzML + spectral library)
                  │
          ┌───────┼───────┐
          │       │       │
    [OpenSwath  [OpenMS  [OpenMS                        Level 5
     targeted    DIA-NN   FeatureFinderMRM
     extraction] wrapper] (untargeted)]
          │       │       │
    [pyprophet [pyprophet │                             Level 6
     score]     score]    │
          │       │       │
          └───────┼───────┘
                  │
          [CONVERGENCE 2] ◄──────────────────           Level 7
          (multi-engine FDR-controlled results)
          [OpenMS FeatureLinker]
                  │
          ┌───────┼───────┐
          │       │       │
    [MapAligner [Normalizer [QC stats]                  Level 8
     (RT align)] (TIC/      │
                  median)]  │
          │       │         │
          └───────┼─────────┘
                  │
          [CONVERGENCE 3] ◄──────────────────           Level 9
          (aligned + normalized + QC)
                  │
          [MSstats/python differential]
          [CONVERGENCE 4] ◄── design matrix             Level 10
          [report]
```

**Longest path:** ThermoRaw(DDA) -> FeatureFinder -> library build -> convergence1 -> OpenSwath -> pyprophet -> convergence2 -> MapAligner -> convergence3 -> MSstats -> report (depth 10)

**Tool List (9 tools):**
| Tool | Conda package | Purpose |
|------|--------------|---------|
| ThermoRawFileParser | `thermorawfileparser` | RAW to mzML conversion |
| OpenMS | `openms` | Feature finding, library building, alignment |
| OpenSwath | `openms` (included) | Targeted DIA extraction |
| pyprophet | `pyprophet` | FDR scoring (semi-supervised) |
| MSstats | `bioconductor-msstats` | Differential protein analysis |
| samtools | `samtools` | Index operations |
| Python/pandas | `pandas` | Data manipulation, reporting |
| R/ggplot2 | `r-ggplot2` | Volcano plots |
| multiqc | `multiqc` | QC aggregation |

**Data Source:** PRIDE PXD014414 (HeLa DIA benchmark, ~500 MB subset)
- OpenSwath tutorial data: http://openswath.org/en/latest/docs/tutorial.html
- Alternative: PXD004873 (LFQbench DIA)
- Galaxy DIA tutorial: Zenodo 4544489

**Data Size:** ~500 MB (3-4 DIA mzML files + 1 DDA library file + FASTA database)
**Estimated Runtime:** ~2.5h on 8 CPUs

**Key Domain-Specific Traps:**
1. **Spectral library is a separate branch**: DDA runs must be processed independently to build the library; agents often skip this and try library-free mode only
2. **pyprophet two-level FDR**: Must apply peptide-level FDR first, then protein-level; applying only one level gives inflated false discoveries
3. **RT alignment**: Retention time alignment (iRT) must happen BEFORE quantification; without it, DIA windows map to wrong peptides
4. **Decoy generation**: OpenSwath requires decoy transitions in the library; must run `OpenSwathDecoyGenerator` -- missing decoys = no FDR control
5. **MSstats input format**: Requires specific column format (ProteinName, PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge, IsotopeLabelType, Condition, BioReplicate, Run, Intensity); column name mismatches cause silent errors

**Why the DAG is genuinely complex:**
- DDA library branch runs in parallel with DIA sample processing, converging at extraction
- Three parallel extraction/scoring approaches
- pyprophet adds hidden depth within two branches
- Alignment + normalization + QC triple branch after scoring convergence
- MSstats requires design matrix input (cross-branch dependency from experimental metadata)

---

### TASK 14: Haplotype Phasing and Imputation

**Name:** `haplotype-phasing`
**Description:** Phase genotypes and impute missing variants from genotype array data using a reference panel, with quality assessment and haplotype-based association testing.

**DAG Diagram (depth=10, convergence=4, tools=10):**
```
 study.vcf.gz        reference_panel.vcf.gz
       │                     │
 [bcftools norm] ────── [bcftools view               Level 1
                         (subset chr)]
       │                     │
 [bcftools +fixploidy]  [bcftools +fixref             Level 2
       │                 (strand align)]
       │                     │
 ┌─────┴─────┐               │
 │           │               │
[bcftools  [plink2           │                        Level 3
 stats      --freq            │
 (pre-QC)]  --missing]       │
 │           │               │
 └─────┬─────┘               │
       │                     │
 [plink2 QC filters] ────── │                        Level 4
 (--geno 0.02 --maf 0.01)   │
       │                     │
       └──────────┬──────────┘
                  │
          [CONVERGENCE 1] ◄──────────────────        Level 5
          (QC'd study + reference panel aligned)
          [bcftools isec (shared sites)]
                  │
          ┌───────┼───────┐
          │       │       │
    [Eagle2     [SHAPEIT4 [bcftools                   Level 6
     phasing]    phasing]  concat
                           (scaffold)]
          │       │       │
          └───────┼───────┘
                  │
          [CONVERGENCE 2] ◄──────────────────        Level 7
          (phased haplotypes, best of 2 methods)
          [bcftools +compare phasing concordance]
                  │
          [Minimac4 imputation] ◄── ref panel         Level 8
                  │
          ┌───────┼───────────┐
          │       │           │
    [bcftools   [plink2     [python                   Level 9
     filter      --r2        imputation
     (info>0.3)] (LD check)] accuracy
                              (if truth)]
          │       │           │
          └───────┼───────────┘
                  │
          [CONVERGENCE 3] ◄──────────────────
          (filtered + LD + accuracy)
                  │
          [plink2 --glm association]
          [CONVERGENCE 4] ◄── pheno + PCs             Level 10
          [report]
```

**Longest path:** norm -> fixploidy -> freq -> QC -> convergence1 -> Eagle2 -> convergence2 -> Minimac4 -> filter -> convergence3 -> association -> report (depth 10)

**Tool List (10 tools):**
| Tool | Conda package | Purpose |
|------|--------------|---------|
| bcftools | `bcftools` | VCF manipulation, normalization, stats |
| PLINK2 | `plink2` | QC, frequency, association |
| Eagle2 | `eagle` | Statistical phasing |
| SHAPEIT4 | `shapeit4` | Statistical phasing |
| Minimac4 | `minimac4` | Genotype imputation |
| Beagle | `beagle` | Alternative phasing/imputation |
| samtools | `samtools` | Index operations |
| Python/pandas | `pandas` | Accuracy, reporting |
| R/qqman | `r-qqman` | Manhattan plots |
| multiqc | `multiqc` | QC aggregation |

**Data Source:** 1000 Genomes Phase 3 chr22
- Reference panel: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
- Study genotypes: HapMap3 chip data or simulated from 1KG
- Michigan Imputation Server test data

**Data Size:** ~300 MB (study VCF + reference panel VCF for chr22 + genetic map)
**Estimated Runtime:** ~2h on 8 CPUs

**Key Domain-Specific Traps:**
1. **Strand alignment**: Study and reference must be on same strand; strand flips cause systematic errors but no crashes -- only detected by checking A/T and C/G SNPs
2. **Chromosome naming**: `chr22` vs `22` mismatch between study and reference silently produces 0 imputed variants
3. **Genetic map**: Eagle and SHAPEIT require genetic map (cM positions); omitting it degrades phasing accuracy by 10-20%
4. **Imputation info score**: Must filter by info score (R-squared > 0.3); unfiltered imputed genotypes contain millions of low-quality variants
5. **Phasing before imputation**: Imputation quality degrades severely without pre-phasing; this is a mandatory ordering constraint

**Why the DAG is genuinely complex:**
- Dual-input (study genotypes + reference panel) with independent preprocessing
- Pre-QC stats and filtering have cross-dependency (freq must finish before filter thresholds)
- Dual phasing methods (Eagle + SHAPEIT) converge with comparison
- Post-imputation three-way branch (quality filter + LD check + accuracy)
- Imputation depends on both phased study data AND reference panel (second convergence of the dual inputs)
- Association testing requires phased/imputed genotypes + phenotypes + PCA (triple input convergence)

---

### TASK 15: Metatranscriptomics

**Name:** `metatranscriptomics`
**Description:** Analyze metatranscriptomic data to profile active microbial community composition and functional gene expression, with rRNA removal and dual taxonomic/functional assignment.

**DAG Diagram (depth=10, convergence=4, tools=11):**
```
 R1.fastq.gz     R2.fastq.gz
     │                │
 [fastp QC] ────── [fastp QC]                     Level 1
     │                │
     └───────┬────────┘
             │
 [SortMeRNA rRNA removal] ◄── rRNA databases      Level 2
             │
     ┌───────┴───────────────┐
     │                       │
 (non-rRNA reads)        (rRNA reads)
     │                       │
     │               [python rRNA                   Level 3
     │                community profile]
     │                       │
 ┌───┼──────────┐            │
 │   │          │            │
[STAR       [MetaPhlAn4     │                       Level 4
 align]      (taxonomy)]    │
 │   │          │            │
 │  [HUMAnN     │            │                      Level 5
 │   (func.     │            │
 │    profiling)]│            │
 │       │      │            │
 │   ┌───┴───┐  │            │
 │   │       │  │            │
 │ [gene    [pathway │       │                      Level 6
 │  families] abund]  │       │
 │   │       │  │     │       │
 │   └───┬───┘  │     │       │
 │       │      │     │       │
 └───────┼──────┘     │       │
         │            │       │
 [CONVERGENCE 1] ◄───┘       │                     Level 7
 (alignment + func + taxonomy)
         │                    │
 ┌───────┼───────────┐       │
 │       │           │       │
[HUMAnN [MetaPhlAn  [python  │                      Level 8
 regroup  merge      expr    │
 (GO/EC)] (species   filter) │
         table)]             │
 │       │           │       │
 └───────┼───────────┘       │
         │                   │
 [CONVERGENCE 2] ◄──────────┘                      Level 9
 (func groups + species + expr + rRNA profile)
         │
 ┌───────┼───────────┐
 │       │           │
[python  [python    [python                         Level 9
 diversity active    func
 (Shannon/ species   enrichment]
  Chao1)]  ranking]
 │       │           │
 └───────┼───────────┘
         │
 [CONVERGENCE 3] ◄──────────────────
 (diversity + active species + enrichment)
         │
 [CONVERGENCE 4] ◄── QC + rRNA%                    Level 10
 [python report]
```

**Longest path:** fastp -> SortMeRNA -> STAR -> HUMAnN -> pathway abundance -> convergence1 -> HUMAnN regroup -> convergence2 -> enrichment -> convergence3 -> report (depth 10)

**Tool List (11 tools):**
| Tool | Conda package | Purpose |
|------|--------------|---------|
| fastp | `fastp` | Read QC |
| SortMeRNA | `sortmerna` | rRNA removal |
| STAR | `star` | Genome alignment |
| MetaPhlAn | `metaphlan` | Taxonomic profiling |
| HUMAnN | `humann` | Functional profiling |
| samtools | `samtools` | BAM operations |
| bedtools | `bedtools` | Coverage operations |
| Python/scipy | `scipy` | Diversity, enrichment |
| Python/pandas | `pandas` | Data manipulation |
| multiqc | `multiqc` | QC aggregation |
| R/vegan | `r-vegan` | Ecological diversity |

**Data Source:** HMP2 metatranscriptomics (subset)
- IBDMDB: https://ibdmdb.org/ (SRR5947006)
- Galaxy metatranscriptomics tutorial: Zenodo 4776250
- nf-core does not have a dedicated metatranscriptomics pipeline

**Data Size:** ~400 MB (PE FASTQ + rRNA databases + MetaPhlAn + HUMAnN databases subset)
**Estimated Runtime:** ~2.5h on 8 CPUs

**Key Domain-Specific Traps:**
1. **rRNA removal is mandatory**: 80-95% of metatranscriptome reads are rRNA; skipping SortMeRNA makes all downstream quantification meaningless
2. **rRNA as community profile**: The rRNA reads themselves provide 16S-like community composition -- agents must process BOTH fractions
3. **HUMAnN requires MetaPhlAn**: HUMAnN internally runs MetaPhlAn for community-stratified functional profiling; must provide MetaPhlAn database path
4. **Normalization**: Gene family abundances must be normalized to copies-per-million (CPM); raw RPKs are not comparable across samples
5. **DNA contamination**: Metatranscriptomics can have DNA carryover; without DNase treatment info, must note this caveat

**Why the DAG is genuinely complex:**
- rRNA/non-rRNA split creates fundamental two-track processing (both informative)
- non-rRNA track branches into alignment vs taxonomy vs function
- HUMAnN depends on MetaPhlAn taxonomy (hidden dependency within functional branch)
- Convergence requires functional groups + species tables + expression + rRNA community
- Three-way diversity/ranking/enrichment branch from merged data
- Cross-branch: rRNA percentage feeds QC AND community comparison with MetaPhlAn taxonomy

---

## MEDIUM-LOW PRIORITY (Specialized): Tasks 16-20

---

### TASK 16: Quantitative Proteomics (DDA Label-Free, Dual Search Engine)

**NOTE:** `radseq-popgen` already exists in the benchmark. Replacing with DDA-LFQ proteomics which is uncovered.

**Name:** `dda-lfq-proteomics`
**Description:** Quantify proteins from data-dependent acquisition label-free mass spectrometry using dual search engines (Comet + MS-GF+), FDR control with Percolator, and differential protein abundance analysis.

**DAG Diagram (depth=10, convergence=4, tools=10):**
```
 sample1.raw  sample2.raw  sample3.raw  database.fasta
     │            │            │              │
 [ThermoRaw   [ThermoRaw   [ThermoRaw   [DecoyDB       Level 1
  FileParser]   FileParser]   FileParser]  Generator]
     │            │            │              │
 [OpenMS PeakPickerHiRes] ── per file         │         Level 2
     │            │            │              │
     └────────────┼────────────┘              │
                  │                           │
     ┌────────────┴────────────┐              │
     │                         │              │
 [Comet search]          [MS-GF+ search] ◄───┘         Level 3
     │                         │
 [OpenMS IDFileConverter] [OpenMS IDFileConverter]      Level 4
     │                         │
     └────────────┬────────────┘
                  │
          [CONVERGENCE 1] ◄──────────────────          Level 5
          [OpenMS ConsensusID / PeptideIndexer]
                  │
          [OpenMS PSMFeatureExtractor]                  Level 6
                  │
          [Percolator FDR control]                      Level 7
                  │
          ┌───────┼───────────┐
          │       │           │
    [OpenMS    [OpenMS      [OpenMS                    Level 8
     ProteomicsLFQ FeatureFinderIdentification
     (quant)]    (intensity)]  IDFilter
                              (1% FDR)]
          │       │           │
          └───────┼───────────┘
                  │
          [CONVERGENCE 2] ◄──────────────────          Level 9
          (quantified + filtered protein groups)
                  │
          ┌───────┼───────────┐
          │       │           │
    [MSstats  [python      [python                     Level 9
     diff      protein      QC stats
     analysis] coverage]    (ID rates)]
          │       │           │
          └───────┼───────────┘
                  │
          [CONVERGENCE 3] ◄──────────────────
                  │
          [CONVERGENCE 4] ◄── QC per-sample            Level 10
          [report]
```

**Longest path:** ThermoRaw -> PeakPicker -> Comet -> IDFileConverter -> ConsensusID -> PSMFeatureExtractor -> Percolator -> ProteomicsLFQ -> MSstats -> report (depth 10)

**Tool List (10 tools):**
| Tool | Conda package | Purpose |
|------|--------------|---------|
| ThermoRawFileParser | `thermorawfileparser` | RAW to mzML conversion |
| OpenMS | `openms` | Peak picking, ID merging, quantification |
| Comet | `comet-ms` | Database search engine 1 |
| MS-GF+ | `msgf_plus` | Database search engine 2 |
| Percolator | `percolator` | FDR control (semi-supervised) |
| MSstats | `bioconductor-msstats` | Differential protein analysis |
| Python/pandas | `pandas` | Data manipulation |
| R/ggplot2 | `r-ggplot2` | Volcano, coverage plots |
| samtools | `samtools` | FASTA operations |
| multiqc | `multiqc` | QC aggregation |

**Data Source:** PRIDE PXD000001 (TMT but can use LFQ) or PXD028735
- nf-core/quantms test data
- Galaxy proteomics tutorial: Zenodo 1489208
- Alternative: PXD004682 (LFQ benchmark, UPS1 spike-in)

**Data Size:** ~400 MB (3-4 RAW files + FASTA database)
**Estimated Runtime:** ~2.5h on 8 CPUs

**Key Domain-Specific Traps:**
1. **Dual search engine merging**: Comet and MS-GF+ use different scoring functions; must use ConsensusID or PeptideIndexer to properly merge before FDR control
2. **Target-decoy approach**: Must generate reversed decoy database BEFORE search; searching without decoys = no FDR control
3. **Percolator features**: PSMFeatureExtractor must run BEFORE Percolator; Percolator without proper features gives suboptimal FDR estimation
4. **LFQ normalization**: Label-free quantification requires intensity normalization across runs; unnormalized data shows run-to-run technical variation, not biology
5. **Missing values**: LFQ data has many missing values (50-70%); must use imputation (MinDet, KNN) or handle in MSstats model

**Why the DAG is genuinely complex:**
- Dual search engine fan-out with ConsensusID convergence (genuine multi-method agreement)
- Percolator creates depth within FDR control branch
- Post-FDR triple branch (quantification, feature intensity, filtering)
- MSstats requires both quantified proteins AND experimental design (cross-input dependency)
- Each sample is independently processed then all converge at quantification (N-to-1 pattern)

---

### TASK 17: Pharmacogenomics (CYP2D6 etc.)

**Name:** `pharmacogenomics`
**Description:** Determine pharmacogenomic star alleles from WGS data, focusing on CYP2D6 and other CPIC-defined genes, with drug interaction reporting.

**DAG Diagram (depth=9, convergence=4, tools=10):**
```
 sample.R1.fq.gz    sample.R2.fq.gz
       │                  │
   [fastp QC] ──────── [fastp QC]                  Level 1
       │                  │
       └────────┬─────────┘
                │
        [bwa-mem2 align]                            Level 2
        (full genome)
                │
        [samtools sort + index]                     Level 3
                │
        [picard MarkDuplicates]                     Level 4
                │
        ┌───────┼───────────────────┐
        │       │                   │
   [gatk HC   [samtools           [mosdepth          Level 5
    (CYP region  view               (CYP2D6
     intervals)] (extract           region
        │        CYP2D6)]           coverage)]
        │       │                   │
   [bcftools  [samtools             │
    norm]      depth per-base]      │
        │       │                   │
        │  [CONVERGENCE 1] ◄────────┘               Level 6
        │  (CYP2D6 BAM + coverage)
        │       │
        │  [python CYP2D6 CN estimation]             Level 7
        │  (read depth ratio vs CYP2D7)
        │       │
        └───────┴───────┐
                        │
                [CONVERGENCE 2] ◄────────            Level 8
                (VCF + CN + coverage)
                        │
        ┌───────────────┼───────────────┐
        │               │               │
  [python star-    [bcftools         [python          Level 8
   allele calling   annotate         drug-gene
   (diplotype)]     (PharmVar        interaction
        │            database)]       lookup]
        │               │               │
        └───────────────┼───────────────┘
                        │
                [CONVERGENCE 3] ◄────────
                (diplotype + annotation + drugs)
                        │
                ┌───────┼───────┐
                │       │       │
          [python    [python   [python               Level 9
           metabolizer activity phenotype
           status]    score]    prediction]
                │       │       │
                └───────┼───────┘
                        │
                [CONVERGENCE 4] ◄── coverage QC
                [clinical report]
```

**Longest path:** fastp -> bwa -> sort -> markdup -> samtools view CYP2D6 -> depth -> convergence1 -> CN estimation -> convergence2 -> star-allele -> convergence3 -> metabolizer status -> report (depth 9 + sub-steps)

**Tool List (10 tools):**
| Tool | Conda package | Purpose |
|------|--------------|---------|
| fastp | `fastp` | Read QC |
| bwa-mem2 | `bwa-mem2` | Alignment |
| samtools | `samtools` | BAM operations, depth |
| picard | `picard` | Duplicate marking |
| GATK4 | `gatk4` | Variant calling (targeted) |
| bcftools | `bcftools` | VCF operations, annotation |
| mosdepth | `mosdepth` | Coverage analysis |
| bedtools | `bedtools` | Region extraction |
| Python/pandas | `pandas` | Star allele calling, reporting |
| Python/scipy | `scipy` | CN estimation statistics |

**Data Source:** GIAB HG002 WGS (CYP2D6 region, chr22:42126000-42132000)
- 1000 Genomes high-coverage WGS: ftp.1000genomes.ebi.ac.uk
- PharmVar database: https://www.pharmvar.org/
- CPIC guidelines: https://cpicpgx.org/

**Data Size:** ~500 MB (WGS FASTQ chr22 + reference + PharmVar VCF + CPIC tables)
**Estimated Runtime:** ~2h on 8 CPUs

**Key Domain-Specific Traps:**
1. **CYP2D6/CYP2D7 homology**: CYP2D6 has a highly homologous pseudogene (CYP2D7) nearby; reads mismap between them. Must use region-specific extraction and careful MAPQ filtering
2. **Copy number variation**: CYP2D6 has common deletions (*5) and duplications (*1xN, *2xN); pure SNP-based calling misses these entirely
3. **Star allele nomenclature**: PGx uses star alleles (*1, *2, *4, etc.) not rsIDs; must translate VCF variants to PharmVar star allele definitions
4. **Diplotype to phenotype**: *1/*4 = extensive/poor metabolizer heterozygote; must consult CPIC translation tables for activity scoring
5. **Phase matters**: CYP2D6 *4 + *10 on same chromosome (cis) vs different chromosomes (trans) give different metabolizer predictions

**Why the DAG is genuinely complex:**
- Three parallel extractions from same BAM (variant calling, region extraction, coverage)
- Copy number estimation requires BOTH CYP2D6 depth AND CYP2D7 depth (cross-gene comparison)
- Star allele calling requires BOTH VCF (SNP-based) AND copy number (structural)
- Drug interaction lookup is a separate branch requiring curated database
- Metabolizer status requires diplotype + activity score + phenotype prediction (triple convergence)
- Pharmacogenomics-specific terminology creates domain knowledge barrier

---

### TASK 18: Single-Cell ATAC-seq

**Name:** `scatac-seq`
**Description:** Analyze single-cell ATAC-seq data for chromatin accessibility profiling at single-cell resolution, including cell clustering, peak calling, motif enrichment, and gene activity scoring.

**DAG Diagram (depth=10, convergence=4, tools=11):**
```
 fragments.tsv.gz (or FASTQ + barcode)
          │
    [python/snapatac2 import] ◄── barcodes         Level 1
          │
    ┌─────┴──────────┐
    │                │
 [snapatac2        [python                          Level 2
  qc_metrics]       barcode stats]
  (TSS enrich,
   unique frags)
    │                │
    └────────┬───────┘
             │
     [snapatac2 filter cells] ◄── QC thresholds     Level 3
             │
     [snapatac2 add_tile_matrix]                     Level 4
     (genome-wide bins)
             │
     ┌───────┼───────────┐
     │       │           │
 [snapatac2 [snapatac2  [snapatac2                   Level 5
  select_    spectral     scrublet
  features]  (dim red)]   (doublets)]
     │       │           │
     └───────┼───────────┘
             │
     [CONVERGENCE 1] ◄──────────────────             Level 6
     (features + embedding + clean cells)
             │
     [snapatac2 leiden clustering]                   Level 7
             │
     ┌───────┼───────────────────┐
     │       │                   │
 [MACS2    [snapatac2          [snapatac2            Level 8
  callpeak  gene_matrix]        export to
  per-       (gene activity     AnnData)]
  cluster]   scoring)]
     │       │                   │
     │  [CONVERGENCE 2] ◄────────┘
     │  (clusters + gene activity + AnnData)
     │       │
     └───────┴───────┐
                     │
             [CONVERGENCE 3] ◄──────────────         Level 9
             (peaks + clusters + gene activity)
                     │
             ┌───────┼───────────┐
             │       │           │
       [HOMER    [deeptools   [python                Level 9
        findMotifs computeMatrix  diff
        Enrichment (peak         accessibility
        (per-       heatmap)]    per cluster)]
        cluster)]
             │       │           │
             └───────┼───────────┘
                     │
             [CONVERGENCE 4] ◄── QC stats            Level 10
             [python report]
```

**Longest path:** import -> QC -> filter -> tile_matrix -> spectral -> convergence1 -> leiden -> MACS2 per-cluster -> convergence3 -> HOMER motif -> convergence4 -> report (depth 10)

**Tool List (11 tools):**
| Tool | Conda package | Purpose |
|------|--------------|---------|
| SnapATAC2 | `snapatac2` | Core scATAC analysis |
| MACS2 | `macs2` | Peak calling per cluster |
| HOMER | `homer` | Motif enrichment |
| deeptools | `deeptools` | Heatmaps |
| samtools | `samtools` | BAM operations |
| bedtools | `bedtools` | Interval operations |
| Python/scanpy | `scanpy` | Visualization (UMAP) |
| Python/anndata | `anndata` | Data structure |
| Python/pandas | `pandas` | Data manipulation |
| Python/matplotlib | `matplotlib` | Plotting |
| multiqc | `multiqc` | QC aggregation |

**Data Source:** 10x Genomics PBMC 5k scATAC-seq
- https://www.10xgenomics.com/resources/datasets/5-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-next-gem-v-1-1-1-1-standard-2-0-0
- Fragment file (~200 MB) + peaks BED
- Alternative: ArchR tutorial dataset

**Data Size:** ~300 MB (fragments.tsv.gz + reference + peaks BED + gene annotation)
**Estimated Runtime:** ~2h on 8 CPUs

**Key Domain-Specific Traps:**
1. **TSS enrichment threshold**: Must compute TSS enrichment score per cell and filter (typical cutoff: TSS >= 4); low TSS enrichment cells are dead/debris
2. **Doublet detection**: scATAC doublets are harder to detect than scRNA; Scrublet adapted for ATAC uses accessibility profiles
3. **Per-cluster peak calling**: Must call peaks PER CLUSTER (not bulk); bulk peaks miss rare cell type-specific regulatory elements
4. **Gene activity scoring**: Gene activity != gene expression; computed from accessibility in gene body + promoter, requires correct TSS annotation
5. **Fragment file format**: Assumes pre-aligned fragments.tsv.gz (10x format: chr, start, end, barcode, count); agents must not try to re-align

**Why the DAG is genuinely complex:**
- Feature selection + dimensionality reduction + doublet detection triple branch
- Per-cluster peak calling requires clustering results (dependency chain)
- Gene activity scoring requires both clustering AND genomic annotation
- Motif enrichment + heatmap + differential accessibility triple final branch
- Cross-branch: QC metrics from early filtering feed into final report
- Iterative potential: clustering quality might require re-embedding (not modeled but domain-specific)

---

### TASK 19: Genome Scaffolding with Long Reads

**Name:** `genome-scaffolding`
**Description:** Scaffold a fragmented short-read assembly using long reads through multiple scaffolding approaches, compare results, and assess improvement.

**DAG Diagram (depth=10, convergence=4, tools=10):**
```
 contigs.fasta (short-read assembly)    long_reads.fastq.gz
         │                                    │
     [QUAST initial assessment]          [NanoPlot QC]          Level 1
         │                                    │
         └──────────────┬─────────────────────┘
                        │
                [CONVERGENCE 1] ◄──────────────                Level 2
                (contigs + long reads available)
                        │
        ┌───────────────┼───────────────┐
        │               │               │
   [minimap2 align  [minimap2 align [ntLink                    Level 3
    + RagTag          + LINKS         scaffolding]
    scaffold]         scaffolding]
        │               │               │
   [RagTag         [abyss-scaffold     │                       Level 4
    patch]           (LINKS output)]    │
        │               │               │
        └───────────────┼───────────────┘
                        │
                [CONVERGENCE 2] ◄──────────────                Level 5
                (3 scaffold sets)
                [python select best / merge]
                        │
                ┌───────┼───────────┐
                │       │           │
          [QUAST     [minimap2   [BUSCO                        Level 6
           scaffold   re-align   (completeness)]
           assessment] reads]
                │       │           │
                └───────┼───────────┘
                        │
                [CONVERGENCE 3] ◄──────────────                Level 7
                (QUAST + mapping + BUSCO)
                        │
                ┌───────┼───────────┐
                │       │           │
          [bedtools [python       [python                      Level 8
           getfasta  N-gap         scaffold
           (gap      analysis]     comparison
           flanks)]                vs initial]
                │       │           │
                └───────┼───────────┘
                        │
                [CONVERGENCE 4] ◄──────────────                Level 9
                (gaps + comparison + QC)
                        │
                [python report]                                 Level 10
```

**Longest path:** QUAST initial -> convergence1 -> minimap2 + RagTag scaffold -> RagTag patch -> convergence2 -> QUAST scaffold -> convergence3 -> gap analysis -> convergence4 -> report (depth 10)

**Tool List (10 tools):**
| Tool | Conda package | Purpose |
|------|--------------|---------|
| minimap2 | `minimap2` | Long-read alignment |
| RagTag | `ragtag` | Reference-guided scaffolding |
| LINKS | `links` | Long-read scaffolding |
| ntLink | `ntlink` | Minimizer-based scaffolding |
| ABySS | `abyss` | Scaffold module |
| QUAST | `quast` | Assembly assessment |
| BUSCO | `busco` | Gene completeness |
| samtools | `samtools` | BAM/FASTA operations |
| bedtools | `bedtools` | Gap analysis |
| Python/pandas | `pandas` | Reporting |

**Data Source:** E. coli K12 short-read assembly + ONT long reads
- Short-read assembly: SPAdes output from SRR1770413 (already assembled)
- Long reads: SRR15058234 (ONT E. coli, subset to 50x)
- Alternative: nf-core test assembly + simulated long reads

**Data Size:** ~300 MB (short-read contigs + long reads + reference for comparison)
**Estimated Runtime:** ~1.5h on 8 CPUs

**Key Domain-Specific Traps:**
1. **Scaffolding vs gap-filling**: RagTag scaffold joins contigs, RagTag patch fills gaps -- these are different operations that must be run sequentially
2. **minimap2 preset**: Must use `-ax map-ont` for ONT, `-ax map-hifi` for HiFi; wrong preset = poor alignment = poor scaffolding
3. **ntLink k-mer size**: Default k=32 works for bacterial; larger genomes need k=40+. Wrong k = fragmented scaffolds or chimeric joins
4. **LINKS parameters**: `-d` (distance) and `-k` (k-mer) are critical; agents often use defaults that are inappropriate for the read length distribution
5. **Misassembly introduction**: Scaffolders can CREATE misassemblies; must check with QUAST `--mis-size` and compare misassembly counts before/after

**Why the DAG is genuinely complex:**
- Initial assessment must complete before scaffolding (baseline metrics needed for comparison)
- Three parallel scaffolding approaches with different algorithmic strategies
- Some scaffolders (RagTag) have internal 2-step chains (scaffold -> patch)
- Triple QC after scaffolding (QUAST + mapping + BUSCO)
- Gap analysis branch requires scaffold FASTA + original contigs (cross-branch dependency)
- Final comparison needs BOTH initial and scaffold QUAST results

---

### TASK 20: eDNA Environmental Metabarcoding

**Name:** `edna-metabarcoding`
**Description:** Analyze environmental DNA metabarcoding data for aquatic biodiversity assessment, from raw amplicon reads through species identification and community ecology analysis.

**DAG Diagram (depth=10, convergence=4, tools=10):**
```
 R1.fastq.gz    R2.fastq.gz (multiplexed)
     │               │
     └───────┬───────┘
             │
 [cutadapt demultiplex] ◄── barcodes              Level 1
             │
 [cutadapt primer removal]                         Level 2
             │
 ┌───────────┼───────────┐
 │           │           │
[vsearch   [vsearch    [fastqc                     Level 3
 --merge    --filter     QC]
 pairs]     (quality)]
 │           │           │
 └─────┬─────┘           │
       │                 │
 [vsearch dereplicate]   │                         Level 4
       │                 │
 ┌─────┼─────┐           │
 │     │     │           │
[vsearch [swarm  [UNOISE3                          Level 5
 --cluster clustering] (denoise)]
 97%]
 │     │     │
 └─────┴──┬──┘
          │
  [CONVERGENCE 1] ◄──────────────────              Level 6
  (OTUs/ASVs consensus)
  [vsearch --uchime chimera removal]
          │
  ┌───────┼───────────┐
  │       │           │
[BLAST  [vsearch    [obitools                      Level 7
 (NCBI   --usearch   ecotag]
 nt)]    global
         (custom DB)]
  │       │           │
  └───────┴─────┬─────┘
                │
        [CONVERGENCE 2] ◄──────────────────        Level 8
        (taxonomy from multiple methods)
        [LCA consensus]
                │
        ┌───────┼───────────┐
        │       │           │
  [python    [python      [python                  Level 9
   species    diversity    occupancy
   list +     (alpha +     modeling
   IUCN       beta)]       (detection
   status]                  prob)]
        │       │           │
        └───────┼───────────┘
                │
        [CONVERGENCE 3] ◄──────────────────
                │
        [CONVERGENCE 4] ◄── QC + negative controls Level 10
        [python report]
```

**Longest path:** cutadapt demux -> primer removal -> vsearch merge -> dereplicate -> vsearch cluster -> convergence1 -> chimera removal -> BLAST -> convergence2 -> diversity -> convergence3 -> report (depth 10)

**Tool List (10 tools):**
| Tool | Conda package | Purpose |
|------|--------------|---------|
| cutadapt | `cutadapt` | Demultiplexing, primer trimming |
| vsearch | `vsearch` | Merging, filtering, clustering, chimera removal |
| SWARM | `swarm` | Alternative OTU clustering |
| BLAST | `blast` | Taxonomic assignment (NCBI) |
| OBITools | `obitools` | eDNA-specific taxonomy (ecotag) |
| FastQC | `fastqc` | Read quality |
| samtools | `samtools` | Sequence operations |
| Python/pandas | `pandas` | Data manipulation |
| Python/scipy | `scipy` | Diversity, occupancy modeling |
| R/vegan | `r-vegan` | Community ecology |

**Data Source:** NCBI eDNA from aquatic monitoring
- Zenodo 4741491 (fish eDNA metabarcoding, 12S/COI)
- Riaz et al. 2011 dataset
- Alternative: DRYAD doi:10.5061/dryad.5x69p8d0h

**Data Size:** ~200 MB (multiplexed amplicon FASTQ + reference DB + negative controls)
**Estimated Runtime:** ~1h on 8 CPUs

**Key Domain-Specific Traps:**
1. **Primer removal is not adapter trimming**: eDNA primers (12S MiFish, COI Leray) have specific sequences that must be removed with exact matching; generic adapter trimming leaves primer artifacts
2. **Negative control processing**: eDNA studies MUST include field blanks and extraction blanks; species appearing in negatives must be flagged or removed
3. **OTU vs ASV choice**: 97% OTU clustering vs denoising (UNOISE3) gives different resolution; fish eDNA typically uses OTUs because intraspecific variation is low
4. **Chimera removal after clustering**: Chimera detection must happen on dereplicated sequences (abundance-aware); running on raw reads misses abundance-skewed chimeras
5. **LCA (Lowest Common Ancestor)**: When BLAST returns multiple hits, must use LCA to assign taxonomy at the appropriate level, not just take top hit
6. **Detection probability**: eDNA detection is stochastic; occupancy modeling requires replicate PCRs per sample

**Why the DAG is genuinely complex:**
- Merge -> filter -> dereplicate -> cluster is a 4-step sequential chain
- Three parallel clustering methods (vsearch OTU, SWARM, UNOISE3) converge
- Three parallel taxonomy methods (BLAST, vsearch global, ecotag) converge with LCA
- Post-taxonomy: species list + diversity + occupancy modeling triple branch
- Negative controls create a cross-branch dependency that affects taxonomy filtering
- OTU/ASV table feeds into multiple downstream analyses (cross-branch fan-out)

---

## LOWER PRIORITY (Niche but Important): Tasks 21-25

---

### TASK 21: Viral Phylodynamics (Molecular Clock)

**Name:** `viral-phylodynamics`
**Description:** Perform molecular clock analysis on viral sequences to estimate evolutionary rates, divergence times, and epidemic dynamics using Bayesian methods.

**DAG Diagram (depth=10, convergence=4, tools=10):**
```
 sequences.fasta     metadata.tsv (dates, locations)
       │                   │
 [mafft align] ──── [python parse dates]            Level 1
       │                   │
 ┌─────┼───────┐           │
 │     │       │           │
[trimal [IQ-TREE [python   │                        Level 2
 gap     model   clock     │
 filter] test]   signal    │
 │     │       (root-to-   │
 │     │       tip)]       │
 │     │       │           │
 └─────┴───┬───┘           │
           │               │
   [CONVERGENCE 1] ◄───────┘                        Level 3
   (cleaned alignment + model + dates)
           │
   [IQ-TREE ML tree]                                Level 4
   (best model, bootstrap)
           │
   [TreeTime molecular clock]                       Level 5
   (--clock-filter)
           │
   ┌───────┼───────────────┐
   │       │               │
[TreeTime [TreeTime     [TreeTime                    Level 6
 mugration  skyline       ancestral
 (geo)]     (Neeff)]      (reconstruct)]
   │       │               │
   └───────┴───────┬───────┘
                   │
           [CONVERGENCE 2] ◄────────                Level 7
           (migration + skyline + ancestral)
                   │
   ┌───────────────┼───────────────┐
   │               │               │
[augur         [python          [python              Level 8
 export         transmission    rate
 (Auspice       chain           estimation
  JSON)]        reconstruction]  report]
   │               │               │
   └───────────────┼───────────────┘
                   │
           [CONVERGENCE 3] ◄────────                Level 9
           (Auspice + chains + rates)
                   │
           [CONVERGENCE 4] ◄── root-to-tip QC       Level 10
           [python report]
```

**Longest path:** mafft -> trimal -> convergence1 -> IQ-TREE -> TreeTime clock -> TreeTime skyline -> convergence2 -> transmission chain -> convergence3 -> report (depth 10)

**Tool List (10 tools):**
| Tool | Conda package | Purpose |
|------|--------------|---------|
| MAFFT | `mafft` | Multiple sequence alignment |
| trimAl | `trimal` | Alignment trimming |
| IQ-TREE | `iqtree` | ML phylogenetics, model selection |
| TreeTime | `treetime` | Molecular clock, skyline |
| Augur | `augur` | Nextstrain pipeline tools |
| BEAST2 | `beast2` | Bayesian phylodynamics (optional) |
| Python/pandas | `pandas` | Data manipulation |
| Python/matplotlib | `matplotlib` | Plotting |
| samtools | `samtools` | FASTA operations |
| bcftools | `bcftools` | Variant operations |

**Data Source:** NCBI Virus SARS-CoV-2 sequences (100-200 sequences)
- Nextstrain Zika tutorial: https://docs.nextstrain.org/en/latest/tutorials/zika.html
- Nextstrain SARS-CoV-2 example data
- Alternative: HIV-1 subtype B sequences from Los Alamos

**Data Size:** ~50 MB (FASTA sequences + metadata TSV)
**Estimated Runtime:** ~1.5h on 8 CPUs (BEAST2 is compute-intensive, IQ-TREE moderate)

**Key Domain-Specific Traps:**
1. **Date format**: TreeTime requires decimal dates (2020.5) not calendar dates (2020-07-01); wrong format = clock analysis failure
2. **Root-to-tip regression**: Must check temporal signal (R-squared > 0.2) BEFORE running molecular clock; no signal = clock model is inappropriate
3. **Clock model selection**: Strict vs relaxed clock; SARS-CoV-2 uses strict clock, HIV uses relaxed -- wrong choice invalidates all rate estimates
4. **Outgroup rooting**: Tree must be rooted correctly for clock analysis; midpoint rooting often gives wrong root for viruses with recombination
5. **Recombination**: Must screen for recombination before tree building; recombinant sequences create artifactually long branches

**Why the DAG is genuinely complex:**
- Alignment cleaning has triple branch (gap filter, model test, clock signal) all needed for ML tree
- ML tree -> molecular clock is sequential dependency
- Three parallel TreeTime analyses (migration, skyline, ancestral reconstruction)
- Post-clock triple branch (Auspice export, transmission chains, rate estimation)
- Root-to-tip regression from early stages feeds into final QC assessment (long cross-branch dependency)
- Metadata (dates, locations) feeds into multiple stages (clock, mugration, augur export)

---

### TASK 22: Microsatellite Instability Detection

**Name:** `msi-detection`
**Description:** Detect microsatellite instability from paired tumor-normal WGS/WES data using multiple MSI callers and correlate with mutation burden and mismatch repair gene status.

**DAG Diagram (depth=10, convergence=4, tools=10):**
```
 tumor.bam              normal.bam
     │                      │
 [samtools index] ──── [samtools index]             Level 1
     │                      │
     └──────────┬───────────┘
                │
        [CONVERGENCE 1] ◄──────────────────        Level 2
        (T+N indexed BAMs)
                │
     ┌──────────┼──────────────────┐
     │          │                  │
 [msisensor-  [msisensor2      [mantis               Level 3
  pro scan     (paired          (kmer-based
  + msi]       analysis)]       MSI scoring)]
     │          │                  │
     └──────────┴──────────┬───────┘
                           │
                   [CONVERGENCE 2] ◄──────────      Level 4
                   (3-caller MSI results)
                   [python consensus scoring]
                           │
                ┌──────────┼──────────┐
                │          │          │
          [gatk HC     [mosdepth   [bcftools          Level 5
           (MMR gene    (MSI loci   view
           regions)]    coverage)]  (MMR gene
                │          │        mutations)]
                │          │          │
          [bcftools    │          │                    Level 6
           filter      │          │
           (PASS)]     │          │
                │      │          │
                └──────┼──────────┘
                       │
               [CONVERGENCE 3] ◄──────────────      Level 7
               (MSI + MMR mutations + coverage)
                       │
               ┌───────┼───────────┐
               │       │           │
         [python    [python      [python              Level 8
          TMB        MSI-H       MMR gene
          estimation classification report]
          (from VCF)] (scoring)]
               │       │           │
               └───────┼───────────┘
                       │
               [CONVERGENCE 4] ◄── T/N QC            Level 9
                       │
               [python clinical report]               Level 10
```

**Longest path:** index -> convergence1 -> msisensor-pro -> convergence2 -> gatk HC(MMR) -> bcftools filter -> convergence3 -> TMB -> convergence4 -> report (depth 10)

**Tool List (10 tools):**
| Tool | Conda package | Purpose |
|------|--------------|---------|
| msisensor-pro | `msisensor-pro` | MSI detection (distribution-based) |
| msisensor2 | `msisensor2` | MSI detection (paired analysis) |
| MANTIS | `mantis-msi` | MSI detection (kmer-based) |
| GATK4 | `gatk4` | MMR gene variant calling |
| samtools | `samtools` | BAM operations |
| bcftools | `bcftools` | VCF operations |
| mosdepth | `mosdepth` | Coverage at MSI loci |
| bedtools | `bedtools` | Interval operations |
| Python/pandas | `pandas` | Scoring, reporting |
| multiqc | `multiqc` | QC aggregation |

**Data Source:** TCGA MSI-H sample or SEQC2 reference
- CPTAC MSI reference: https://cptac-data-portal.georgetown.edu/
- nf-core/sarek includes MSIsensor
- Alternative: cell line with known MSI status (HCT116 = MSI-H)

**Data Size:** ~500 MB (T+N BAM subset for target regions + reference + MSI loci BED)
**Estimated Runtime:** ~2h on 8 CPUs

**Key Domain-Specific Traps:**
1. **MSI loci BED file**: msisensor-pro requires a homopolymer/microsatellite loci file generated by `msisensor-pro scan`; using wrong reference genome version gives zero loci
2. **MSI-H vs MSS threshold**: msisensor-pro uses 3.5% (WES) or 10% (WGS) threshold; wrong threshold for data type = wrong classification
3. **MMR genes**: Must specifically check MLH1, MSH2, MSH6, PMS2; agents may not know which genes define mismatch repair deficiency
4. **Tumor purity**: Low tumor purity dilutes MSI signal; must either correct for purity or flag low-purity samples
5. **TMB correlation**: MSI-H tumors typically have TMB > 10 mut/Mb; discordant MSI/TMB should be flagged as suspicious

**Why the DAG is genuinely complex:**
- Triple MSI caller with consensus scoring (convergence 2)
- Parallel MMR gene analysis (variant calling, coverage, mutation filtering) creates nested diamond
- MSI classification + TMB estimation + MMR gene status must all converge for clinical interpretation
- Cross-branch: tumor purity estimation from initial BAM stats affects MSI thresholding
- Clinical decision tree: MSI-H requires BOTH high MSI score AND/OR MMR gene mutations

---

### TASK 23: Repeat Element Analysis (TE Annotation)

**Name:** `repeat-element-annotation`
**Description:** Annotate and classify transposable elements in a genome assembly, estimate TE activity, and analyze TE-gene interactions.

**DAG Diagram (depth=10, convergence=4, tools=10):**
```
 genome.fasta
       │
 ┌─────┼───────────────────┐
 │     │                   │
[RepeatModeler  [RepeatMasker    [python                Level 1-3
 (de novo TE     (known TE       genome stats
  library        library)]       (GC, size)]
  building)]         │
 │(slow: 2h)     [parse RM                             Level 4
 │                output]
 │                   │
 └───────┬───────────┘
         │
 [CONVERGENCE 1] ◄──────────────────────               Level 5
 (de novo + known TE libraries merged)
 [cat + cd-hit dedup libraries]
         │
 [RepeatMasker (comprehensive)]                         Level 6
 (merged library)
         │
 ┌───────┼───────────────────┐
 │       │                   │
[perl   [bedtools          [python                      Level 7
 parseRM  intersect          TE age
 landscape] (TE vs genes)]   distribution
 │       │                   (Kimura)]
 │       │                   │
 └───────┴─────────┬─────────┘
                   │
           [CONVERGENCE 2] ◄────────                    Level 8
           (landscape + gene proximity + age)
                   │
           ┌───────┼───────────┐
           │       │           │
     [python    [python      [python                    Level 9
      TE class   gene         solo-LTR
      summary    disruption   detection]
      (SINE/     analysis]
       LINE/
       LTR/DNA)]
           │       │           │
           └───────┼───────────┘
                   │
           [CONVERGENCE 3] ◄────────
                   │
           [CONVERGENCE 4] ◄── genome stats             Level 10
           [python report]
```

**Longest path:** RepeatModeler (de novo, slowest) -> convergence1 -> RepeatMasker comprehensive -> parseRM landscape -> convergence2 -> TE class summary -> convergence3 -> report (depth 10)

**Tool List (10 tools):**
| Tool | Conda package | Purpose |
|------|--------------|---------|
| RepeatModeler | `repeatmodeler` | De novo TE library construction |
| RepeatMasker | `repeatmasker` | TE annotation and masking |
| cd-hit | `cd-hit` | Library deduplication |
| bedtools | `bedtools` | TE-gene intersection |
| samtools | `samtools` | FASTA operations |
| BLAST | `blast` | Sequence homology |
| Python/pandas | `pandas` | Data manipulation |
| Python/matplotlib | `matplotlib` | Landscape plots |
| Perl | `perl` | parseRM scripts |
| multiqc | `multiqc` | QC aggregation |

**Data Source:** Drosophila melanogaster genome (dm6, ~140 MB)
- NCBI: GCF_000001215.4
- Alternative: Arabidopsis thaliana TAIR10 (~120 MB)
- Dfam library: https://www.dfam.org/releases/current/families/

**Data Size:** ~200 MB (genome FASTA + Dfam/RepBase library)
**Estimated Runtime:** ~3h on 8 CPUs (RepeatModeler is compute-intensive)

**Key Domain-Specific Traps:**
1. **RepeatModeler runtime**: De novo TE library building takes 1-3h even for small genomes; agents may timeout or skip this crucial step
2. **Library merging**: Must merge de novo and known libraries, then deduplicate with cd-hit; using only known library misses species-specific TEs
3. **Kimura distance**: TE age estimation uses Kimura 2-parameter distance from consensus; requires parsing RepeatMasker `.align` output (not `.out`)
4. **TE classification**: RepeatMasker uses hierarchical classification (DNA/TcMar-Tc1, SINE/Alu, LINE/L1, LTR/Gypsy); must correctly parse the class/family hierarchy
5. **Simple repeats vs TEs**: RepeatMasker annotates both; simple repeats (microsatellites) must be filtered out for TE analysis

**Why the DAG is genuinely complex:**
- RepeatModeler (slow, de novo) runs in parallel with RepeatMasker (fast, known library)
- Library merge convergence is mandatory before comprehensive RepeatMasker run
- Post-annotation triple branch (landscape, gene proximity, age distribution)
- TE class + gene disruption + solo-LTR analysis triple final branch
- RepeatModeler runtime creates asymmetric diamond (one branch 10x slower)
- Cross-dependency: genome stats from initial analysis feed into normalization and final report

---

### TASK 24: Multi-omics Integration (RNA + ATAC)

**Name:** `multiomics-rna-atac`
**Description:** Integrate matched RNA-seq and ATAC-seq data from the same samples to identify regulatory circuits linking chromatin accessibility to gene expression.

**DAG Diagram (depth=10, convergence=5, tools=12):**
```
 RNA_R1.fq.gz  RNA_R2.fq.gz    ATAC_R1.fq.gz  ATAC_R2.fq.gz
     │              │               │               │
 [fastp] ────── [fastp]         [fastp] ────── [fastp]         Level 1
     │              │               │               │
     └──────┬───────┘               └───────┬───────┘
            │                               │
     [STAR align]                    [bowtie2 align]            Level 2
            │                               │
     [samtools sort]                 [samtools sort]             Level 3
            │                               │
     ┌──────┼──────┐                 ┌──────┼──────┐
     │      │      │                 │      │      │
 [Salmon  [StringTie [featureCounts [picard [MACS2  [deeptools  Level 4
  quant]   assemble]  (gene counts)] dedup]  peaks]  bamCov]
     │      │      │                 │      │      │
     └──┬───┘      │                 └──┬───┘      │
        │          │                    │          │
 [CONVERGENCE 1]   │             [CONVERGENCE 2]   │            Level 5
 (expression:       │             (accessibility:   │
  Salmon+StringTie) │              dedup+peaks)     │
        │          │                    │          │
        └──────────┼────────────────────┘          │
                   │                               │
           [CONVERGENCE 3] ◄───────────────────────┘            Level 6
           (RNA expression + ATAC peaks + coverage)
                   │
           ┌───────┼───────────────┐
           │       │               │
     [bedtools  [HOMER           [python                        Level 7
      closest    findMotifs       correlation
      (peak→     Enrichment       (peak score
       gene)]    (peaks)]         vs expression)]
           │       │               │
           └───────┼───────────────┘
                   │
           [CONVERGENCE 4] ◄───────────────────                Level 8
           (peak-gene links + motifs + correlation)
                   │
           ┌───────┼───────────┐
           │       │           │
     [DESeq2    [python      [python                            Level 9
      diff       TF-target    regulatory
      expression  network]     circuit
      + diff                   scoring]
      accessibility]
           │       │           │
           └───────┼───────────┘
                   │
           [CONVERGENCE 5] ◄── QC from both tracks              Level 10
           [python integrative report]
```

**Longest path:** ATAC fastp -> bowtie2 -> sort -> picard dedup -> MACS2 peaks -> convergence2 -> convergence3 -> bedtools peak-gene -> convergence4 -> regulatory circuit -> convergence5 -> report (depth 10)

**Tool List (12 tools):**
| Tool | Conda package | Purpose |
|------|--------------|---------|
| fastp | `fastp` | Read QC (both modalities) |
| STAR | `star` | RNA-seq alignment |
| bowtie2 | `bowtie2` | ATAC-seq alignment |
| samtools | `samtools` | BAM operations |
| Salmon | `salmon` | Transcript quantification |
| StringTie | `stringtie` | Transcript assembly |
| MACS2 | `macs2` | ATAC peak calling |
| picard | `picard` | Duplicate removal |
| deeptools | `deeptools` | Signal tracks |
| HOMER | `homer` | Motif analysis |
| DESeq2 | `bioconductor-deseq2` | Differential analysis |
| bedtools | `bedtools` | Peak-gene assignment |

**Data Source:** ENCODE matched RNA-seq + ATAC-seq (K562 or GM12878)
- RNA-seq: ENCSR000AEL (K562)
- ATAC-seq: ENCSR868FGM (K562)
- Subset to chr22 for both modalities

**Data Size:** ~600 MB (RNA FASTQ + ATAC FASTQ + genome + annotation)
**Estimated Runtime:** ~3h on 8 CPUs

**Key Domain-Specific Traps:**
1. **ATAC-seq Tn5 shift**: Must apply +4/-5 bp shift to ATAC-seq reads for Tn5 insertion site correction; without this, peak summits are offset by 4-5 bp
2. **Peak-to-gene assignment**: Distance-based (nearest TSS) is simplistic; should use window approach (e.g., 100 kb upstream, 10 kb downstream) and consider TAD boundaries
3. **Normalization across modalities**: RNA counts and ATAC peak scores have different distributions; must normalize each separately before correlation
4. **ATAC-seq mitochondrial reads**: Must remove chrM reads (typically 30-60% of ATAC reads); forgetting this dilutes true signal
5. **Paired analysis**: Must ensure RNA and ATAC data are from same biological samples; mismatched samples give random correlation
6. **MACS2 parameters for ATAC**: Must use `--nomodel --shift -100 --extsize 200` or `--nolambda` for ATAC-seq; ChIP-seq defaults are wrong

**Why the DAG is genuinely complex:**
- 5 convergence points (maximum complexity)
- Dual-modality input with completely independent preprocessing tracks
- RNA track has internal Salmon/StringTie convergence
- ATAC track has internal dedup/peaks convergence
- Cross-modality convergence is the key integration point
- Three parallel integration analyses (peak-gene, motifs, correlation)
- Differential expression + differential accessibility + regulatory circuits triple final branch
- QC from both modalities feeds into final report

---

### TASK 25: Methylation Array Analysis (Illumina EPIC)

**Name:** `methylation-array-epic`
**Description:** Analyze Illumina EPIC methylation array data for differential methylation analysis, including normalization, batch correction, DMR detection, and pathway enrichment.

**DAG Diagram (depth=10, convergence=4, tools=10):**
```
 idat_files/ (Red + Green per sample)    sample_sheet.csv
       │                                      │
 [sesame/minfi readIDat] ◄────────────────────┘    Level 1
       │
 [sesame pOOBAH detection p-values]                  Level 2
       │
 ┌─────┼──────────────┐
 │     │              │
[sesame [sesame      [python                         Level 3
 QC     noob          sample QC
 (failed normalize]   (sex check,
 probes)]              genotype)]
 │     │              │
 └─────┼──────────────┘
       │
 [CONVERGENCE 1] ◄──────────────────                Level 4
 (QC'd + normalized beta values)
       │
 ┌─────┼──────────────────┐
 │     │                  │
[limma [python            [sesame                    Level 5
 combat/ PCA               age prediction
 batch   (batch            (epigenetic
 correct] check)]          clock)]
 │     │                  │
 └─────┼──────────────────┘
       │
 [CONVERGENCE 2] ◄──────────────────                Level 6
 (batch-corrected + PCA + age estimates)
       │
 [limma differential methylation]                    Level 7
 (CpG-level DMP)
       │
 ┌─────┼──────────────────────┐
 │     │                      │
[DMRcate [missMethyl        [python                  Level 8
 (DMR     GO/KEGG            volcano +
  regions)] enrichment]       manhattan]
 │     │                      │
 └─────┼──────────────────────┘
       │
 [CONVERGENCE 3] ◄──────────────────                Level 9
 (DMRs + pathways + plots)
       │
 [CONVERGENCE 4] ◄── QC + sample info               Level 10
 [python report]
```

**Longest path:** readIDat -> pOOBAH -> noob normalize -> convergence1 -> limma batch -> convergence2 -> limma DMP -> DMRcate -> convergence3 -> convergence4 -> report (depth 10)

**Tool List (10 tools):**
| Tool | Conda package | Purpose |
|------|--------------|---------|
| SeSAMe | `bioconductor-sesame` | IDAT reading, normalization, QC |
| minfi | `bioconductor-minfi` | Alternative IDAT processing |
| limma | `bioconductor-limma` | Batch correction, differential methylation |
| DMRcate | `bioconductor-dmrcate` | Differentially methylated regions |
| missMethyl | `bioconductor-missmethyl` | GO/KEGG enrichment |
| R/ggplot2 | `r-ggplot2` | Visualization |
| Python/pandas | `pandas` | Data manipulation |
| Python/matplotlib | `matplotlib` | Plots |
| samtools | `samtools` | General operations |
| multiqc | `multiqc` | QC aggregation |

**Data Source:** GEO methylation array datasets
- GSE168160 (EPIC array, matched cancer/normal, ~500 MB)
- Alternative: GSE149609 (aging study)
- FlowSorted.Blood.EPIC reference: Bioconductor package
- NOTE: IDAT files can be obtained from GEO via supplementary files

**Data Size:** ~400 MB (IDAT files for 8-12 samples + annotation + sample sheet)
**Estimated Runtime:** ~1.5h on 8 CPUs

**Key Domain-Specific Traps:**
1. **Probe type bias**: EPIC arrays have Type I and Type II probes with different dynamic ranges; must use noob or SWAN normalization to correct, not simple quantile normalization
2. **Detection p-value filtering**: Must remove probes with pOOBAH/detection p > 0.05; including failed probes adds noise that looks like differential methylation
3. **Cross-reactive probes**: ~5% of EPIC probes cross-hybridize; must filter using published lists (Chen et al. 2013, Pidsley et al. 2016)
4. **Sex chromosome handling**: Must either remove sex chromosome probes or account for sex differences; mixing sexes without correction creates spurious DMPs
5. **Batch effects**: Array data is notoriously batch-sensitive; ComBat correction requires careful design to not remove biological signal
6. **Enrichment bias**: Standard GO enrichment is biased for methylation arrays because gene length correlates with number of CpG probes; missMethyl corrects for this

**Why the DAG is genuinely complex:**
- Triple QC branch (failed probes, normalization, sample QC) converges before batch correction
- Batch correction + PCA + epigenetic clock triple branch (second diamond)
- DMP -> DMR + pathway + visualization triple branch (third diamond)
- Cross-branch: sample QC (sex, genotype) from early stages feeds into statistical model design
- Epigenetic clock is a unique cross-cutting analysis that requires normalized but NOT batch-corrected data
- Probe annotation (manifest) is a hidden dependency at multiple stages

---

## Summary Table

| # | Task ID | Domain | Depth | Conv | Tools | Data Size | Runtime | Data Source |
|---|---------|--------|-------|------|-------|-----------|---------|-------------|
| 1 | `hic-3d-conformation` | 3D Genomics | 10 | 4 | 10 | 400 MB | 2.5h | 4DN/Rao 2014 |
| 2 | `circrna-detection` | Non-coding RNA | 9 | 4 | 11 | 400 MB | 2h | nf-core/circrna |
| 3 | `cnv-detection-wes` | Cancer Genomics | 10 | 4 | 10 | 500 MB | 3h | GIAB/SEQC2 |
| 4 | `immune-repertoire` | Immunology | 12 | 4 | 10 | 300 MB | 1.5h | nf-core/airrflow |
| 5 | `gwas-association` | Pop Genetics | 10 | 4 | 10 | 200 MB | 1.5h | 1000G Phase 3 |
| 6 | `germline-wes-gatk` | Clinical Genomics | 12 | 4 | 12 | 600 MB | 3h | GIAB NA12878 |
| 7 | `structural-variant-multi` | Clinical Genomics | 10 | 4 | 12 | 700 MB | 3.5h | GIAB HG002 |
| 8 | `clinical-wgs-interpretation` | Clinical WGS | 12 | 5 | 14 | 800 MB | 3.5h | GIAB HG001 |
| 9 | `rna-editing-detection` | Transcriptomics | 10 | 4 | 10 | 500 MB | 2.5h | ENCODE K562 |
| 10 | `nascent-transcription` | Transcription | 9 | 4 | 10 | 300 MB | 1.5h | ENCODE PRO-seq |
| 11 | `longread-rna-isoform` | Long-read RNA | 10 | 4 | 11 | 400 MB | 2h | SG-NEx |
| 12 | `mag-recovery` | Metagenomics | 12 | 5 | 14 | 600 MB | 3h | CAMI2/nf-core/mag |
| 13 | `dia-proteomics` | Proteomics | 10 | 4 | 9 | 500 MB | 2.5h | PRIDE PXD014414 |
| 14 | `haplotype-phasing` | Pop Genetics | 10 | 4 | 10 | 300 MB | 2h | 1000G Phase 3 |
| 15 | `metatranscriptomics` | Metagenomics | 10 | 4 | 11 | 400 MB | 2.5h | HMP2 |
| 16 | `dda-lfq-proteomics` | Proteomics | 10 | 4 | 10 | 400 MB | 2.5h | PRIDE PXD000001 |
| 17 | `pharmacogenomics` | Clinical PGx | 9 | 4 | 10 | 500 MB | 2h | GIAB HG002 |
| 18 | `scatac-seq` | Single-cell | 10 | 4 | 11 | 300 MB | 2h | 10x PBMC 5k |
| 19 | `genome-scaffolding` | Assembly | 10 | 4 | 10 | 300 MB | 1.5h | E. coli K12 |
| 20 | `edna-metabarcoding` | Ecology | 10 | 4 | 10 | 200 MB | 1h | Zenodo eDNA |
| 21 | `viral-phylodynamics` | Phylogenetics | 10 | 4 | 10 | 50 MB | 1.5h | Nextstrain |
| 22 | `msi-detection` | Cancer Genomics | 10 | 4 | 10 | 500 MB | 2h | TCGA/SEQC2 |
| 23 | `repeat-element-annotation` | Genomics | 10 | 4 | 10 | 200 MB | 3h | Drosophila dm6 |
| 24 | `multiomics-rna-atac` | Multi-omics | 10 | 5 | 12 | 600 MB | 3h | ENCODE K562 |
| 25 | `methylation-array-epic` | Epigenomics | 10 | 4 | 10 | 400 MB | 1.5h | GEO EPIC |

**Complexity statistics:**
- Average depth: 10.2 (range: 9-12)
- Average convergence points: 4.1 (range: 4-5)
- Average tools: 10.6 (range: 9-14)
- All tools conda-installable (verified via micromamba search bioconda+conda-forge)
- All data < 1 GB, all runtimes < 4h on 8 CPUs
- 125 domain-specific traps total (5 per task)

**Overlap analysis with existing 44 tasks:**
- `hicar-chromatin` exists (HiCAR multi-omic) -- Task 1 (`hic-3d-conformation`) covers STANDARD Hi-C with completely different tools (cooler/cooltools vs pairtools+MACS2). No overlap.
- `sv-detection` exists (bacterial MRSA SVs, depth 6, Delly only) -- Task 7 (`structural-variant-multi`) covers HUMAN WGS with 4 callers + SURVIVOR merge. Different scale and complexity.
- `radseq-popgen` exists -- replaced Task 16 with `dda-lfq-proteomics` (dual search engine label-free quantitative proteomics).
- `somatic-germline-dual` exists -- Tasks 3/6/8 cover different aspects (CNV detection, germline-only WES, clinical WGS interpretation).
- All 25 tasks cover genuinely distinct domains from the existing benchmark.
