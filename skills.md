# Population Genetics Tool Skills

**Scope:** population genetics / ancient DNA analysis.  
**Input assumption:** raw genotype data are mainly `EIGENSTRAT`, `PACKEDANCESTRYMAP`, or `PLINK bed/bim/fam` format.

## 1. Genotype Data Handling & QC

- Tool: `PLINK 1.9/2.0`, `EIGENSOFT/CONVERTF`
- Capability:
  - Filter SNPs and individuals by missingness, MAF, chromosome, LD, and sample/population list.
  - Convert between PLINK and EIGENSTRAT/PACKEDANCESTRYMAP formats.
  - Prepare merged ancient/modern genotype panels for PCA, ADMIXTURE, ADMIXTOOLS, TreeMix, and qpAdm workflows.

## 2. PCA & Population Structure Screening

- Tool: `smartpca`, `PLINK --pca`, `ADMIXTOOLS2`
- Capability:
  - Run PCA on selected reference populations.
  - Project ancient or low-coverage individuals onto modern PCA space.
  - Detect outliers, batch effects, broad clines, and population structure.

## 3. Global Ancestry Clustering

- Tool: `ADMIXTURE`, `STRUCTURE`, `fastSTRUCTURE`, `NGSadmix`
- Capability:
  - Estimate global ancestry components under unsupervised or supervised models.
  - Run multiple K values, cross-validation, and replicate clustering.
  - Export Q matrices for ancestry barplots.

## 4. f-statistics, qpWave, qpAdm & qpGraph

- Tool: `ADMIXTOOLS 1`, `ADMIXTOOLS2`, `admixr`
- Capability:
  - Compute `f2`, `f3`, `f4`, D-statistics, outgroup-f3, and f4-ratio tests.
  - Run `qpWave`, `qpAdm`, and `qpGraph` using EIGENSTRAT or PLINK-derived panels.
  - Test clade structure, ancestry proportions, rotating outgroups, and admixture graph fit.

## 5. Population Trees & Migration Graphs

- Tool: `TreeMix`, `qpGraph`, `MiqoGraph`
- Capability:
  - Infer population split graphs with migration edges.
  - Compare alternative graph topologies.
  - Generate graph hypotheses for downstream qpAdm/qpGraph validation.

## 6. Admixture Dating

- Tool: `DATES`, `ALDER`
- Capability:
  - Estimate admixture dates from ancestry-LD decay.
  - Test target/source combinations.
  - Report fitted curves and jackknife standard errors.

## 7. Haplotype Sharing & Fine-scale Structure

- Tool: `ChromoPainter`, `fineSTRUCTURE`, `GLOBETROTTER`, `MOSAIC`
- Capability:
  - Paint recipient genomes as mosaics of donor haplotypes.
  - Build coancestry matrices and fine-scale population clusters.
  - Infer recent admixture sources and dates from haplotype sharing.

**Note:** These tools usually require phased haplotype input, so PLINK/EIGENSTRAT data often need conversion or reprocessing.

## 8. Local Ancestry & Archaic Ancestry Tracts

- Tool: `RFMix`, `admixfrog`, `FLARE`
- Capability:
  - Infer local ancestry along chromosomes.
  - Detect source-specific, introgressed, or archaic ancestry fragments.
  - Handle low-coverage or contaminated ancient DNA when using specialized tools such as `admixfrog`.

## 9. Kinship, ROH & IBD

- Tool: `PLINK`, `KING`, `READ`, `lcMLkin`, `hapROH`, `ancIBD`
- Capability:
  - Detect close relatives, duplicates, and related individuals.
  - Estimate runs of homozygosity for inbreeding/endogamy analysis.
  - Detect long IBD segments between ancient individuals.

## 10. Demographic Modeling & Simulation

- Tool: `fastsimcoal2`, `dadi`, `moments`, `momi2`, `msprime/tskit`
- Capability:
  - Model divergence, migration, bottlenecks, growth, and admixture.
  - Simulate genotype data under explicit demographic scenarios.
  - Compare observed summary statistics with simulated or expected statistics.

**Note:** `PSMC`, `MSMC`, and `SMC++` are important SMC-based methods, but they usually require WGS/VCF-level data rather than sparse SNP-array-style EIGENSTRAT/PLINK panels.

## 11. Selection & Polygenic Adaptation

- Tool: `selscan`, `PLINK`, `BayEnv/BayPass`, `custom PBS/FST scripts`
- Capability:
  - Compute EHH, iHS, XP-EHH, nSL, FST, PBS, and allele-frequency differentiation.
  - Scan candidate adaptive loci or population-specific drift.
  - Summarize allele-frequency or polygenic-score changes through time.

## 12. Visualization & Reporting

- Tool: `R/ggplot2`, `patchwork`, `pheatmap`, `ComplexHeatmap`, `pophelper`, `pong`, `CLUMPAK`, `ggtree`, `ape`, `igraph`, `Graphviz`, `Plotly`, `Shiny`, `Streamlit`, `QGIS`, `sf`, `geopandas`, `leaflet`
- Capability:
  - Draw PCA plots, ADMIXTURE/STRUCTURE barplots, f4 heatmaps, qpAdm result grids, TreeMix/qpGraph diagrams, ROH/IBD tracks, and geographic ancestry maps.
  - Create interactive dashboards for filtering populations, models, SNP cutoffs, Z-scores, p-values, and ancestry components.
  - Export publication-ready PDF/SVG/PNG figures.

## References / Basis

- User-provided reading list: `Reading list for population genetic analysis methods` by Choongwon Jeong, 20 August 2020.
- Additional commonly used tools added for modern workflows: `ADMIXTOOLS2`, `DATES`, `admixfrog`, `FLARE`, `hapROH`, `ancIBD`, `msprime/tskit`, and broad visualization tools.
