# AGENTS.md

## Purpose

This agent supports a population genetics / ancient DNA research workflow.

Detailed software capabilities are maintained separately in `skills.md`. Do not duplicate tool-by-tool descriptions here. This file defines the agent’s working principles inside the project repository.

## Primary Responsibilities

1. Search, verify, and summarize research papers.
2. Analyze genotype data in EIGENSTRAT, PACKEDANCESTRYMAP, or PLINK format.
3. Visualize population genetic analysis results.
4. Organize research notes and interpretations into `wiki/`.
5. Track progress, decisions, and reproducible changes with Git.

## Directory Structure

- `raw/`: original source materials and unmodified input files.
- `wiki/`: paper summaries, method notes, analysis notes, and interpretation logs.
- `output/`: reports, tables, figures, summaries, and action plans.
- `scripts/`: reusable scripts for analysis, plotting, conversion, and QC.
- `config/`: sample lists, population lists, parameters, plotting settings, and manifests.
- `logs/`: run logs, software versions, errors, and troubleshooting notes.
- `archive/`: deprecated or superseded outputs kept for traceability.

## Source Preservation

- Do not overwrite or modify original files in `raw/`.
- Record corrections to original sources in separate notes.
- Preserve metadata such as source, access date, DOI, URL, software version, command options, and reference panel.
- Keep original genotype panels, metadata, and sample annotations unchanged.
- Save derived files with explicit suffixes such as `.filtered`, `.autosome`, `.tv`, `.ldpruned`, `.projected`, or date tags.

## Paper Search & Summarization

- Prefer peer-reviewed papers, official documentation, method papers, and clearly labeled preprints.
- Verify title, authors, year, journal/preprint server, DOI, and URL before citing.
- Do not invent references, links, DOIs, quotations, or paper claims.
- Mark unverifiable claims as `needs verification`.
- Summaries should include:
  - research question,
  - data and sample set,
  - methods/tools,
  - main findings,
  - limitations,
  - relevance to the project,
  - useful figure/table/supplement locations.
- Distinguish authors’ claims, documentation claims, later interpretation, local interpretation, and speculation.

## Genotype Data Analysis

- Assume input is usually EIGENSTRAT, PACKEDANCESTRYMAP, or PLINK unless stated otherwise.
- Do not assume BAM, VCF, sequence-level, or phased haplotype data are available.
- Check whether a requested method is compatible with the available genotype format.
- Record reproducibility details:
  - input prefix and format,
  - sample/population lists,
  - SNP and chromosome filters,
  - missingness filters,
  - transition/transversion filtering,
  - pseudo-haploid treatment,
  - SNP/sample/population counts,
  - software version,
  - commands or scripts.
- For f-statistics, report statistic form, population order, estimate, SE, Z-score, and SNP count.
- For qpAdm/qpWave, report target, sources, left/right sets, fixed or rotating outgroups, p-value, proportions, SEs, SNP count, and model status.
- Do not treat PCA, ADMIXTURE, TreeMix, or qpAdm as final historical conclusions without checking assumptions and alternative models.

## Visualization

- Make figures reproducible from saved scripts and config files.
- Keep population order, labels, colors, shapes, and plot limits explicit.
- Prefer publication-ready `PDF`, `SVG`, or high-resolution `PNG`.
- Save plot input tables when possible.
- Include relevant context: statistic definition, comparison direction, Z-score threshold, SNP cutoff, error bar definition, and source/target/outgroup labels.
- Use filenames that include analysis type, SNP panel, date, and major filters.

## Wiki Writing

- Record key claims with supporting evidence locations.
- Mark unsupported or uncertain content as `needs verification`.
- Document conflicting sources instead of forcing one conclusion.
- Separate paper summaries, method notes, analysis notes, and interpretation notes.
- Analysis notes should include purpose, inputs, commands/scripts, parameters, outputs, interpretation, caveats, and next steps.
- Keep entries concise but traceable.

## Git & Progress Tracking

- Commit at meaningful checkpoints.
- Use commit messages that state what changed and why.
- Do not commit large raw datasets unless explicitly intended.
- Do not commit credentials, tokens, keys, or sensitive personal data.
- Keep scripts, configs, logs, and outputs clearly separated.
- Record major analysis decisions in `wiki/`, `logs/`, or commit messages.
- Check working tree status before major changes.
- Document changed assumptions or parameters when modifying scripts.

## Response Principles

- Be direct about uncertainty.
- Do not present unsupported claims as facts.
- Distinguish source-backed information, conventions, inference, and speculation.
- Ask for missing details when the request is ambiguous.
- Provide code body rather than executing code unless execution is explicitly requested.
- Consult relevant `wiki/` material before reports or action plans when available.
- Do not use sensitive project materials for external searches or transmission.
- If a method is unsuitable for the available input format, state the limitation and suggest a compatible alternative.
