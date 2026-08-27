# Thesis figure scripts

The R scripts that produced the figures in the dissertation. One script per
figure, except for the last, which produces two.

Each script is self contained: it reads the result tables produced earlier in
the analysis, draws the figure, and writes a PNG. None of them recompute any
statistics. The differential testing, annotation and simulation work all happens
upstream, in the reditR package and the analysis scripts, and these scripts only
present the results of it.

## What each script produces

| Script | Figure | Dataset |
|---|---|---|
| `fig01_diffedit_summary_endothelial2.R` | Differential editing testing summary | Endothelial, ADAR1 and ADAR2 knockdown |
| `fig02_top_genes_endothelial2.R` | Genes with the most differentially edited sites | Endothelial |
| `fig03_differential_editing_timecourse_mouse.R` | Change in editing across the dehydration time course | Mouse kidney endothelium |
| `fig04_editing_ratio_compartment_timepoint_mouse.R` | Editing level by compartment and timepoint | Mouse kidney endothelium |
| `fig05_per_gene_delta_endothelial2.R` | Per gene change in editing ratio | Endothelial |
| `fig06_parametric_sweep_reditr.R` | Test calibration and power against between sample variance | Simulated |
| `fig07_variant_consequences_endothelial2.R` | Where editing sites fall within genes | Endothelial |
| `fig08_09_qc_simulation_checks.R` | Two figures: simulator parameter recovery, and dispersion of simulated against real data | Simulated and all observed |

Figures 3 and 4 are companions. Figure 4 shows the editing level itself, so the
three compartments can be compared with each other. Figure 3 shows the change
from each compartment's own control, so the compartments can be compared with
where they started.

Figures 8 and 9 are the check on Figure 6. Figure 6 compares the three
statistical tests on simulated data, and that comparison only means something if
the simulated data is realistic, which is what those two figures test.

## Requirements

R with `data.table`, `ggplot2` and `patchwork`. Figure 1 additionally needs
`ggforce` for the circles in its middle panel, and figures 8 and 9 need the
`reditR` package itself, because they generate simulated data using it.

## Running them

Each script sets its own input and output directories near the top, as absolute
paths on the machine the analysis was run on. Those paths need changing before
the scripts will run anywhere else.

Run any of them with:

```
Rscript fig01_diffedit_summary_endothelial2.R
```

The PNG is written to the output directory set in the script, and a short
summary of the numbers behind the figure is printed to the console. That printed
output is worth reading, since several of the figures deliberately keep numbers
off the image and report them separately.

## A note on the comments

The comments explain why each figure is drawn the way it is, not just what the
code does. Several of the figures involved choices that change how a reader
interprets them: which denominator a count is shown against, whether an axis is
logarithmic, whether a gap in the data is drawn as a gap or quietly bridged.
Those choices are recorded in the scripts alongside the code that implements
them, so the reasoning is not lost.
