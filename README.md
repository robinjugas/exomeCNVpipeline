# exomeCNVpipeline
exomeCNVpipeline

tri mody - config["analysis_mode"]:

paired - somatic_paired - pro vetsinu nastroju bude stejne jako reference_cohort -
with normals cohorts - testxcontrol_cohort 1xX
baseline - versus baselina created from test samples


podle modu - seznam calleru
plus respektovat ty vybrane


sample status 3 typy:
status = "tumor|normal|control"
tumor & normal -  pokud somatic_paired
normal & control - pokud reference_cohort



somatic_paired:
ADD donor to json config
HLAVNE U CNVkitu osetrit, jinde je to 1xX


exomeDpeth / cnMOPS / panelcnMOPS - pouze kohorta vs vzorek, neumi paired mode


?
jaka verze GTF slouzila k BEDu
https://cnvkit.readthedocs.io/en/stable/germline.html   CNVKIT GERMLINE!!!!!

With no control samples

Alternatively, you can create a “flat” reference of neutral copy number (i.e. log2 0.0) for each probe from the target and antitarget interval files. This still computes the GC content of each region if the reference genome is given.

cnvkit.py reference -o FlatReference.cnn -f ucsc.hg19.fa -t targets.bed -a antitargets.bed

Possible uses for a flat reference include:

    Extract copy number information from one or a small number of tumor samples when no suitable reference or set of normal samples is available. The copy number calls will not be quite as accurate, but large-scale CNVs should still be visible.
    Create a “dummy” reference to use as input to the batch command to process a set of normal samples. Then, create a “real” reference from the resulting *.targetcoverage.cnn and *.antitargetcoverage.cnn files, and re-run batch on a set of tumor samples using this updated reference.
    Evaluate whether a given paired or pooled reference is suitable for an analysis by repeating the CNVkit analysis with a flat reference and comparing the CNAs found with both the original and flat reference for the same samples.

<!-- https://cnvkit.readthedocs.io/en/stable/heterogeneity.html -->

conda activate /home/rj/4TB/CEITEC/snakemake_conda/00d90c2430ea88c11fcf578383e3ad6b_