## Characterize HIV drug resistance mutations

Call mutations in pooled HIV samples, reporting drug resistance mutations across
a wide range of variation frequencies. The goal is to generate a set of
drug resistance mutations across multiple frequencies for identifying resistance
mutations in pooled sequencing inputs.

Uses deconvoluted BAM inputs from [bcbio](https://bcbio-nextgen.readthedocs.io/en/latest/)
or another pipeline with contaminating human reads removed. Run
[hivmmer](https://github.com/kantorlab/hivmmer) to call variants, producing AAVF
output. Add drug resistance mutations from the Stanford HIVdb and export to CSV
for uptake into R or Python tools for post analysis.

### Install

Create an isolated conda environment with hivmmer and required Python supporting
packages:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p env
env/bin/conda install -y -c kantorlab hivmmer

# https://github.com/winhiv/aavf-tools
env/bin/conda install -y -c conda-forge lxml
env/bin/pip install git+https://github.com/winhiv/PyAAVF.git@master#egg=PyAAVF-0
env/bin/pip install git+https://github.com/phac-nml/asipython.git@master#egg=asipython-0
env/bin/pip install git+https://github.com/winhiv/aavf-tools.git@master#egg=aavf-tools-0
```
On the Harvard Odyssey cluster, this is pre-installed at:
```
/net/hsphs10/srv/export/hsphs10/share_root/hsphfs1/share_root/chb/projects/novitsky_hiv/hivmmer
```

### Run

1. Retrieve inputs for hivmmer and drug resistance annotation, a POL reference and
   Standard HIVdb XML file:
```
mkdir -p input
cd input
wget https://raw.githubusercontent.com/kantorlab/hivmmer/master/test/pol.hmm
wget -O hivdb.xml https://raw.githubusercontent.com/winhiv/aavf-tools/master/aavf_resistance/data/sample.xml
```
2. Prepare a CSV file with project and sample names to process. Input sample
   BAMs get looked up from a bcbio output directory as `PROJECT/PROJECT/final/SAMPLE/SAMPLE*.bam`:
```
FC04187,GEN00128326
FC04187,GEN00128306
FC04229,GEN00132449
```
3. Run the analysis feeding the input sample file, base directory with projects
   and BAMs, and the POL reference HIVdb XML file:
```
export PATH=/net/hsphs10/srv/export/hsphs10/share_root/hsphfs1/share_root/chb/projects/novitsky_hiv/hivmmer/env/bin:$PATH
INDIR=/n/regal/novitsky_lab/vnovitsky
python hiv_var_resist.py input/samples.csv $INDIR input/pol.hmm input/hivdb.xml
```
Individual output AAVF files for each sample will be in 
`work/SAMPLE/SAMPLE-hivmmer.hmmsearch2.annotated.aavf` and a summary CSV file
with drug resistance calls from all samples is in `samples-calls.csv`.
