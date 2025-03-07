# motif-mark

## A Tool for Motif Annotation

### Tool Takes In:

- FASTA file
- Text file containg motifs to be annotated

### To Run:

1. Create conda environment from provided yaml file
     `conda env create -f env/my_pycairo.yml`
2. Activate the conda environment
     `conda activate my_pycairo`
3. Run the script
     `python motif-mark-oop.py -f <FASTA file> -m <Motif file>`  