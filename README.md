# motif-mark

## A Tool for Motif Annotation

### Tool Takes In

- FASTA file
- Text file containg motifs to be annotated

### Current Functionality

1. Takes in input files 
2. Finds exons, introns, and motifs within the provided FASTA sequences
3. Creates an image containing a visual representation of the annotated FASTA sequences  

### To Run

1. Create conda environment from provided yaml file

     `conda env create -f env/my_pycairo.yml`
   
2. Activate the conda environment

     `conda activate my_pycairo`
   
3. Run the script

     `python motif-mark-oop.py -f <FASTA file> -m <Motif file>`  
