# Written Programs  
A repository for useful programs used in in bioinformatic analyses. Brief descriptions for each below.  

## Python Scripts  

`combine_seqs.py` - This is used for concatenating completed orthogroup sequences in a single copy ortholog analysis. It will generate one file with all alignments for each species so to be used in `raxml`, `iqtree` etc etc.  
```
./combine_seqs.py --help
usage: combine_seqs.py [-h] --expected_file EXPECTED_FILE --output_file
                       OUTPUT_FILE --ortholog_files ORTHOLOG_FILES
                       [ORTHOLOG_FILES ...] [--line_length LINE_LENGTH]

Concatenate ortholog FASTA files across species and validate sequence lengths.

optional arguments:
  -h, --help            show this help message and exit
  --expected_file EXPECTED_FILE
                        File with expected species list (lines like >species)
  --output_file OUTPUT_FILE
                        Path to write concatenated FASTA
  --ortholog_files ORTHOLOG_FILES [ORTHOLOG_FILES ...]
                        Ortholog FASTA files to concatenate
  --line_length LINE_LENGTH
                        Wrap FASTA sequences at this line length (0 = no
                        wrapping)
```
`complete_missing_species.py` - This takes in the full list of species in a single copy ortholog analysis, and completes the entries if you are using SCOs in >x% of species. It will calculate the length of the alignment in each SCO, and then complete any missing species using `---` as the padding sequences.  
```
./complete_missing_species.py --help
usage: complete_missing_species.py [-h] --species_file SPECIES_FILE
                                   --input_dir INPUT_DIR
                                   [--input_ext INPUT_EXT] --output_dir
                                   OUTPUT_DIR [--output_ext OUTPUT_EXT]
                                   [--line_length LINE_LENGTH]
                                   [--missing_report MISSING_REPORT]

Fill missing species in orthogroup FASTA files with padded sequences.

optional arguments:
  -h, --help            show this help message and exit
  --species_file SPECIES_FILE
                        File with expected species list (lines like >species)
  --input_dir INPUT_DIR
                        Directory with cleaned orthogroup files
  --input_ext INPUT_EXT
                        Input file extension to select (default: .fa)
  --output_dir OUTPUT_DIR
                        Directory to write completed orthogroup files
  --output_ext OUTPUT_EXT
                        Output file extension (default: .fasta)
  --line_length LINE_LENGTH
                        Max characters per sequence line (default: 60)
  --missing_report MISSING_REPORT
                        Optional path to write a table of missing species per
                        orthogroup
```

`iqtree_mod_select.py` - this takes the outputs from generated genetrees from an `IQtree3` run and makes a tsv file with the model you want.  
```
./iqtree_mod_select.py --help
usage: test.py [-h] -i INPUT -o OUTPUT [-crit {AIC,BIC,AICc}]

Extract best-fit models from IQ-TREE logs.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Root directory containing SCO directories.
  -o OUTPUT, --output OUTPUT
                        Output TSV file.
  -crit {AIC,BIC,AICc}, --criterion {AIC,BIC,AICc}
                        Criterion for model selection (default: AIC)
```
`parse_antismash.py` - parses all the knownclusterblast results for a number of samples for all regions. You will need a basic python biological environment. It will note which programs need to be loaded so you may need to activate a `conda` environmnet with `python` biological programs.  
```
./parse_antismash.py --help
usage: parse_antismash.py [-h] --input_root INPUT_ROOT --output_root OUTPUT_ROOT [--cluster_subdir CLUSTER_SUBDIR] [--samples SAMPLES [SAMPLES ...]]

Parse ClusterBlast/antismash .txt files into long tables

options:
  -h, --help            show this help message and exit
  --input_root INPUT_ROOT
                        Path to parent directory that contains sample directories (e.g. /path/to/barcodes)
  --output_root OUTPUT_ROOT
                        Where to write per-sample results (a subdir per sample will be created)
  --cluster_subdir CLUSTER_SUBDIR
                        Name of subdirectory inside each sample that contains the .txt files (default: knownclusterblastdirectory)
  --samples SAMPLES [SAMPLES ...]
                        Optional: list of sample directory names, paths, or glob patterns (e.g. "barcode*", "sample01", "/abs/path/sampleX"). If omitted, all subdirectories of --input_root are processed.
```

`summarise_proteinortho.py` - Will take the results of a `proteinortho` run and generate `.txt` files with summaries for single copy orthologs (SCOs) between all species, and SCO in > xx% of species from 50% upwards in increments of 5. The outputs can then be used to then extract the SCOs from the proteomes for downstream analyses. 

```
usage: summarise_proteinortho.py [-h] [-e EXTENSION] [-o OUTPUT_DIR] tsv_file fasta_dir

Process ProteinOrtho output to extract SCOs at different percentage thresholds.

positional arguments:
  tsv_file              ProteinOrtho .tsv file
  fasta_dir             Directory with input FASTA files

options:
  -h, --help            show this help message and exit
  -e EXTENSION, --extension EXTENSION
                        FASTA file extension (default: .faa)
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Directory for output SCO files
```


## Bash Scripts  

`busco_stats.sh` - Script to loop through all files from a BUSCO analysis and generate a long form data frame that can be used in R visualisation.  

`summary_stats.sh` - Script to loop through raw reads, trimmed reads, and contaminant removed reads to get number of reads. This is slow and I can definetly improve it to make it quicker.  


## Perl (ewwwwwwww)

`extract_salmon_aln.pl` - getting alignment rates from all salmon directories in a directory.  Requires `bioperl` conda environment activated. 
```
perl \
extract_salmon_aln.pl \
/scratch/alpine/beyo2625/apal_genet/align \ ## directory with all salmon align directories in it
./mapping_summary.tsv ## output file to write the results to. 
```


