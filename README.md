# Written Programs  
A repository for useful programs used in in bioinformatic analyses. Brief descriptions will be below

`concat.sh` - This is used for concatenating completed orthogorup sequences in a single copy orthology analysis when using >x% species in the analysis. 
```
./combine.pl \
expected.txt \ ## list of species identifiers used in analysis
concat_sco_101.fasta \ ## output file name
sco_comp_g90/SCOgreat101.txt.OrthoGroup* ## completed orthogroup sequences in a directory
```
`complete_ogs.sh` - This takes in the full list of species in a single copy ortholog analysis, and completes the entries if you are using SCOs in >x % species
```
./add_missing_species.sh \
-s /path/to/species_list.txt \ ## list of all species being used
-i /path/to/orthogroups \ ## path to the extracted orthogroups 
-o /path/to/output_dir ## path to directory to write the completed orthogroups to 
```

`busco_stats.sh` - Script to loop through all files from a BUSCO analysis and generate a long form data frame that can be used in R visualisation.  

`summary_stats.sh` - Script to loop through raw reads, trimmed reads, and contaminant removed reads to get number of reads. This is slow and I can definetly improve it to make it quicker. 
