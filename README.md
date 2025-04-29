# Written Programs  
A repository for useful programs used in in bioinformatic analyses. Brief descriptions will be below

`concat.sh` - This is used for concatenating completed orthogorup sequences in a single copy orthology analysis when using >x% species in the analysis. 
```
./combine.pl \
expected.txt \ ## list of species identifiers used in analysis
concat_sco_101.fasta \ ## output file name
sco_comp_g90/SCOgreat101.txt.OrthoGroup* ## completed orthogroup sequences in a directory
```
