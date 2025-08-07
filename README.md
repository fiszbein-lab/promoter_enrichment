Code for determining the enriched gate of a promoter in 
Uriostegui-Arcos, et al. 2025.
- - - -
### Dependencies
- python 3.12.3
- biopython 1.84
- scipy 1.14.1
- statsmodels 0.14.3

### Description
`normalize.py` takes the read count for each promoter in each
replicate in each gate and divides that value over the depth of the replicate;
it then further divides the depth-normalized value over the percentage 
representation of the promoter in the enrichment library, recorded in
`out/enrichment-count.tsv`, which is derived from bulk amplicon-seq of the
integrated promotors and is thus used as a proxy for variation in integration
efficiencies irrespective of gate. 

`enrichment.py` requires the output of `normalize.py`. It takes the mean of the 
replicate values for each promoter in each gate, performs Chi-squared tests, 
and uses Bonferroni correction to adjust the *p*-values. If the adjusted
*p*-value is less than 0.05, the enriched gate is determined using a percentage 
threshhold greater than 50%, e.g., if the cut-off is 75%, the gate responsible
for >75% of the normalized count values is designated enriched gate, whereas if 
no gate meets that threshhold, no enriched gate is designated. 



