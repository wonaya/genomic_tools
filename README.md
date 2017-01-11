# genomic_tools
various python scripts

1. fai to bed<br />
```python fai_to_bed.py [fai]```
2. chipcompare-0.1.py<br />
```python chipcompare-0.1.py --test_1 input1.broadPeak --test_2 input2.broadPeak --Ncoord nonmappable.bed --fai genome.fai --output output_prefix```
3. chiprnacompare-0.1.py<br />
```python chiprnacompare-0.1py --chip chip.broadPeak --rnaseq genes.fpkm_tracking --gff3 maize.gff3 --macs_score 50```
4. coverage_calculate.py<br />
```python coverage_calculate.py [bam] [fai] [-P] [read length]```
5. bedgz_to_summary_bisulfite.py<br />
```python bedgz_to_summary_bisulfite.py <geno name eg. B73 looks for file named B73_methratio.txt from $CWD>```
6. count_cg_sites.py<br />
```python count_cg_sites.py [fasta]``` gives bed file with count of CG, CHG, CHH sites of 100bp windows
