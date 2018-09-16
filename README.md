# StayingAlive
GWAS via Cox proportional hazards model

# How to run

```
python stayinalive.py run_coxph \
--vcf-path genotypes.gz \
--phenotype-path phenotypes.csv \
--output-path out.csv \
--id-column iid \
--covar-columns [sex,pc1,pc2] \
--event-column event_status \
--time-column event_time \
--chrom 6 \ # optional
--start 31000000 \  # optional
--end 32000000 # optional
```
