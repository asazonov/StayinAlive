# StayingAlive
GWAS via Cox proportional hazards model

#How to run

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
--start 32590771-500 \  # optional
--end 32590771+500 # optional
```
