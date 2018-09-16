from cyvcf2 import VCF
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
import fire

class StayinAlive(object):

    def run_coxph(self, vcf_path, phenotype_path, output_path, id_column,
                  covar_columns, event_column, time_column,
                  chrom=None, start=None, end=None):

        def parse_phenotypes(path, id_column, covar_columns,
                             event_column, time_column):
            df = pd.read_csv(path)
            df = df.set_index(id_column)
            columns_to_keep = covar_columns + [event_column, time_column]
            df = df[columns_to_keep]
            return df

        def gt_types_convert(gt_types):
            gt_types = gt_types.astype(float)
            # gt_types is array of 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT
            gt_types[gt_types == 2] = np.nan  # missing
            gt_types[gt_types == 3] = 2  # hom_alt
            return gt_types

        df_pheno = parse_phenotypes(phenotype_path,
                                    id_column, covar_columns,
                                    event_column, time_column)

        cph = CoxPHFitter()
        vcf = VCF(vcf_path)
        sample_list = vcf.samples
        # TODO: handle partial information (e.g. chr only)
        region = ''
        if (chrom) and (start) and (end):
            region = f'{chrom}:{start}-{end}'
        f = open(output_path, "a")
        hdr = 'rsid,chr,start,end,ref,alt,maf,hr,lower_ci,upper_ci,se,z,p'
        f.write(hdr + '\n')

        for v in vcf(region):
            ref = v.REF
            alt = ''.join(v.ALT)  # assuming biallelics, which is a sin
            genotypes_numeric = gt_types_convert(v.gt_types)

            df_geno = pd.DataFrame({'id': sample_list,
                                   'genotype': genotypes_numeric})
            df_geno = df_geno.set_index('id')
            df = pd.merge(df_geno, df_pheno, left_index=True, right_index=True)
            df = df.dropna()
            try:
                fit = cph.fit(df, event_col=event_column,
                              duration_col=time_column)
                genotype_fit = fit.summary.loc['genotype']
                f.write(','.join(str(x) for x in
                        [
                        v.ID,
                        v.CHROM,
                        v.start,
                        v.end,
                        ref,
                        alt,
                        v.aaf,
                        genotype_fit["exp(coef)"],
                        genotype_fit["lower 0.95"],
                        genotype_fit["upper 0.95"],
                        genotype_fit["se(coef)"],
                        genotype_fit["z"],
                        genotype_fit["p"]
                        ]
                ) + '\n')
            except ValueError as e:
                print(f'{v.start} {v.end} {v.ID} failed, {str(e)}')
        f.close()


if __name__ == '__main__':
    fire.Fire(StayinAlive)
