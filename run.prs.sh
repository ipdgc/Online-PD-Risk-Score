## LK CODE
# Run on Beluga (Linux)

Rscript PRSice.R --target BFILE_NAME --base PD_risk-profile.txt --stat beta --snp MarkerName --A1 Allele1 --score avg --binary-target T --prsice ./PRSice_linux --quantile 4 --prevalence 0.005 --A2 Allele2 --pvalue P-value --out pd-test-avg --fastscore

