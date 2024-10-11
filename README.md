# snakepipes_ATAC-seq
## about author

> author: [马天宇 | Tianyu Ma]
>
> email: tianyu2697@gmail.com
>
## doc
`snakepipes_ATAC-seq` is a standard snakemake pipeline for ATAC-seq sequencing data
## env
```shell
conda env create -f conda_env.yml
```
## run
```shell
# run Jupyter notebook to abtain the config
# run this cmd
# or
# open notebook and run all cells
runipy step.01.GetFileName.ipynb
# dry run for test
snakemake -pr -j 10 -s step.02.Snakefile.smk.py -n
# run calculation
snakemake -pr -j 10 -s step.02.Snakefile.smk.py
```
## project structure