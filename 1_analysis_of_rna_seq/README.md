# Seriola dumerili

Sex:

    Ma: Male
    Fe: Female
    Mi: Mix

Tissue:

    Br: Brain
    Ey: Eye
    Gi: Gill
    Go: Gonad
    He: Heart
    Ki: Kidney
    Li: Liver
    Mu: Muscle
    In: Intestines
    Pi: Pituitary
    Sp: Spleen
    St: Stomach

To run the workflow:

    for fn in 1_SnakePrepare.smk 2_SnakeMapping.smk 3_SnakeExpression.smk 4_SnakeDenovoMapping.smk 5_SnakeTracks.smk 6_SnakeAssembly.smk; do
        snakemake --scheduler greedy -j --ri -k --rerun-triggers mtime --use-conda --conda-not-block-search-path-envvars -s $fn
    done
