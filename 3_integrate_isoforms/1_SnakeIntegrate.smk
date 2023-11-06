#!/usr/bin/env snakemake
include: "0_SnakeCommon.smk"
outdir = "results/integrate"

rule all:
    input:
        outdir + "/ngs_tgs",
        outdir + "/ngs_tgs.bed",
        outdir + "/ngs_tgs.gtf",
    

rule integrate:
    input:
        ngs_gtf = GTFS["ngs"],
        tgs_gtf = GTFS["tgs"]
    output:
        ngs_bed = outdir + "/ngs.bed",
        tgs_bed = outdir + "/tgs.bed",
        txt = outdir + "/ngs_tgs.filelist",
        out = directory(outdir + "/ngs_tgs")
    log:
        outdir + "/ngs_tgs.log"
    params:
        prefix = outdir + "/ngs_tgs/ngs_tgs"
    shell:
        """
        set +u; source activate TAMA
        python ~/software/tama/tama_go/format_converter/tama_format_gtf_to_bed12_stringtie.py {input.ngs_gtf} {output.ngs_bed}
        python ~/software/tama/tama_go/format_converter/tama_format_gtf_to_bed12_stringtie.py {input.tgs_gtf} {output.tgs_bed}
        echo -e "{output.ngs_bed}\\tno_cap\\t2,1,2\\tNGS" >> {output.txt}
        echo -e "{output.tgs_bed}\\tno_cap\\t1,1,1\\tTGS" >> {output.txt}
        mkdir {output.out}
        python ~/software/tama/tama_merge.py -f {output.txt} -p {params.prefix} &> {log}
        """

rule post_integrate:
    input:
        beddir = rules.integrate.output.out
    output:
        bed = outdir + "/ngs_tgs.bed",
        bed_gz = outdir + "/ngs_tgs.sorted.bed.gz",
        gtf = outdir + "/ngs_tgs.gtf",
        gtf_gz = outdir + "/ngs_tgs.sorted.gtf.gz"
    shell:
        """
        cp {input.beddir}/ngs_tgs.bed {output.bed}
        sort -k1,1 -k2,2n {output.bed} | bgzip -c > {output.bed_gz}
        tabix -p bed {output.bed_gz}
        set +u; source activate TAMA
        python ~/software/tama/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py {output.bed} {output.gtf}
        sort -k1,1 -k4,4n {output.gtf} | bgzip -c > {output.gtf_gz}
        tabix -p gff {output.gtf_gz}
        """