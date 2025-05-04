home_dir = "/lustre/pulkit.h/snakemake_local/non_coding_rna_snakemake"
input_dir= "new_data_output"
output_dir="test_output"
logs_dir="test_logs"
benchmark_dir = "test_benchmark"
rule all:
    input:
        make_input = f'{home_dir}/input/{input_dir}/',
        cscape_out = expand(f'{home_dir}/{output_dir}/cscape/output_cscape_grch38.txt',home_dir= home_dir),
        candra_out= expand(f"{home_dir}/{output_dir}/candra/output_candra.txt",home_dir= home_dir),
        funseq_out = expand(f"{home_dir}/{output_dir}/funseq/funseq2_output",home_dir= home_dir),
        annovar_out= expand(f"{home_dir}/{output_dir}/annovar/annovar_output-remove.hg19_multianno.vcf",home_dir= home_dir),
        cadd_out = expand(f"{home_dir}/{output_dir}/cadd",home_dir= home_dir),
        fathmm_out= expand(f"{home_dir}/{output_dir}/fathmm/fathmm_out.txt",home_dir= home_dir),
        sift_out= expand(f"{home_dir}/{output_dir}/sift/output_sift",home_dir= home_dir),

###make the input files for different tools
rule makeinput:
    input: 
        vcf = f'{home_dir}/input/new_data/Breast-AdenoCa_pacwg_hg19_10000.vcf',
    output:
        input_dir = directory(f"{home_dir}/input/{input_dir}"),
        cscape=f"{home_dir}/input/{input_dir}/cscape_input.txt",
        candra=f"{home_dir}/input/{input_dir}/candra_input.txt",
        funseq2=f"{home_dir}/input/{input_dir}/funseq2_input.txt",
        annovar=f"{home_dir}/input/{input_dir}/annovar_input.vcf",
        cadd=f"{home_dir}/input/{input_dir}/cadd_input.vcf",
        fathmm=f"{home_dir}/input/{input_dir}/fathmm_input.txt",
        sift=f"{home_dir}/input/{input_dir}/sift_input.vcf"
    shell:
        """
        mkdir -p {output.input_dir}
        ./fileconversion.sh {input} {output.input_dir}
        """
####require headers#####
rule spliceai:
    input:
        #vcf = f'{home_dir}/input/input.vcf',
        vcf = f'{home_dir}/input/{input_dir}/spliceai_input.vcf',
        ref = f'{home_dir}/input/reference/gatk_hg38_reference.fa',
        annot= f'{home_dir}/input/annotations/grch38.txt'
    output:
        vcf = f'{home_dir}/{output_dir}/spliceai/output_spliceai/variants.vcf'
    benchmark:
        f"{home_dir}/{benchmark_dir}/spliceai_benchmark.txt"
    log:
        f"{home_dir}/{logs_dir}/spliceai_logs.txt"
    shell:
        """
        spliceai -I {input.vcf} -O {output.vcf} -R {input.ref} -A {input.annot} 2> {log}
        """
#####only snp s and comma separated files with chrom,pos,ref,alt without any header####
rule cscape:
    input:
        #vcf = f'{home_dir}/cscape/cscape_input.txt'
        vcf = f'{home_dir}/input/{input_dir}/cscape_input.txt'
    output:
        csv = f'{home_dir}/{output_dir}/cscape/output_cscape_grch38.txt'
    benchmark:
        f"{home_dir}/{benchmark_dir}/cscape_benchmark.txt"
    log:
        f"{home_dir}/{logs_dir}/cscape_logs.txt"
    shell:
        """
        cd /lustre/pulkit.h/snakemake_local/non_coding_rna_snakemake/cscape
        pwd
        python3 /lustre/pulkit.h/snakemake_local/non_coding_rna_snakemake/cscape/cscape_query.py {input.vcf} -n /lustre/pulkit.h/snakemake_local/non_coding_rna_snakemake/cscape/cscape_noncoding.vcf.gz -o {output.csv} -c /lustre/pulkit.h/snakemake_local/non_coding_rna_snakemake/cscape/cscape_coding.vcf.gz 2> {log}
        """
#####only chrom pos ref alt, tab delimited files####
rule candra:
    input:
        vcf= f"{home_dir}/input/{input_dir}/candra_input.txt",
        script= f"{home_dir}/candra/CanDrA.v1.0/open_candra.pl",
    output:
        tsv= f"{home_dir}/{output_dir}/candra/output_candra.txt",
    benchmark:
        f"{home_dir}/{benchmark_dir}/candra_benchmark.txt"
    log:
        f"{home_dir}/{logs_dir}/candra_logs.txt"
    shell:
        """
        cd /lustre/pulkit.h/snakemake_local/non_coding_rna_snakemake/candra/CanDrA.v1.0
        pwd
        cat input.txt
        perl {input.script} BRCA {input.vcf} > {output.tsv} 2> {log}
        """


#funseq2 specific vcf needed just add ##fileformat=VCFv4.0 
#can give user annotation db , add multiple dbs in the directory and mention the path to the directory in the comamnd
#can use funseq2 dbs ex breast cancer db
#output folder path is needed in the command
rule funseq:
    input:
        #vcf = f"{home_dir}/funseq2-1.0/GRCH38_TEST.vcf",
        vcf = f"{home_dir}/input/{input_dir}/funseq2_input.txt",
        script = f"{home_dir}/funseq2-1.0/run.sh",
        db_brc= f"{home_dir}/funseq2-1.0/data/cancer_recurrence/Breast.recur",
        db_user= f"{home_dir}/funseq2-1.0/data/user_annotations"

    output: 
        dir = directory(f"{home_dir}/{output_dir}/funseq/funseq2_output")
    benchmark:
        f"{home_dir}/{benchmark_dir}/funseq_benchmark.txt"
    log:
        f"{home_dir}/{logs_dir}/funseq_logs.txt"
    shell:
        """
        cd /lustre/pulkit.h/snakemake_local/non_coding_rna_snakemake/funseq2-1.0
        {input.script} -f {input.vcf} -maf 1 -m 1 -len 20 -inf vcf -outf vcf -nc -o {output.dir} -s 1.5  -ua {input.db_user} -db {input.db_brc} 2> {log}
        """
####only acctept #chrom, id, pos ref, alt, qual, filter and tab delimited files####        
rule annovar:
    input:
        #vcf= "/lustre/pulkit.h/snakemake_local/annovar/vep_output_filtered_morfee.vcf",
        vcf = f"{home_dir}/input/{input_dir}/annovar_input.vcf",
        humandb= "/lustre/pulkit.h/snakemake_local/annovar/humandb/",
        db_gene= "/lustre/pulkit.h/snakemake_local/annovar/humandb/morfee.lof_metrics.by_gene.txt"
    output:
        out= f"{home_dir}/{output_dir}/annovar/annovar_output-remove.hg19_multianno.vcf",
    params:
        out= f"{home_dir}/{output_dir}/annovar/annovar_output"
    benchmark:
        f"{home_dir}/{benchmark_dir}/annovar_benchmark.txt"
    log:
        f"{home_dir}/{logs_dir}/annovar_logs.txt"
    shell:
        """
        cd /lustre/pulkit.h/snakemake_local/annovar
        ./table_annovar.pl {input.vcf} {input.humandb} -buildver hg19 -out {params.out}-remove -protocol refGene,gwasCatalog,avsnp150,clinvar_20190305,gnomad211_genome,dbnsfp35a -operation gx,r,f,f,f,f -nastring . -vcfinput -xreffile {input.db_gene} 2> {log}
        """
#### tab delimited files with vcf extension####
rule cadd:
    input:
        #vcf = f"{home_dir}/cadd/test/input.vcf",
        vcf = f"{home_dir}/input/{input_dir}/cadd_input.vcf",
    output:
        out = directory(f"{home_dir}/{output_dir}/cadd")    
    benchmark:
        f"{home_dir}/{benchmark_dir}/cadd_benchmark.txt"
    log:
        f"{home_dir}/{logs_dir}/cadd_logs.txt"
    params: 
        genome_ver= "GRCh37",
        cores = 20
    shell:
        """
        cd /lustre/pulkit.h/snakemake_local/non_coding_rna_snakemake/cadd
        ./CADD.sh -a -c {params.cores} -g {params.genome_ver} -o {output.out} {input.vcf} 2> {log}
        """

#### column heading CHROM,POS,REF,ALT. only accept 4 cols and commaseparated files####
rule fathmm:
    input:
        #vcf= f"{home_dir}/fathmm/test.txt",
        vcf = f"{home_dir}/input/{input_dir}/fathmm_input.txt",
        db = f"{home_dir}/fathmm/fathmm-MKL_Current.tab.gz",
        script = f"{home_dir}/fathmm/fathmm-MKL.py"
    output:
        out= f"{home_dir}/{output_dir}/fathmm/fathmm_out.txt"    
    benchmark:
        f"{home_dir}/{benchmark_dir}/fathmm_benchmark.txt"
    log:
        f"{home_dir}/{logs_dir}/fathmm_logs.txt"
    shell:
        """
        python {input.script} {input.vcf} {output.out} {input.db} 2> {log}
        """
#### 8 col min required; no extra lines at eof; # on heading line #####
rule sift:
    input:
        jar= f"{home_dir}/sift/SIFT4G_Annotator.jar",
        #vcf = f"{home_dir}/sift/GRCH38_TEST.vcf",
        vcf = f"{home_dir}/input/{input_dir}/sift_input.vcf",
        db = f"{home_dir}/sift/databse/GRCh37.74/",
    output:
        out= directory(f"{home_dir}/{output_dir}/sift/output_sift"),
    benchmark:
        f"{home_dir}/{benchmark_dir}/sift_benchmark.txt"
    log:
        f"{home_dir}/{logs_dir}/sift_logs.txt"
    shell:
        """
        java -jar {input.jar} -c -i {input.vcf} -d {input.db} -r {output.out} 2> {log}
        """