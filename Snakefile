shell.executable("/bin/bash")
#shell.prefix("source $HOME/.bashrc; ")

import os

configfile: "config.yaml"
TEMP=config['temp_dir']
MAGDIR=config['mag_dir']
IDS, = glob_wildcards(os.path.join(MAGDIR,"{id}.fa"))
KEGG_DB = config['kegg_db']
KEGG_DIR = config['kegg_dir']
DIAMOND_OUT = os.path.join(KEGG_DIR, 'diamond_db', os.path.split(KEGG_DB)[-1]).strip('.pep')
BLASTPREFIX = os.path.join(KEGG_DIR, 'blastdb', os.path.split(KEGG_DB)[-1]).strip('.pep')
BLASTDB = expand(BLASTPREFIX+'.00.{ext}', ext =["phr", "pin", "psq"])
BUSCO_GROUP = config['busco_group']
rule all:
	 input: expand("metaeuk/{magid}.faa", magid=IDS), expand("metaeuk/{magid}.faa.headersMap.tsv", magid=IDS),
            expand('genemark/{magid}/output/gmhmm.mod', magid=IDS)
            expand('maker2/{magid}/{magid}.all.maker.genemark.proteins.fasta',magid=IDS), expand('maker2/{magid}/{magid}.gff3', magid=IDS), expand('pfam/{magid}.pfam', magid=IDS), expand('kegg/diamond/{magid}.blastp.tab', magid=IDS), 
            expand(os.path.join('busco',BUSCO_GROUP,'{magid}','run_'+BUSCO_GROUP,"full_table.tsv"), magid=IDS)
localrules: run_make_dir_maker, static_file, configure_busco, download_lineage, download_busco 

rule metaeuk:
	input: mag=os.path.join(MAGDIR,'{magid}.fa'), db=config['metaeuk_db']
	output: 
		faa='metaeuk/{magid}.faa',
		gff='metaeuk/{magid}.faa.headersMap.tsv'
	params: t = os.path.join(TEMP,'{magid}-tmp'),other = '--protein 1 --threads 16' 
	conda: "envs/metaeuk.yaml"
	shell:'''
          mkdir -p {params.t}
          metaeuk easy-predict {input.mag} {input.db} {output.faa} {params.t} {params.other}
          '''

rule train_genemark:
    input: mag=os.path.join(MAGDIR,'{magid}.fa')
    output: 'genemark/{magid}/output/gmhmm.mod'
    params: folder_name='genemark/{magid}', tmpdir='/vortexfs1/scratch/halexander/genemarktmp/{magid}', min_contig=10000
    conda: 'envs/genemaker.yaml'
    shell:'''  
          mkdir -p {params.tmpdir}
          export TMPDIR={params.tmpdir}
          mkdir -p {params.folder_name}
          cd {params.folder_name}
          /vortexfs1/home/halexander/bin/gm_et_linux_64/gmes_petap.pl --ES --min_contig {params.min_contig} --sequence ../../{input.mag}  
          '''

rule run_make_dir_maker:
    input: 'genemark/{magid}/output/gmhmm.mod'
    output: os.path.join('maker2','{magid}','maker_opts.ctl')
    params: name='{magid}',outdir=os.path.join('maker2','{magid}')
    shell:'''
          mkdir -p {params.outdir}
          cp static/maker_opts.ctl {output}
          cd {params.outdir}
          sed -e 's/XXX/{params.name}/g' maker_opts.ctl > tmp
          mv tmp maker_opts.ctl
          '''

rule run_maker: 
    input: 
        mag=os.path.join(MAGDIR,'{magid}.fa'), 
        #gmes='genemark/{magid}/output/gmhmm.mod',
        conf=os.path.join('maker2','{magid}','maker_opts.ctl')
    output: 'maker2/{magid}/{magid}.maker.output/{magid}_master_datastore_index.log'
    params: maker_dir=os.path.join('maker2','{magid}'), tmpdir='/vortexfs1/scratch/halexander/maker2tmp/{magid}'
    conda: 'envs/maker2.yaml'
    shell:'''
          mkdir -p {params.tmpdir}
          export TMPDIR={params.tmpdir}
          cd {params.maker_dir}
          maker -g ../../{input.mag} -c 6 maker_opts.ctl 
          '''

rule combine_maker:
    input:'maker2/{magid}/{magid}.maker.output/{magid}_master_datastore_index.log'
    output: 'maker2/{magid}/{magid}.all.maker.genemark.proteins.fasta','maker2/{magid}/{magid}.all.maker.genemark.transcripts.fasta','maker2/{magid}/{magid}.all.maker.proteins.fasta','maker2/{magid}/{magid}.all.maker.transcripts.fasta'
    params: name='maker2/{magid}/{magid}'
    conda: 'envs/maker2.yaml'
    shell:''' 
          fasta_merge -d {input} -o {params.name}
          '''

rule combine_maker_gff3:
    input:'maker2/{magid}/{magid}.maker.output/{magid}_master_datastore_index.log'
    output: 'maker2/{magid}/{magid}.gff3'
    conda: 'envs/maker2.yaml'
    shell:''' 
          gff3_merge -d {input} -o {output}
          '''
rule pfam:
    input: 'maker2/{magid}/{magid}.all.maker.proteins.fasta'
    output: 'pfam/{magid}.pfam'
    threads: 8
    params:
        pfam=config["pfam_dir"]
    conda: "envs/pfam.yaml"
    shell:'''
          pfam_scan.pl -outfile {output} -as -cpu {threads} -fasta {input} -dir {params.pfam}
          '''
rule compute_diamond_index:
    input: KEGG_DB
    output: DIAMOND_OUT+'.dmnd'
    conda:
        "envs/diamond.yaml"
    params:
        db = DIAMOND_OUT
    shell:
        """
        diamond makedb --in {input} --db {params.db}
        """

rule diamond_map:
    input:
        dmnd = DIAMOND_OUT + '.dmnd',
        pep = 'maker2/{magid}/{magid}.all.maker.proteins.fasta'
    output: 'kegg/diamond/{magid}.blastp.tab'
    conda:
        "envs/diamond.yaml"
    params:
        other="--outfmt 6 --sensitive --block-size 2.0 -p 8"
    shell:
        """
        diamond blastp --db {input.dmnd} -q {input.pep} -o {output} {params.other} 
        """

rule static_file:
    output: os.path.join('static','config.ini')
    conda: 
        'envs/busco.yaml'
    shell:
        '''
        busco_configurator.py $BUSCO_CONFIG_FILE {output}
        
        '''

rule download_busco:
    output: os.path.join('busco_downloads', 'file_versions.tsv')
    params: url='https://busco-data.ezlab.org/v4/data/', 
    shell:
        '''
        wget {params.url}file_versions.tsv
        cp file_versions.tsv {output}
        '''
rule download_lineage: 
    input: os.path.join('busco_downloads', 'file_versions.tsv')
    output: directory(os.path.join('busco_downloads', 'lineages', BUSCO_GROUP))
    params: url='https://busco-data.ezlab.org/v4/data/', bg=BUSCO_GROUP
    shell: 
        '''
        date=$(grep ^{params.bg} {input}| cut -f2)
        wget {params.url}lineages/{params.bg}.$date.tar.gz
        tar -xvzf {params.bg}.$date.tar.gz
        rm {params.bg}.$date.tar.gz
        mv {params.bg} {output}
        '''

rule configure_busco:
    input:  fastafile = 'maker2/{magid}/{magid}.all.maker.proteins.fasta', conf = os.path.join('static', 'config.ini')
    output: os.path.join('busco','conf-files','{magid}.ini')
    params: directory = os.path.join('busco',BUSCO_GROUP), name='{magid}'
    shell:
        '''
        cp {input.conf} {output}
        sed -i 's%;out = BUSCO_run%out = {params.name}%' {output} 
        sed -i 's%;out_path = /path/to/output_folder%out_path = {params.directory}%' {output}
        sed -i 's%;download_path = ./busco_downloads/%download_path = ./busco_downloads/%' {output}
        '''

rule busco:  
    input:
        fastafile = 'maker2/{magid}/{magid}.all.maker.proteins.fasta', config = os.path.join('busco','conf-files','{magid}.ini'), busco_data = directory(os.path.join('busco_downloads', 'lineages', BUSCO_GROUP)) 

    output: 
        os.path.join('busco',BUSCO_GROUP,'{magid}','run_'+BUSCO_GROUP, "full_table.tsv")
    params: 
        busco_level = BUSCO_GROUP,
        CPUs = 8
    conda:
        'envs/busco.yaml'
    shell:
        '''
        busco -i {input.fastafile} -f -l {params.busco_level} -m protein --cpu {params.CPUs} --config {input.config} --offline 
        '''
