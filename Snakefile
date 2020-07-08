shell.executable("/bin/bash")
#shell.prefix("source $HOME/.bashrc; ")

import os

configfile: "config.yaml"
TEMP=config['temp_dir']
IDS, = glob_wildcards("euk-mags/{id}.fa")

rule all:
	 input: expand("output/{magid}.faa", magid=IDS), expand("output/{magid}.faa.headersMap.tsv", magid=IDS),
            #expand('genemark/{magid}/output/gmhmm.mod', magid=IDS)
            expand('maker2/{magid}/{magid}.all.maker.genemark.proteins.fasta',magid=IDS)
localrules: run_make_dir_maker

rule metaeuk:
	input: mag='euk-mags/{magid}.fa',
            db=config['metaeuk_db']
	output: 
		faa='output/{magid}.faa',
		gff='output/{magid}.faa.headersMap.tsv'
	params: t = os.path.join(TEMP,'{magid}-tmp'),other = '--protein 1 --threads 16' 
	conda: "envs/metaeuk.yaml"
	shell:'''
          mkdir -p {params.t}
          metaeuk easy-predict {input.mag} {input.db} {output.faa} {params.t} {params.other}
          '''

rule train_genemark:
    input: mag='euk-mags/{magid}.fa'
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
        mag='euk-mags/{magid}.fa',
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
 

    
        
