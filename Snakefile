shell.executable("/bin/bash")
#shell.prefix("source $HOME/.bashrc; ")

import os

configfile: "config.yaml"
TEMP=config['temp_dir']
IDS, = glob_wildcards("euk-mags/{id}.fa")

rule all:
	 input: expand("output/{magid}.faa", magid=IDS), expand("output/{magid}.faa.headersMap.tsv", magid=IDS)

rule metaeuk:
	input: mag='euk-mags/{magid}.fa',
            db=config['metaeuk_db']
	output: 
		faa='output/{magid}.faa',
		gff='output/{magid}.faa.headersMap.tsv'
	params: t = os.path.join(TEMP,'{magid}-tmp'),other = '--protein 1 --threads 16' 
	conda: "metaeuk.yaml"
	shell:'''
          mkdir -p {params.t}
          metaeuk easy-predict {input.mag} {input.db} {output.faa} {params.t} {params.other}
          '''
