#! /usr/bin/env python
import pandas as pd
import glob
import os
import argparse

def main():
    p = argparse.ArgumentParser()
    p.add_argument('-d', '--directory')
    p.add_argument('-o', '--outfile')
    args = p.parse_args()
    outdf = pd.DataFrame(columns = ['Complete', 'Duplicated', 'Fragmented', 'Missing'])
    for n, tab in enumerate(glob.glob(os.path.join(args.directory, 'eukaryota_odb10', '*', 'run_eukaryota_odb10', 'full_table.tsv'))):
        name = tab.split('/')[-3]  
        df = pd.read_csv(tab, sep='\t', header=2)
        gp = df.groupby(['# Busco id','Status']).count().reset_index().groupby('Status').count().iloc[:,0] 
        outdf.loc[name]=gp
    outdf = outdf.fillna(0)
    outdf['Complete_ALL']=outdf.Complete + outdf.Duplicated
    outdf['Complete_Duplciated']=outdf.Duplicated
    outdf['Complete_Single'] = outdf.Complete
    outdf = outdf.drop(['Duplicated','Complete'], axis=1)
    outdf.to_csv(args.outfile, sep='\t')

if __name__ == '__main__':
    main()
