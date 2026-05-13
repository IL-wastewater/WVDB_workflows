#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse
from Bio import Entrez, SeqIO
# Set your email (required by NCBI)
Entrez.email = ""

def load_checkv(checkv_file):
    checkv_df = pd.read_csv(checkv_file, sep='\t')
    checkv_df.rename(columns={'contig_id':'contig'}, inplace=True)

    return checkv_df


def load_genomad(genomad_file):
    '''
    load geNomad virus summary and extract taxonomic levels to columns
    '''
    genomad_tax_df = pd.read_csv(genomad_file, sep='\t') # 20539
    genomad_tax_df = genomad_tax_df.rename(columns={'seq_name':'contig'})
    genomad_tax_df[['Domain', 'Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family']] = genomad_tax_df.taxonomy.str.split(';', expand=True)
    genomad_tax_df[['Domain', 'Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family']] = genomad_tax_df[['Domain', 'Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family']].replace('', np.nan)

    return genomad_tax_df


def load_rdrpcatch(rdrpcatch_file):
    '''
    load RdRpCATCH taxonomy and extract taxonomic levels to columns
    '''
    rdrp_df = pd.read_csv(rdrpcatch_file, sep='\t')
    rdrp_df = rdrp_df.rename(columns={'Contig_name':'contig'})
    rdrp_df['MMseqs_Taxonomy_2bLCA'] = rdrp_df['MMseqs_Taxonomy_2bLCA'].str.replace(' ', '_') # get rid of spaces in taxonomic names
    rdrp_df['Domain'] = rdrp_df.MMseqs_Taxonomy_2bLCA.str.extract(r'd_(\w+);?')
    rdrp_df['Realm'] = rdrp_df.MMseqs_Taxonomy_2bLCA.str.extract(r'd_\w+;-_(\w+);?')
    rdrp_df['Kingdom'] = rdrp_df.MMseqs_Taxonomy_2bLCA.str.extract(r';k_(\w+);?')
    rdrp_df['Phylum'] = rdrp_df.MMseqs_Taxonomy_2bLCA.str.extract(r';p_(\w+);?')
    rdrp_df['Order'] = rdrp_df.MMseqs_Taxonomy_2bLCA.str.extract(r';o_(\w+);?')
    rdrp_df['Class'] = rdrp_df.MMseqs_Taxonomy_2bLCA.str.extract(r';c_(\w+);?')
    rdrp_df['Family'] = rdrp_df.MMseqs_Taxonomy_2bLCA.str.extract(r';f_(\w+);?')
    rdrp_df['Genus'] = rdrp_df.MMseqs_Taxonomy_2bLCA.str.extract(r';g_(\w+);?')
    rdrp_df['Species'] = rdrp_df.MMseqs_Taxonomy_2bLCA.str.extract(r';s_([\w\d-]+);?')
    rdrp_df['Lowest'] = rdrp_df.MMseqs_Taxonomy_2bLCA.str.extract(r'_([\w\d-]+)$')

    rdrp_tax_df = rdrp_df[['contig','Domain', 'Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Lowest', 'MMseqs_Taxonomy_2bLCA', 'Best_hit_bitscore']].copy()
    #rdrp_tax_df = rdrp_tax_df.replace(np.nan, None)
    
    return rdrp_tax_df

def get_consensus_tax(rdrpcatch_df, genomad_df, checkv_df, outfile):
    '''make orders.csv file for RNAVirHost
    rnavirhost predict -i {fasta} --taxa {contig2order.csv} -o {outdir}
    note that outdir cannot already exist and contig2order requires specific headers (see renaming in function)
    '''
    orders_df = pd.merge(rdrpcatch_df[['contig', 'Order']], genomad_df[['contig', 'Order']], on='contig', how='outer', suffixes=('_R', '_G'))
    # add checkV list of contigs since some contigs were not hit by geNomad or RdRpCATCH
    orders_df = orders_df.merge(checkv_df[['contig']], how='outer')
    orders_df = orders_df.replace(np.nan, 'Unclassified') # consistent naming of unclassified
    # prioritize RdRpCATCH order, if available
    orders_df['Order_consensus'] = orders_df['Order_R']
    # use geNomad order if RdRpcatch didn't classify
    orders_df.loc[orders_df.Order_R == 'Unclassified', 'Order_consensus'] = orders_df.Order_G
    orders_df = orders_df[['contig', 'Order_consensus']]
    out_df = orders_df.rename(columns={'contig':'', 'Order_consensus':'y|virus order'})
    
    out_df.to_csv(outfile, index=False)
    
    return orders_df


def load_rnavirhost(rnavirhost):
    rnavirhost_df = pd.read_csv(rnavirhost, skiprows=1, names=['contig', 'Order', 'host_domain', 'host_chordata', 'evidence'])
    rnavirhost_df = rnavirhost_df.replace(np.nan, 'unknown')
    rnavirhost_df.rename(columns={'host_domain':'rnavirhost_domain', 'host_chordata':'rnavirhost_chordata'}, inplace=True)
    rnavirhost_df = rnavirhost_df[['contig', 'rnavirhost_domain', 'rnavirhost_chordata']]
    return rnavirhost_df
    

def process_blastn(tsv: str, db: str):
    '''
    starting with blast result ani file output from nayfach script (tsv) against a specific db (name of db)
    determine whether hits were genus or species level
    if there were multiple hits, select just the best one
    '''
    df = pd.read_csv(tsv, sep='\t')

    # filter for high identity and high coverage
    df['hit_level'] = 'none' 
    df.loc[(df.pid >=70) & (df.qcov >=85), 'hit_level'] = 'genus'
    df.loc[(df.pid >=95) & (df.qcov >=85), 'hit_level'] = 'species'
    df['idXcov'] = df.pid * df.qcov
    # if there are multiple hits, take the best one (groupby and idxmax)
    idx_max_values = df.groupby('qname')['idXcov'].idxmax()
    df = df.loc[idx_max_values]
    df.rename(columns={'qname':'contig',
                       'tname':f'tname_{db}',
                       'pid':f'pid_{db}',
                       'qcov':f'qcov_{db}',
                       'hit_level': f'hit_level_{db}'}, 
              inplace=True)
    df.drop(columns=['num_alns','tcov', 'idXcov'], inplace=True)

    return df


def batch_extract_hosts_lineages(genbank_ids):
    '''for a list of genbank ids, use a batch Entrez job to get 
    taxonomy, host, and isolation source if available'''
    
    ids_str = ",".join(genbank_ids)
    handle = Entrez.efetch(db="nucleotide", id=ids_str, rettype="gb", retmode="text")
    records = list(SeqIO.parse(handle, "genbank"))
    handle.close()
    results = []
    for record in records:
        # Extract hosts and isolation source
        hosts = []
        isolation_source = None
        for feature in record.features:
            # note: not all records have host and isolation source
            if feature.type == "source":
                if "host" in feature.qualifiers:
                    hosts.extend(feature.qualifiers["host"])
                if "isolation_source" in feature.qualifiers:
                    isolation_source = feature.qualifiers.get("isolation_source", [None])[0]
        
        hosts = ','.join(hosts) if hosts else None
        isolation_source = isolation_source if isolation_source else None
        
        # Extract lineage
        lineage = record.annotations.get("taxonomy", None)
        lineage = ';'.join(lineage) # make a string from the list
        # combine all and save
        results.append((record.id, hosts, isolation_source, lineage))
    return results


def get_host_lineage_blasthit(corent, target_col='tname', batch_size=100):
    '''
    get host and lineage of best blast hit, starting with an ani tsv from blastn output
    ani tsv generated by nayfach ani script
    return an updated ani tsv with host, isolation source, and lineage of the hit as listed in the genbank entry
    also write these new columns to their own tsv
    '''
    genbank_id_list = list(corent[target_col].unique()) # get dereplicated list of genbank_ids hit by blast nt
    blastn_results = []
    for i in range(0, len(genbank_id_list), batch_size): # NCBI recommends <=500 per request
        print(f'querying entrez with batch {i}')
        
        batch_ids = genbank_id_list[i:i+batch_size]
        blastn_results.extend(batch_extract_hosts_lineages(batch_ids))
    
    # Convert results to DataFrame and merge with original
    blastn_results_df = pd.DataFrame(blastn_results, columns=[target_col, 'hosts_ntBlastHit', 'isolation_source_ntBlastHit', 'lineage_ntBlastHit'])
    blastn_results_df['lowest_tax_blasthit'] = blastn_results_df.lineage_ntBlastHit.str.extract(r'.*;(.*$)')
    blastn_results_df = corent[['contig', target_col]].merge(blastn_results_df, on=target_col, how='left')
    
    return blastn_results_df


def main():
    parser = argparse.ArgumentParser(
        description="Parse input arguments for annotation processing."
    )

    parser.add_argument(
        "--annotation_dir",
        required=True,
        help="Path to the annotation directory"
    )

    parser.add_argument(
        "--out_dir",
        required=True,
        help="Path to the output directory"
    )

    parser.add_argument(
        "--ictv_simplified",
        required=True,
        help="Path to the ictv_simplified file"
    )

    args = parser.parse_args()

    # find files within the annotation_dir
    checkv = f'{args.annotation_dir}/votus_final_quality_summary.tsv'
    genomad = f'{args.annotation_dir}/votus_final_virus_summary_genomad.tsv'
    rdrpcatch = f'{args.annotation_dir}/votus_final_rdrpcatch_output_annotated.tsv'
    rnavirhost = f'{args.annotation_dir}/votus_final_rnavirhost.csv'
    core_nt_out = f'{args.annotation_dir}/votus_final-vs-corentJuly25.ani.tsv'
    imgvr_out = f'{args.annotation_dir}/votus_final-vs-IMGVRv5.ani.tsv'
    host_lineage_out = f'{args.out_dir}/votus_final-vs-corentJuly25_host_tax_info.tsv'
    taxhost_outfile = f'{args.out_dir}/votus_final_full_annotation.tsv'

    # load files
    checkv_df = load_checkv(checkv)
    genomad_tax_df = load_genomad(genomad)
    rdrp_tax_df = load_rdrpcatch(rdrpcatch)
    corent = process_blastn(core_nt_out, 'nt')
    imgvr = process_blastn(imgvr_out, 'imgvr')

    # pull info from NCBI entrez using accessions of best BLAST hits
    for_entrez = corent[corent.hit_level_nt.isin(['genus', 'species'])].copy()
    blastn_results_df = get_host_lineage_blasthit(for_entrez, target_col='tname_nt')
    blastn_results_df.to_csv(host_lineage_out, sep='\t', index=False)

    # checks and reporting:
    provirus_count = len(genomad_tax_df[genomad_tax_df.contig.str.contains('provirus')])
    high_kmerfreq_count = len(checkv_df[checkv_df.kmer_freq>1.2])
    if provirus_count > 0:
        contigs = genomad_tax_df[genomad_tax_df.contig.str.contains('provirus')].contig.to_list()
        print(f'geNomad reports {provirus_count} proviruses. Consider trimming before continuing: {contigs}')
    if high_kmerfreq_count > 0:
        contigs = checkv_df[checkv_df.kmer_freq>1.2].contig.to_list()
        print(f'checkV reports {high_kmerfreq_count} contigs with kmer_freq > 1.2. Consider trimming before continuing: {contigs}')
    
    if len(rdrp_tax_df) != len(rdrp_tax_df.contig.unique()):
        print('RdRpCATCH output has multiple RdRp genes per contig')

    print(f'Annotation counts\ncheckV: {len(checkv_df)}\nRdRpCATCH: {len(rdrp_tax_df)}\ngeNomad: {len(genomad_tax_df)}')

    # generate consensus taxonomy for RNAVirHost
    consensus_out='/Users/kantor4/data/trilab_ldrd/virus_db_analysis/clusters_v3/order_consensus.csv'
    orders_df = get_consensus_tax(rdrp_tax_df, genomad_tax_df, checkv_df, consensus_out)

    # after running rnavirhost predict using the consensus_orders.csv created above
    rnavirhost_df = load_rnavirhost(rnavirhost)
    print(f'RNAVirHost: {len(rnavirhost_df)}')


    # merge RdRpCatch, geNomad, checkV, RNAVirHost, annotations and best blast hit info
    taxhost_df = pd.merge(genomad_tax_df[['contig', 'Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family']],
                        rdrp_tax_df[['contig', 'Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Lowest']], 
                        on='contig', how='outer', suffixes=['_gNd', '_RdRp'])
    taxhost_df = taxhost_df.rename(columns={'Genus':'Genus_RdRp', 'Species':'Species_RdRp', 'Lowest':'Lowest_RdRp'}) # add suffix to these columns too
    taxhost_df = taxhost_df.merge(checkv_df, on='contig', how='outer')
    taxhost_df = taxhost_df.merge(rnavirhost_df, on='contig', how='outer')
    taxhost_df.loc[taxhost_df.rnavirhost_domain.isna(), 'rnavirhost_domain'] = 'Unknown' # fill unknowns
    taxhost_df = taxhost_df.merge(orders_df, on='contig', how='left')
    taxhost_df = taxhost_df.merge(corent, on='contig', how='left')
    taxhost_df = taxhost_df.merge(blastn_results_df, on=['contig', 'tname_nt'], how='left')
    taxhost_df = taxhost_df.merge(imgvr, how='left', on='contig')

    # get consensus family and then add ICTV info by family
    taxhost_df['Family_consensus'] = taxhost_df['Family_RdRp']
    taxhost_df.loc[taxhost_df.Family_RdRp.isna(), 'Family_consensus'] = taxhost_df['Family_gNd']
    taxhost_df.loc[taxhost_df.Family_consensus.isna(), 'Family_consensus'] = 'Unclassified'

    ictv_fam_df = pd.read_csv(args.ictv_simplified, sep='\t') # family level table from ICTV
    taxhost_df = taxhost_df.merge(ictv_fam_df, how='left', left_on='Family_consensus', right_on='Family_ICTV')
    taxhost_df.loc[taxhost_df.host_ICTV_simple.isna(), 'host_ICTV_simple'] = 'Unknown'

    # save final table
    taxhost_df.to_csv(taxhost_outfile, sep='\t', index=False)


if __name__ == '__main__':
    main()