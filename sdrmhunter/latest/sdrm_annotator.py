#!python3

import logging
import re
import subprocess

from functools import reduce

import os
import pandas as pd
import numpy as np
from Bio import SeqIO, Seq
from Bio.Alphabet import generic_dna

ALIGNED_GENE_SEQ_COL = 'alignedGeneSequences'
INPUT_SEQUENCE_COL = 'inputSequence'
SDRMS_COL = 'SDRMs'
GENE_COL = 'gene'


def extract_drms_from_sierra_json(json):
    """
    Extracts DRM-related information from json files obtained with [sierra](https://hivdb.stanford.edu/DR/webservices/)
    :param json: str, filepath to the sierra output json file
    :return: str, filepath for the output tab-delimited table containing DRM information:
    input_sequence_id as index, subtypeText and various DRMs as columns
    """
    df = pd.read_json(json)
    df[INPUT_SEQUENCE_COL] = df[INPUT_SEQUENCE_COL].apply(lambda d: d['header'])
    df.set_index(keys=INPUT_SEQUENCE_COL, drop=True, inplace=True)

    # split a list of genes inside the alignedGeneSequences column into several rows
    # then split the dictionaries inside s rows into several columns
    gene_df = df[ALIGNED_GENE_SEQ_COL].apply(pd.Series, 1).stack().apply(pd.Series)
    gene_df[GENE_COL] = gene_df[GENE_COL].apply(lambda d: d['name']).astype('category', ordered=True)
    gene_df[SDRMS_COL] = gene_df[SDRMS_COL].apply(lambda l: [d['text'] for d in l])
    gene_df.index = gene_df.index.droplevel(-1)

    # Put all the DRMs together and make them columns
    gene_df[SDRMS_COL] = \
        gene_df.apply(lambda row: {'%s:%s' % (row[GENE_COL], m): True for m in row[SDRMS_COL]}, axis=1)

    def join_dicts(ds):
        return reduce(lambda d1, d2: {**d1, **d2}, ds, {})

    gene_df = gene_df.groupby(gene_df.index)[SDRMS_COL].apply(list).apply(join_dicts).apply(pd.Series)

    lbls_to_drop = [ALIGNED_GENE_SEQ_COL]
    return df.drop(labels=lbls_to_drop, axis=1).join(gene_df)


def join_dicts(ds):
    return reduce(lambda d1, d2: {**d1, **d2}, ds, {})


class SierraException(Exception):

    def __init__(self, message):
        self.message = message


def extract_sdrms(fasta, tab, bunch_size=1000):
    """
    Extracts SDRM-related information for the sequences in fasta file
    using [sierra](https://hivdb.stanford.edu/DR/webservices/)
    :param fasta: str, filepath for the input sequences in fasta format
    :return: str, filepath for the output tab-delimited table containing DRM information:
    input_sequence_id as index, subtypeText and various DRMs as columns
    """
    sequences = SeqIO.parse(fasta, "fasta", alphabet=generic_dna)

    def next_n_sequences(seq_iterable, n=bunch_size, skip_ids=None):
        for _ in seq_iterable:
            if skip_ids and str(_.id) in skip_ids:
                continue
            yield _
            n -= 1
            if n == 0:
                break

    i = 0
    master_df = None

    temp_tab = '{}.temp.tab'.format(tab)
    if os.path.exists(temp_tab):
        master_df = pd.read_table(temp_tab, index_col=0, header=0)
        skip_ids = set(master_df.index.map(str))
    else:
        skip_ids = set()

    gql = '{}.gql'.format(fasta)
    create_gql_file(gql)

    while True:
        seq_bunch = [SeqIO.SeqRecord(id=seq.id, description=seq.description, seq=Seq.Seq(str(seq.seq).replace('_', '-')))
                     for seq in next_n_sequences(sequences, skip_ids=skip_ids)]
        if not seq_bunch:
            break
        fasta_i = '{}_{}.fa'.format(fasta, i)
        count = SeqIO.write(seq_bunch, fasta_i, "fasta")
        logging.info("Extracted %d sequences to analyse with sierrapy" % count)
        json_i = '{}.json'.format(fasta_i)

        success = False
        for i in range(100):
            try:
                subprocess.Popen(['sierrapy', 'fasta', fasta_i, '-o', json_i, '-q', gql], stdout=subprocess.PIPE)\
                    .communicate()
                logging.info('Successfully loaded sierra annotations.')
                success = True
                break
            except:
                logging.warning('Failed to load sierra annotations, retrying...')
                pass
        if not success:
            raise SierraException('Stopping as could not extract SDRM annotations after 100 tries, '
                                  'please make sure that you have internet access.')

        df = extract_drms_from_sierra_json(json_i)
        if master_df is None:
            master_df = df
        else:
            master_df = pd.concat([master_df, df], axis=0, ignore_index=False)
        master_df.to_csv(temp_tab, sep='\t')

        remove_file_if_you_can(fasta_i)
        remove_file_if_you_can(json_i)
        i += 1
    master_df.to_csv(tab, sep='\t')
    remove_file_if_you_can(temp_tab)
    remove_file_if_you_can(gql)


def prettify_sdrms(tab):
    """
    Prettify mutation file, by replacing True/False by resistant/sensitive,
    as well as searching for multiple possible letter outcome mutations
    and replacing the corresponding real ones with unknown state,
    i.e. M184IV would make M184I and M184V unknown

    :param tab: SDRM file
    :return: prettified SDRM file
    """
    df = pd.read_table(tab, header=0, index_col=0)
    mutations = [_ for _ in df.columns if _.startswith('RT:') or _.startswith('PR:')]
    for mutation in mutations:
        df[mutation] = df[mutation].fillna(False).astype(bool).map(
            {True: 'resistant', False: 'sensitive', None: 'sensitive'}).astype(str)

    # search for multiple possible letter outcome mutations and replace the corresponding real ones with unknown state,
    # i.e. M184IV would make M184I and M184V unknown
    for mutation in mutations:
        to_letters = re.search(r'[A-Z]+$', mutation)[0]
        if len(to_letters) == 1:
            continue
        start = re.sub(r'[A-Z]+$', '', mutation)
        for mut in ('{}{}'.format(start, letter) for letter in to_letters):
            if mut in df.columns:
                df[mut] = np.where(df[mutation] == 'resistant',
                                   (np.where(df[mut] == 'resistant', df[mut], None)), df[mut])
        df.drop(axis=1, labels=[mutation], inplace=True)
    df.to_csv(tab, sep='\t')


def remove_file_if_you_can(filepath):
    try:
        os.remove(filepath)
    except:
        logging.error('Could not remove the file: {}'.format(filepath))
        pass


def create_gql_file(gql):
    """
    Creates a SDRM gql file for [sierra](https://hivdb.stanford.edu/DR/webservices/).
    :param filepath: str, filepath to write the gql file to.
    :return: void
    """
    with open(gql, 'w+') as fp:
        fp.write('''
inputSequence {
    header
},
subtypeText,
alignedGeneSequences {
    gene {
        name 
    },
    SDRMs:mutations(filterOptions:[SDRM]) {
        text
    }
}''')


if '__main__' == __name__:
    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt="%Y-%m-%d %H:%M:%S",
                        filename=None)

    import argparse
    parser = argparse.ArgumentParser(description="Extracts SDRM annotations for a fasta alignment.")
    parser.add_argument('--fasta', required=True,
                        type=str, help="the input sequence alignment in fasta format.")
    parser.add_argument('--output', required=False, default=None,
                        type=str, help="the output SDRM annotation file in tab-delimited format.")
    params = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt="%Y-%m-%d %H:%M:%S",
                        filename=None)

    output = params.output if params.output else '{}.drms.tab'.format(params.fasta)

    extract_sdrms(params.fasta, output, 1000)
    prettify_sdrms(output)

