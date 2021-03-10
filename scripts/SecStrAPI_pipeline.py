'''
This Python3 script runs the full pipeline for secondary structure annotation of a protein family, 
including: downloading the list of family members defined by CATH and Pfam, downloading their structures, 
selecting a non-redundant set, annotation, multiple sequence alignment for individual SSEs, 
formatting into SecStrAPI format, formatting into TSV format for further analyses.

Before running, make your copy of SecStrAPI_master_settings.json and modify it 
to set your family of interest, data directory, options, paths to your template annotation etc.
Unwanted steps of the pipeline can be commented out in the main function.

Example usage:
    python3  SecStrAPI_master.py  SecStrAPI_master_settings.json 
'''

import argparse
from typing import Dict, Any
import os
from os import path
import shutil
import json
from collections import namedtuple
import glob

import lib
from lib import RedirectIO
import SecStrAnnotator_batch
from secstrapi_data_preparation import domains_from_pdbeapi
from secstrapi_data_preparation import domain_lists_to_SecStrAPI_format
from secstrapi_data_preparation import get_taxids
from secstrapi_data_preparation import classify_taxids
from secstrapi_data_preparation import select_best_domain_per_uniprot
from secstrapi_data_preparation import select_by_taxon_group
from secstrapi_data_preparation import simplify_domain_list
from secstrapi_data_preparation import extract_pdb_domain_list
from secstrapi_data_preparation import download_from_pdbe
from secstrapi_data_preparation import collect_annotations
from secstrapi_data_preparation import extract_sequences
from secstrapi_data_preparation import align_sequences
from secstrapi_data_preparation import add_reference_residues
from secstrapi_data_preparation import select_sse_fields
from secstrapi_data_preparation import divide_annotations_by_pdb
from secstrapi_data_preparation import annotation_json_to_tsv
from secstrapi_data_preparation import annotation_json_to_bulges_tsv
from secstrapi_data_preparation import get_full_sequences

Settings = namedtuple('Settings', [
    'data_dir', 
    'cath_family_id', 
    'pfam_family_id', 
    'api_version', 
    'structure_cache_dir',
    'template_domain',
    'template_annotation_file',
    'template_structure_file',
    'n_threads',
    'secstrannotator_dll',
    'secstrannotator_options',
    'sses_for_generic_numbering'
    ])

CHECKPOINT_FILE = 'checkpoints.txt'

#  MAIN  #####################################################################################

def parse_args() -> Dict[str, Any]:
    '''Parse command line arguments.'''
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('settings_file', help='JSON file with settings', type=str)
    parser.add_argument('--resume', help=f'Read checkpoints from {CHECKPOINT_FILE} and continue the pipeline from the last completed task', action='store_true')
    args = parser.parse_args()
    return vars(args)


def main(settings_file: str, resume: bool = False) -> None:
    '''Run the full pipeline for secondary structure annotation of a protein family'''
    
    with open(settings_file, 'r', encoding=lib.DEFAULT_ENCODING) as f:
        settings_dict = json.load(f)
        settings = Settings(**settings_dict)

    if settings.secstrannotator_dll is not None:
        secstrannotator_dll = path.join(path.dirname(path.abspath(__file__)), settings.secstrannotator_dll)
    else:
        secstrannotator_dll = path.join(path.dirname(path.abspath(__file__)), '..', 'SecStrAnnotator.dll')

    print(f'Results will be in {settings.data_dir}')
    os.makedirs(settings.data_dir, exist_ok=True)
    os.chdir(settings.data_dir)
    os.makedirs('structures', exist_ok=True)
    template_id = settings.template_domain.split(',')[0]
    shutil.copy(settings.template_annotation_file, path.join('structures', f'template_{template_id}-template.sses.json'))
    shutil.copy(settings.template_structure_file, path.join('structures', f'template_{template_id}.cif'))
    
    pipeline = lib.Pipeline(checkpoint_file=CHECKPOINT_FILE)

    # Get domains from CATH and Pfam
    pipeline.add_task(None, print, '\n=== Get domains from CATH and Pfam ===')
    pipeline.add_task('get domains - CATH', domains_from_pdbeapi.main, settings.cath_family_id, join_domains_in_chain=True, stdout='set_cath.simple.json')
    pipeline.add_task('get domains - Pfam', domains_from_pdbeapi.main, settings.pfam_family_id, join_domains_in_chain=True, stdout='set_pfam.simple.json')

    # Merge domain lists and format them in SecStrAPI format (also adds UniProt refs)
    pipeline.add_task(None, print, '\n=== Merge domain lists and format them in SecStrAPI format (also adds UniProt refs) ===')
    pipeline.add_task('merge domain lists', domain_lists_to_SecStrAPI_format.main,
        ['CATH', settings.cath_family_id, 'set_cath.simple.json', 'Pfam', settings.pfam_family_id, 'set_pfam.simple.json'],
        api_version=settings.api_version,
        stdout='set_ALL.json')

    # Get NCBI taxons and groups (Euka/Bact/Arch/Viru)
    pipeline.add_task(None, print, '\n=== Get NCBI taxons and groups (Euka/Bact/Arch/Viru) ===')
    pipeline.add_task('download NCBI taxonomy', lib.get_from_ftp, 'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz', 'taxdump.tar.gz')
    pipeline.add_task('unzip NCBI taxonomy', lib.extract_from_tar, 'taxdump.tar.gz', 'nodes.dmp', 'ncbi_taxonomy_nodes.dmp')
    pipeline.add_task('get taxids for Set-ALL', get_taxids.main, 'set_ALL.json', stdout='domains_taxons.tsv')
    pipeline.add_task('classify domains by superkingdom', classify_taxids.main, 'ncbi_taxonomy_nodes.dmp', 'domains_taxons.tsv', stdout='domains_taxons_groups.tsv')

    # Select nonredundant set (with best quality)
    pipeline.add_task(None, print, '\n=== Select nonredundant set (with best quality) ===')
    pipeline.add_task('select Set-NR', select_best_domain_per_uniprot.main, 'set_ALL.json', stdout='set_NR.json')
    pipeline.add_task('select Set-NR-Bact', select_by_taxon_group.main, 'set_NR.json', 'domains_taxons_groups.tsv', 'Bact', stdout='set_NR_Bact.json')
    pipeline.add_task('select Set-NR-Euka', select_by_taxon_group.main, 'set_NR.json', 'domains_taxons_groups.tsv', 'Euka', stdout='set_NR_Euka.json')

    # Simplify domain lists (for SecStrAnnotator)
    pipeline.add_task(None, print, '\n=== Simplify domain lists (for SecStrAnnotator) ===')
    pipeline.add_task('simplify domain list - Set-ALL', simplify_domain_list.main, 'set_ALL.json', stdout='set_ALL.simple.json')
    pipeline.add_task('simplify domain list Set-NR', simplify_domain_list.main, 'set_NR.json', stdout='set_NR.simple.json')
    pipeline.add_task('create domain list for SecStrAPI', extract_pdb_domain_list.main, 'set_ALL.json', stdout='AnnotationList.json')

    # Download CIF files
    pipeline.add_task(None, print, '\n=== Download CIF files ===')
    pipeline.add_task('download CIF files', download_from_pdbe.main, 'set_ALL.simple.json', 'structures', structure_format='cif', no_gzip=True, cache=settings.structure_cache_dir)

    # Annotate
    pipeline.add_task(None, print, '\n=== Annotate ===')
    pipeline.add_task('run SecStrAnnotator', SecStrAnnotator_batch.main, 'structures', f'template_{settings.template_domain}', 'set_ALL.simple.json', 
        threads=int(settings.n_threads), dll=secstrannotator_dll, options=settings.secstrannotator_options)
    pipeline.add_task('remove accidentally downloaded files', lambda: [os.remove(filename) for filename in glob.glob('*.cif')])  # Remove files accidentally downloaded by PyMOL

    # Collect annotations and put them to SecStrAPI format
    pipeline.add_task(None, print, '\n=== Collect annotations and put them to SecStrAPI format ===')
    pipeline.add_task('collect annotations into annotations_ALL.json', collect_annotations.main, 'set_ALL.json', 'structures', stdout='annotations_ALL.json')
    pipeline.add_task('collect annotations into annotations_NR.json', collect_annotations.main, 'set_NR.json', 'structures', stdout='annotations_NR.json')
    pipeline.add_task('collect annotations into annotations_NR_Bact.json', collect_annotations.main, 'set_NR_Bact.json', 'structures', stdout='annotations_NR_Bact.json')
    pipeline.add_task('collect annotations into annotations_NR_Euka.json', collect_annotations.main, 'set_NR_Euka.json', 'structures', stdout='annotations_NR_Euka.json')
    pipeline.add_task('extract sequences - Set-ALL', extract_sequences.main, 'annotations_ALL.json', 'sequences_ALL')
    pipeline.add_task('extract sequences - Set-NR', extract_sequences.main, 'annotations_NR.json', 'sequences_NR')
    pipeline.add_task('extract sequences - Set-NR-Bact', extract_sequences.main, 'annotations_NR_Bact.json', 'sequences_NR_Bact')
    pipeline.add_task('extract sequences - Set-NR-Euka', extract_sequences.main, 'annotations_NR_Euka.json', 'sequences_NR_Euka')

    # Perform no-gap sequence alignment and create sequence logos (from Set-NR)
    pipeline.add_task(None, print, '\n=== Perform no-gap sequence alignment and create sequence logos (from Set-NR) ===')
    pipeline.add_task('align sequences - Set-NR', align_sequences.main, 'annotations_NR.json', alignments_dir='aligments_NR', trees_dir='trees_NR', logos_dir='logos_NR', matrices_dir='alignment_matrices_NR', labels_for_matrices=settings.sses_for_generic_numbering, stdout='logo_statistics_NR.tsv')
    pipeline.add_task('copy matrix file', shutil.copy, path.join('alignment_matrices_NR', 'ALL.json'), 'alignment_matrices_NR.json')
    pipeline.add_task('align sequences - Set-NR-Bact', align_sequences.main, 'annotations_NR_Bact.json', alignments_dir='aligments_NR_Bact', trees_dir='trees_NR_Bact', logos_dir='logos_NR_Bact', stdout='logo_statistics_NR_Bact.tsv')
    pipeline.add_task('align sequences - Set-NR-Euka', align_sequences.main, 'annotations_NR_Euka.json', alignments_dir='aligments_NR_Euka', trees_dir='trees_NR_Euka', logos_dir='logos_NR_Euka', stdout='logo_statistics_NR_Euka.tsv')

    # Realign sequences from Set-ALL to the alignment from Set-NR and add reference residue information
    pipeline.add_task(None, print, '\n=== Realign sequences from Set-ALL to the alignment from Set-NR and add reference residue information ===')
    pipeline.add_task('add reference residues - Set-NR', add_reference_residues.main, 'annotations_NR.json', 'aligments_NR', labels=settings.sses_for_generic_numbering, label2auth_dir='structures', stdout='annotations_with_reference_residues_NR.json')
    pipeline.add_task('add reference residues - Set-ALL', add_reference_residues.main, 'annotations_ALL.json', 'aligments_NR', labels=settings.sses_for_generic_numbering, label2auth_dir='structures', stdout='annotations_with_reference_residues_ALL.json')

    # Divide annotations into per-PDB files
    pipeline.add_task(None, print, '\n=== Divide annotations into per-PDB files ===')
    pipeline.add_task('remove unnecessary fields from annotations', select_sse_fields.main, 'annotations_with_reference_residues_ALL.json', stdout='annotations_with_reference_residues_ALL-selected_fields.json')
    pipeline.add_task('divide annotations 1-PDB-per-file', divide_annotations_by_pdb.main, 'annotations_with_reference_residues_ALL-selected_fields.json', 'annotations_ALL', min_dir='annotations_ALL_min')

    # Prepare TSV tables for analyses
    pipeline.add_task(None, print, '\n=== Prepare TSV tables for analyses ===')
    pipeline.add_task('create TSV table with annotations - Set-NR', annotation_json_to_tsv.main, 'annotations_with_reference_residues_NR.json', add_missing_sses=True, stdout='annotations_with_reference_residues_NR.tsv')
    pipeline.add_task('create TSV table with annotations - Set-ALL', annotation_json_to_tsv.main, 'annotations_with_reference_residues_ALL.json', add_missing_sses=True, stdout='annotations_with_reference_residues_ALL.tsv')
    pipeline.add_task('create TSV table with beta-bulges - Set-NR', annotation_json_to_bulges_tsv.main, 'annotations_with_reference_residues_NR.json', stdout='beta_bulges_NR.tsv')
    pipeline.add_task('create TSV table with beta-bulges - Set-ALL', annotation_json_to_bulges_tsv.main, 'annotations_with_reference_residues_ALL.json', stdout='beta_bulges_ALL.tsv')

    # Get full sequences
    pipeline.add_task(None, print, '\n=== Get full sequences ===')
    pipeline.add_task('get full FASTA sequences', get_full_sequences.main, 'set_NR.json', stdout='set_NR.fasta')

    if resume:
        print('\nResuming pipeline')
        pipeline.resume()
    else:        
        pipeline.start()


if __name__ == '__main__':
    args = parse_args()
    main(**args)
