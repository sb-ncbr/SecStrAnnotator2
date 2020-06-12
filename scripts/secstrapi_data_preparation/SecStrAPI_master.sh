#!/bin/bash

#####################################################################################################################
# SETTINGS

TODAY="20200323"
API_VERSION="1.0"
N_THREADS="8"

FAMILY="CytochromesP450"
CATH_FAMILY_ID="1.10.630.10"
PFAM_FAMILY_ID="PF00067"

DATA_DIR="$HOME/Workspace/C#/SecStrAnnot2_data/SecStrAPI/$FAMILY-$TODAY"
STRUCTURE_CACHE_DIR="$DATA_DIR/../cached_structures/"  # Optional

TEMPLATE="2NNJ,A,:"  # Uppercase to avoid collisions with query protein
TEMPLATE_ANNOTATION_FILE="$DATA_DIR/../2NNJ-template-strict.sses.json"
TEMPLATE_STRUCTURE_FILE="$DATA_DIR/../2NNJ-canonical_rotation.cif"
SECSTRANNOTATOR_OPTIONS="--soft  --label2auth  --maxmetric 25,0.5,0.5  --verbose"
ALIGNED_SSE_LABELS="B,C,E,H,I,J,K,L,J',K',K'',1-1,1-2,1-3,1-4,1-5,2-2,3-1,3-2"  # List of SSEs for which a generic residue numbering should be established (comma-separated list or "all")

#####################################################################################################################
# SCRIPT PATHS

SCRIPT_DIR=`dirname $0`  # Location of this script
DOMAINS_FROM_PDBEAPI="$SCRIPT_DIR/domains_from_pdbeapi.py"
DOMAINS_TO_SECSTRAPI_FORMAT="$SCRIPT_DIR/domain_lists_to_SecStrAPI_format.py"
GET_TAXIDS="$SCRIPT_DIR/get_taxids.py"
CLASSIFY_TAXIDS="$SCRIPT_DIR/classify_taxids.py"
SELECT_BEST_DOMAINS="$SCRIPT_DIR/select_best_domain_per_uniprot.py"
SELECT_DOMAINS_BY_TAXON_GROUP="$SCRIPT_DIR/select_by_taxon_group.py"
SIMPLIFY_DOMAIN_LIST="$SCRIPT_DIR/simplify_domain_list.py"
EXTRACT_PDB_DOMAIN_LIST="$SCRIPT_DIR/extract_pdb_domain_list.py"
DOWNLOAD_DOMAINS="$SCRIPT_DIR/download_from_pdbe.py"
SECSTRANNOTATOR_BATCH="$SCRIPT_DIR/SecStrAnnotator_batch.py"
SECSTRANNOTATOR_DLL="$SCRIPT_DIR/SecStrAnnotator.dll"
COLLECT_ANNOTATIONS="$SCRIPT_DIR/collect_annotations.py"
SELECT_SSE_FIELDS="$SCRIPT_DIR/select_sse_fields.py"
DIVIDE_ANNOTATIONS="$SCRIPT_DIR/divide_annotations_by_pdb.py"
EXTRACT_SEQUENCES="$SCRIPT_DIR/extract_sequences.py"
ALIGN_SEQUENCES="$SCRIPT_DIR/align_sequences.py"
ADD_REFERENCE_RESIDUES="$SCRIPT_DIR/add_reference_residues.py"
GET_FULL_SEQUENCES="$SCRIPT_DIR/get_full_sequences.py"
JSON_TO_TSV="$SCRIPT_DIR/annotation_json_to_tsv.py"
JSON_TO_BULGES_TSV="$SCRIPT_DIR/annotation_json_to_bulges_tsv.py"

#####################################################################################################################
# MAIN PIPELINE

echo "Results in $DATA_DIR"

echo; echo '=== Get domains from CATH and Pfam ==='
mkdir  -p  $DATA_DIR
python3  $DOMAINS_FROM_PDBEAPI  --join_domains_in_chain  $CATH_FAMILY_ID  >  $DATA_DIR/set_cath_$TODAY.simple.json 
python3  $DOMAINS_FROM_PDBEAPI  --join_domains_in_chain  $PFAM_FAMILY_ID  >  $DATA_DIR/set_pfam_$TODAY.simple.json 

echo; echo '=== Merge domain lists and format them in SecStrAPI format (also adds UniProt refs) ==='
python3  $DOMAINS_TO_SECSTRAPI_FORMAT  --api_version $API_VERSION  \
    CATH  $CATH_FAMILY_ID  $DATA_DIR/set_cath_$TODAY.simple.json  \
    Pfam  $PFAM_FAMILY_ID  $DATA_DIR/set_pfam_$TODAY.simple.json  \
    >  $DATA_DIR/set_ALL_$TODAY.json

echo; echo '=== Get NCBI taxons and groups (Euka/Bact/Arch/Viru) ==='
wget  ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz  -O $DATA_DIR/taxdump.tar.gz
tar  xzf  $DATA_DIR/taxdump.tar.gz  nodes.dmp  -O  > $DATA_DIR/ncbi_taxonomy_nodes.dmp
python3  $GET_TAXIDS  $DATA_DIR/set_ALL_$TODAY.json  >  $DATA_DIR/domain_taxons.tsv
python3  $CLASSIFY_TAXIDS  $DATA_DIR/ncbi_taxonomy_nodes.dmp  $DATA_DIR/domain_taxons.tsv  >  $DATA_DIR/domain_taxons_groups.tsv

echo; echo '=== Select nonredundant set (with best quality) ==='
python3  $SELECT_BEST_DOMAINS  $DATA_DIR/set_ALL_$TODAY.json  >  $DATA_DIR/set_NR_$TODAY.json
python3  $SELECT_DOMAINS_BY_TAXON_GROUP  $DATA_DIR/set_NR_$TODAY.json  $DATA_DIR/domain_taxons_groups.tsv  Bact  >  $DATA_DIR/set_NR_Bact_$TODAY.json
python3  $SELECT_DOMAINS_BY_TAXON_GROUP  $DATA_DIR/set_NR_$TODAY.json  $DATA_DIR/domain_taxons_groups.tsv  Euka  >  $DATA_DIR/set_NR_Euka_$TODAY.json

echo; echo '=== Simplify domain lists (for SecStrAnnotator) ==='
python3  $SIMPLIFY_DOMAIN_LIST  $DATA_DIR/set_ALL_$TODAY.json  >  $DATA_DIR/set_ALL_$TODAY.simple.json
python3  $SIMPLIFY_DOMAIN_LIST  $DATA_DIR/set_NR_$TODAY.json  >  $DATA_DIR/set_NR_$TODAY.simple.json
python3  $EXTRACT_PDB_DOMAIN_LIST  $DATA_DIR/set_ALL_$TODAY.json  >  $DATA_DIR/AnnotationList.json

echo; echo '=== Download CIF files ==='
python3  $DOWNLOAD_DOMAINS  $DATA_DIR/set_ALL_$TODAY.simple.json  $DATA_DIR/structures/  --format cif  --no_gzip  --cache $STRUCTURE_CACHE_DIR

echo; echo '=== Annotate ==='
cp  $TEMPLATE_ANNOTATION_FILE  $DATA_DIR/structures/${TEMPLATE:0:4}-template.sses.json
cp  $TEMPLATE_STRUCTURE_FILE  $DATA_DIR/structures/${TEMPLATE:0:4}.cif
python3  $SECSTRANNOTATOR_BATCH  --dll $SECSTRANNOTATOR_DLL \
    --threads $N_THREADS  --options " $SECSTRANNOTATOR_OPTIONS " \
    $DATA_DIR/structures  $TEMPLATE  $DATA_DIR/set_ALL_$TODAY.simple.json 

echo; echo '=== Collect annotations and put them to SecStrAPI format ==='
python3  $COLLECT_ANNOTATIONS  $DATA_DIR/set_ALL_$TODAY.json  $DATA_DIR/structures/  >  $DATA_DIR/annotations_ALL.json
python3  $COLLECT_ANNOTATIONS  $DATA_DIR/set_NR_$TODAY.json  $DATA_DIR/structures/  >  $DATA_DIR/annotations_NR.json
python3  $COLLECT_ANNOTATIONS  $DATA_DIR/set_NR_Bact_$TODAY.json  $DATA_DIR/structures/  >  $DATA_DIR/annotations_NR_Bact.json
python3  $COLLECT_ANNOTATIONS  $DATA_DIR/set_NR_Euka_$TODAY.json  $DATA_DIR/structures/  >  $DATA_DIR/annotations_NR_Euka.json
python3  $EXTRACT_SEQUENCES  $DATA_DIR/annotations_ALL.json  $DATA_DIR/sequences_ALL/
python3  $EXTRACT_SEQUENCES  $DATA_DIR/annotations_NR.json  $DATA_DIR/sequences_NR/
python3  $EXTRACT_SEQUENCES  $DATA_DIR/annotations_NR_Bact.json  $DATA_DIR/sequences_NR_Bact/
python3  $EXTRACT_SEQUENCES  $DATA_DIR/annotations_NR_Euka.json  $DATA_DIR/sequences_NR_Euka/

echo; echo '=== Perform no-gap sequence alignment and create sequence logos (from Set-NR) ==='
python3  $ALIGN_SEQUENCES  $DATA_DIR/annotations_NR.json  --alignments_dir $DATA_DIR/aligments_NR/  --trees_dir $DATA_DIR/trees_NR/  --logos_dir $DATA_DIR/logos_NR/
python3  $ALIGN_SEQUENCES  $DATA_DIR/annotations_NR_Bact.json  --alignments_dir $DATA_DIR/aligments_NR_Bact/  --trees_dir $DATA_DIR/trees_NR_Bact/  --logos_dir $DATA_DIR/logos_NR_Bact/
python3  $ALIGN_SEQUENCES  $DATA_DIR/annotations_NR_Euka.json  --alignments_dir $DATA_DIR/aligments_NR_Euka/  --trees_dir $DATA_DIR/trees_NR_Euka/  --logos_dir $DATA_DIR/logos_NR_Euka/

echo; echo '=== Realign sequences from Set-ALL to the alignment from Set-NR and add reference residue information ==='
python3  $ADD_REFERENCE_RESIDUES  $DATA_DIR/annotations_NR.json  $DATA_DIR/aligments_NR/  --labels $ALIGNED_SSE_LABELS  --label2auth_dir $DATA_DIR/structures/  >  $DATA_DIR/annotations_with_reference_residues_NR.json
python3  $ADD_REFERENCE_RESIDUES  $DATA_DIR/annotations_ALL.json  $DATA_DIR/aligments_NR/  --labels $ALIGNED_SSE_LABELS  --label2auth_dir $DATA_DIR/structures/  >  $DATA_DIR/annotations_with_reference_residues_ALL.json

echo; echo '=== Divide annotations into per-PDB files ==='
python3  $SELECT_SSE_FIELDS  $DATA_DIR/annotations_with_reference_residues_ALL.json  >  $DATA_DIR/annotations_with_reference_residues_ALL-selected_fields.json
python3  $DIVIDE_ANNOTATIONS  $DATA_DIR/annotations_with_reference_residues_ALL-selected_fields.json  $DATA_DIR/annotations_ALL/  --min_dir $DATA_DIR/annotations_ALL_min/

echo; echo '=== Prepare TSV tables for analyses ==='
python3  $JSON_TO_TSV  --add_missing_sses  $DATA_DIR/annotations_with_reference_residues_NR.json  >  $DATA_DIR/annotations_with_reference_residues_NR.tsv
python3  $JSON_TO_TSV  --add_missing_sses  $DATA_DIR/annotations_with_reference_residues_ALL.json  >  $DATA_DIR/annotations_with_reference_residues_ALL.tsv
python3  $JSON_TO_BULGES_TSV  $DATA_DIR/annotations_with_reference_residues_NR.json  >  $DATA_DIR/beta_bulges_NR.tsv
python3  $JSON_TO_BULGES_TSV  $DATA_DIR/annotations_with_reference_residues_ALL.json  >  $DATA_DIR/beta_bulges_ALL.tsv

echo; echo '=== Get full sequences ==='
python3  $GET_FULL_SEQUENCES  $DATA_DIR/set_NR_$TODAY.json  >  $DATA_DIR/set_NR_$TODAY.fasta
