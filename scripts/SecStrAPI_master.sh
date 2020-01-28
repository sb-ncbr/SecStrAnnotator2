#!/bin/bash
TODAY="20200128"
API_VERSION="1.0"
N_THREADS="8"

DATA_DIR="/home/adam/Workspace/C#/SecStrAnnot2_data/SecStrAPI/testing_$TODAY"
SCRIPT_DIR="/home/adam/Workspace/C#/SecStrAnnot2/scripts"

DOMAINS_FROM_PDBEAPI="/home/adam/Workspace/Python/Ubertemplate/Ubertemplate/domains_from_pdbeapi.py"
DOMAINS_TO_SECSTRAPI_FORMAT="$SCRIPT_DIR/domain_lists_to_SecStrAPI_format.py"
SELECT_BEST_DOMAINS="$SCRIPT_DIR/select_best_domain_per_uniprot.py"
SELECT_DOMAINS_BY_TAXON_GROUP="$SCRIPT_DIR/select_by_taxon_group.py"
SIMPLIFY_DOMAIN_LIST="$SCRIPT_DIR/simplify_domain_list.py"
DOWNLOAD_DOMAINS="$SCRIPT_DIR/download_from_pdbe.py"
SECSTRANNOTATOR_BATCH="$SCRIPT_DIR/SecStrAnnotator2_batch.py"
SECSTRANNOTATOR_DLL="/home/adam/Workspace/C#/SecStrAnnot2/bin/Debug/netcoreapp2.0/SecStrAnnotator2.dll"
COLLECT_ANNOTATIONS="$SCRIPT_DIR/collect_annotations.py"
DIVIDE_ANNOTATIONS="$SCRIPT_DIR/divide_annotations_by_pdb.py"
EXTRACT_SEQUENCES="$SCRIPT_DIR/extract_sequences.py"
ALIGN_SEQUENCES="$SCRIPT_DIR/align_sequences.py"
ADD_REFERENCE_RESIDUES="$SCRIPT_DIR/add_reference_residues.py"

TEMPLATE="2nnj,A,:"
TEMPLATE_ANNOTATION_FILE="$DATA_DIR/../2nnj-template-strict.sses.json"
SECSTRANNOTATOR_OPTIONS="--soft  --label2auth  --maxmetric 25,0.5,0.5  --verbose"  # Use with --verbose for CYP Anatomy analyses, without --verbose for SecStrAPI
ALIGNED_SSE_LABELS="A,B,C,D,E,H,I,J,K,L"
# ALIGNED_SSE_LABELS="all"
# ALIGNED_SSE_LABELS="A,B,C,D,E,F,G,H,I,J,K,L"
# ALIGNED_SSE_LABELS="A',B',B'',B''',F',G',J',K',K'',K''',L'"
# ALIGNED_SSE_LABELS="1.1,1.2,1.3,1.4,1.5,2.1,2.2,3.1,3.2,3.3,4.1,4.2,4.3"

#####################################################################################################################

echo '=== Get domains from CATH and Pfam ==='
mkdir $DATA_DIR
python3  $DOMAINS_FROM_PDBEAPI  --numbering label  --allow_null_domain_name  --join_domains_in_chain  1.10.630.10  >  $DATA_DIR/cyps_cath_$TODAY.simple.json 
# Downloading https://www.ebi.ac.uk/pdbe/api/mappings/1.10.630.10
# Found 1386 domains in 728 PDB entries.
python3  $DOMAINS_FROM_PDBEAPI  --numbering label  --allow_null_domain_name  --join_domains_in_chain  PF00067  >  $DATA_DIR/cyps_pfam_$TODAY.simple.json 
# Downloading https://www.ebi.ac.uk/pdbe/api/mappings/PF00067
# Found 1732 domains in 926 PDB entries.

echo '=== Merge domain lists and format them in SecStrAPI format (also adds UniProt refs) ==='
python3  $DOMAINS_TO_SECSTRAPI_FORMAT  --api_version $API_VERSION  \
    CATH  1.10.630.10  $DATA_DIR/cyps_cath_$TODAY.simple.json  \
    Pfam  PF00067  $DATA_DIR/cyps_pfam_$TODAY.simple.json  \
    >  $DATA_DIR/cyps_ALL_$TODAY.json

echo '=== Get NCBI taxons and groups (Euka/Bact/Arch/Viru) ==='
wget  ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz  -O $DATA_DIR/taxdump.tar.gz
tar  xzf  $DATA_DIR/taxdump.tar.gz  nodes.dmp  -O  > $DATA_DIR/ncbi_taxonomy_nodes.dmp
python3  $SCRIPT_DIR/get_taxids.py  $DATA_DIR/cyps_ALL_$TODAY.json  >  $DATA_DIR/domain_taxons.tsv
python3  $SCRIPT_DIR/classify_taxids.py  $DATA_DIR/ncbi_taxonomy_nodes.dmp  $DATA_DIR/domain_taxons.tsv  >  $DATA_DIR/domain_taxons_groups.tsv

echo '=== Select nonredundant set (with best quality) ==='
python3  $SELECT_BEST_DOMAINS  $DATA_DIR/cyps_ALL_$TODAY.json  >  $DATA_DIR/cyps_NR_$TODAY.json
# python3  $SELECT_DOMAINS_BY_TAXON_GROUP  $DATA_DIR/cyps_NR_$TODAY.json  $DATA_DIR/domain_taxons_groups.tsv  Bact  >  $DATA_DIR/cyps_NR_Bact_$TODAY.json
# python3  $SELECT_DOMAINS_BY_TAXON_GROUP  $DATA_DIR/cyps_NR_$TODAY.json  $DATA_DIR/domain_taxons_groups.tsv  Euka  >  $DATA_DIR/cyps_NR_Euka_$TODAY.json

echo '=== Simplify domain lists (for SecStrAnnotator) ==='
python3  $SIMPLIFY_DOMAIN_LIST  $DATA_DIR/cyps_ALL_$TODAY.json  >  $DATA_DIR/cyps_ALL_$TODAY.simple.json
python3  $SIMPLIFY_DOMAIN_LIST  $DATA_DIR/cyps_NR_$TODAY.json  >  $DATA_DIR/cyps_NR_$TODAY.simple.json
python3 extract_pdb_domain_list.py  $DATA_DIR/cyps_ALL_$TODAY.json  >  $DATA_DIR/AnnotationList.json

echo '=== Download CIF files ==='
python3  $DOWNLOAD_DOMAINS  $DATA_DIR/cyps_ALL_$TODAY.simple.json  $DATA_DIR/structures/  --format cif  --no_gzip  --cache $DATA_DIR/../cached_structures/
# Downloaded 916 PDB entries, failed to download 0 PDB entries

echo '=== Annotate ==='
cp  $TEMPLATE_ANNOTATION_FILE  $DATA_DIR/structures/${TEMPLATE:0:4}-template.sses.json
python3  $SECSTRANNOTATOR_BATCH  --dll $SECSTRANNOTATOR_DLL \
    --threads $N_THREADS  --options " $SECSTRANNOTATOR_OPTIONS " \
    $DATA_DIR/structures  $TEMPLATE  $DATA_DIR/cyps_ALL_$TODAY.simple.json 

echo '=== Collect annotations and put them to SecStrAPI format ==='
python3  $COLLECT_ANNOTATIONS  $DATA_DIR/cyps_ALL_$TODAY.json  $DATA_DIR/structures/  >  $DATA_DIR/annotations_ALL.json
python3  $COLLECT_ANNOTATIONS  $DATA_DIR/cyps_NR_$TODAY.json  $DATA_DIR/structures/  >  $DATA_DIR/annotations_NR.json
# python3  $COLLECT_ANNOTATIONS  $DATA_DIR/cyps_NR_Bact_$TODAY.json  $DATA_DIR/structures/  >  $DATA_DIR/annotations_NR_Bact.json
# python3  $COLLECT_ANNOTATIONS  $DATA_DIR/cyps_NR_Euka_$TODAY.json  $DATA_DIR/structures/  >  $DATA_DIR/annotations_NR_Euka.json
python3  $EXTRACT_SEQUENCES  $DATA_DIR/annotations_ALL.json  $DATA_DIR/sequences_ALL/
python3  $EXTRACT_SEQUENCES  $DATA_DIR/annotations_NR.json  $DATA_DIR/sequences_NR/
# python3  $EXTRACT_SEQUENCES  $DATA_DIR/annotations_NR_Bact.json  $DATA_DIR/sequences_NR_Bact/
# python3  $EXTRACT_SEQUENCES  $DATA_DIR/annotations_NR_Euka.json  $DATA_DIR/sequences_NR_Euka/

echo '=== Perform no-gap sequence alignment and create sequence logos (from Set-NR) ==='
python3  $ALIGN_SEQUENCES  $DATA_DIR/annotations_NR.json  --alignments_dir $DATA_DIR/aligments_NR/  --trees_dir $DATA_DIR/trees_NR/  --logos_dir $DATA_DIR/logos_NR/
# python3  $ALIGN_SEQUENCES  $DATA_DIR/annotations_NR_Bact.json  --alignments_dir $DATA_DIR/aligments_NR_Bact/  --trees_dir $DATA_DIR/trees_NR_Bact/  --logos_dir $DATA_DIR/logos_NR_Bact/
# python3  $ALIGN_SEQUENCES  $DATA_DIR/annotations_NR_Euka.json  --alignments_dir $DATA_DIR/aligments_NR_Euka/  --trees_dir $DATA_DIR/trees_NR_Euka/  --logos_dir $DATA_DIR/logos_NR_Euka/

echo '=== Realign sequences from Set-ALL to the alignment from Set-NR and add reference residue information ==='
python3  $ADD_REFERENCE_RESIDUES  $DATA_DIR/annotations_ALL.json  $DATA_DIR/aligments_NR/  --labels $ALIGNED_SSE_LABELS  --label2auth_dir $DATA_DIR/structures/  >  $DATA_DIR/annotations_with_reference_residues_ALL.json
python3  $ADD_REFERENCE_RESIDUES  $DATA_DIR/annotations_NR.json  $DATA_DIR/aligments_NR/  --labels $ALIGNED_SSE_LABELS  --label2auth_dir $DATA_DIR/structures/  >  $DATA_DIR/annotations_with_reference_residues_NR.json

echo '=== Divide annotations into per-PDB files ==='
python3  $DIVIDE_ANNOTATIONS  $DATA_DIR/annotations_with_reference_residues_ALL.json  $DATA_DIR/annotations_ALL/  --min_dir $DATA_DIR/annotations_ALL_min/
