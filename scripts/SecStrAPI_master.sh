TODAY="20190607"
API_VERSION="1.0"
N_THREADS="8"

DATA_DIR="/home/adam/Workspace/C#/SecStrAnnot2_data/SecStrAPI/testing_$TODAY"
SCRIPT_DIR="/home/adam/Workspace/C#/SecStrAnnot2/scripts"

DOMAINS_FROM_PDBEAPI="/home/adam/Workspace/Python/Ubertemplate/Ubertemplate/domains_from_pdbeapi.py"
DOMAINS_TO_SECSTRAPI_FORMAT="$SCRIPT_DIR/domain_lists_to_SecStrAPI_format.py"
SELECT_BEST_DOMAINS="$SCRIPT_DIR/select_best_domain_per_uniprot.py"
SIMPLIFY_DOMAIN_LIST="$SCRIPT_DIR/simplify_domain_list.py"
DOWNLOAD_DOMAINS="$SCRIPT_DIR/download_from_pdbe.py"
SECSTRANNOTATOR_BATCH="$SCRIPT_DIR/SecStrAnnotator2_batch.py"
SECSTRANNOTATOR_DLL="/home/adam/Workspace/C#/SecStrAnnot2/bin/Debug/netcoreapp2.0/SecStrAnnot2.dll"
COLLECT_ANNOTATIONS="$SCRIPT_DIR/collect_annotations.py"
DIVIDE_ANNOTATIONS="$SCRIPT_DIR/divide_annotations_by_pdb.py"
EXTRACT_SEQUENCES="$SCRIPT_DIR/extract_sequences.py"
ALIGN_SEQUENCES="$SCRIPT_DIR/align_sequences.py"
ADD_PIVOT_RESIDUES="$SCRIPT_DIR/add_pivot_residues.py"

TEMPLATE="1og2,A,:"
TEMPLATE_ANNOTATION_FILE="$DATA_DIR/../1og2-template.sses.json"
SECSTRANNOTATOR_OPTIONS="--soft --label2auth"
ALIGNED_SSE_LABELS="A,B,C,D,E,H,I,J,K,L"
# ALIGNED_SSE_LABELS="all"
# ALIGNED_SSE_LABELS="A,B,C,D,E,F,G,H,I,J,K,L"
# ALIGNED_SSE_LABELS="A',B',B'',B''',F',G',J',K',K'',K''',L'"
# ALIGNED_SSE_LABELS="1.1,1.2,1.3,1.4,1.5,2.1,2.2,3.1,3.2,3.3,4.1,4.2,4.3"

#####################################################################################################################

# # Get domains from CATH and Pfam
# mkdir $DATA_DIR
# python3  $DOMAINS_FROM_PDBEAPI  --numbering label  --allow_null_domain_name  --join_domains_in_chain  1.10.630.10  >  $DATA_DIR/cyps_cath_$TODAY.simple.json 
# # Downloading https://www.ebi.ac.uk/pdbe/api/mappings/1.10.630.10
# # Found 728 PDB entries.
# python3  $DOMAINS_FROM_PDBEAPI  --numbering label  --allow_null_domain_name  --join_domains_in_chain  PF00067  >  $DATA_DIR/cyps_pfam_$TODAY.simple.json 
# # Downloading https://www.ebi.ac.uk/pdbe/api/mappings/PF00067
# # Found 899 PDB entries.

# # Merge domain lists and format them in SecStrAPI format (also adds UniProt refs)
# python3  $DOMAINS_TO_SECSTRAPI_FORMAT  --api_version 1.0  \
#     CATH  1.10.630.10  $DATA_DIR/cyps_cath_$TODAY.simple.json  \
#     Pfam  PF00067  $DATA_DIR/cyps_pfam_$TODAY.simple.json  \
#     >  $DATA_DIR/cyps_all_$TODAY.json

# # Get NCBI taxons and groups (Euka/Bact/Arch/Viru)
# wget  ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz  -O $DATA_DIR/taxdump.tar.gz
# tar  xzf  $DATA_DIR/taxdump.tar.gz  nodes.dmp  -O  > $DATA_DIR/ncbi_taxonomy_nodes.dmp
# python3  $SCRIPT_DIR/get_taxids.py  $DATA_DIR/cyps_all_$TODAY.json  >  $DATA_DIR/domain_taxons.tsv
# python3  $SCRIPT_DIR/classify_taxids.py  $DATA_DIR/ncbi_taxonomy_nodes.dmp  $DATA_DIR/domain_taxons.tsv  >  $DATA_DIR/domain_taxons_groups.tsv

# # Select nonredundant set (with best quality)
# python3  $SELECT_BEST_DOMAINS  $DATA_DIR/cyps_all_$TODAY.json  >  $DATA_DIR/cyps_best_$TODAY.json

# # Simplify domain lists (for SecStrAnnotator)
# python3  $SIMPLIFY_DOMAIN_LIST  $DATA_DIR/cyps_all_$TODAY.json  >  $DATA_DIR/cyps_all_$TODAY.simple.json
# python3  $SIMPLIFY_DOMAIN_LIST  $DATA_DIR/cyps_best_$TODAY.json  >  $DATA_DIR/cyps_best_$TODAY.simple.json

# # Download CIF files
# python3  $DOWNLOAD_DOMAINS  $DATA_DIR/cyps_all_$TODAY.simple.json  $DATA_DIR/structures/  --format cif  --no_gzip  --cache $DATA_DIR/../cached_structures/
# # Downloaded 916 PDB entries, failed to download 0 PDB entries

# Annotate
cp  $TEMPLATE_ANNOTATION_FILE  $DATA_DIR/structures/
python3  $SECSTRANNOTATOR_BATCH  --dll $SECSTRANNOTATOR_DLL \
    --threads $N_THREADS  --options " $SECSTRANNOTATOR_OPTIONS " \
    $DATA_DIR/structures  $TEMPLATE  $DATA_DIR/cyps_all_$TODAY.simple.json 

# Collect annotations and put them to SecStrAPI format
python3  $COLLECT_ANNOTATIONS  $DATA_DIR/cyps_all_$TODAY.json  $DATA_DIR/structures/  >  $DATA_DIR/annotations_all.json
python3  $COLLECT_ANNOTATIONS  $DATA_DIR/cyps_best_$TODAY.json  $DATA_DIR/structures/  >  $DATA_DIR/annotations_best.json
python3  $EXTRACT_SEQUENCES  $DATA_DIR/annotations_all.json  $DATA_DIR/sequences_all/
python3  $EXTRACT_SEQUENCES  $DATA_DIR/annotations_best.json  $DATA_DIR/sequences_best/

# Perform no-gap sequence alignment and create sequence logos (from Set-NR)
python3  $ALIGN_SEQUENCES  $DATA_DIR/annotations_best.json  --alignments_dir $DATA_DIR/aligments_best/  --trees_dir $DATA_DIR/trees_best/  --logos_dir $DATA_DIR/logos_best/

# Realign sequences from Set-ALL to the alignment from Set-NR and add pivot residue information
python3  $ADD_PIVOT_RESIDUES  $DATA_DIR/annotations_all.json  $DATA_DIR/aligments_best/  --labels $ALIGNED_SSE_LABELS  --label2auth_dir $DATA_DIR/structures/  >  $DATA_DIR/annotations_with_pivots_all.json
python3  $ADD_PIVOT_RESIDUES  $DATA_DIR/annotations_best.json  $DATA_DIR/aligments_best/  --labels $ALIGNED_SSE_LABELS  --label2auth_dir $DATA_DIR/structures/  >  $DATA_DIR/annotations_with_pivots_best.json

# Divide annotations into per-PDB files
python3  $DIVIDE_ANNOTATIONS  $DATA_DIR/annotations_with_pivots_all.json  $DATA_DIR/annotations_all/
