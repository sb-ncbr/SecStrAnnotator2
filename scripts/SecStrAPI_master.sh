DATA_DIR="/home/adam/Workspace/C#/SecStrAnnot2_data/SecStrAPI/testing_20190314"
# STRUCTURES_DIR="$DATA_DIR/structures"
# ANNOTATIONS_DIR="$DATA_DIR/annotations"
# BEST_ANNOTATIONS_DIR="$DATA_DIR/annotations_best"
# SEQUENCES_DIR="$DATA_DIR/sequences"
# BEST_SEQUENCES_DIR="$DATA_DIR/sequences_best"
# ALIGNMENTS_DIR="$DATA_DIR/aligments"

SCRIPT_DIR="/home/adam/Workspace/C#/SecStrAnnot2/scripts"
DOMAINS_FROM_PDBEAPI="/home/adam/Workspace/Python/Ubertemplate/Ubertemplate/domains_from_pdbeapi.py"
DOMAINS_TO_SECSTRAPI_FORMAT="$SCRIPT_DIR/domain_lists_to_SecStrAPI_format.py"
SELECT_BEST_DOMAINS="$SCRIPT_DIR/select_best_domain_per_uniprot.py"
SIMPLIFY_DOMAIN_LIST="$SCRIPT_DIR/simplify_domain_list.py"
DOWNLOAD_DOMAINS="$SCRIPT_DIR/download_from_pdbe.py"
SECSTRANNOTATOR_BATCH="$SCRIPT_DIR/SecStrAnnotator2_batch.py"
SECSTRANNOTATOR_DLL="/home/adam/Workspace/C#/SecStrAnnot2/bin/Debug/netcoreapp2.0/SecStrAnnot2.dll"
COLLECT_ANNOTATIONS="$SCRIPT_DIR/collect_annotations.py"
ALIGN_SEQUENCES="$SCRIPT_DIR/align_sequences.py"

TEMPLATE="1og2,A,:"
TEMPLATE_ANNOTATION_FILE="$DATA_DIR/1og2-template.sses.json"
SECSTRANNOTATOR_OPTIONS="--soft"

TODAY="20190314"
API_VERSION="1.0"
N_THREADS="8"


# Get domains from CATH and Pfam
python3  $DOMAINS_FROM_PDBEAPI  --numbering label  --allow_null_domain_name  --join_domains_in_chain  1.10.630.10  >  $DATA_DIR/cyps_cath_$TODAY.simple.json 
# Downloading https://www.ebi.ac.uk/pdbe/api/mappings/1.10.630.10
# Found 728 PDB entries.
python3  $DOMAINS_FROM_PDBEAPI  --numbering label  --allow_null_domain_name  --join_domains_in_chain  PF00067  >  $DATA_DIR/cyps_pfam_$TODAY.simple.json 
# Downloading https://www.ebi.ac.uk/pdbe/api/mappings/PF00067
# Found 883 PDB entries.

# Merge domain lists and format them in SecStrAPI format
python3  $DOMAINS_TO_SECSTRAPI_FORMAT  --api_version 1.0  \
    CATH  1.10.630.10  $DATA_DIR/cyps_cath_$TODAY.simple.json  \
    Pfam  PF00067  $DATA_DIR/cyps_pfam_$TODAY.simple.json  \
    >  $DATA_DIR/cyps_all_$TODAY.json

# Select nonredundant set (with best quality)
python3  $SELECT_BEST_DOMAINS  $DATA_DIR/cyps_all_$TODAY.json  >  $DATA_DIR/cyps_best_$TODAY.json

# Simplify domain lists (for SecStrAnnotator)
python3  $SIMPLIFY_DOMAIN_LIST  $DATA_DIR/cyps_all_$TODAY.json  >  $DATA_DIR/cyps_all_$TODAY.simple.json
python3  $SIMPLIFY_DOMAIN_LIST  $DATA_DIR/cyps_best_$TODAY.json  >  $DATA_DIR/cyps_best_$TODAY.simple.json

# Download CIF files
python3  $DOWNLOAD_DOMAINS  --format cif  --no_gzip  $DATA_DIR/cyps_all_$TODAY.simple.json  $DATA_DIR/structures/

# Annotate
cp  $TEMPLATE_ANNOTATION_FILE  $DATA_DIR/structures/
python3  $SECSTRANNOTATOR_BATCH  --dll $SECSTRANNOTATOR_DLL \
    --threads $N_THREADS  --options " $SECSTRANNOTATOR_OPTIONS " \
    $DATA_DIR/structures  $TEMPLATE  $DATA_DIR/cyps_all_$TODAY.simple.json 

# Collect annotations and put them to SecStrAPI format
python3  $COLLECT_ANNOTATIONS  $DATA_DIR/cyps_all_$TODAY.json  $DATA_DIR/structures/  $DATA_DIR/annotations/  --sequences_directory $DATA_DIR/sequences/
python3  $COLLECT_ANNOTATIONS  $DATA_DIR/cyps_best_$TODAY.json  $DATA_DIR/structures/  $DATA_DIR/annotations_best/  --sequences_directory $DATA_DIR/sequences_best/

# Perform no-gap sequence alignment and create sequence logos
python3  $ALIGN_SEQUENCES  $DATA_DIR/annotations_best/all.json  --alignments_dir $DATA_DIR/aligments_best/  --trees_dir $DATA_DIR/trees_best/  --logos_dir $DATA_DIR/logos_best/


# #TODO make separate script for extracting sequences
# #TODO perform no-gap alignment, realign and add pivot residues to the annotations for SecStrAPI (find out how to obtain their auth_ numbers!)
