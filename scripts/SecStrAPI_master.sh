DATA_DIR='/home/adam/Workspace/C#/SecStrAnnot2_data/SecStrAPI'
STRUCTURES_DIR=$DATA_DIR/structures
ANNOTATIONS_DIR=$DATA_DIR/annotations
SEQUENCES_DIR=$DATA_DIR/sequences
ALIGNMENTS_DIR=$DATA_DIR/aligments
DOMAINS_FROM_PDBEAPI='/home/adam/Workspace/Python/Ubertemplate/Ubertemplate/domains_from_pdbeapi.py'
MERGE_DOMAINS='/home/adam/Workspace/C#/SecStrAnnot2/scripts/merge_domain_lists.py'
DOWNLOAD_DOMAINS='/home/adam/Workspace/C#/SecStrAnnot2/scripts/download_from_pdbe.py'
SECSTRANNOTATOR_BATCH='/home/adam/Workspace/C#/SecStrAnnot2/scripts/SecStrAnnotator2_batch.py'
SECSTRANNOTATOR_DLL='/home/adam/Workspace/C#/SecStrAnnot2/bin/Debug/netcoreapp2.0/SecStrAnnot2.dll'
COLLECT_ANNOTATIONS='/home/adam/Workspace/C#/SecStrAnnot2/scripts/collect_annotations.py'
TEMPLATE='1og2,A,:'
TEMPLATE_ANNOTATION_FILE=$DATA_DIR/1og2-template.sses.json

TODAY='20190313'
API_VERSION='1.0'
N_THREADS='8'


# # Get domains from CATH and Pfam
# python3 $DOMAINS_FROM_PDBEAPI --numbering label --allow_null_domain_name --join_domains_in_chain 1.10.630.10  >  $DATA_DIR/cyps_cath_$TODAY.json 
# # Downloading https://www.ebi.ac.uk/pdbe/api/mappings/1.10.630.10
# # Found 728 PDB entries.
# python3 $DOMAINS_FROM_PDBEAPI --numbering label --allow_null_domain_name --join_domains_in_chain PF00067  >  $DATA_DIR/cyps_pfam_$TODAY.json 
# # Downloading https://www.ebi.ac.uk/pdbe/api/mappings/PF00067
# # Found 883 PDB entries.

# # Merge domain lists
# python3 $MERGE_DOMAINS  --api_version 1.0 \
#     CATH 1.10.630.10 $DATA_DIR/cyps_cath_$TODAY.json  \
#     Pfam PF00067 $DATA_DIR/cyps_pfam_$TODAY.json  \
#     --simple_output $DATA_DIR/cyps_merged_simple_$TODAY.json \
#     >  $DATA_DIR/cyps_merged_$TODAY.json

# # Download CIF files
# python3 $DOWNLOAD_DOMAINS  $DATA_DIR/cyps_merged_simple_$TODAY.json  $STRUCTURES_DIR

# # Annotate
# cp $TEMPLATE_ANNOTATION_FILE $STRUCTURES_DIR/
# python3 $SECSTRANNOTATOR_BATCH --dll $SECSTRANNOTATOR_DLL \
#     $STRUCTURES_DIR  $TEMPLATE  $DATA_DIR/cyps_merged_simple_$TODAY.json \
#     --threads $N_THREADS  --options ''

# Collect annotations and put them to SecStrAPI format
python3 $COLLECT_ANNOTATIONS  $DATA_DIR/cyps_merged_$TODAY.json  $STRUCTURES_DIR  $ANNOTATIONS_DIR  --sequences_directory $SEQUENCES_DIR

#TODO select Set-NR, perform no-gap alignment, realign and add pivot residues to the annotations for SecStrAPI (find out how to obtain their auth_ numbers!)
