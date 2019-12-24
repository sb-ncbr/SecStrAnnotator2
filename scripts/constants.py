#  CONSTANTS  ##############################################################################

# SecStrAPI format details
DOMAINS_IN_DICT = True

# Field names in SecStrAPI format - general and domain specification
API_VERSION = 'api_version'
ANNOTATIONS = 'annotations'
PDB = 'pdb'
CHAIN = 'chain_id'
AUTH_CHAIN = 'auth_chain_id'
RANGES = 'ranges'
AUTH_RANGES = 'auth_ranges'
UNIPROT_ID = 'uniprot_id'
UNIPROT_NAME = 'uniprot_name'
MAPPINGS = 'domain_mappings'
DOMAIN_NAME = 'domain'
SOURCE = 'source'
FAMILY_ID = 'family'

# Field names in SecStrAPI format - SSE annotation (from SecStrAnnotator format)
SSES = 'secondary_structure_elements'
NESTED_SSES = 'nested_sses'
CONNECTIVITY = 'beta_connectivity'

LABEL = 'label'
SEQUENCE = 'sequence'
COMMENT = 'comment'

CHAIN_ID = 'chain_id'
START = 'start'
END = 'end'
TYPE = 'type'

AUTH_CHAIN_ID = 'auth_chain_id'
AUTH_START = 'auth_start'
AUTH_START_INS = 'auth_start_ins_code'
AUTH_END = 'auth_end'
AUTH_END_INS = 'auth_end_ins_code'

PIVOT_RESIDUE = 'reference_residue'
AUTH_PIVOT_RESIDUE = 'auth_reference_residue'
AUTH_PIVOT_RESIDUE_INS_CODE = 'auth_reference_residue_ins_code'