#  CONSTANTS  ##############################################################################

# SecStrAPI format details
DOMAINS_IN_DICT = False

# Field names in SecStrAPI format - general and domain specification
API_VERSION = 'api_version'
ANNOTATIONS = 'annotations'
PDB = 'pdb'
CHAIN = 'chain'
RANGES = 'ranges'
UNIPROT_ID = 'uniprot_id'
UNIPROT_NAME = 'uniprot_name'
MAPPINGS = 'domain_mappings'
DOMAIN_NAME = 'domain'
SOURCE = 'source'
FAMILY_ID = 'family'

# Field names in SecStrAPI format - SSE annotation (from SecStrAnnotator format)
SSES = 'secondary_structure_elements'
CONNECTIVITY = 'beta_connectivity'
COMMENT = 'comment'
LABEL = 'label'
SEQUENCE = 'sequence'
START = 'start'
END = 'end'
PIVOT_RESIDUE = 'pivot_residue'