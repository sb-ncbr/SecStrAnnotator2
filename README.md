# SecStrAnnotator

**SecStrAnnotator Suite** is a collection of tools for annotation of secondary structure elements (SSEs), namely helices and beta-strands, in protein structures. Its core part is **SecStrAnnotator** – a C#/.NET Core program implementing the annotation algorithm. Additional scripts in Python and R provide related functionality, such as batch processing, analysis and visualization of the annotation results, format conversions etc.

More information is available at <https://webchem.ncbr.muni.cz/Wiki/SecStrAnnotator>.

If you found this software helpful, please cite:

- Midlik A, Hutařová Vařeková I, Hutař J, Moturu TR, Navrátilová V, Koča J, Berka K, Svobodová Vařeková R (2019) Automated family-wide annotation of secondary structure elements. In: Kister AE (ed) Protein supersecondary structures: Methods and protocols, Humana Press. ISBN 978-1-4939-9160-0.

## Requirements

For SecStrAnnotator:

- .NET Core Runtime 3.0 or 3.1
- PyMOL
- SecStrAnnotator_config.json – configuration file (modify according to your system)

For additional scripts:

- Python3 (>=3.6)
- bash
- R, rstudio

## Main parts

### `SecStrAnnotator.dll`

The core program, annotates a query protein domain (e.g. chain A in 1tqn) according to provided template annotation (e.g. chain A in 2nnj). Can also be used without a template to detect secondary structure elements without annotating (`--onlyssa`).

Keep in mind that the chains and residues are numbered according to the label_* numbering scheme in mmCIF file format (i.e. chain identifier is `label_asym_id`, residue number is `label_seq_id`).

Example usage:

    dotnet SecStrAnnotator.dll --help
    dotnet SecStrAnnotator.dll --onlyssa examples/ 1tqn,A  # Detect
    dotnet SecStrAnnotator.dll examples/ 2nnj,A 1tqn,A     # Detect and annotate

### `scripts/SecStrAnnotator_batch.py`

Runs `SecStrAnnotator.dll` on multiple query protein domains in one batch.

Example usage:

    python3 SecStrAnnotator_batch.py --help
    python3 SecStrAnnotator_batch.py examples/ 2nnj,A examples/cyp_family_sample.json --threads 4

### `scripts/secstrapi_data_preparation/`

A pipeline for annotating the whole protein family, including: downloading the list of family members defined by CATH and Pfam, downloading their structures, selecting a non-redundant set, annotation, multiple sequence alignment for individual SSEs, formatting into SecStrAPI format, formatting into TSV format for further analyses.

Example usage:

    bash scripts/secstrapi_data_preparation/SecStrAPI_master.sh

Before running, modify the SETTINGS section of the `SecStrAPI_master.sh` to set your family of interest, data directory, options, paths to your template annotation (TEMPLATE_ANNOTATION_FILE, TEMPLATE_STRUCTURE_FILE). Unwanted steps of the pipeline can be commented out in the MAIN PIPELINE section.

### `scripts/R_sec_str_anatomy_analysis/`

Statistical analysis of the annotation results on the whole protein family.

Example usage:

- Launch `rstudio` from this directory
- Set the path to your annotation data (DATADIR) in `sec_str_anatomy.R`
- Modify the family-specific settings in `sec_str_anatomy_settings.R`
- Run `sec_str_anatomy.R` line by line

### `scripts/secstrapi_plugin/`

PyMOL plugin for visualization of secondary structure annotations directly on the 3D structure. By default the annotations are downloaded from SecStrAPI (<https://webchem.ncbr.muni.cz/API/SecStr>).

Example usage:

- Launch PyMOL
- Install the plugin via Plugin Manager or run `run .../secstrapi_plugin.py`
- Run `annotate_sec_str 1tqnA`
