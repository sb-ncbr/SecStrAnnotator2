# SecStrAnnotator

**SecStrAnnotator Suite** is a collection of tools for annotation of secondary structure elements (SSEs), namely helices and beta-strands, in protein structures. Its core part is **SecStrAnnotator** – a C#/.NET Core program implementing the annotation algorithm. Additional scripts in Python and R provide related functionality, such as batch processing, analysis and visualization of the annotation results, format conversions etc.

More information is available at <https://webchem.ncbr.muni.cz/Wiki/SecStrAnnotator>.

If you found this software helpful, please cite:

- Midlik A, Hutařová Vařeková I, Hutař J, Moturu TR, Navrátilová V, Koča J, Berka K, Svobodová Vařeková R (2019) Automated family-wide annotation of secondary structure elements. In: Kister AE (ed) Protein supersecondary structures: Methods and protocols, Humana Press. ISBN 978-1-4939-9160-0.

## Requirements

For SecStrAnnotator:

- .NET runtime 6.0 (for SecStrAnnotator<=2.2: .NET Core runtime 3.0 or 3.1)
- PyMOL
- SecStrAnnotator_config.json – configuration file (modify according to your system)

For additional scripts:

- Python3 (>=3.6)
- Python3 packages listed in `scripts/requirements.txt` (install by `pip3 install -r scripts/requirements.txt`)
- R, rstudio
- R packages `dplyr`, `tidyr`, `readr`, `ggplot2`

## Main parts

### `SecStrAnnotator.dll`

The core program, annotates a query protein domain (e.g. chain A in 1tqn) according to provided template annotation (e.g. chain A in 2nnj). Can also be used without a template to detect secondary structure elements without annotating (`--onlyssa`).

Keep in mind that the chains and residues are numbered according to the label_* numbering scheme in mmCIF file format (i.e. chain identifier is `label_asym_id`, residue number is `label_seq_id`).

Before running SecStrAnnotator, make sure that `SecStrAnnotator_config.json` is set correctly (most importantly, that `"PymolExecutable"` points to the PyMOL executable in your system, such as `pymol` or `C:/ProgramData/PyMOL/Scripts/pymol.exe`).

Example usage:

    dotnet SecStrAnnotator.dll --help
    dotnet SecStrAnnotator.dll --onlyssa examples/ 1tqn,A  # Detect
    dotnet SecStrAnnotator.dll examples/ 2nnj,A 1tqn,A     # Detect and annotate

### `scripts/SecStrAnnotator_batch.py`

Runs `SecStrAnnotator.dll` on multiple query protein domains in one batch.

Example usage:

    python3 SecStrAnnotator_batch.py --help
    python3 SecStrAnnotator_batch.py examples/ 2nnj,A examples/cyp_family_sample.json --threads 4

### `scripts/SecStrAPI_pipeline.py`

A pipeline for annotating the whole protein family and preparing data for analysis and for publishing via SecStrAPI (includes downloading the list of family members defined by CATH and Pfam, downloading their structures, selecting a non-redundant set Set-NR, annotation, multiple sequence alignment for individual SSEs, formatting into SecStrAPI format, formatting into TSV format for further analyses etc.)

Example usage:

    python3 scripts/SecStrAPI_pipeline.py scripts/SecStrAPI_pipeline_settings.json --resume

Before running, modify the settings in `SecStrAPI_pipeline_settings.json` to set the following:

- `"data_dir"` - directory for results
- `"cath_family_id"` - family identifier in CATH database (`null` = don't use CATH)
- `"pfam_family_id"` - family identifier in Pfam database (`null` = don't use Pfam)
- `"api_version"` - SecStrAPI version information to put into annotation files for SecStrAPI
- `"structure_cache_dir"` - directory with pre-downloaded structures (named \<PDB_ID\>.cif); structures not found there will be downloaded (`null` = don't use any cache directory, download all the structures)
- `"template_domain"` - specification of the template domain in format `"<PDB_ID>,<CHAIN>,<RESIDUE_RANGES>"`
- `"template_annotation_file"` - file with prepared annotation of the template domain in SecStrAnnotator format (must exist)
- `"template_structure_file"` - file with template domain structure (must exist)
- `"n_threads"` - number of CPU threads to run SecStrAnnotator in
- `"secstrannotator_dll"` - path to SecStrAnnotator executable (DLL), absolute or relative to `scripts/` directory (`null` = default `../SecStrAnnotator.dll`)
- `"secstrannotator_options"` - any options to pass to SecStrAnnotator
- `"sses_for_generic_numbering"` - comma-separated list of SSEs for which a generic numbering scheme should be established (`"all"` = establish for all SSEs)

Each completed step of the pipeline is logged into `checkpoints.txt` in the directory with results. If any step fails, the pipeline can be resumed from the last completed step by using option `--resume`.

Individual steps of the pipeline can be skipped by commenting out the respective lines in the `main` function in `SecStrAPI_pipeline.py`. Individual steps can also be performed by directly executing the scripts in `scripts/secstrapi_data_preparation`.

### `scripts/R_sec_str_anatomy_analysis/`

Statistical analysis of the annotation results on the whole protein family.

Example usage:

- Launch `rstudio` from directory `scripts/R_sec_str_anatomy_analysis/`
- Set the path to your annotation data (`DATADIR`) in `sec_str_anatomy.R`
- Modify the family-specific settings in `sec_str_anatomy_settings.R`
- Run `sec_str_anatomy.R` line by line

### `scripts/secstrapi_plugin/`

PyMOL plugin for visualization of secondary structure annotations directly on the 3D structure. By default the annotations are downloaded from SecStrAPI (<https://webchem.ncbr.muni.cz/API/SecStr>).

Example usage:

- Launch PyMOL
- Install the plugin via Plugin Manager or run `run .../secstrapi_plugin.py`
- Run `annotate_sec_str 1tqnA`
