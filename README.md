# bioiain
[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/Thecnopapa/bioiain)

By Iain Visa (iainvisa@gmail.com)


![Bioiain logo](src/bioiain/bioiain_logo.png)

Toolbox for structural analysis of proteins.

> [!WARNING]
> **WIP** EVERYTHING IS UNDER DEVELOPMENT!

Many features are still not commented/documented or even mentioned. Feel free to explore and/or use any functions.

Can be downloaded from the [PiPy repository](https://pypi.org/project/bioiain/) but note the used version as any function might change during development.

`pip install bioiain`

`pip install bioiain[ml]` for machine learning dependencies (pytorch).

`pip isntall bioiain[aleph]` to use ALEPH dependant functions.

> [!NOTE]
> `pip install bioiain[ml,aleph]` will install all dependencies at once.

> [!WARNING]
> Some dependencies must be installed through conda and are not automatically installed.
> Although not needed for most modules, the can be installed through the following commands as needed:
> 
> `conda install -c conda-forge pymol-open-source` (>=2.5.0 tested)
> 
> `conda install -c bioconda clustalw` (==2.1 tested)
> 
> `conda install -c bioconda mmseqs2`
> 
> `conda install -c bioconda foldseek` (>=8 tested)


If you were to use this and find any issue I'll be happy to fix it :D


# INFO

Module wiki can be found at [wiki.bioiain.com](https://wiki.bioiain.com) .

Relevant python code can be found in the `src` folder within their respective folders.
The `test` folder is for active development, and it's contents are **constantly** deleted/modified, and are not included in the package.

> [!NOTE]
> Preset workflows are being developed, including the [projectDimer](https://gitHub.com/thecnopapa/projectdimer) workflow. (WIP)


**Protein Framework**

Originally based on Biopython's hierarchy, but no longer dependant on it. Classes for structures, chains, residues, and atoms are included for manipulation and analysis of protein models. 
 
Unlike Biopython, residues and atoms do not share the base entity framework as they behave in significantly different ways. Also respective classes for nucleotides, ligands and water are included and under development.

Bioiain is designed to be expandable, and custom classes are encouraged to match your purpose. Two custom classes (`FragmentedStructure` and `Fragment`) are included within the `aleph` module to deal with protein fragments.

Includes general-purpose tools and pipelines for importing, processing, saving, and exporting structures in mmCIF format (but PDB exports are slightly supported)

> [!IMPORTANT]
> PDB Parsing is not supported yet, so the input so far must be mmCIF

>[!NOTE]
> WIP: Allow PDB parsing, dealing with structures with several models, cast data from respective Biopython objects.


**Symmetries**

This framework is designed to work with all the information available in crystallographic structures, therefore symmetry is considered by default when available.


**Machine Learning**

Still at a very early stage, Bioiain includes a PyTorch-based ML framework to simplify the development and training of ML models using on structural data.

It is intended for the development of new model architectures, embeddings, protein languages, classification, and other ML based approaches to structural biology.

> [!NOTE]
> It will probably not outperform any other library, but you won't have to worry about file organisation, logging, or blowing your RAM by loading all your data into memeory.

This includes a `BaseModel` class with all the utilities commonly used during train/test/eval/inference of models. This allows setting up a model with custom layers and/or sub-models without dealing with losses, batching, logging, saving, loading, or exporting.

> [!NOTE]
> This is designed with medium complexity models in mind, such as Variational Auto-Encoders, therefore might be over-complicated for a simple classifier, and I have not tried to set up a transformer, but could be complicated.

This framework works best with the provided `EmbeddingDataset` class, which doubles as a dataset and dataloader. And can be easily set  up wilenusing (or not) the provided Protein Framework.

The framework relies on the provided `Embedding` class and subclasses. Which allow the generation and integration of custom-made embeddings, at protein or residue level.

This includes (local) integrated logging using Tensorboard.

>[!NOTE]
> WIP: Still on the first iteration, will be reworked at some point.

> [!IMPORTANT]
> Integration with external datasets is still not dealt with and will likely require some processing.


**ALEPH**

Characteristic vectors are a powerful abstraction of protein structure, and can be calculated with [ALEPH](https://doi.org/10.1107/S2059798320001679), through direct integration within the FragmentedStructure and Fragment classes included in the Protein Framework.


**Tools**

Utility functions to use and parse some external tools are included. For now this includes:
 - [DSSP](https://doi.org/10.1002/bip.360221211)
 - [PISA](https://doi.org/10.1016/j.jmb.2007.05.022)


**Utilities**

Additionally, a large set of utilities is included, from logging, to common mathematical operations.

The `utilities.files.StructureDataset` facilitates working with large numbers of protein structures, and can fetch your structures of interest straight from the PDB.

**Visualization**

For structural visualisation, a custom PyMol scripting framework is included, replacing heavy sessions with generative commands.

Some common matplotlib utilities are also included. Including 3D visualisation of structures.

**Databases**

Utility functions download, parse, and query some online databases are included. For now this includes:
 - [Plinder](https://plinder-org.github.io/) (protein ligand interactions dataset)

> [!NOTE]
> Currently under development (separate repo): [UniProt](https://www.uniprot.org/) and [COSMIC](https://cancer.sanger.ac.uk/cosmic/) databases

