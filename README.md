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

If you were to use this and find any issue I'll be happy to fix it :D


# INFO

Relevant python code can be found in the `src` folder within their relevant folders.
The `test` is for development use, and it's contents will be **probably** deleted/modified at some point, and are not included in the package.

> [!NOTE]
> Preset workflows are being developed, including the [projectDimer](Https://GitHub.com/Thecnopapa/projectDimer) workflow. (WIP)


**Protein Framework**

Originally based on Biopython's hierarchy, but no longer dependant on it. Classes for structures and chains are included for manipulation and analysis of protein models. Designed to be expandable, custom classes are encouraged to match each purpose.
 
Unlike Biopython, residues and atoms do not share the base entity framework as they behave in significantly different ways. Also respective classes for nucleotides, ligands and water are included.

Includes general-purpose tools and pipelines for importing, processing,saving, and exporting structures in mmCIF format (but PDB is still slightly supported)

> [!IMPORTANT]
> PDB Parsing is not supported yet, so the input so far must be mmCIF

>[!NOTE]
> WIP: Allow PDB parsing, dealing with structures with several models, cast data from respective Biopython objects.


**Symmetries**

This framework is designed to work with all the information available in crystallographic structures, therefore symmetry is considered when available.


**Machine Learning**

Still at a very early stage, Bioiain includes a PyTorch-based ML framework to simplify the development and training of ML models, focused on structural data.

This includes a base model with all the utilities commonly used during train/test/eval/inference of models.

Also a dataset/embedding framework is also set up with integrations with the Protein Framework.

This includes integrated logging using Tensorboard.


**ALEPH**

Characteristic vectors are a powerful abstraction of protein structure, and can be calculated with [ALEPH](https://doi.org/10.1107/S2059798320001679), through direct integration within the FragmentedStructure and Fragment classes included in the Protein Framework.


**Tools**

Utility functions to use and parse some external tools are included. For now this includes:
 - DSSP
 - PISA


**Utilities**

Additionally, a large set of utilities is included, from logging, to common mathematical operations.


**Visualization**

For structural visualisation, a custom PyMol scripting framework is included, replacing heavy sessions with generative commands.

Some common matplotlib utilities are also included.

**Databases**

Utility functions download, parse, and query some online databases are included. For now this includes:
 - [Plinder](https://plinder-org.github.io/) (protein ligand interactions dataset)

> [!NOTE]
> Currently under development (separate repo): [UniProt](https://www.uniprot.org/) and [COSMIC](https://cancer.sanger.ac.uk/cosmic/) databases

