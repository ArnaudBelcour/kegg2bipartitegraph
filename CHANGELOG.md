# Changelog

# Kegg2bipartitegraph v0.0.2 (2023-11-29)

## Add

* GitHub Actions to publish on PyPI.

## Fix

* Wrong URL in setup.py.

# Kegg2bipartitegraph v0.0.1 (2023-11-29)

First release of `kegg2bipartitegraph` with reference database created from KEGG version `108`.

## Add

* Creation of the subcommand `reference` to create the reference database.
* Creation of reconstruction subcommand from files: `reconstruct_from_esmecata` from esmecata results, `reconstruct_from_eggnog` from eggnog-mapper results and `reconstruct_from_kofamkoala` from Kofam-Koala.
* Creation of reconstruction subcommand from KEGG organisms `reconstruct_from_organism`.
* The goal of `kegg2bipartitegraph` is to create metabolic networks similar to the ones created in the article made by [Adèle Weber Zendrera et al. (2021)](https://www.nature.com/articles/s41598-021-91486-8). To this end, it uses an intern database created from KEGG that tries to apply the same rules that were used in this article (but at this moment it is not a complete match). The metabolic networks are created in sbml and graphml format. There are also files to show KEGG pathways and modules identified (they are also present in the sbml file).
* There is a modification of the scope code proposed by [Adèle Weber Zendrera et al. (2021)](https://www.nature.com/articles/s41598-021-91486-8) to include it in the tools but it is not used.
* Another idea of `kegg2bipartitegraph` is to create a reproducible environment meaning it contains an intern database that is used for the prediction. One a new version of KEGG is released, a new intern database is created and the older one is archived and can be used instead.