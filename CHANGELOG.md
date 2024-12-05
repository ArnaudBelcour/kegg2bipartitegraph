# Changelog

# Kegg2bipartitegraph v0.2.2 (2024-12-05)

## Fix

* issue in pyproject.toml leading to missing files in archive.

# Kegg2bipartitegraph v0.2.1 (2024-12-05)

## Add

* `module_class.tsv` creation to `k2bg organism`.

# Kegg2bipartitegraph v0.2.0 (2024-12-05)

## Add

* `scope` subcommand to compute scope using graphml files and seed files.
* Script (`seeds.py`) to create several seed files used in the [article of Webber-Zendrera et al. (2021)](https://doi.org/10.1038/s41598-021-91486-8) when creating reference model. Add the seed files to the archive.
* Warning to k2bg organism if there is a difference between the version of KEGG in the package and the version of KEGG online.

## Fix

* Issue with GitHub Actions.

## Modify

* Update reference model to KEGG 112.
* Reference model now contains reactions to remove (either glycan or ubiquitous), they are removed when creating specific models.

## Remove

* Data folder containing seed files as they are now inside the package.

# Kegg2bipartitegraph v0.1.2 (2024-11-26)

## Fix

* Issue with GitHub Actions.

# Kegg2bipartitegraph v0.1.1 (2024-11-26)

## Fix

* An issue with gene setAssociation when gene ID is only numeric.
* Test with organism.

# Kegg2bipartitegraph v0.1.0 (2024-01-09)

WARNING: Rename command `kegg2bipartitegraph` into `k2bg`.

## Add

* `reconstruct_from_genbank`: using EC numbers from GenBank files to create metabolic graphs. Add also associated test.
* Mapping between Gene Ontology terms and EC numbers using Gene Ontology `ec2go` into a new file `ec_to_gos.tsv`. Used in `reconstruct_from_esmecata`, `reconstruct_from_eggnog` and `reconstruct_from_genbank`.
* Hierarchy for metabolites, pathways and modules. This is used to identify classes of modules in organisms.
* Function to handle gene compatibility to SBML.

## Fix

* Issue with gene IDs not compatible with SBML format.

## Modify

* Rename command `kegg2bipartitegraph` into `k2bg`.
* Several optimisation to increase speed on large dataset.
* Update license year.
* Update reference model to KEGG 109.

# Kegg2bipartitegraph v0.0.3 (2023-12-11)

## Add

* `reconstruct_from_picrust`: using KEGG Orthologs and EC numbers from picrust2 to create metabolic graphs for each sequences.

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