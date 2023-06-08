# kegg2bipartitegraph

kegg2bipartitegraph is a Python package to create KEGG graphs. The main idea of this package is to create metabolic graphs from KEGG database according to the ones used in the article [Weber Zendrera et al. (2021)](https://www.nature.com/articles/s41598-021-91486-8). In this article, the authors creates the metabolic networks from the organism of KEGG (accessible in this [github repository](https://github.com/AWebZen/FunctionalPrediction5000species)). Using annotation (at this moment only from [EsMeCaTa](https://github.com/AuReMe/esmecata)), kegg2bipartitegraph maps the EC and reconstruct metabolic networks associated with the organism.

## Installation

This package can be installed as a python repository with:

```pip install -e . ```

If it becomes mature enough, a pip package could be created.

## Usage

It is divided in two steps:

- `kegg2bipartitegraph reference` an optional ones that creates the reference data (especially the universal reference metabolic graphs). By default, these data are precomputed and available within the package located in `kegg2bipartitegraph/data/kegg_model`.

- `kegg2bipartitegraph esmecata` which takes as input the annotation output folder from EsMeCaTa and reconstruct the metabolic networks associated with each taxon.

## Output

The `kegg2bipartitegraph reference` is to be used only if you want to update the KEGG reference data. First, delete the data contain in `kegg2bipartitegraph/data/kegg_model`, then use this command to download all the required data. This step is long, it is advised to not use it.

It will create 4 files:

- `kegg_model.sbml`: a universal graph containing most of the reactions contained in KEGG database. Such as in the graph made by [Weber Zendrera et al. (2021)](https://www.nature.com/articles/s41598-021-91486-8), 14 cofactors have been removed (H2O, ATP, ADP, NAD+, NADH, NADP+, NADPH, CO2, ammonia, sulfate, thioredoxin, phosphate, pyrophosphate (PPi), and H+). Also the stoechiometry is simplified as these metabolic networks are created in order to be used in topological analysis. **So they are supposed to be used with other metabolic methods (such as Constraint-Based Modelling)**.

- several mapping files to go from annotation (especially EC number) to KEGG reactions: `kegg_compound_name.tsv`, `kegg_mapping.tsv` and `kegg_pathways.tsv`.

The `kegg2bipartitegraph esmecata` command will reconstruct draft metabolic networks by mapping the annotation from EsMeCaTa with the metabolic graphs contained in kegg2bipartitegraph.

Then for each taxon of EsMeCaTa, it will create multiple files:

- a sbml file containing the metabolic network that can be used with topological analysis methods (such as [MeneTools](https://github.com/cfrioux/MeneTools), [MiSCoTo](https://github.com/cfrioux/miscoto) or [Metage2Metabo](https://github.com/AuReMe/metage2metabo)).

- a tv file indicating the pathways contained in the metabolic networks, their completness ratio and the associated reactions.


