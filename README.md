# kegg2bipartitegraph

kegg2bipartitegraph is a Python package to create KEGG graphs. The main idea of this package is to create metabolic graphs from KEGG database according to the ones used in the article [Weber Zendrera et al. (2021)](https://www.nature.com/articles/s41598-021-91486-8). In this article, the authors creates the metabolic networks from the organism of KEGG (accessible in this [github repository](https://github.com/AWebZen/FunctionalPrediction5000species)). Using annotation (from EsMeCaTa, eggnog-mapper, KofamKOALA) or a KEGG organism ID, **kegg2bipartitegraph** maps the EC and reconstruct metabolic networks associated with the organism following the proposal of this article.

## Installation

This package can be installed as a python repository with:

```pip install -e . ```

If it becomes mature enough, a pip package could be created.

## Usage

It is divided in different parts:

- `kegg2bipartitegraph reference` an optional ones that creates the reference data (especially the universal reference metabolic graphs). By default, these data are precomputed and available within the package located in `kegg2bipartitegraph/data/kegg_model`.

- subcommand to reconstruct metabolic graphs from different inputs:
    - `kegg2bipartitegraph reconstruct_from_esmecata` takes as input the annotation output folder from [EsMeCaTa](https://github.com/AuReMe/esmecata) and reconstruct the metabolic networks associated with each taxon.
    - `kegg2bipartitegraph reconstruct_from_eggnog` takes as input the annotation output file from [eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper) to map the EC to KEGG reactions.
    - `kegg2bipartitegraph reconstruct_from_kofamkoala` takes as input the result from [KofamKOALA](https://www.genome.jp/tools/kofamkoala/).
    - `kegg2bipartitegraph reconstruct_from_organism` takes as input an organism ID from KEGG (such as `hsa` for human or `eco` for *Escherichia coli*).

## Reference model

The `kegg2bipartitegraph reference` is to be used only if you want to update the KEGG reference data. First, delete the data contain in `kegg2bipartitegraph/data/kegg_model`, then use this command to download all the required data. This step is long, it is advised to not use it.

It will create 4 files:

- `kegg_model.sbml`: a universal graph containing most of the reactions contained in KEGG database. Such as in the graph made by [Weber Zendrera et al. (2021)](https://www.nature.com/articles/s41598-021-91486-8), 14 cofactors have been removed (H2O, ATP, ADP, NAD+, NADH, NADP+, NADPH, CO2, ammonia, sulfate, thioredoxin, phosphate, pyrophosphate (PPi), and H+). Also the stoechiometry is simplified as these metabolic networks are created in order to be used in topological analysis. **So they are not supposed to be used with other methods (such as Constraint-Based Modelling)**.

- several mapping files to go from annotation (especially EC number) to KEGG reactions: `kegg_compound_name.tsv`, `kegg_mapping.tsv` and `kegg_pathways.tsv`.

## Output files of other command

The other subcommands will reconstruct draft metabolic networks by mapping the annotation with the metabolic graphs contained in kegg2bipartitegraph.

Then it will create multiple files:

- a sbml file containing the metabolic network that can be used with topological analysis methods (such as [MeneTools](https://github.com/cfrioux/MeneTools), [MiSCoTo](https://github.com/cfrioux/miscoto) or [Metage2Metabo](https://github.com/AuReMe/metage2metabo)).

- a graphml file containing the metabolic network as a bipartite graph. At this moment, it is not used, but I am currently adaptating the scope method of [Weber Zendrera et al. (2021)](https://www.nature.com/articles/s41598-021-91486-8) to automatise its use with this package.

- tsv files indicating the pathways/modules contained in the metabolic networks, their completness ratio and the associated reactions.

- a tsv file showing KO information if the option has been used.

- several statistics/metadata/log files.
