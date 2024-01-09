[![PyPI version](https://img.shields.io/pypi/v/kegg2bipartitegraph.svg)](https://pypi.org/project/kegg2bipartitegraph/) [![GitHub license](https://img.shields.io/github/license/AuReMe/metage2metabo.svg)](https://github.com/AuReMe/metage2metabo/blob/master/LICENSE) [![KEGG version](https://img.shields.io/badge/KEGG-109-brightgreen)](https://www.genome.jp/kegg/docs/upd_all.html)

# kegg2bipartitegraph

kegg2bipartitegraph is a Python package to create KEGG graphs. The main idea of this package is to create metabolic graphs from KEGG database according to the ones used in the article [Weber Zendrera et al. (2021)](https://www.nature.com/articles/s41598-021-91486-8). In this article, the authors creates the metabolic networks from the organism of KEGG (accessible in this [github repository](https://github.com/AWebZen/FunctionalPrediction5000species)). Using annotation (from EsMeCaTa, eggnog-mapper, KofamKOALA or GenBank files) or a KEGG organism ID, **kegg2bipartitegraph** maps the EC (and Gene Ontology Terms) to KEGG reactions and reconstruct metabolic networks associated with the organism following the proposal of the article [Weber Zendrera et al. (2021)](https://www.nature.com/articles/s41598-021-91486-8).

## Installation

This package can be installed using pip:

```pip install kegg2bipartitegraph```

## Usage

Kegg2bipartitegraph can be called by using the command `k2bg`. It is divided in different parts:

- `k2bg reference` an optional ones that creates the reference data (especially the universal reference metabolic graphs). By default, these data are precomputed and available within the package located in `kegg2bipartitegraph/data/kegg_model`. It is not meant to be used by the user as it takes a long time to process, only used the precomputed reference files or archived files.

- subcommands to reconstruct metabolic graphs from different inputs:
    - `k2bg reconstruct_from_esmecata` takes as input the annotation output folder from [EsMeCaTa](https://github.com/AuReMe/esmecata) and reconstruct the metabolic networks associated with each taxon.
    - `k2bg reconstruct_from_eggnog` takes as input the annotation output file from [eggnog-mapper](https://github.com/eggnogdb/eggnog-mapper) to map the EC to KEGG reactions.
    - `k2bg reconstruct_from_kofamkoala` takes as input the result from [KofamKOALA](https://www.genome.jp/tools/kofamkoala/).
    - `k2bg reconstruct_from_picrust` takes as input the result folder from [picrust2](https://github.com/picrust/picrust2).
    - `k2bg reconstruct_from_organism` takes as input an organism ID from KEGG (such as `hsa` for human or `eco` for *Escherichia coli*). You can find the list of the accessile organisms in [KEGG website](https://www.genome.jp/kegg/catalog/org_list.html).
    - `k2bg reconstruct_from_genbank` takes as input a folder containing GenBank files.

## Online / Offline requirements

Multiple subcommands can be used to reconstruct draft networks. Some of them required an internet connection to work, you can see which ones in the following table:

| Subcommands  | Online  | Offline  |
|---|---|---|
| reconstruct_from_esmecata  | (Mapping of KOs)  | X (without mapping KOs)  |
|  reconstruct_from_eggnog |   |  X |
|  reconstruct_from_kofamkoala |   | X  |
|  reconstruct_from_picrust |   | X  |
|  reconstruct_from_genbank |   | X  |
|  reconstruct_from_organism | X  |   |
|  reference | X  |   |

## Reference model

The `k2bg reference` is to be used only if you want to update the KEGG reference data. First, delete the data contain in `kegg2bipartitegraph/data/kegg_model`, then use this command to download all the required data. This step is long, it is advised to not use it.

It will create several files:

- `kegg_model.sbml`: a universal graph containing most of the reactions contained in KEGG database. Such as in the graph made by [Weber Zendrera et al. (2021)](https://www.nature.com/articles/s41598-021-91486-8), 14 cofactors have been removed (H2O, ATP, ADP, NAD+, NADH, NADP+, NADPH, CO2, ammonia, sulfate, thioredoxin, phosphate, pyrophosphate (PPi), and H+). Also the stoechiometry is simplified as these metabolic networks are created in order to be used in topological analysis. **So they are not supposed to be used with other methods (such as Constraint-Based Modelling)**.

- `kegg_model.graphml`: the metabolic bipartite graph associated with the `kegg_model.sbml` file.

- several mapping files to go from annotation (especially EC number) to KEGG reactions: `kegg_compound_name.tsv`, `kegg_mapping.tsv` and `kegg_pathways.tsv`. Also a file to use KEGG hierarchy for pathway/module/metabolite: `kegg_hierarchy.json`. A mapping file `ec_to_gos.tsv` to convert Gene Ontology terms to EC number using the [go2ec file](https://www.ebi.ac.uk/GOA/EC2GO) provided by the Gene Ontology Consortium.

- `kegg_metadata.json`: a metadata file showing the metadata for the creation of the reference files for kegg2bipartitegraph.

## Output files of other command

The other subcommands will reconstruct draft metabolic networks by mapping the annotation with the metabolic graphs contained in kegg2bipartitegraph.

Then it will create multiple files:

- a sbml file containing the metabolic network that can be used with topological analysis methods (such as [MeneTools](https://github.com/cfrioux/MeneTools), [MiSCoTo](https://github.com/cfrioux/miscoto) or [Metage2Metabo](https://github.com/AuReMe/metage2metabo)).

- a graphml file containing the metabolic network as a bipartite graph. At this moment, it is not used, but I am currently adaptating the scope method of [Weber Zendrera et al. (2021)](https://www.nature.com/articles/s41598-021-91486-8) to automatise its use with this package.

- tsv files indicating the pathways/modules contained in the metabolic networks, their completness ratio and the associated reactions.

- a tsv file showing KO information if the option has been used.

- a `module_class.tsv` showing the absence/presence of generic module classes in the organism.

- several statistics/metadata/log files.

## Citation

At this moment, there are no articles for kegg2bipartitegraph, if you use it and want to cite it, you can cite this GitHub.

Also, please cite the following article:

- the article made by Adèle Weber Zendrera et al. (2021) that proposed this method:

Weber Zendrera, A., Sokolovska, N. & Soula, H.A. Functional prediction of environmental variables using metabolic networks. Scientific Reports  11, 12192 (2021). https://doi.org/10.1038/s41598-021-91486-8

- the `KEGG database`:

Kanehisa, M.,  Furumichi, M., Sato, Y., Kawashima, M., Ishiguro-Watanabe, M. KEGG for taxonomy-based analysis of pathways and genomes, Nucleic Acids Research, Volume 51, Issue D1, Pages D587–D592 (2023). https://doi.org/10.1093/nar/gkac963

Kanehisa, M., Goto, S. KEGG: Kyoto Encyclopedia of Genes and Genomes, Nucleic Acids Research, Volume 28, Issue 1, Pages 27–30 (2000). https://doi.org/10.1093/nar/28.1.27

Kanehisa, M. Toward understanding the origin and evolution of cellular organisms. Protein Science. 28: 1947–1951 (2019). https://doi.org/10.1002/pro.3715

- `bioservices` for the query on KEGG:

Cokelaer, T., Pultz, D., Harder, L., M., Serra-Musach, J., Saez-Rodriguez, J., BioServices: a common Python package to access biological Web Services programmatically, Bioinformatics, Volume 29, Issue 24, Pages 3241–3242 (2013). https://doi.org/10.1093/bioinformatics/btt547

- `libsbml` for the handling of the SBML:

Bornstein B. J., Keating S. M., Jouraku, A., Hucka, M., LibSBML: an API Library for SBML, Bioinformatics, Volume 24, Issue 6, Pages 880–881 (2008). https://doi.org/10.1093/bioinformatics/btn051

- `networkx` for the creation of the graphml:

Hagberg A. A., Schult D. A., Swart P. J. Exploring Network Structure, Dynamics, and Function using NetworkX, in: Varoquaux, G., Vaught, T., Millman, J. (Eds.), . Presented at the Proceedings of the Python in Science Conference (SciPy) 2008. 11–15. http://conference.scipy.org/proceedings/SciPy2008/paper_2/


If you have used the subcommand `reconstruct_from_genbank`, please also cite:

- `Biopython` for GenBank parsing:

Cock, P.J.A., Antao, T., Chang, J.T., Chapman, B.A., Cox, C.J., Dalke, A., Friedberg, I., Hamelryck, T., Kauff, F., Wilczynski, B., de Hoon, M.J.L. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 2009, 25, 1422–1423 https://doi.org/10.1093/bioinformatics/btp163.