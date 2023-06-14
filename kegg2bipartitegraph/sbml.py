# Copyright (C) 2021-2023 Arnaud Belcour - Inria Dyliss
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>

import libsbml

def create_sbml_from_kegg_reactions(kegg_sbml_model_path, taxon_reactions):
    """Create a SBML model from the KEGG reference model and reference found for the taxon.

    Args:
        kegg_sbml_model_path (str): KEGG code for an organism
        taxon_reactions (dict): reaction ID as key and associated genes as value
    """
    # Read the reference KEGG sbml file.
    # Use it to create the organism sbml file.
    reader = libsbml.SBMLReader()
    kegg_document = reader.readSBML(kegg_sbml_model_path)
    kegg_model = kegg_document.getModel()

    # Remove all the gene products from the model.
    model_fbc = kegg_model.getPlugin('fbc')
    model_fbc.setStrict(True)

    remove_gene_products = []
    for gene_product in model_fbc.getListOfGeneProducts():
        remove_gene_products.append(gene_product.id)
    for gene_product in remove_gene_products:
        model_fbc.removeGeneProduct(gene_product)

    # Keep only the reactions find during the reconstruction process.
    genes = []
    remove_reactions = []
    kept_metabolites = []
    for reaction in kegg_model.getListOfReactions():
        if reaction.id in taxon_reactions:
            # Add the gene associated with the organism.
            r_fbc: "libsbml.FbcReactionPlugin" = reaction.getPlugin("fbc")
            gpr_association = r_fbc.createGeneProductAssociation()

            for gene in taxon_reactions[reaction.id]:
                if gene not in genes:
                    gene_prod = model_fbc.createGeneProduct()
                    gene_prod.setId(gene), 'add gene %s' %gene
                    gene_prod.setName(gene)
                    gene_prod.setLabel(gene)
                    genes.append(gene)
            gpr = ' or '.join(taxon_reactions[reaction.id])
            gpr_association.setAssociation(gpr, True, True)
            kept_metabolites.extend([i.species for i in reaction.getListOfReactants()])
            kept_metabolites.extend([i.species for i in reaction.getListOfProducts()])
        else:
            remove_reactions.append(reaction.id)

    # Remove reactions not found in organism.
    for reaction_id in remove_reactions:
        kegg_model.removeReaction(reaction_id)
    # Remove metabolites not found in organism.
    remove_metabolites = set([m.id for m in kegg_model.getListOfSpecies()]) - set(kept_metabolites)
    for metabolite_id in remove_metabolites:
        kegg_model.removeSpecies(metabolite_id)

    return kegg_document, kegg_model

