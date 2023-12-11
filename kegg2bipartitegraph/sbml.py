# Copyright (C) 2021-2023 Arnaud Belcour - Inria, Univ Rennes, CNRS, IRISA Dyliss
# Univ. Grenoble Alpes, Inria, Microcosme
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

def create_sbml_from_kegg_reactions(kegg_sbml_model_path, taxon_reactions, taxon_pathways, taxon_modules):
    """Create a SBML model from the KEGG reference model and reference found for the taxon.

    Args:
        kegg_sbml_model_path (str): KEGG SBML reference file
        taxon_reactions (dict): reaction ID as key and associated genes as value
        taxon_pathways (list): list of pathways in organism
        taxon_modules (lsit): list of modules in organism

    Returns:
        kegg_document: libsbml document
        kegg_model: libsbml model
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

    # Keep only groups associated with pathway/module present in the organisms.
    model_groups = kegg_model.getPlugin("groups")
    group_to_delete = []
    members_to_delete = {}

    # If condition required for compatibility with kegg archive 106.
    if model_groups is not None:
        for group in model_groups.getListOfGroups():
            group_id = group.id
            if group_id not in taxon_pathways and group_id not in taxon_modules:
                group_to_delete.append(group_id)
            else:
                for rxn_member in group.getListOfMembers():
                    if rxn_member.id not in taxon_reactions:
                        if group_id not in members_to_delete:
                            members_to_delete[group_id] = [rxn_member.id]
                        else:
                            members_to_delete[group_id].append(rxn_member.id)

        # Remove group not present in modules/pathways of organism.
        for group in group_to_delete:
            model_groups.removeGroup(group)

        # Remove reaction member of group not present in reaction list of organism.
        for group in members_to_delete:
            group_to_modify = model_groups.getGroup(group)
            for member in members_to_delete[group]:
                group_to_modify.removeMember(member)

    return kegg_document, kegg_model

