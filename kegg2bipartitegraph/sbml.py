# Copyright (C) 2021-2024 Arnaud Belcour - Inria, Univ Rennes, CNRS, IRISA Dyliss
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

SBML_CHARACTER_TO_REPLACE = ['-', '|', '/', '(', ')', "'", '=', '#', '*', '.',
                        ':', '!', '+', '[', ']', ',', ' ']

def gene_to_sbml(gene):
    """ Replace incorrect character for SBML in gene by unicode value.

    Args:
        gene (str): gene with uncorrect character

    Returns:
        gene (str): gene with unicode value
    """
    for character in SBML_CHARACTER_TO_REPLACE:
        gene = gene.replace(character, "__" + str(ord(character)) + "__")
    return gene


def libsbml_check(value, message):
    """If 'value' is None, prints an error message constructed using
    'message' and then exits with status code 1.  If 'value' is an integer,
    it assumes it is a libSBML return status code.  If the code value is
    LIBSBML_OPERATION_SUCCESS, returns without further action; if it is not,
    prints an error message constructed using 'message' along with text from
    libSBML explaining the meaning of the code, and exits with status code 1.

    Args:
        value (int): return value by a libsbml function or class
        message (str): associated string message
    """
    if value == None:
        raise SystemExit('LibSBML returned a null value trying to ' + message + '.')
    elif type(value) is int:
        if value == libsbml.LIBSBML_OPERATION_SUCCESS:
            return
        else:
            err_msg = 'Error encountered trying to ' + message + '.' \
                 + 'LibSBML returned error code ' + str(value) + ': "' \
                 + libsbml.OperationReturnValue_toString(value).strip() + '"'
            raise TypeError(err_msg)
    else:
        return


def initiate_sbml_model(model_name):
    """Initiate a libsbml model.

    Args:
        model_name (str): Name of the model
        taxon_reactions (dict): reaction ID as key and associated genes as value
        taxon_pathways (list): list of pathways in organism
        taxon_modules (lsit): list of modules in organism

    Returns:
        document: libsbml document
        model: libsbml model
        model_fbc: FBC of libsbml model
        model_groups: groups of libsbml model
    """
    # Set libsbml namespace.
    sbml_ns = libsbml.SBMLNamespaces(3, 1, 'groups', 1)  # SBML L3V1 with groups
    sbml_ns.addPackageNamespace("fbc", 2)  # fbc-v2

    # Create document
    document = libsbml.SBMLDocument(sbml_ns)
    model = document.createModel(model_name)

    document.enablePackage(libsbml.FbcExtension.getXmlnsL3V1V2(), 'fbc', True)
    document.setPackageRequired("fbc", False)
    model_fbc = model.getPlugin('fbc')
    model_fbc.setStrict(True)
    model_groups = model.getPlugin("groups")

    # Set units.
    libsbml_check(model,                              'create model')
    libsbml_check(model.setTimeUnits("second"),       'set model-wide time units')
    libsbml_check(model.setExtentUnits("mole"),       'set model units of extent')
    libsbml_check(model.setSubstanceUnits('mole'),    'set model substance units')

    math_ast = libsbml.parseL3Formula('FLUX_VALUE')
    libsbml_check(math_ast, 'create AST for rate expression')

    # Set compartments.
    compart = model.createCompartment()
    libsbml_check(compart,'create compartment')
    libsbml_check(compart.setId('c'),'set compartment id c')
    libsbml_check(compart.setSize(1),'set size for compartment id c')
    libsbml_check(compart.setConstant(True),'set constant for compartment id c')
    libsbml_check(compart.setName("cytosol"),'set compartment name cytosol')

    # Set default bound values.
    default_lb = model.createParameter()
    default_lb.setId('default_lower_bound')
    default_lb.setValue(-1000)
    default_lb.setConstant(True)

    default_ub = model.createParameter()
    default_ub.setId('default_upper_bound')
    default_ub.setValue(1000)
    default_ub.setConstant(True)

    zero_bound = model.createParameter()
    zero_bound.setId('default_zero_bound')
    zero_bound.setValue(0)
    zero_bound.setConstant(True)

    return document, model, model_fbc, model_groups


def create_sbml_from_kegg_reactions(org_name, reference_reactions, reference_species, reference_groups, taxon_reactions, taxon_pathways, taxon_modules):
    """Create a SBML model from the KEGG reference model and reference found for the taxon.

    Args:
        org_name (str): name of the corresponding organism
        reference_reactions (dict): dictionary of reference model reactions with reaction ID as key and libsbml reaction object as value
        reference_species (list): list of libsbml species object
        reference_groups (list): list of lisbml group object
        taxon_reactions (dict): reaction ID as key and associated genes as value
        taxon_pathways (list): list of pathways in organism
        taxon_modules (lsit): list of modules in organism

    Returns:
        kegg_document: libsbml document
        kegg_model: libsbml model
    """
    # Create KEGG model.
    document, model, model_fbc, model_groups = initiate_sbml_model(org_name)

    # Add reactions find during the reconstruction process.
    already_added_genes = []
    kept_metabolites = []

    # Extract reaction found both in organism and in reference.
    reaction_to_adds = set(list(taxon_reactions.keys())).intersection(set(list(reference_reactions.keys())))
    for reaction_id in reaction_to_adds:
        reaction = reference_reactions[reaction_id]
        ref_r_fbc: "libsbml.FbcReactionPlugin" = reaction.getPlugin("fbc")
        # Remove reference Gene Product Association.
        ref_r_fbc.unsetGeneProductAssociation()

        # Add genes to organism model. 
        for gene in taxon_reactions[reaction_id]:
            if gene not in already_added_genes:
                characters_in_genes = set([char for gene in taxon_reactions[reaction_id] for char in list(gene)])
                if len(set(SBML_CHARACTER_TO_REPLACE).intersection(characters_in_genes)) > 0:
                    gene = gene_to_sbml(gene)
                else:
                    gene = gene
                if gene.isnumeric():
                    gene = 'g' + gene
                gene_prod = model_fbc.createGeneProduct()
                gene_prod.setId(gene), 'add gene %s' %gene
                gene_prod.setName(gene)
                gene_prod.setLabel(gene)
                already_added_genes.append(gene)

        # Add new GPR association to reference reactions.
        gpr_association = ref_r_fbc.createGeneProductAssociation()
        if len(taxon_reactions[reaction_id]) > 0:
            # Handle incorrect character in genes names.
            characters_in_genes = set([char for gene in taxon_reactions[reaction_id] for char in list(gene)])
            if len(set(SBML_CHARACTER_TO_REPLACE).intersection(characters_in_genes)) > 0:
                genes = [gene_to_sbml(gene) for gene in taxon_reactions[reaction_id]]
            else:
                genes = [gene for gene in taxon_reactions[reaction_id]]
            if any([gene.isnumeric() for gene in genes]):
                genes = ['g'+gene for gene in taxon_reactions[reaction_id]]
            libsbml_check(gpr_association.setAssociation(' or '.join(genes), True, True), "set gpr: ")

        # Add reaction to organism model.
        org_reaction = reaction.clone()
        model.addReaction(org_reaction)

        # Get list of metabolites to add in organism model.
        kept_metabolites.extend([i.species for i in reaction.getListOfReactants()])
        kept_metabolites.extend([i.species for i in reaction.getListOfProducts()])

    # Add metabolites from reaction to organism model.
    kept_metabolites = set(kept_metabolites)
    for species in reference_species:
        if species.id in kept_metabolites:
            model.addSpecies(species)

    # Add groups in metabolic networks.
    # If condition required for compatibility with kegg archive 106.
    if reference_groups is not None:
        for group in reference_groups.getListOfGroups():
            group_id = group.id
            if group_id in taxon_pathways or group_id in taxon_modules:
                org_group = model_groups.createGroup()
                org_group.setId(group.id)
                org_group.setName(group.name)
                org_group.setKind("partonomy")
                for rxn_member in group.getListOfMembers():
                    if rxn_member.id in taxon_reactions:
                        member = org_group.createMember()
                        member.setId(rxn_member.id)
                        member.setIdRef(rxn_member.id)

    return document, model

