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

import csv
import os

import libsbml
import kegg2bipartitegraph

from kegg2bipartitegraph.sbml import libsbml_check, initiate_sbml_model

kegg2bipartitegraph_path = kegg2bipartitegraph.__path__[0]

# Seed created in Weber Zendrera et al. (2021): https://doi.org/10.1038/s41598-021-91486-8
# from https://github.com/AWebZen/FunctionalPrediction5000species/blob/2ff7c4fd4092a8b565da9d3fa2f4b557d9954bf0/utils_general.py

# Letort c et al, 2001 : https://mediadb.systemsbiology.net/defined_media/media/322/
MESO_MEDIUM = [ 
"C00568", #4-Aminobenzoate
"C00147", #Adenine
"C00041", #Alanine
"C01342", "C00158", #Ammonium citrate
"C00062", #Arginine
"C00072", #Ascorbate
"C00152", #Asparagine
"C00049", #Aspartate
"C00120", #Biotin
"C08130", "C00076", "C00698", #Calcium chloride anhydrous, Ca, Cl
"C00864", #Calcium pantothenate
"C00175", #Cobalt chloride: Cb
"C00070", "C00059", #Cupric sulfate: Cu, SO42-
"C00097", #Cysteine
"C00740", "C00065", #DL-Serine
"C06420", "C00082", #DL-Tyrosine
"C00855", "C00073", #DL-methionine
"C14818", "C14819",  #Ferrous chloride
"C00504", #Folate
"C00064", #Glutamine
"C00037", #Glycine
"C00242", #Guanine
"C00135", #Histidine
"C00294", "C00081", "C00104",#Inosine, Inosine TP, IDP
"C00407", #Isoleucine
"C00025", #L-Glutamate
"C00243", #Lactose
"C00123", #Leucine
"C00725", #Lipoate
"C00047", #Lysine
"C00305", #Magnesium chloride
"C00034", "C19610", "C19611", #Mn, Mn2+, Mn3+ :Manganese sulfate
"C00253", #Nicotinate
"C00295", #Orotate
"C00079", #Phenylalanine
"C13197", #Potassium dibasic phosphate
"C00238", "C00009", #Potassium dihydrogen phosphate
"C00148", #Proline
"C00534", #Pyridoxamine HCl
"C00314", #Pyridoxine HCl
"C00255", #Riboflavin
"C01330", "C00033", #Sodium acetate
"C00378", #Thiamine HCl
"C00188", #Threonine
"C00214", #Thymidine
"C00078", #Tryptophan
"C00106", #Uracil
"C00183", #Valine
"C05776", #Vitamin B12
"C00385", #Xanthine
"C00038", #Zinc sulfate
]

#Maria-Paz Cortes https://www.frontiersin.org/articles/10.3389/fmicb.2017.02462/full
PSYCHROMED = [ 
    "C00407", # L-Isoleucine
    "C00123", # L-Leucine
    "C00183", # L-Valine
    "C00073", # L-Met
    "C00135", # L-His
    "C00062", # L-Arg
    "C00097", # L-Cys
    "C00037", # Gly
    "C00041", # L-Ala
    "C00049", # L-Asp
    "C00025", # L-Glu
    "C00064", # L-Gln
    "C00079", # L-Phe
    "C00148", # L-Pro
    "C00065", # L-Ser
    "C00188", # L-Threonine
    "C00082", # L-Tyrosine
    "C00314", # Pyridoxine
    "C00378", # Thiamine
    "C00864", # D-pantothenate
    "C14819", # Fe3+
    "C00122", # Fumarate
    "C00149", "C00497", # (L/D)-Malate
    "C00042", # Succinate
    "C00026", # Ketoglutarate (2-Oxoglutarate; Oxoglutaric acid; 2-Ketoglutaric acid; alpha-Ketoglutaric acid)
    "C00031", # D-Glu
    "C00124", # D-Gal
    "C00051", # Glutathione
    "C00116", # Glycerol
    "C01342", # Ammonium
    "C00147", # Adenine
    ]

#https://mediadb.systemsbiology.net/defined_media/media/382/ Rinker kd et al, 2000 Thermotoga maritima
HYPERTHERM_MED = [
"C00568", # 4-Aminobenzoate
"C01342", # Ammonium chloride
"C00120", # Biotin
"C12486", # Boric acid
"C08130", "C00698", # Calcium chloride anhydrous
"C00076", "C00864", # Calcium pantothenate
"C00504", # Folate
"C00725", # Lipoate
"C07755", "C00305", # Magnesium chloride
"C00208", # Maltose
"C00253", # Nicotinate
"C00238", # Potassium chloride: K
# "C13197", # Potassium dibasic phosphate
"C08219", "C01382", # Potassium iodide: K iodide, iodine
"C00314", # Pyridoxine HCl
"C11178", # Resazurin
"C00255", # Riboflavin
"C01330", "C01324", # Sodium bromide
# "", # Sodium chloride
"C00059", # Sodium sulfate
# "", # Sodium sulfide
"C20679", # Sodium tungstate
"C13884", # Strontium chloride
"C00378", # Thiamine HCl
"C05776", # Vitamin B12
]

#https://mediadb.systemsbiology.net/defined_media/media/227/ Suzuki et al, 2001 Hydrogenobacter thermophilus TK-6
THERM_MED = [
"C01342", "C00698", # Ammonium chloride
"C12486", # Boric acid
"C00076", # Calcium chloride anhydrous
 # Cupric chloride: Cu2+, Cl-
"C00070", "C00059",  # Cupric sulfate
"C07755", "C00305", # Magnesium sulfate
"C00034", "C19610", "C19611", #Mn, Mn2+, Mn3+ :Manganese sulfate
"C00150", #Molybdenum trioxide
"C14818", # Ferrous sulfate
"C00238", "C00009", "C13197", # Potassium dibasic phosphate
"C08219", # Potassium dihydrogen phosphate
"C01330", # Sodium chloride
"C00038", #Zinc sulfate
]

def write_sbml_species_file(metabolites, compounds, output_sbml):
    # Initiate libsbml model.
    document, model, model_fbc, model_groups = initiate_sbml_model('seed')
    for metabolite in metabolites:
        metabolite_id = metabolite.id
        s = model.createSpecies()
        libsbml_check(s, 'create species')
        libsbml_check(s.setId(metabolite_id), 'set species id %s' %metabolite_id)
        libsbml_check(s.setMetaId(metabolite_id), 'set species meta id %s' %metabolite_id)
        libsbml_check(s.setBoundaryCondition(False), 'set boundaryCondition to False')
        libsbml_check(s.setHasOnlySubstanceUnits(False), 'set setHasOnlySubstanceUnits to False')
        libsbml_check(s.setConstant(False), 'set setConstant to False')
        libsbml_check(s.setInitialAmount(0.0), 'set initAmount')
        libsbml_check(s.setName(compounds[metabolite_id]), 'set species Name {0}'.format(compounds[metabolite_id]))
        libsbml_check(s.setCompartment('c'), 'set species compartment c')
    libsbml.writeSBMLToFile(document, output_sbml)

def write_txt_seed_file(metabolites, txt_file):
    with open(txt_file, 'w') as open_output_file:
        for metabolite in metabolites:
            open_output_file.write(metabolite+'\n')

def create_seeds_file(output_folder):
    kegg_model_path = os.path.join(kegg2bipartitegraph_path, 'data', 'kegg_model',)
    kegg_model_sbml = os.path.join(kegg_model_path, 'kegg_model.sbml')
    reader = libsbml.SBMLReader()
    kegg_document = reader.readSBML(kegg_model_sbml)
    reference_kegg_model = kegg_document.getModel()
    reference_species = reference_kegg_model.getListOfSpecies()

    kegg_compound_file_path = os.path.join(kegg_model_path, 'kegg_compound_name.tsv')
    compounds = {}
    with open(kegg_compound_file_path, 'r') as output_file:
        csvreader = csv.reader(output_file, delimiter='\t')
        next(csvreader)
        for line in csvreader:
            compounds[line[0]] = line[1]

    # Create seeds file from KEGG model and list of metabolites.
    meso_medium_seeds = [metabolite for metabolite in reference_species if metabolite.id in MESO_MEDIUM]
    write_sbml_species_file(meso_medium_seeds, compounds, os.path.join(output_folder, 'seeds_meso_medium.sbml'))
    write_txt_seed_file(MESO_MEDIUM, os.path.join(output_folder, 'seeds_meso_medium.txt'))

    psychromed_medium_seeds = [metabolite for metabolite in reference_species if metabolite.id in PSYCHROMED]
    write_sbml_species_file(psychromed_medium_seeds, compounds, os.path.join(output_folder, 'seeds_psychromed.sbml'))
    write_txt_seed_file(PSYCHROMED, os.path.join(output_folder, 'seeds_psychromed.txt'))

    hyperterm_med_seeds = [metabolite for metabolite in reference_species if metabolite.id in HYPERTHERM_MED]
    write_sbml_species_file(hyperterm_med_seeds, compounds, os.path.join(output_folder, 'seed_hyperterm_med.sbml'))
    write_txt_seed_file(HYPERTHERM_MED, os.path.join(output_folder, 'seed_hyperterm_med.txt'))

    therm_med_seeds = [metabolite for metabolite in reference_species if metabolite.id in THERM_MED]
    write_sbml_species_file(therm_med_seeds, compounds, os.path.join(output_folder, 'seeds_therm_med.sbml'))
    write_txt_seed_file(THERM_MED, os.path.join(output_folder, 'seeds_therm_med.txt'))

    union_all_seed_ids = set(MESO_MEDIUM).union(set(PSYCHROMED)).union(set(HYPERTHERM_MED)).union(set(PSYCHROMED))
    union_all_seeds = [metabolite for metabolite in reference_species if metabolite.id in union_all_seed_ids]
    write_sbml_species_file(union_all_seeds, compounds, os.path.join(output_folder, 'seed_union_all.sbml'))
    write_txt_seed_file(union_all_seed_ids, os.path.join(output_folder, 'seed_union_all.txt'))