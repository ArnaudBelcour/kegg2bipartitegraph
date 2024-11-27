import os
from cobra.io import write_sbml_model, read_sbml_model
from cobra import Model

import kegg2bipartitegraph

kegg2bipartitegraph_path = kegg2bipartitegraph.__path__[0]
# Seed created in Weber Zendrera et al. (2021): https://doi.org/10.1038/s41598-021-91486-8
# from https://github.com/AWebZen/FunctionalPrediction5000species/blob/2ff7c4fd4092a8b565da9d3fa2f4b557d9954bf0/utils_general.py

MESO_MEDIUM = [ #letort c et al, 2001 : https://mediadb.systemsbiology.net/defined_media/media/322/
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


PSYCHROMED = [ #Maria-Paz Cortes https://www.frontiersin.org/articles/10.3389/fmicb.2017.02462/full
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


HYPERTHERM_MED = [ #https://mediadb.systemsbiology.net/defined_media/media/382/ Rinker kd et al, 2000 Thermotoga maritima
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


THERM_MED = [ #https://mediadb.systemsbiology.net/defined_media/media/227/ Suzuki et al, 2001 Hydrogenobacter thermophilus TK-6
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

def write_txt_seed_file(metabolites, txt_file):
    with open(txt_file, 'w') as open_output_file:
        for metabolite in metabolites:
            open_output_file.write(metabolite)


kegg_model_sbml = os.path.join(kegg2bipartitegraph_path, 'data', 'kegg_model', 'kegg_model.sbml')
kegg_model = read_sbml_model(kegg_model_sbml)

# Create seeds file from KEGG model and list of metabolites.
meso_medium_seeds = [metabolite for metabolite in kegg_model.metabolites if metabolite.id in MESO_MEDIUM]
species_model = Model('Seeds_Meso_Medium')
species_model.add_metabolites(meso_medium_seeds)
write_sbml_model(species_model, 'seeds_meso_medium.sbml')
write_txt_seed_file(MESO_MEDIUM, 'seeds_meso_medium.txt')

psychromed_seeds = [metabolite for metabolite in kegg_model.metabolites if metabolite.id in PSYCHROMED]
species_model = Model('Seeds_Psychromed')
species_model.add_metabolites(psychromed_seeds)
write_sbml_model(species_model, 'seeds_psychromed.sbml')
write_txt_seed_file(PSYCHROMED, 'seeds_psychromed.txt')

hyperterm_med_seeds = [metabolite for metabolite in kegg_model.metabolites if metabolite.id in HYPERTHERM_MED]
species_model = Model('Seeds_Hypertherm_Med')
species_model.add_metabolites(hyperterm_med_seeds)
write_sbml_model(species_model, 'seed_hyperterm_med.sbml')
write_txt_seed_file(HYPERTHERM_MED, 'seed_hyperterm_med.txt')

therm_med_seeds = [metabolite for metabolite in kegg_model.metabolites if metabolite.id in THERM_MED]
species_model = Model('Therm_Med')
species_model.add_metabolites(therm_med_seeds)
write_sbml_model(species_model, 'seeds_therm_med.sbml')
write_txt_seed_file(THERM_MED, 'seeds_therm_med.txt')

union_all_seeds = set(MESO_MEDIUM).union(set(PSYCHROMED)).union(set(HYPERTHERM_MED)).union(set(PSYCHROMED))
kegg_model_union_all_seeds = [metabolite for metabolite in kegg_model.metabolites if metabolite.id in union_all_seeds]
species_model = Model('Union_All_Seeds')
species_model.add_metabolites(kegg_model_union_all_seeds)
write_sbml_model(species_model, 'union_all_seeds.sbml')
write_txt_seed_file(union_all_seeds, 'union_all_seeds.txt')