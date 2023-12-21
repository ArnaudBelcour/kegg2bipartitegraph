import shutil
import subprocess
import libsbml

from kegg2bipartitegraph.genbank import create_gbff_network

def test_create_gbff_network():
    create_gbff_network('genbank', 'test_out')

    reader = libsbml.SBMLReader()
    sbml_document = reader.readSBML('test_out/sbml/betaox.sbml')
    sbml_model = sbml_document.getModel()

    expected_reactions = ['R00238', 'R00390', 'R00391', 'R00829', 'R00927', 'R01177', 'R01280', 'R01778', 'R01975', 'R01976', 'R02685',
                          'R03026', 'R03045', 'R03224', 'R03276', 'R03778', 'R03858', 'R03991', 'R04100', 'R04137', 'R04170', 'R04203',
                          'R04204', 'R04224', 'R04737', 'R04738', 'R04739', 'R04740', 'R04741', 'R04742', 'R04743', 'R04744', 'R04745',
                          'R04746', 'R04747', 'R04748', 'R04749', 'R04756', 'R05305', 'R05506', 'R05576', 'R05586', 'R05595', 'R06411',
                          'R06412', 'R06941', 'R06942', 'R07314', 'R07889', 'R07890', 'R07891', 'R07893', 'R07894', 'R07895', 'R07897',
                          'R07898', 'R07899', 'R07937', 'R07953', 'R08091', 'R08093', 'R08094', 'R08095', 'R00389', 'R01176']
    expected_metabolites = ['C00010', 'C00020', 'C00024', 'C00040', 'C00091', 'C00100', 'C00136', 'C00154', 'C00249', 'C00264', 'C00332',
                            'C00512', 'C00527', 'C00638', 'C00640', 'C00658', 'C00877', 'C00894', 'C01086', 'C01122', 'C01144', 'C01832',
                            'C01944', 'C02232', 'C02593', 'C02843', 'C02944', 'C03069', 'C03221', 'C03344', 'C03345', 'C03460', 'C03561',
                            'C04405', 'C05067', 'C05116', 'C05258', 'C05259', 'C05260', 'C05261', 'C05262', 'C05263', 'C05264', 'C05265',
                            'C05266', 'C05267', 'C05268', 'C05269', 'C05270', 'C05271', 'C05272', 'C05273', 'C05274', 'C05275', 'C05276',
                            'C05279', 'C05280', 'C05668', 'C05998', 'C06000', 'C06714', 'C06715', 'C07118', 'C11945', 'C11946', 'C11947',
                            'C14144', 'C14145', 'C16169', 'C16173', 'C16328', 'C16329', 'C16330', 'C16331', 'C16332', 'C16333', 'C16334',
                            'C16335', 'C16336', 'C16337', 'C16338', 'C16339', 'C16376', 'C16389', 'C16466', 'C16468', 'C16469', 'C16470',
                            'C16471', 'C00246', 'C21925', 'C21926']
    expected_modules_pathways_ids = ['M00013', 'M00032', 'M00085', 'M00086', 'M00087', 'M00088', 'M00095', 'M00113', 'M00373', 'M00374',
                                     'M00375', 'M00376', 'M00849', 'M00861', 'M00878', 'M00957', 'map00061', 'map00062', 'map00071', 'map00280',
                                     'map00310', 'map00360', 'map00362', 'map00380', 'map00410', 'map00592', 'map00620', 'map00627', 'map00630',
                                     'map00640', 'map00642', 'map00650', 'map00720', 'map00900', 'map00907', 'map00930', 'map01040', 'map01100',
                                     'map01110', 'map01120', 'map01200', 'map01212', 'map01220']

    found_reactions = [reaction.id for reaction in sbml_model.getListOfReactions()]
    found_metabolites = [reaction.id for reaction in sbml_model.getListOfSpecies()]

    model_groups = sbml_model.getPlugin("groups")
    modules_pathways_ids = [group.id for group in model_groups.getListOfGroups()]

    assert sorted(found_reactions) == sorted(expected_reactions)
    assert sorted(found_metabolites) == sorted(expected_metabolites)
    assert sorted(modules_pathways_ids) == sorted(expected_modules_pathways_ids)

    shutil.rmtree('test_out')

def test_create_gbff_network_cli():
    subprocess.call(['k2bg', 'reconstruct_from_genbank', '-i', 'genbank', '-o', 'test_out'])
    reader = libsbml.SBMLReader()
    sbml_document = reader.readSBML('test_out/sbml/betaox.sbml')
    sbml_model = sbml_document.getModel()

    expected_reactions = ['R00238', 'R00390', 'R00391', 'R00829', 'R00927', 'R01177', 'R01280', 'R01778', 'R01975', 'R01976', 'R02685',
                          'R03026', 'R03045', 'R03224', 'R03276', 'R03778', 'R03858', 'R03991', 'R04100', 'R04137', 'R04170', 'R04203',
                          'R04204', 'R04224', 'R04737', 'R04738', 'R04739', 'R04740', 'R04741', 'R04742', 'R04743', 'R04744', 'R04745',
                          'R04746', 'R04747', 'R04748', 'R04749', 'R04756', 'R05305', 'R05506', 'R05576', 'R05586', 'R05595', 'R06411',
                          'R06412', 'R06941', 'R06942', 'R07314', 'R07889', 'R07890', 'R07891', 'R07893', 'R07894', 'R07895', 'R07897',
                          'R07898', 'R07899', 'R07937', 'R07953', 'R08091', 'R08093', 'R08094', 'R08095', 'R00389', 'R01176']
    expected_metabolites = ['C00010', 'C00020', 'C00024', 'C00040', 'C00091', 'C00100', 'C00136', 'C00154', 'C00249', 'C00264', 'C00332',
                            'C00512', 'C00527', 'C00638', 'C00640', 'C00658', 'C00877', 'C00894', 'C01086', 'C01122', 'C01144', 'C01832',
                            'C01944', 'C02232', 'C02593', 'C02843', 'C02944', 'C03069', 'C03221', 'C03344', 'C03345', 'C03460', 'C03561',
                            'C04405', 'C05067', 'C05116', 'C05258', 'C05259', 'C05260', 'C05261', 'C05262', 'C05263', 'C05264', 'C05265',
                            'C05266', 'C05267', 'C05268', 'C05269', 'C05270', 'C05271', 'C05272', 'C05273', 'C05274', 'C05275', 'C05276',
                            'C05279', 'C05280', 'C05668', 'C05998', 'C06000', 'C06714', 'C06715', 'C07118', 'C11945', 'C11946', 'C11947',
                            'C14144', 'C14145', 'C16169', 'C16173', 'C16328', 'C16329', 'C16330', 'C16331', 'C16332', 'C16333', 'C16334',
                            'C16335', 'C16336', 'C16337', 'C16338', 'C16339', 'C16376', 'C16389', 'C16466', 'C16468', 'C16469', 'C16470',
                            'C16471', 'C00246', 'C21925', 'C21926']
    expected_modules_pathways_ids = ['M00013', 'M00032', 'M00085', 'M00086', 'M00087', 'M00088', 'M00095', 'M00113', 'M00373', 'M00374',
                                     'M00375', 'M00376', 'M00849', 'M00861', 'M00878', 'M00957', 'map00061', 'map00062', 'map00071', 'map00280',
                                     'map00310', 'map00360', 'map00362', 'map00380', 'map00410', 'map00592', 'map00620', 'map00627', 'map00630',
                                     'map00640', 'map00642', 'map00650', 'map00720', 'map00900', 'map00907', 'map00930', 'map01040', 'map01100',
                                     'map01110', 'map01120', 'map01200', 'map01212', 'map01220']

    found_reactions = [reaction.id for reaction in sbml_model.getListOfReactions()]
    found_metabolites = [reaction.id for reaction in sbml_model.getListOfSpecies()]

    model_groups = sbml_model.getPlugin("groups")
    modules_pathways_ids = [group.id for group in model_groups.getListOfGroups()]

    assert sorted(found_reactions) == sorted(expected_reactions)
    assert sorted(found_metabolites) == sorted(expected_metabolites)
    assert sorted(modules_pathways_ids) == sorted(expected_modules_pathways_ids)

    shutil.rmtree('test_out')
