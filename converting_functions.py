import requests
from PIL import Image
from io import BytesIO
import rdkit.Chem as Chem
from rdkit.Chem import Descriptors
from libchebipy import ChebiEntity
from molmass import Formula
from pubchempy import *
def convert_name2inchikey(name):
    r = requests.get('http://cts.fiehnlab.ucdavis.edu/rest/convert/Chemical%20Name/InChIKey/'+ name)
    if r.text.find('"results":[]') == -1:
        inchikey = r.text[r.text.find('result')+11:r.text.find('result')+11 +27]
        return(inchikey)
    else:
        return('Not Found')
def convert_name2chebi(name):
    r = requests.get('http://cts.fiehnlab.ucdavis.edu/rest/convert/Chemical%20Name/ChEBI/'+ name)
    if r.text.find('CHEBI:') == -1:
        return('Not Found')
    else:
        Chebi = r.text[r.text.find('CHEBI:')+6:r.text.find('"',r.text.find('CHEBI:'))]
    return(Chebi)
def convert_inchikey2name(inchikey):
    r = requests.get('https://cts.fiehnlab.ucdavis.edu/rest/convert/InChIKey/Chemical%20Name/' + inchikey)
    name = r.text[r.text.find('result')+10:r.text.find(',',r.text.find('result')+10)].replace('"','')
    return(name)
def convert_inchikey2smiles(inchikey):
    r = requests.get('https://cactus.nci.nih.gov/chemical/structure/InChIKey='+ inchikey +'/smiles')
    if r.text.find('\n') != -1:
        smiles = r.text[:r.text.find('\n')]
    else:
        smiles = r.text
    return(smiles)
def convert_inchikey2formula(inchikey):
    r = requests.get('https://cactus.nci.nih.gov/chemical/structure/InChIKey='+ inchikey +'/Formula')
    if r.text.find('\n') != -1:
        smiles = r.text[:r.text.find('\n')]
    else:
        smiles = r.text
    return(smiles)
def convert_inchi2pubchem(inchi):
    cpd_id = str(get_compounds(inchi,'inchi')[0])
    cpd_id = cpd_id[cpd_id.find('(')+1:cpd_id.find(')')]
    return(cpd_id)
def get_TPSA(smiles):
    m = Chem.MolFromSmiles(smiles)
    tpsa = Descriptors.TPSA(m)
    return(tpsa)
def convert_name2smiles(name):
    inchikey = convert_name2inchikey(name)
    smiles = convert_inchikey2smiles(inchikey)
    return(smiles)
def convert_Chebi2formula(ChebiID):
    CE = ChebiEntity(ChebiID)
    formula = CE.get_formula()
    return(formula)
def AddHs_2smiles(smiles):
    m = Chem.MolFromSmiles(smiles)
    m2 =Chem.AddHs(m)
    new_smiles = Chem.MolToSmiles(m2)
    return(new_smiles)
def inchikey_2_png(inchikey):
    r = requests.get('http://cactus.nci.nih.gov/chemical/structure/InChIKey='+ inchikey+'/image')
    i = Image.open(BytesIO(r.content))
    i.show()
def get_monoisotopic_mass_smiles(smiles):
    mass = Descriptors.ExactMolWt(Chem.MolFromSmiles(smiles))
    return(mass)
def get_monoisotopic_mass_formula(formula):
    f= Formula(formula)
    monoisotope_mass = f.isotope.mass
    return(monoisotope_mass)
def name_2_Monoisotopic_mass_smiles(name):
    smiles = convert_name2smiles(name)
    if smiles != '<h1>Page not found (404)</h1>':
        mass = Descriptors.ExactMolWt(Chem.MolFromSmiles(smiles))
        return(mass)
    else:
        return('Not Found')
def name_2_Monoisotopic_mass_chebi(name):
    chebi = convert_name2chebi(name)
    if chebi != 'Not Found':
        formula = convert_Chebi2formula(chebi)
        mass = get_monoisotopic_mass_formula(formula)
        return(mass)
    else:
        return('Not Found')
def convert_name2formula(name):
    chebi = convert_name2chebi(name)
    if chebi != 'Not Found':
        formula = convert_Chebi2formula(chebi)
        return(formula)
    else:
        return('Not Found')
def create_name_formula_mass_tsv_names(names,outputfile):
    file = open(outputfile,'w')
    file.write('Name\tFormula\tExact Mass\n')
    for name in names:
        formula = 'Not Found'
        mass = "Not Found"
        formula = convert_name2formula(name)
        if formula != 'Not Found':
            mass = get_monoisotopic_mass_formula(formula)
        else:
            inchikey = convert_name2inchikey(name)
            if inchikey != 'Not Found':
                formula = convert_inchikey2formula(inchikey)
                mass = get_monoisotopic_mass_formula(formula)
            else:
                pass
        
        file.write(name + '\t' + formula + '\t' + mass + '\n')
    file.close()
    print(outputfile + ' was created!')
def create_name_formula_mass_tsv_inchis(inchis,outputfile):
    file = open(outputfile,'w')
    file.write('Name\tFormula\tExact Mass\n')
    for inchi in inchis:
        name = ''
        formula = 'Not Found'
        mass = "Not Found"
        pubchemid = convert_inchi2pubchem(inchi)
        c  = Compound.from_cid(pubchemid)
        if len(c.synonyms) > 0:            
            name = c.synonyms[0]
        else:
            name = c.iupac_name
        formula = c.molecular_formula
        mass = get_monoisotopic_mass_formula(formula)
        file.write(name + '\t' + formula + '\t' + str(mass) + '\n')
    file.close()
    print(outputfile + ' was created!')
