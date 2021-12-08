import requests
import json
import os
import periodictable as pt
from rdkit import Chem
from rdkit.Chem import Draw

periodic_table = pt.core.default_table()

TIMEOUT = 5
MOLECULE_ID_FILE = 'data/molecule_ids'
JSON_RAW_PATH = 'data/json_raw/'
JSON_PATH = 'data/json/'
XYZ_RAW_PATH = 'data/xyz_raw/'
XYZ_PATH = 'data/xyz/'
IMG_PATH = 'data/img/'

'''--------------------------------------------------------
                   LOADIND IDS
--------------------------------------------------------'''
mol_ids = []
with open(MOLECULE_ID_FILE, 'r') as file:
    for line in file:
        mol_ids.append(line.strip())

'''--------------------------------------------------------
         DOWNLOADING DATA FROM CHEMSPIDER
--------------------------------------------------------'''
'''for name in mol_ids.keys():
    mol_id = mol_ids[name]['chemspider_id']
    img_url = 'https://www.chemspider.com/ImagesHandler.ashx?id={0}&w=100&h=100'.format(mol_id)
    try:
        r = requests.get(img_url, timeout=TIMEOUT)
    except Exception as e:
        print('Error: Request failed for SpiderChem ID {0}. Continuing.\nURL:{1}'.format(mol_id, img_url))
        print(e)
        continue
    if r.status_code == 200:
        with open('data/img/{0}.png'.format(mol_id), 'wb') as img_file:
            img_file.write(r.content)
    else:
        print('Error: could not download image for SpiderChem ID {0} from URL:\n{1}'.format(mol_id, img_url))

    mol_url = 'https://www.chemspider.com/FilesHandler.ashx?type=str&striph=yes&id={0}'.format(mol_id)
    try:
        r = requests.get(mol_url, timeout=TIMEOUT)
    except Exception as e:
        print('Error: Request failed for SpiderChem ID {0}. Continuing.\nURL:{1}'.format(mol_id, mol_url))
        print(e)
        continue
    if r.status_code == 200:
        with open('data/mol/{0}.mol'.format(mol_id), 'wb') as mol_file:
            mol_file.write(r.content)
    else:
        print('Error: could not download molecule for SpiderChem ID {0} from URL:\n{1}'.format(mol_id, mol_url))'''

'''--------------------------------------------------------
         DOWNLOADING DATA FROM PUBCHEM
--------------------------------------------------------'''

for mol_id in mol_ids:
    file_name = '{0}.json'.format(mol_id)
    if not file_name in os.listdir(JSON_RAW_PATH):
        json_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{0}/JSON/?response_type=save'.format(mol_id)
        try:
            r = requests.get(json_url)
        except Exception as e:
            print('Error: Request failed for PubChem ID {0}. Continuing.\nURL:{1}'.format(mol_id, json_url))
            print(e)
            continue
        if r.status_code == 200:
            with open(JSON_RAW_PATH + file_name, 'wb') as json_file:
                json_file.write(r.content)
        else:
            print('Error: could not download JSON for PubChem ID {0} from URL:\n{1}'.format(mol_id, json_url))
    else:
        print(JSON_RAW_PATH + file_name + ' already exists. Skipping ...')

    if not file_name in os.listdir(XYZ_RAW_PATH):
        json_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/{0}/record/JSON/?record_type=3d&response_type=save'.format(mol_id)
        try:
            r = requests.get(json_url)
        except Exception as e:
            print('Error: Request failed for PubChem ID {0}. Continuing.\nURL:{1}'.format(mol_id, json_url))
            print(e)
            continue
        if r.status_code == 200:
            with open(XYZ_RAW_PATH + file_name, 'wb') as json_file:
                json_file.write(r.content)
        else:
            print('Error: could not download XYZ JSON for PubChem ID {0} from URL:\n{1}'.format(mol_id, json_url))
    else:
        print(XYZ_RAW_PATH + file_name + ' already exists. Skipping ...')

'''--------------------------------------------------------
               GENERATING PROCESSED JSON
--------------------------------------------------------'''
def processJson(mol_id):
    pubchem_url = 'https://pubchem.ncbi.nlm.nih.gov/compound/{0}'.format(mol_id)
    name = None
    smiles = None
    formula = None
    mol_weight = None
    mol_weight_unit = None
    atom_count = None
    bond_count = None
    xyz = None

    file_name = '{0}.json'.format(mol_id)
    img_file = '{0}.png'.format(mol_id)

    with open(JSON_RAW_PATH + file_name, 'r') as json_raw:
        json_raw = json.load(json_raw)
        for i in json_raw:
            if i == 'Record':
                name = json_raw[i]['RecordTitle']
                for j in json_raw[i]:
                    if j == 'Section':
                        for k in json_raw[i][j]:
                            if k['TOCHeading'] == 'Names and Identifiers':
                                for l in k['Section']:
                                    if l['TOCHeading'] == 'Computed Descriptors':
                                        for m in l['Section']:
                                            if m['TOCHeading'] == 'Canonical SMILES':
                                                smiles = m['Information'][0]['Value']['StringWithMarkup'][0]['String']
                                    if l['TOCHeading'] == 'Molecular Formula':
                                        formula = l['Information'][0]['Value']['StringWithMarkup'][0]['String']
                            if k['TOCHeading'] == 'Chemical and Physical Properties':
                                for l in k['Section']:
                                    if l['TOCHeading'] == 'Computed Properties':
                                        for m in l['Section']:
                                            if m['TOCHeading'] == 'Molecular Weight':
                                                mol_weight = m['Information'][0]['Value']['StringWithMarkup'][0]['String']
                                                mol_weight_unit = m['Information'][0]['Value']['Unit']

    with open(XYZ_RAW_PATH + file_name, 'r') as xyz_raw:
        xyz_raw = json.load(xyz_raw)
        element_numbers = xyz_raw['PC_Compounds'][0]['atoms']['element']
        x = xyz_raw['PC_Compounds'][0]['coords'][0]['conformers'][0]['x']
        y = xyz_raw['PC_Compounds'][0]['coords'][0]['conformers'][0]['y']
        z = xyz_raw['PC_Compounds'][0]['coords'][0]['conformers'][0]['z']
        atom_count = len(element_numbers)
        bond_count = len(xyz_raw['PC_Compounds'][0]['bonds']['aid1'])

        xyz_string = str(atom_count) + '\n\n'
        for i, e in enumerate(element_numbers):
            xyz_string += str(periodic_table[e]) + '\t' + str(x[i]) + '\t' + str(y[i]) + '\t' + str(z[i]) + '\n'
        xyz = xyz_string

    json_processed = {
        'pubchem_id': mol_id,
        'pubchem_url': pubchem_url,
        'name': name,
        'formula': formula,
        'smiles': smiles,
        'mol_weight': str(mol_weight) + ' ' + mol_weight_unit,
        'atom_count': atom_count,
        'bond_count': bond_count,
        'xyz': xyz
    }
    with open(JSON_PATH + file_name, 'w', encoding='utf-8') as json_processed_file:
        json.dump(json_processed, json_processed_file, ensure_ascii=False, indent=4)

    mol = Chem.MolFromSmiles(smiles)
    Draw.MolToFile(mol, IMG_PATH + img_file, size=(100, 100), imageType='PNG')

for mol_id in mol_ids:
    processJson(mol_id)

