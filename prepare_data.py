import requests
import os
import utils

TIMEOUT = 10
MOLECULE_ID_FILE = 'data/molecule_ids'
JSON_RAW_PATH = 'data/json_raw/'
XYZ_RAW_PATH = 'data/xyz_raw/'

'''--------------------------------------------------------
                   LOADIND IDS
--------------------------------------------------------'''
mol_ids = []
with open(MOLECULE_ID_FILE, 'r') as file:
    for line in file:
        mol_ids.append(line.strip())

'''--------------------------------------------------------
         DOWNLOADING DATA FROM PUBCHEM
--------------------------------------------------------'''

for mol_id in mol_ids:
    file_name = '{0}.json'.format(mol_id)
    if not file_name in os.listdir(JSON_RAW_PATH):
        json_url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{0}/JSON/?response_type=save'.format(mol_id)
        try:
            r = requests.get(json_url, timeout=TIMEOUT)
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
            r = requests.get(json_url, timeout=TIMEOUT)
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

for mol_id in mol_ids:
    utils.createJson(mol_id)