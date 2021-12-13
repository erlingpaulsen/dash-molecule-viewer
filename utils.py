import io
import base64
import json
import re
import pandas as pd
from dash_bio_utils import xyz_reader
from dash import dcc
from dash import html
import dash_bootstrap_components as dbc
import periodictable as pt
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image

periodic_table = pt.core.default_table()
tryptamine = Chem.MolFromSmiles('C1=CC=C2C(=C1)C(=CN2)CCN')
lysergamide = Chem.MolFromSmiles('CN1CC(C=C2C1CC3=CNC4=CC=CC2=C34)C(=O)N')
phenetylamine = Chem.MolFromSmiles('C1=CC=C(C=C1)CCN')

def generateMolImg(mol_id, smiles, template):
    # Align molecule to template and generate 2D skeletal image from molecule SMILES

    mol = Chem.MolFromSmiles(smiles)

    # Align molecule to the base class template molecule
    AllChem.Compute2DCoords(template)
    AllChem.Compute2DCoords(mol)
    AllChem.GenerateDepictionMatching2DStructure(mol, template)

    # Draw molecule
    d = rdMolDraw2D.MolDraw2DCairo(1000, 1000)
    d.drawOptions().minFontSize = 72
    d.drawOptions().maxFontSize = 72
    d.drawOptions().bondLineWidth = 10
    d.drawOptions().setBackgroundColour((1, 1, 1, 1))
    d.DrawMolecule(mol)
    d.FinishDrawing()

    # Convert to RGBA and make treansparent backgroun and adopt for "dark mode"
    img = Image.open(io.BytesIO(d.GetDrawingText()))
    img = img.convert('RGBA')
    data = img.getdata()
    newData = []
    t = 180
    for item in data:
        if item[0] > t and item[1] > t and item[2] > t:
            newData.append((255, 255, 255, 0))
        elif item[0] < t and item[1] < t and item[2] < t:
            newData.append((200, 200, 200, 255))
        else:
            newData.append((item[0] + 80, item[1] + 80, item[2] + 80, item[3]))
    img.putdata(newData)

    # Save as thumbnail
    img.thumbnail((300, 300))
    img.save('data/img/{0}.png'.format(mol_id), 'PNG')

    # Return image and base64 encoded string
    buf = io.BytesIO()
    img.save(buf, format='PNG')
    img_b64 = base64.b64encode(buf.getvalue()).decode('UTF-8')
    return img_b64

def getMolClass(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol.HasSubstructMatch(lysergamide):
        mol_class = 'Lysergamide'
        template = lysergamide
    elif mol.HasSubstructMatch(tryptamine):
        mol_class = 'Tryptamine'
        template = tryptamine
    elif mol.HasSubstructMatch(phenetylamine):
        mol_class = 'Phenetylamine'
        template = phenetylamine
    else:
        mol_class = 'Not found'
        template = None
    return mol_class, template

def subscriptFormula(formula):
    match = re.split('(\d+)', formula)
    return [html.Sub(c) if c.isnumeric() else c for c in match]

def parsePubChemJson(mol_id):
    synonyms = []
    with open('data/json_raw/{0}.json'.format(mol_id), 'r') as json_raw:
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
                                    elif l['TOCHeading'] == 'Molecular Formula':
                                        formula = l['Information'][0]['Value']['StringWithMarkup'][0]['String']
                                    elif l['TOCHeading'] == 'Synonyms':
                                        for m in l['Section']:
                                            if m['TOCHeading'] == 'MeSH Entry Terms' or m['TOCHeading'] == 'Depositor-Supplied Synonyms':
                                                synonyms += [n['String'] for n in m['Information'][0]['Value']['StringWithMarkup']]
                            if k['TOCHeading'] == 'Chemical and Physical Properties':
                                for l in k['Section']:
                                    if l['TOCHeading'] == 'Computed Properties':
                                        for m in l['Section']:
                                            if m['TOCHeading'] == 'Molecular Weight':
                                                mol_weight = m['Information'][0]['Value']['StringWithMarkup'][0][
                                                    'String']
                                                mol_weight_unit = m['Information'][0]['Value']['Unit']
    if len(synonyms) == 0:
        synonyms = ['None']

    return name, synonyms, smiles, formula, mol_weight, mol_weight_unit

def parsePubChemXyz(mol_id):
    with open('data/xyz_raw/{0}.json'.format(mol_id), 'r') as xyz_raw:
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
        xyz = xyz_reader.read_xyz(datapath_or_datastring=xyz_string, is_datafile=False)
    return xyz, atom_count, bond_count

def createJson(mol_id):
    name, synonyms, smiles, formula, mol_weight, mol_weight_unit = parsePubChemJson(mol_id)
    xyz, atom_count, bond_count = parsePubChemXyz(mol_id)
    mol_class, template = getMolClass(smiles)
    img_b64 = generateMolImg(mol_id, smiles, template)

    json_processed = {
        'pubchem_id': mol_id,
        'pubchem_url': 'https://pubchem.ncbi.nlm.nih.gov/compound/{0}'.format(mol_id),
        'name': name,
        'class': mol_class,
        'synonyms': synonyms,
        'formula': formula,
        'smiles': smiles,
        'mol_weight': mol_weight,
        'mol_weight_unit': mol_weight_unit,
        'atom_count': atom_count,
        'bond_count': bond_count,
        'xyz': xyz,
        'img': img_b64
    }

    with open('data/json/{0}.json'.format(mol_id), 'w', encoding='utf-8') as f:
        json.dump(json_processed, f, ensure_ascii=False, indent=4)

def createDataFrame(mol_ids):
    jsons = [json.load(open('data/json/{0}.json'.format(mol_id))) for mol_id in mol_ids]
    data = {
        'pubchem_id': [],
        'pubchem_url': [],
        'pubchem_link': [],
        'name': [],
        'class': [],
        'synonyms': [],
        #'long_name': [],
        #'base_class': [],
        #'sub_class': [],
        'formula': [],
        'smiles': [],
        'atom_count': [],
        'bond_count': [],
        'mol_weight': [],
        'xyz': [],
        'img': [],
        'img_b64': [],
    }
    for j in jsons:
        data['pubchem_id'].append(j['pubchem_id'])
        data['pubchem_url'].append(j['pubchem_url'])
        #data['pubchem_link'].append('<p style="text-align: center;"><a href="' + j['pubchem_url'] + '" target="_blank">' + j['pubchem_id'] + '</a></p>')
        data['pubchem_link'].append(dcc.Link(href=j['pubchem_url'], children=j['pubchem_id'], target='_blank', className='link-info'))
        data['name'].append(j['name'])
        data['class'].append(j['class'])
        data['synonyms'].append(j['synonyms'])
        #data['long_name'].append(j['long_name'])
        #data['base_class'].append(j['base_class'])
        #data['sub_class'].append(j['sub_class'])
        #data['formula'].append(subscriptFormula(j['formula']))
        data['formula'].append(html.P(subscriptFormula(j['formula'])))
        data['smiles'].append(j['smiles'])
        data['atom_count'].append(j['atom_count'])
        data['bond_count'].append(j['bond_count'])
        data['mol_weight'].append(j['mol_weight'] + ' ' + j['mol_weight_unit'])
        data['xyz'].append(j['xyz'])
        #data['img'].append("<img src='data:image/png;base64," + j['img'] + "' class='img-fluid'>")
        data['img'].append(html.Img(className='img-fluid pe-0', src='data:image/png;base64,' + j['img'], style={'maxWidth': '8vw'}))
        data['img_b64'].append(j['img'])

    df = pd.DataFrame(data=data)
    df.sort_values(by=['atom_count', 'bond_count'], inplace=True)
    df.index = range(len(df))
    return df

def getTableColumns():
    return [
        'pubchem_link',
        'name',
        'class',
        'formula',
        'atom_count',
        'bond_count',
        'mol_weight',
        'img'
    ]

def getPrettyName(name):
    if name == 'pubchem_link':
        return 'PubChem ID'
    elif name == 'name':
        return 'Common Name'
    elif name == 'class':
        return 'Base Class'
    elif name == 'formula':
        return 'Formula'
    elif name == 'smiles':
        return 'SMILES'
    elif name == 'atom_count':
        return '# Atoms'
    elif name == 'bond_count':
        return '# Bonds'
    elif name == 'mol_weight':
        return 'Molecular Weight'
    elif name == 'img':
        return '2D Skeletal'
    else:
        return 'NA'

def createTable(df, cols):
    table_header = [
        html.Thead(html.Tr([html.Th(getPrettyName(c)) for c in cols], className='fs-6 bold text-wrap table-primary'))
    ]
    rows = []
    for i in df.index:
        rows.append(
            html.Tr(
                [html.Td(df.loc[i, c], className='text-wrap', style={'maxWidth': '8vw'}) for c in cols],
                id='row-{}'.format(str(i)),
                className='fs-6'
            )
        )
    table_body = [html.Tbody(rows)]
    table = dbc.Table(
        table_header + table_body,
        bordered=True,
        hover=True,
        responsive=True,
        striped=True,
        id='molecule-table'
    )
    return table
