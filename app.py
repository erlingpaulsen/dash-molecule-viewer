import dash
from dash import dcc
from dash import html
from dash.dependencies import Input, Output
import dash_bio as dashbio
from dash_bio_utils import xyz_reader
import dash_bootstrap_components as dbc
from periodictable import *
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdmolfiles import MolToXYZBlock

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.DARKLY])

colors = {
    'background': '#111111',
    'text': 'darkgrey'
}

serotonin = Chem.MolFromMolFile('mol_files/5013.mol')
if serotonin is None:
    print('Error')
else:
    serotonin = Chem.AddHs(serotonin)
    AllChem.EmbedMolecule(serotonin, randomSeed=0xf00d)

data = xyz_reader.read_xyz(datapath_or_datastring=MolToXYZBlock(serotonin), is_datafile=False)

app.layout = html.Div(style={'backgroundColor': colors['background']}, children=[
    dashbio.Speck(
        data=data,
        view={
            'resolution': 400,
            'brightness': 0.45,
            'ao': 0.75,
            'aoRes': 1,
            'outline': 0,
            'atomScale': 0.2,
            'relativeAtomScale': 0.7,
            'bonds': True,
            'zoom': 0.1
        }
    )
])

if __name__ == '__main__':
    app.run_server(debug=True)
