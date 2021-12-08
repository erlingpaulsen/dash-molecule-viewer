import dash
from dash import dcc
from dash import html
from dash.dependencies import Input, Output
import dash_bio as dashbio
from dash_bio_utils import xyz_reader
import dash_bootstrap_components as dbc
import json
import base64

MOLECULE_ID_FILE = 'data/molecule_ids'
JSON_PATH = 'data/json/'
IMG_PATH = 'data/img/'

FULL = 12
HALF = 6
THIRD = 4
FOURTH = 3
SIXTH = 2

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.DARKLY])
server = app.server

'''--------------------------------------------------------
                   LOADIND IDS
--------------------------------------------------------'''
mol_ids = []
with open(MOLECULE_ID_FILE, 'r') as file:
    for line in file:
        mol_ids.append(line.strip())

'''--------------------------------------------------------
          LOADING DATA AND BUILDING LAYOUT
--------------------------------------------------------'''
molecule_list_table_rows = []
molecule_speck = []
molecule_view_button_group = html.Div([
    dbc.RadioItems(
        id="molecule-speck-preset-views",
        className="btn-group",
        inputClassName="btn-check",
        labelClassName="btn btn-outline-primary",
        labelCheckedClassName="active",
        options=[
            {'label': 'Ball and stick', 'value': 'stickball'},
            {'label': 'Default', 'value': 'default'}
        ],
        value='stickball',
    )],
    className="radio-group",
)

for i, mol_id in enumerate(mol_ids):
    json_filename = '{0}.json'.format(mol_id)
    img_filename = '{0}.png'.format(mol_id)
    encoded_image = base64.b64encode(open(IMG_PATH + img_filename, 'rb').read())
    with open(JSON_PATH + json_filename, 'r') as f:
        json_file = json.load(f)
        xyz_data = xyz_reader.read_xyz(datapath_or_datastring=json_file['xyz'], is_datafile=False)
        molecule_list_table_row = html.Tr([
            html.Td(dcc.Link(json_file['pubchem_id'], href=json_file['pubchem_url'], target='_blank')),
            html.Td(json_file['name']),
            html.Td(json_file['formula']),
            html.Td(json_file['atom_count']),
            html.Td(json_file['bond_count']),
            html.Td(json_file['mol_weight']),
            html.Td(html.Img(src='data:image/png;base64,{}'.format(encoded_image.decode())))
        ])
        molecule_list_table_rows.append(molecule_list_table_row)

        if i == 0:
            molecule_speck.append(html.H4(json_file['name'], style={'textAlign': 'center'}))
            molecule_speck.append(molecule_view_button_group)
            molecule_speck.append(dashbio.Speck(
                id='molecule-speck',
                data=xyz_data,
                view={
                    'resolution': 700,
                    'brightness': 0.45,
                    'ao': 0.75,
                    'aoRes': 1,
                    'outline': 0,
                    'atomScale': 0.2,
                    'relativeAtomScale': 0.7,
                    'bonds': True,
                    'zoom': 0.1
                }
            ))

title = dbc.Col(html.H1('Molecular Explorer', style={'textAlign': 'center'}), width=FULL)
sub_titles = [
    dbc.Col(html.H3('List of Molecules', style={'textAlign': 'center'}), width=HALF),
    dbc.Col(html.H3('3D Viewer', style={'textAlign': 'center'}), width=HALF)
]
molecule_list_table_header = [html.Thead(html.Tr([
    html.Th('PubChem ID'),
    html.Th('Name'),
    html.Th('Formula'),
    html.Th('Atom Count'),
    html.Th('Bond Count'),
    html.Th('Molecular Weight'),
    html.Th('2D Skeletal')
]))]
molecule_list_table_body = [html.Tbody(molecule_list_table_rows)]
molecule_list_table = dbc.Table(molecule_list_table_header + molecule_list_table_body, bordered=True, hover=True)

app.layout = dbc.Container(fluid=True, children=[
    dbc.Row(title, justify='center', align='center'),
    dbc.Row(sub_titles, justify='center', align='center'),
    dbc.Row([
        dbc.Col(molecule_speck, width=HALF),
        dbc.Col(molecule_list_table, width=HALF)
    ])
])



'''@app.callback(
    Output(component_id='my-output', component_property='children'),
    Input(component_id='my-input', component_property='value')
)
def update_3d_viewer(input_value):
    return 'Output: {}'.format(input_value)'''


@app.callback(
    Output('molecule-speck', 'presetView'),
    Input('molecule-speck-preset-views', 'value')
)
def update_speck_view(preset_name):
    return preset_name

if __name__ == '__main__':
    app.run_server(debug=True)
