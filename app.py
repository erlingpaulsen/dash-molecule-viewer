import dash
from dash import dcc
from dash import html
from dash.dependencies import Input, Output
import dash_bio as dashbio
import dash_bootstrap_components as dbc
import utils

MOLECULE_ID_FILE = 'data/molecule_ids'

app = dash.Dash(
    __name__,
    external_stylesheets=[dbc.themes.DARKLY],
    meta_tags=[
        {"name": "viewport", "content": "width=device-width, height=device-height, initial-scale=1"}
    ]
)
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
df = utils.createDataFrame(mol_ids)
dt_cols = utils.getTableColumns()
tryptamine_table = utils.createTable(df, dt_cols)

page_title = [
    html.H1('Molecular Explorer', className='text-center')
]

page_info = dbc.Card(
    'Use the tabs to explore different molecules in the base classes Tryptamines, Phenetylamines and Lysergamides. Click to select a molecule from the table to interactively explore its 3D molecular structure and known synonyms.',
    body=True,
    class_name='text-white bg-primary border mt-2 fst-italic overflow-auto',
    style={'maxWidth': '90vw', 'maxHeight': '10vh'}
)

molecule_speck_view_toggle = html.Div([
    dbc.Checklist(
        options=[
            {'label': 'Toggle view', 'value': 1}
        ],
        value=[1],
        id="molecule-speck-preset-views",
        switch=True,
        label_class_name='fs-6'
    )
])

molecule_speck_view = html.Div([
    dashbio.Speck(
        id='molecule-speck',
        view={
            'resolution': 300,
            'zoom': 0.085
        }
    )],
    className='overflow-hidden p-0 m-0'
)

molecule_3d_viewer = [
    dbc.Card(
        dbc.CardBody([
            html.H4('3D Viewer', className='card-title text-center'),
            html.H5(id='molecule-speck-title',  className='card-text text-center fst-italic'),
            molecule_speck_view_toggle,
            molecule_speck_view,
        ]),
        class_name='mt-2 border overflow-hidden',
        style={'maxHeight': '40vh', 'maxWidth': '20vw', 'minWidth': '10vw'}
    ),
    dbc.Card(
        dbc.CardBody([
            html.H4('Synonyms', className='card-title text-center'),
            html.Ul(id='molecule-synonyms', className='card-text text-wrap fs-6', style={'maxWidth': '20vw'})
        ]),
        class_name='mt-2 border overflow-auto',
        style={'maxHeight': '30vh', 'maxWidth': '20vw', 'minWidth': '10vw'}
    )
]

tryptamines = [
    dbc.Card(
        dbc.Row(
            [
                dbc.Col([
                    dbc.CardBody([
                        html.H4('Tryptamines', className='card-title'),
                        html.P('Tryptamines are a group of monoamine alkaloids (indolalkylamines), derived from the amino acid tryptophan, that can be found in natural sources including plants, fungi, microbes, and amphibia', className='card-text'),
                    ])
                ], xs=6, md=6, lg=8),
                dbc.Col([
                    dbc.CardImg(
                        src='data:image/png;base64,' + df[df['name'] == 'Tryptamine']['img_b64'],
                        class_name='img-fluid rounded-end'
                    )
                ], xs=6, md=5, lg=3)
            ],
            justify='start'
        )
    ),
    html.Div([
            tryptamine_table
        ],
    )
]

molecule_table_tabs = dbc.Tabs(
    [
        dbc.Tab(tryptamines, label='Tryptamines'),
        dbc.Tab('test', label='Phenetylamines', disabled=True),
        dbc.Tab('test2', label='Lysergamides', disabled=True),
    ]
)

molecule_table_section = dbc.Card(
    dbc.CardBody([
        html.P('Base Classes', className='card-title text-center h3'),
        molecule_table_tabs
    ]),
    class_name='mt-2 border overflow-auto',
    style={'maxHeight': '85vh', 'maxWidth': '70vw'}
)

app.layout = dbc.Container(fluid=True, children=[
    dbc.Row(
        [
            dbc.Col(page_title, width='auto')
        ],
        justify='center',
        align='center',
        class_name='g-2'
    ),
    dbc.Row(
        [
            dbc.Col(page_info, width='auto')
        ],
        justify='center',
        align='center',
        class_name='g-2'
    ),
    dbc.Row(
        [
            dbc.Col(
                molecule_3d_viewer,
                width='auto',
                class_name='g-2'
            ),
            dbc.Col(
                molecule_table_section,
                width='auto',
                class_name='g-2'
            ),
        ],
        class_name='g-2',
        justify='center'
    )
])

@app.callback(
    Output('molecule-speck-title', 'children'),
    Output('molecule-speck', 'data'),
    Output('molecule-synonyms', 'children'),
    [Input('row-{}'.format(str(i)), 'n_clicks') for i in df.index]
)
def update_speck_data(*clicks):
    ctx = dash.callback_context
    if not ctx.triggered:
        i = 0
    else:
        row_id = ctx.triggered[0]['prop_id'].split('.')[0]
        i = int(row_id.split('-')[-1])
    return df.loc[i, 'name'], df.loc[i, 'xyz'], [html.Li(s) for s in df.loc[i, 'synonyms']]

@app.callback(
    Output('molecule-speck', 'presetView'),
    Input('molecule-speck-preset-views', 'value')
)
def update_speck_view(v):
    if len(v) > 0:
        return 'stickball'
    else:
        return 'default'

if __name__ == '__main__':
    app.run_server(debug=True)

