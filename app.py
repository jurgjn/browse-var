
import ast, gzip, random, os, tempfile, time, sqlite3, urllib.request
import matplotlib, matplotlib.colors, matplotlib.pyplot as plt, seaborn as sns, pandas as pd, streamlit as st, streamlit_ext as ste, st_aggrid, py3Dmol, stmol

st.set_page_config(
    page_title='ClinVar variants mapped to (raw) pockets',
    page_icon='🔬',
    layout='wide',
)
#st.cache_resource.clear()

def uf(x):
    return '{:,}'.format(x)

def select_dataframe_row(df_, selected_row_index, height=200):
    gb = st_aggrid.GridOptionsBuilder.from_dataframe(df_)
    gb.configure_selection(selection_mode='single', use_checkbox=True, pre_selected_rows=[ selected_row_index ])
    gb.configure_grid_options(domLayout='normal')
    #gb.configure_pagination()
    #https://github.com/PablocFonseca/streamlit-aggrid/issues/57
    gb.configure_grid_options(onFirstDataRendered=st_aggrid.JsCode("""
    function(e) { 
        e.api.ensureIndexVisible(%d, 'middle');
    }
    """ % (selected_row_index,)).js_code)
    gridOptions = gb.build()
    gridResponse = st_aggrid.AgGrid(df_,
        gridOptions=gridOptions,
        #update_mode=st_aggrid.GridUpdateMode.SELECTION_CHANGED,
        fit_columns_on_grid_load=True,
        height=height,
        width='100%',
        enable_enterprise_modules=False,
        allow_unsafe_jscode=True,
    )
    if not(len(gridResponse['selected_rows']) > 0): time.sleep(5) # Prevent annoying row-not-selected errors during loading
    if len(gridResponse['selected_rows']) == 0: return None
    return gridResponse['selected_rows'][0]

@st.cache_resource
def read_af2_v4_(af2_id):
    url_ = f'https://alphafold.ebi.ac.uk/files/AF-{af2_id}-F1-model_v4.pdb'
    with urllib.request.urlopen(url_) as url:
        return url.read().decode('utf-8')

@st.cache_resource
def read_clinvar_pocket():
    #df_ = pd.read_csv('data/af2.obabel_hxr.autosite.summary.scaled.swiss_coverage.in_clinvar.tsv', sep='\t').query('resid_swissmodel_coverage == "none"').reset_index(drop=True)
    fp_ = 'data/af2.obabel_hxr.autosite.summary.scaled.swiss_coverage.in_clinvar.sqlite'
    with sqlite3.connect(fp_) as con:
        df_ = pd.read_sql(f'select * from pockets', con)
    return df_

st.write(f'# ClinVar variants mapped to (raw) pockets (n={uf(len(read_clinvar_pocket()))})')

col1, col2 = st.columns([0.8, 0.2])
with col1:
    row_displ_ = select_dataframe_row(read_clinvar_pocket()[['uniprot_id', 'struct_id', 'n_resid', 'pocket_id', 'score', 'mean_pLDDT', 'score_integrated', 'score_integrated_scaled', 'resid_swissmodel_coverage', 'CLNDN', 'CLNSIG', 'Amino_acid_position', 'Amino_acid_change']], selected_row_index=0)
    row_index_ = int(row_displ_['_selectedRowNodeInfo']['nodeId'])
    row_ = read_clinvar_pocket().iloc[row_index_].to_dict()
    af2_id_ = row_['uniprot_id']
    struct_id_ = row_['struct_id']
    pocket_id_ = row_['pocket_id']
    #st.write(f'{uf(len(read_clinvar_pocket()))} entries shown')
    st.write(f'{af2_id_} in [UniProt](https://www.uniprot.org/uniprotkb/{af2_id_}/entry) / [AlphaFill](https://alphafill.eu/model?id={af2_id_})')

    st.write(f'## Structure for {af2_id_}')
    pdb_ = read_af2_v4_(af2_id_)

    #cmap_ = sns.color_palette('viridis', as_cmap=True)
    #if not(saliency_) is None:
    #    colors_pocket = {i + 1: matplotlib.colors.to_hex(cmap_(val_)) for i, val_ in saliency_.items() }
    #else:
    #colors_pocket = {}
    colors_pocket = {resid_: '#d55e00' if resid_ == row_['Amino_acid_position'] else 'gray' for resid_ in range(1, 3000)} #https://github.com/jurgjn/relmapping/blob/master/scripts/yarp/yarp.py
    #st.write(row_['Amino_acid_position'])
    #st.write(colors_pocket)
    #colors_pocket = {i + 1: matplotlib.colors.to_hex(cmap_(val_)) for i, val_ in saliency_.items() }

    xyzview = py3Dmol.view()

    # Add structure
    xyzview.addModel(pdb_, format='pdb')
    xyzview.setStyle({'model': 0}, {
        'cartoon': {
            'colorscheme': {
                'prop': 'resi',
                'map': colors_pocket,
            },
            'arrows': True,
        }
    })
    xyzview.addStyle({'model': 0, 'resi': row_['Amino_acid_position']}, {'stick': {'color': '#d55e00'}}) #https://colab.research.google.com/github/pb3lab/ibm3202/blob/master/tutorials/lab02_molviz.ipynb

    # Add pocket surface
    #fp = os.path.join(os.path.expanduser('~/euler-home/project/22.12_pocketomes'), row_['cl_file'])
    #with open(fp) as fh:
    #    pocket_ = fh.read()#.decode('ascii')
    pocket_ = gzip.decompress(row_['cl_str']).decode()

    xyzview.addModel(pocket_, format='pdb')
    xyzview.setStyle({'model': -1}, {})
    xyzview.addSurface(py3Dmol.VDW, {'opacity': 0.5, 'color': 'pink'}, {'model': -1})

    # Back matter
    xyzview.setBackgroundColor('#eeeeee')
    xyzview.zoomTo()
    stmol.showmol(xyzview, height=600, width=600)
    #st.write(cmap_)

with col2:
    st.write(row_)
