import streamlit as st
import pandas as pd
from sqlalchemy import create_engine
import psycopg2


# Functions ------------------------------------------------ #
@st.cache_resource
def init_connection():
    username = st.secrets["postgres"]["user"]
    password = st.secrets["postgres"]["password"]
    host = st.secrets["postgres"]["host"]
    database = st.secrets["postgres"]["dbname"]
    port = st.secrets["postgres"]["port"]
    return create_engine(f"postgresql+psycopg2://{username}:{password}@{host}:{port}/{database}")

conn = init_connection()
st.session_state["conn"] = conn

@st.cache_data
def convert_df(df):
    return df.to_csv(sep="\t", index=False, header=False).encode('utf-8')

def create_motif_df(df):
    # gene, motif, start, end , strand
    tmp_dict = df.set_index('gene').T.to_dict('list')
    tmp_list = []
    for gene in tmp_dict:
        motifs = ''.join(tmp_dict[gene]).split(";")
        for m in motifs:
            m = m.split(":")
            aux_lst = [gene, m[0], m[1], m[2], m[3]]
            tmp_list.append(aux_lst)
    out_df = pd.DataFrame(tmp_list, columns=['gene', 'motif', 'start', 'end', 'strand'])
    return out_df


def callback():
    st.session_state.button_clicked2 = True

if "button_clicked2" not in st.session_state:
    st.session_state.button_clicked2 = False

# CSS to inject contained in a string
hide_table_row_index = """
            <style>
            thead tr th:first-child {display:none}
            tbody th {display:none}
            </style>
            """

# Inject CSS with Markdown
st.markdown(hide_table_row_index, unsafe_allow_html=True)

# Main --------------------------------------------------- #
# Upload list of genes
st.markdown("# Gene Selection")
option = st.radio("Input gene list by either uploading a txt file or typing in the genes:", 
                      options=('Upload file', 'Type in genes'))

gene_list = []

if option == 'Upload file':
    uploaded_file = st.file_uploader("")
    if uploaded_file is not None:
        for line in uploaded_file:
            line = line.strip()
            line = line.decode('UTF-8')
            line = line.replace('scga7_', '')
            gene_list.append(line)

        genes = tuple(gene_list)
        sql_query = pd.read_sql(f"SELECT * FROM tfbs_planttfdb WHERE gene in {genes}", conn)

        motif_df = create_motif_df(sql_query)
        
        tsv = convert_df(motif_df)
        st.sidebar.markdown("**Download output**")
        st.sidebar.download_button(label="Download data as TSV", data=tsv, mime="text/tsv", file_name="genes_motifs_loc.tsv")
        
        st.table(motif_df)

        st.session_state["genes"] = genes
        st.session_state["my_genes_motifs_df"] = motif_df
        
else:
    genes = st.text_area("Insert a list of genes, one per line:")
    genes = genes.split()
    for gene in genes:
        gene = gene.replace('scga7_', '')
        gene_list.append(gene)
    if st.button('Submit', on_click=callback) or st.session_state.button_clicked2:
        genes = tuple(gene_list)
        sql_query = pd.read_sql(f"SELECT * FROM tfbs_planttfdb WHERE gene in {genes}", conn)

        motif_df = create_motif_df(sql_query)

        st.table(motif_df)

        tsv = convert_df(motif_df)

        st.session_state["genes"] = genes
        st.session_state["my_genes_motifs_df"] = motif_df

        st.sidebar.markdown("**Download output**")
        st.sidebar.download_button(label="Download data as TSV", data=tsv, mime="text/tsv", file_name="genes_motifs_loc.tsv", on_click=callback)