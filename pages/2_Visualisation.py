import streamlit as st
import io
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib import gridspec
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial import distance

sns.set_style('white')


# Functions ------------------------------------------------ #
def motifdf4plotting(df, conn):
    motif_lst = df['motif'].to_list()
    motif_lst = tuple(motif_lst)
    tfbs_gene_association = pd.read_sql(f"SELECT * FROM tair_association WHERE gene_id in {motif_lst}", conn)
    tmp_dict = tfbs_gene_association.set_index('gene_id').T.to_dict('list')
    tmp = df.values.tolist()
    motifs = []
    for item in tmp:
        motif_start = int(item[2])
        motif_end = int(item[3])
        motif_pos = (motif_start + motif_end) // 2
        motif_pos_corrected = - motif_pos
        motif = tmp_dict[item[1]]
        motif = ''.join(motif)
        aux_lst = [item[0], motif, motif_pos_corrected]
        motifs.append(aux_lst)
    
    outdf = pd.DataFrame(motifs, columns=['seq', 'motif', 'motif_location'])  
    return outdf


def reorder_df_by_list(len_df, motif_df, lst, denovo_df=None):
    # reorder len_df
    len_df2 = len_df[len_df['seq'].isin(lst)]
    len_df2 = len_df2.set_index(['seq'])
    len_df3 = len_df2.reindex(lst)
    len_df3 = len_df3.reset_index()
    len_df3['seq_id'] = len_df3.index
    len_reordered = len_df3.set_index(['seq'])

    # reorder motif_df
    motif_reordered = motif_df[motif_df['seq'].isin(lst)]
    motif_reordered['seq_id'] = motif_reordered.seq.map(len_df3.set_index('seq')['seq_id'].to_dict())
    
    # reorder denovo_df if present
    if denovo_df is None:
        return len_reordered, motif_reordered
    else:
        denovo_reordered = denovo_df[denovo_df['seq'].isin(lst)]
        denovo_reordered['seq_id'] = denovo_reordered.seq.map(len_df3.set_index('seq')['seq_id'].to_dict())
        return len_reordered, motif_reordered, denovo_reordered



# Main --------------------------------------------------- #
# Load data
conn = st.session_state['conn']
genes = st.session_state['genes']
len_df = pd.read_sql(f"SELECT * FROM promoter_len WHERE seq in {genes}", conn)
motif_df = st.session_state['my_genes_motifs_df']
motif4plotting_df = motifdf4plotting(motif_df, conn)

# Reorder genes in the dataframes
len_df2, motif_df2 = reorder_df_by_list(len_df, motif4plotting_df, genes)

# Sliders
fig_height = st.sidebar.slider("Set figure height", min_value=5, max_value=20, value=10)
ncolslegend = st.sidebar.slider("Number of columns in legend", min_value=1, max_value=5)

# Plot figure
fig = plt.figure(figsize=(20,fig_height))
gs = gridspec.GridSpec(nrows=1, ncols=3, width_ratios=[1, 1, 1],
                           wspace=0.0, hspace=0.025, top=0.95, bottom=0.05,
                           left=0.3, right=0.845)
# Create subplots
tmp_ax = plt.subplot(gs[0,1]) # space for gene names
tmp_ax.axis('off')
plt.sca(tmp_ax)
plt.xticks([])
plt.yticks([])
sns.despine(left=True, bottom=True, ax=tmp_ax)

## subplots for plotting actual data
ax = plt.subplot(gs[0,0]) # motifs
plt.sca(ax)
plt.xlabel("Distance from TSS (bp)")
plt.ylabel(" ")
sns.despine(left=True, bottom=True, ax=ax)

legend_ax = plt.subplot(gs[0,2]) # legend
legend_ax.axis('off')

# plot motifs
## plot line that represents the promoter sequence
for entry, index in zip(len_df2.index.values, len_df2.seq_id):
    seq_len = len_df2.iloc[index]['len']
    x1 = [-seq_len, seq_len-seq_len]
    y1 = [entry, entry]
    ax.plot(x1, y1, linewidth=0.5, color='black', zorder=1)

max_len = len_df2['len'].max()

## plot motifs
g1 = sns.scatterplot(x='motif_location', y='seq_id', data=motif_df2, hue='motif',
                    s=100, zorder=5, ax=ax, marker='o')
ax.set_xlim([-(max_len+50), 0])
n_genes = len(len_df2['seq_id'].to_list())
y_value = n_genes - 0.5
ax.set_ylim([-0.5, y_value])
ax.get_legend().remove()
ax.yaxis.set_label_position('right')
ax.yaxis.tick_right()

# plot motif legend in another subplot
label_params = ax.get_legend_handles_labels()
legend_ax.legend(*label_params, loc="upper left", ncol=ncolslegend, title="Motifs")

fig.align_labels()
st.pyplot(fig)

# save fig
fn = "motifs_visualisation.png"
img = io.BytesIO()
plt.savefig(img, format='png', dpi=600, bbox_inches='tight')

btn = st.sidebar.download_button(label="Download image", data=img, file_name=fn, mime="image/png")