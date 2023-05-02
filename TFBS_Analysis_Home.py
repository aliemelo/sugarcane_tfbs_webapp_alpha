import streamlit as st
import pandas as pd

st.set_page_config(layout="wide")

# Add a title and intro text
st.title("TFBS Analysis in sugarcane SP80-3280")
st.markdown("This is a web app to explore TFBSs in the promoter sequences of sugarcane variety SP80-3280.")

st.markdown("### How to use this web app?")
st.markdown("1. Provide a list of target genes in the *Gene Selection* page")
st.markdown("2. Visualise the distribution of motifs in the promoter region of the target genes in the *Visualisation* page")

st.markdown("---")
st.markdown("Please note that this webapp is still under development. Future features that will be implemented include:")
st.markdown("- Selection of experimental data")
st.markdown("- Association of experimental data with the target promoter regions")
st.markdown("- Calculation of candidate promoter architectures")
