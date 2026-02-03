import streamlit as st

st.set_page_config(
    page_title="Biofármacos Suite",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 Biofármacos Suite: Ferramentas de Biologia Molecular")
st.markdown("### Bem-vindo à central de ferramentas in silico do Instituto Butantan.")

col1, col2 = st.columns(2)

with col1:
    st.info("### 🧪 1. Digestão In Silico")
    st.write("""
    Simule géis de agarose com precisão física.
    * Suporte a SnapGene (.dna) e FASTA.
    * Visualização realista de plasmídeos (Supercoiled/Nicked).
    * Exportação de relatórios.
    """)

with col2:
    st.info("### 🧬 2. Mapa de Plasmídeo (Beta)")
    st.write("""
    Visualize mapas circulares e lineares.
    * Identificação automática de ORFs.
    * Anotação de enzimas de restrição.
    * Exportação de imagens para publicação.
    """)

st.divider()
st.success("👈 Selecione a ferramenta desejada no menu lateral para começar.")

# Rodapé
st.markdown("""
<div style="position: fixed; bottom: 0; left: 0; width: 100%; text-align: center; padding: 10px; background-color: #0E1117; color: #666; font-size: 12px; border-top: 1px solid #333;">
    <p><b>Elton Ostetti</b> | Laboratório de Biofármacos - Instituto Butantan</p>
</div>
""", unsafe_allow_html=True)
