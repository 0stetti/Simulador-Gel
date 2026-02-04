import streamlit as st
from dna_features_viewer import BiopythonTranslator, CircularGraphicRecord
from Bio import SeqIO
from io import BytesIO, StringIO
import matplotlib.pyplot as plt

st.set_page_config(page_title="Mapa de Plasmídeo", layout="wide", page_icon="🧬")

st.title("🧬 Visualizador de Mapas de Plasmídeo")

# Upload
uploaded_file = st.file_uploader("Carregue seu arquivo (.gb, .dna, .fasta)", type=["gb", "dna", "fasta"])

if uploaded_file:
    # Processamento simples (similar ao do gel)
    nome = uploaded_file.name
    record = None
    
    try:
        if nome.endswith(".dna"):
            record = SeqIO.read(BytesIO(uploaded_file.getvalue()), "snapgene")
        elif nome.endswith(".gb"):
            record = SeqIO.read(StringIO(uploaded_file.getvalue().decode("utf-8")), "genbank")
        elif nome.endswith(".fasta"):
            record = SeqIO.read(StringIO(uploaded_file.getvalue().decode("utf-8")), "fasta")
            # Fasta não tem features, então é sem graça, mas funciona
    except Exception as e:
        st.error(f"Erro ao ler arquivo: {e}")

    if record:
        st.subheader(f"Mapa: {record.id}")
        
        # Opções de Visualização
        col1, col2 = st.columns([1, 3])
        
        with col1:
            st.info("Configurações")
            show_seq = st.toggle("Mostrar Sequência (Texto)", False)
            map_width = st.slider("Tamanho da Imagem", 6, 15, 8)
            
        with col2:
            # A Mágica do dna_features_viewer
            class CustomTranslator(BiopythonTranslator):
                def compute_feature_label(self, feature):
                    if feature.type == 'restriction_site': return None
                    return BiopythonTranslator.compute_feature_label(self, feature)

            translator = CustomTranslator()
            graphic_record = translator.translate_record(record)
            
            # Desenha Circular
            fig, ax = plt.subplots(figsize=(map_width, map_width))
            ax.set_aspect("equal") # Garante que é um círculo
            
            circular_rec = CircularGraphicRecord(graphic_record.sequence_length, graphic_record.features)
            circular_rec.plot(ax=ax, figure_width=map_width)
            
            st.pyplot(fig)
            
        if show_seq:
            st.text_area("Sequência", str(record.seq))

else:
    st.info("Faça upload de um arquivo GenBank (.gb) ou SnapGene (.dna) para ver as anotações.")
