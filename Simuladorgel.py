import streamlit as st
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio.Restriction import RestrictionBatch, Analysis, CommOnly
from io import StringIO
from Bio import SeqIO

# --- CONFIGURAÇÃO INICIAL ---
st.set_page_config(page_title="Simulador de Gel (Full Enzymes)", layout="wide", page_icon="🧬")

# 1. Carrega TODAS as enzimas comerciais do Biopython
# CommOnly é um conjunto que contém apenas enzimas disponíveis comercialmente
TODAS_ENZIMAS = sorted(list(CommOnly))

# Dados de Ladders (Marcadores)
LADDERS = {
    "1kb Plus (Invitrogen)": [100, 200, 300, 400, 500, 650, 850, 1000, 1650, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 12000],
    "100bp Ladder": [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1517, 2017],
    "High Mass": [1000, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 20000, 48500]
}

def processar_fasta(input_data, is_file=False):
    """Lê FASTA de string ou arquivo e retorna (nome, sequencia_string)."""
    try:
        if is_file:
            stringio = StringIO(input_data.getvalue().decode("utf-8"))
            iterator = SeqIO.parse(stringio, "fasta")
        else:
            if ">" in input_data:
                from io import StringIO
                iterator = SeqIO.parse(StringIO(input_data), "fasta")
            else:
                return "Seq Manual", "".join(input_data.split()).upper()
        
        record = next(iterator)
        return record.id, str(record.seq).upper()
    except Exception:
        return "Erro", ""

def calcular_digestao(sequencia, enzimas, eh_circular):
    """Calcula fragmentos para uma sequência."""
    if not sequencia: return []
    
    seq_obj = Seq(sequencia)
    tamanho_total = len(seq_obj)
    
    if not enzimas:
        return [tamanho_total]
    
    rb = RestrictionBatch(enzimas)
    analise = Analysis(rb, seq_obj, linear=not eh_circular)
    
    cortes = analise.full()
    locais = sorted(list(set([local for lista in cortes.values() for local in lista])))
    
    if not locais:
        return [tamanho_total]
        
    fragmentos = []
    if not eh_circular:
        prev = 0
        for cut in locais:
            fragmentos.append(cut - prev)
            prev = cut
        fragmentos.append(tamanho_total - prev)
    else:
        if len(locais) == 1:
            fragmentos.append(tamanho_total)
        else:
            for i in range(len(locais)-1):
                fragmentos.append(locais[i+1] - locais[i])
            fragmentos.append((tamanho_total - locais[-1]) + locais[0])
            
    return sorted(fragmentos, reverse=True)

# --- INTERFACE ---
st.title("🧪 Simulador de Eletroforese Profissional")
st.markdown(f"Base de dados carregada: **{len(TODAS_ENZIMAS)} enzimas comerciais** disponíveis.")

with st.sidebar:
    st.header("Configurações do Gel")
    num_pocos = st.slider("Número de Poços", min_value=1, max_value=20, value=10)
    st.info("Dica: A lista de enzimas agora é pesquisável. Digite o nome para encontrar.")

dados_para_plotar = []
labels_eixo_x = []

st.subheader("Configuração dos Poços")
cols = st.columns(2)

for i in range(num_pocos):
    col_atual = cols[i % 2]
    with col_atual:
        with st.expander(f"Poço {i+1}", expanded=(i==0)):
            tipo_conteudo = st.radio(f"Tipo do Poço {i+1}:", ["Amostra (DNA)", "Ladder (Marcador)"], key=f"tipo_{i}", horizontal=True)
            
            if tipo_conteudo == "Ladder (Marcador)":
                ladder_nome = st.selectbox("Selecione o Ladder:", list(LADDERS.keys()), key=f"lad_{i}")
                dados_para_plotar.append(LADDERS[ladder_nome])
                labels_eixo_x.append("M") 
            else:
                tab_file, tab_text = st.tabs(["Arquivo", "Texto"])
                seq_final = ""
                nome_seq = f"P{i+1}"
                
                with tab_file:
                    arq = st.file_uploader("Upload FASTA", type=['fasta', 'txt'], key=f"up_{i}")
                    if arq:
                        n, s = processar_fasta(arq, is_file=True)
                        nome_seq, seq_final = n, s
                
                with tab_text:
                    txt = st.text_area("Cole Sequência", height=70, key=f"txt_{i}")
                    if txt and not seq_final:
                        n, s = processar_fasta(txt, is_file=False)
                        if n != "Seq Manual": nome_seq = n
                        seq_final = s
                
                c1, c2 = st.columns(2)
                with c1:
                    is_circular = st.checkbox("Circular?", value=True, key=f"circ_{i}")
                with c2:
                    # AQUI ESTÁ A MUDANÇA: Usamos a lista completa TODAS_ENZIMAS
                    # O Streamlit permite digitar para buscar dentro do multiselect
                    enzimas = st.multiselect("Enzimas:", options=TODAS_ENZIMAS, key=f"enz_{i}")
                
                if seq_final:
                    try:
                        bandas = calcular_digestao(seq_final, enzimas, is_circular)
                        dados_para_plotar.append(bandas)
                        labels_eixo_x.append(str(i+1)) 
                    except Exception as e:
                        dados_para_plotar.append([])
                        labels_eixo_x.append("Erro")
                else:
                    dados_para_plotar.append([])
                    labels_eixo_x.append(str(i+1))

st.divider()
st.subheader("Resultado da Eletroforese")

if any(dados_para_plotar):
    # --- VISUALIZAÇÃO ESTILO REFERÊNCIA ---
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.set_facecolor('#222222')
    fig.patch.set_facecolor('#222222')
    
    for spine in ax.spines.values():
        spine.set_visible(False)

    for idx_poco, bandas in enumerate(dados_para_plotar):
        x_pos = idx_poco + 1
        eh_ladder = labels_eixo_x[idx_poco] == "M"
        cor_banda = 'white'
        
        for tamanho in bandas:
            ax.hlines(y=tamanho, xmin=x_pos-0.3, xmax=x_pos+0.3, colors=cor_banda, linewidth=3)
            
            if eh_ladder:
                label_texto = f"{tamanho/1000:.1f}" if tamanho >= 1000 else f"{tamanho}"
                ax.text(x_pos-0.5, tamanho, label_texto, color='white', fontsize=9, ha='right', va='center')
                ax.hlines(y=tamanho, xmin=x_pos-0.5, xmax=x_pos-0.3, colors='white', linewidth=0.5)

    ax.set_yscale('log')
    # Ajuste fino da escala para parecer com a imagem de referência (bandas maiores em cima)
    ax.set_ylim(20000, 100) 
    
    ax.set_xlim(0, num_pocos + 1)
    ax.set_xticks(range(1, num_pocos + 1))
    ax.set_xticklabels(labels_eixo_x, color='white', fontsize=12, weight='bold')
    
    ax.set_ylabel("Kb", color='white', fontsize=12, weight='bold', rotation=0, ha='right')
    ax.yaxis.set_label_coords(-0.05, 0.95)
    
    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.tick_params(axis='x', colors='white')
    ax.grid(False)
    st.pyplot(fig)
else:
    st.warning("Preencha pelo menos um poço.")
