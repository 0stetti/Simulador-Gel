import streamlit as st
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio.Restriction import RestrictionBatch, Analysis, CommOnly
from io import StringIO
from Bio import SeqIO

# --- CONFIGURAÇÃO INICIAL ---
st.set_page_config(page_title="Simulador de Gel (Massa Proporcional)", layout="wide", page_icon="🧬")

# Carrega TODAS as enzimas comerciais
TODAS_ENZIMAS = sorted(list(CommOnly))

# Dados de Ladders (Marcadores)
# Nota: O 1kb Plus oficial da Invitrogen geralmente vai de 100bp a 12000bp.
# As bandas de referência (mais fortes) costumam ser 1000, 1650 (as vezes) e 3000.
LADDERS = {
    "1kb Plus DNA Ladder": [100, 200, 300, 400, 500, 650, 850, 1000, 1650, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 12000],
    "1kb DNA Ladder (Genérico)": [250, 500, 750, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 8000, 10000],
    "100bp DNA Ladder": [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1517, 2017],
    "High Mass": [1000, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 20000, 48500]
}

def processar_fasta(input_data, is_file=False):
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
st.title("🧪 Simulador de Eletroforese (Física Real)")
st.markdown("Simulação com conservação de massa: fragmentos maiores são mais brilhantes.")

with st.sidebar:
    st.header("Configurações")
    num_pocos = st.slider("Número de Poços", 1, 15, 10)
    st.divider()
    inverter_cores = st.toggle("Inverter Cores (Modo Impressão)", value=False)
    st.caption("Nota: No modo escuro, bandas grossas simulam maior massa de DNA.")

dados_para_plotar = []
labels_eixo_x = []
ladder_names_used = []

cols = st.columns(2)

for i in range(num_pocos):
    col_atual = cols[i % 2]
    with col_atual:
        with st.expander(f"Poço {i+1}", expanded=(i==0)):
            tipo = st.radio(f"Conteúdo {i+1}:", ["Amostra", "Ladder"], key=f"t_{i}", horizontal=True)
            
            if tipo == "Ladder":
                lad = st.selectbox("Ladder:", list(LADDERS.keys()), key=f"l_{i}")
                dados_para_plotar.append(LADDERS[lad])
                labels_eixo_x.append("M")
                ladder_names_used.append(lad)
            else:
                ladder_names_used.append(None)
                tab_f, tab_t = st.tabs(["Arquivo", "Texto"])
                seq, nome = "", f"{i+1}"
                
                with tab_f:
                    up = st.file_uploader("FASTA", key=f"u_{i}")
                    if up: 
                        nome, seq = processar_fasta(up, True)
                with tab_t:
                    txt = st.text_area("Seq", height=70, key=f"tx_{i}")
                    if txt and not seq: 
                        nome_t, seq_t = processar_fasta(txt, False)
                        if nome_t != "Seq Manual": nome = nome_t
                        seq = seq_t
                
                c1, c2 = st.columns(2)
                circ = c1.checkbox("Circular?", True, key=f"c_{i}")
                enz = c2.multiselect("Enzimas", TODAS_ENZIMAS, key=f"e_{i}")
                
                if seq:
                    try:
                        res = calcular_digestao(seq, enz, circ)
                        dados_para_plotar.append(res)
                        labels_eixo_x.append(str(i+1))
                    except:
                        dados_para_plotar.append([])
                        labels_eixo_x.append("Erro")
                else:
                    dados_para_plotar.append([])
                    labels_eixo_x.append(str(i+1))

st.divider()

if any(dados_para_plotar):
    # Cores
    bg = 'white' if inverter_cores else '#1e1e1e'
    band_color = 'black' if inverter_cores else 'white'
    text_color = 'black' if inverter_cores else 'white'

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_facecolor(bg)
    fig.patch.set_facecolor(bg)
    
    for spine in ax.spines.values(): spine.set_visible(False)

    for i, bandas in enumerate(dados_para_plotar):
        x = i + 1
        eh_ladder = (labels_eixo_x[i] == "M")
        ladder_name = ladder_names_used[i]
        
        # Se for amostra, precisamos calcular a massa total para a proporção
        massa_total = sum(bandas) if bandas and not eh_ladder else 1
        
        for tam in bandas:
            # LÓGICA DE ESPESSURA (INTENSIDADE)
            lw = 2.0
            alpha = 0.8
            
            if eh_ladder:
                # Lógica para Ladders (Bandas de Referência Fixas)
                # 1kb Plus e Generico costumam ter 3000 e 1000 fortes
                if tam in [3000, 1000, 500]: 
                    lw = 4.0
                    alpha = 1.0
                elif tam >= 5000:
                    lw = 2.5
                    alpha = 0.9
                else:
                    lw = 1.5
                    alpha = 0.7
            else:
                # LÓGICA DE AMOSTRA (Baseada na Massa)
                # Fração da massa total que este fragmento representa
                fracao_massa = tam / massa_total
                
                # O fragmento maior brilha mais. 
                # Base minima 1.0, Maximo adiciona até 4.0
                lw = 1.5 + (4.5 * fracao_massa) 
                
                # A transparência também ajuda no efeito visual
                # Fragmentos muito pequenos ficam mais "fantasmas"
                alpha = 0.6 + (0.4 * fracao_massa)

            ax.hlines(y=tam, xmin=x-0.35, xmax=x+0.35, colors=band_color, linewidth=lw, alpha=alpha)
            
            # Rótulos do Ladder (Apenas numeros inteiros "pb")
            if eh_ladder:
                ax.text(x-0.5, tam, f"{tam}", color=text_color, fontsize=8, ha='right', va='center')
                # Linha guia fina
                ax.hlines(y=tam, xmin=x-0.5, xmax=x-0.35, colors=text_color, linewidth=0.5, alpha=0.3)

    # Eixos
    ax.set_yscale('log')
    ax.set_ylim(15000, 80) # De 80bp a 15kb
    
    ax.set_xticks(range(1, num_pocos + 1))
    ax.set_xticklabels(labels_eixo_x, color=text_color, fontsize=11, weight='bold')
    
    # Label "pb" no topo
    ax.set_ylabel("pb", color=text_color, fontsize=12, weight='bold', rotation=0, ha='right')
    ax.yaxis.set_label_coords(-0.06, 0.96)
    
    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.tick_params(axis='x', colors=text_color)
    
    st.pyplot(fig)
else:
    st.info("Adicione amostras ou ladders para visualizar.")
