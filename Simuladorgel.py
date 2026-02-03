import streamlit as st
import plotly.graph_objects as go
import math # <--- CORREÇÃO: Importado no início
from Bio.Seq import Seq
from Bio.Restriction import RestrictionBatch, Analysis, CommOnly
from io import StringIO
from Bio import SeqIO

# --- CONFIGURAÇÃO INICIAL ---
st.set_page_config(page_title="Simulador de Gel Interativo", layout="wide", page_icon="🧬")

# Carrega TODAS as enzimas comerciais
TODAS_ENZIMAS = sorted(list(CommOnly))

# Dados de Ladders
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
    
    if not enzimas: return [tamanho_total]
    
    rb = RestrictionBatch(enzimas)
    analise = Analysis(rb, seq_obj, linear=not eh_circular)
    cortes = analise.full()
    locais = sorted(list(set([local for lista in cortes.values() for local in lista])))
    
    if not locais: return [tamanho_total]
        
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
st.title("🧪 Simulador de Eletroforese Interativo")
st.markdown("Passe o mouse sobre as bandas para ver o tamanho exato em pares de base (pb).")

with st.sidebar:
    st.header("Configurações")
    num_pocos = st.slider("Número de Poços", 1, 15, 10)
    st.divider()
    inverter_cores = st.toggle("Inverter Cores (Modo Impressão)", value=False)
    st.caption("Nota: A espessura da banda indica a massa relativa de DNA.")

dados_para_plotar = []
labels_eixo_x = []
nomes_ladders = [] 
detalhes_hover = [] 

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
                nomes_ladders.append(lad)
                detalhes_hover.append(lad)
            else:
                nomes_ladders.append(None)
                tab_f, tab_t = st.tabs(["Arquivo", "Texto"])
                seq, nome = "", f"{i+1}"
                
                with tab_f:
                    up = st.file_uploader("FASTA", key=f"u_{i}")
                    if up: nome, seq = processar_fasta(up, True)
                with tab_t:
                    txt = st.text_area("Seq", height=70, key=f"tx_{i}")
                    if txt and not seq: 
                        nome_t, seq_t = processar_fasta(txt, False)
                        if nome_t != "Seq Manual": nome = nome_t
                        seq = seq_t
                
                c1, c2 = st.columns(2)
                circ = c1.checkbox("Circular?", True, key=f"c_{i}")
                enz = c2.multiselect("Enzimas", TODAS_ENZIMAS, key=f"e_{i}")
                
                info_texto = f"Circular: {circ}<br>Enzimas: {', '.join(enz) if enz else 'Nenhuma'}"
                detalhes_hover.append(info_texto)

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
    bg_color = 'white' if inverter_cores else '#1e1e1e'
    line_color = 'black' if inverter_cores else 'white'
    text_color = 'black' if inverter_cores else 'white'
    
    fig = go.Figure()

    for i, bandas in enumerate(dados_para_plotar):
        x_center = i + 1
        eh_ladder = (labels_eixo_x[i] == "M")
        ladder_name = nomes_ladders[i]
        
        massa_total = sum(bandas) if bandas and not eh_ladder else 1
        
        for tam in bandas:
            # --- LÓGICA DE ESPESSURA E OPACIDADE ---
            width = 2
            opacity = 0.8
            
            if eh_ladder:
                # Destaque para bandas de referência
                if tam in [3000, 1000, 500]: 
                    width = 4
                    opacity = 1.0
                elif tam >= 5000:
                    width = 3
                    opacity = 0.9
                else:
                    width = 1.5
                    opacity = 0.7
            else:
                # Física: Intensidade baseada na massa
                fracao_massa = tam / massa_total
                width = 2 + (6 * fracao_massa) 
                opacity = 0.5 + (0.5 * fracao_massa)

            # --- DESENHO DA BANDA ---
            fig.add_trace(go.Scatter(
                x=[x_center - 0.35, x_center + 0.35],
                y=[tam, tam],
                mode='lines',
                line=dict(color=line_color, width=width),
                opacity=opacity,
                name=f"Poço {labels_eixo_x[i]}",
                showlegend=False,
                hoverinfo='text',
                hovertext=f"<b>Tamanho:</b> {tam} pb<br><b>Poço:</b> {labels_eixo_x[i]}<br>{detalhes_hover[i]}"
            ))

            # --- RÓTULOS LATERAIS DO LADDER ---
            if eh_ladder:
                fig.add_annotation(
                    x=x_center - 0.5,
                    y=tam, # CORREÇÃO: Valor direto, Plotly lida com o log
                    text=f"{tam}",
                    showarrow=False,
                    xanchor="right",
                    font=dict(color=text_color, size=9),
                    yshift=0
                )
                # Linha guia fina
                fig.add_shape(
                    type="line",
                    x0=x_center - 0.5, x1=x_center - 0.35,
                    y0=tam, y1=tam,
                    line=dict(color=text_color, width=0.5),
                    opacity=0.3
                )

    # --- CONFIGURAÇÃO DOS EIXOS E LAYOUT ---
    fig.update_layout(
        plot_bgcolor=bg_color,
        paper_bgcolor=bg_color,
        height=600,
        margin=dict(t=30, b=30, l=60, r=30),
        xaxis=dict(
            tickmode='array',
            tickvals=list(range(1, num_pocos + 1)),
            ticktext=labels_eixo_x,
            tickfont=dict(color=text_color, size=14, family='Arial Black'),
            showgrid=False,
            zeroline=False,
            range=[0.5, num_pocos + 0.5]
        ),
        yaxis=dict(
            type='log',
            # Definimos o range explicitamente invertido: log(max) em cima, log(min) em baixo?
            # Plotly log axis é meio chato. Geralmente [log10(max), log10(min)] inverte.
            range=[math.log10(20000), math.log10(80)], 
            showgrid=False,
            zeroline=False,
            title="pb",
            titlefont=dict(color=text_color, size=14),
            tickfont=dict(color=text_color),
            showticklabels=False 
        )
    )
    
    # Label "pb" no topo
    fig.add_annotation(
        x=-0.05, y=1, xref="paper", yref="paper",
        text="pb", showarrow=False,
        font=dict(color=text_color, size=14, family="Arial Black")
    )

    st.plotly_chart(fig, use_container_width=True)

else:
    st.info("Adicione amostras para gerar o gel.")
