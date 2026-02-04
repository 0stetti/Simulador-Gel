import streamlit as st
import plotly.graph_objects as go
import pandas as pd
import math
from Bio.Seq import Seq
from Bio.Restriction import RestrictionBatch, Analysis, CommOnly
from io import StringIO, BytesIO
from Bio import SeqIO

# --- 1. CONFIGURAÇÃO DA PÁGINA ---
st.set_page_config(
    page_title="BioSpark Studio",
    layout="wide",
    page_icon="🧬",
    initial_sidebar_state="expanded"
)

# --- 2. ESTÉTICA BOHRIUM TURQUESA (CSS HACK) ---
st.markdown("""
<style>
    /* Importar fonte Inter (Padrão moderno) */
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&display=swap');

    /* Fundo Geral com Degradê Turquesa Intenso */
    .stApp {
        background: linear-gradient(180deg, #E0F7FA 0%, #FFFFFF 100%);
        font-family: 'Inter', sans-serif;
    }

    /* Sidebar - Branco Puro e Clean */
    section[data-testid="stSidebar"] {
        background-color: #FFFFFF;
        border-right: 1px solid #E5E7EB;
        box-shadow: 4px 0 15px rgba(0,0,0,0.02);
    }

    /* Títulos e Textos */
    h1, h2, h3 {
        color: #111827 !important;
        font-weight: 700 !important;
        letter-spacing: -0.025em;
    }
    p, label {
        color: #374151;
    }

    /* Transformar os Expanders em "Cards" Flutuantes */
    .stExpander {
        background-color: #FFFFFF;
        border-radius: 12px !important;
        border: 1px solid #E5E7EB !important;
        box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.05), 0 2px 4px -1px rgba(0, 0, 0, 0.03) !important;
        margin-bottom: 1rem;
        transition: all 0.2s ease-in-out;
    }

    .stExpander:hover {
        box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.08) !important;
        border-color: #4F46E5 !important; /* Destaque Roxo ao passar o mouse */
    }

    /* Conteúdo interno do Expander */
    div[data-testid="stExpanderDetails"] {
        background-color: #FFFFFF;
        border-radius: 0 0 12px 12px;
    }

    /* Botões Modernos */
    .stButton > button {
        background-color: #4F46E5; /* Roxo Bohrium */
        color: white;
        border-radius: 8px;
        font-weight: 600;
        border: none;
        padding: 0.5rem 1rem;
        transition: background-color 0.2s;
    }
    .stButton > button:hover {
        background-color: #4338CA;
        color: white;
    }

    /* Inputs e Selectboxes */
    .stTextInput > div > div > input, .stSelectbox > div > div > div {
        background-color: #F9FAFB;
        border-radius: 8px;
        border: 1px solid #D1D5DB;
        color: #111827;
    }

    /* Rodapé */
    .footer {
        width: 100%;
        text-align: center;
        padding-top: 20px;
        padding-bottom: 20px;
        font-size: 12px;
        color: #9CA3AF;
        border-top: 1px solid #E5E7EB;
        margin-top: 30px;
    }
</style>
""", unsafe_allow_html=True)

# --- 3. LÓGICA DO SISTEMA (MANTIDA IGUAL) ---

# Converte enzimas para string
TODAS_ENZIMAS = sorted([str(e) for e in CommOnly])

# Dados de Ladders
LADDERS = {
    "1kb Plus DNA Ladder": [100, 200, 300, 400, 500, 650, 850, 1000, 1650, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 12000],
    "1kb DNA Ladder (Genérico)": [250, 500, 750, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 8000, 10000],
    "100bp DNA Ladder": [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1517, 2017],
    "High Mass": [1000, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 20000, 48500]
}

def processar_upload(input_data):
    try:
        nome_arquivo = input_data.name
        nome_sugerido = nome_arquivo.rsplit('.', 1)[0]

        if nome_arquivo.lower().endswith('.dna'):
            try:
                bytes_io = BytesIO(input_data.getvalue())
                record = SeqIO.read(bytes_io, "snapgene")
                return nome_sugerido, str(record.seq).upper()
            except Exception as e:
                return "Erro", f"Erro .dna: {str(e)}"

        bytes_data = input_data.getvalue()
        try:
            conteudo = bytes_data.decode("utf-8")
        except UnicodeDecodeError:
            conteudo = bytes_data.decode("latin-1")

        if ">" in conteudo:
            try:
                iterator = SeqIO.parse(StringIO(conteudo), "fasta")
                record = next(iterator)
                return record.id if record.id else nome_sugerido, str(record.seq).upper()
            except:
                pass

        linhas = conteudo.splitlines()
        seq_limpa = ""
        for linha in linhas:
            linha = linha.strip()
            if not linha or linha.startswith(">") or linha.startswith(";"): continue
            seq_limpa += linha

        seq_final = "".join(seq_limpa.split()).upper()
        if any(c not in "ATGCNRYKMSWBDHV" for c in seq_final[:100]):
             return "Erro", "Arquivo inválido."

        return nome_sugerido, seq_final

    except Exception as e:
        return "Erro", str(e)

def processar_texto_manual(texto):
    try:
        if ">" in texto:
            iterator = SeqIO.parse(StringIO(texto), "fasta")
            record = next(iterator)
            return record.id, str(record.seq).upper()
        else:
            return "Seq Manual", "".join(texto.split()).upper()
    except:
        return "Erro", ""

def calcular_digestao(sequencia, enzimas, eh_circular):
    if not sequencia or sequencia.startswith("Erro"): return []
    sequencia = "".join([c for c in sequencia if c in "ATGCMRWSYKVHDBN"])
    if not sequencia: return []

    seq_obj = Seq(sequencia)
    tamanho_total = len(seq_obj)

    if eh_circular and not enzimas:
        return [
            (tamanho_total * 1.4, "Nicked (Relaxed)", tamanho_total),
            (tamanho_total * 0.7, "Supercoiled", tamanho_total)
        ]

    if not enzimas:
        return [(tamanho_total, "Linear", tamanho_total)]

    rb = RestrictionBatch(enzimas)
    analise = Analysis(rb, seq_obj, linear=not eh_circular)
    cortes = analise.full()
    locais = sorted(list(set([local for lista in cortes.values() for local in lista])))

    if not locais:
        tipo = "Circ. (Sítio Ausente)" if eh_circular else "Lin. (Não Cortado)"
        return [(tamanho_total, tipo, tamanho_total)]

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

    return [(frag, "Fragmento", frag) for frag in sorted(fragmentos, reverse=True)]

# --- BARRA LATERAL (VISUAL CLEAN) ---
with st.sidebar:
    # Título Moderno e Futurista
    st.markdown("""
    <div style="text-align: center; margin-bottom: 20px;">
        <h1 style="font-family: 'Inter', sans-serif; font-weight: 700; color: #4F46E5; font-size: 28px; letter-spacing: 1px;">
            🧬 BioSpark
        </h1>
    </div>
    """, unsafe_allow_html=True)
    st.markdown("---")

    st.caption("CONFIGURAÇÕES GERAIS")

    # 1. Número de Poços
    num_pocos = st.slider("Número de Poços", 1, 15, 3)

    # 2. Concentração de Agarose
    agarose = st.slider("Agarose (%)", 0.5, 2.0, 1.0, 0.1)

    st.caption("VISUALIZAÇÃO")
    # 3. Estilo Visual (Ordem corrigida)
    estilo_gel = st.selectbox(
        "Tema do Gel",
        ["Profissional (Dark P&B)", "Publicação (Light P&B)", "Neon (Verde/Laranja)"]
    )

    st.markdown("---")

    # --- GUIA DE USO ---
    with st.expander("❓ Guia Rápido"):
        st.markdown("""
        **1. Setup:** Ajuste poços e agarose.
        **2. Amostras:**
        * Faça upload de `.dna` ou `.fasta`.
        * Para plasmídeos, marque **Circular**.
        * Escolha as enzimas.
        **3. Resultado:** O gel é gerado automaticamente.
        """)

    st.markdown("---")
    # --- IDIOMA E CRÉDITOS ---
    st.caption("PREFERÊNCIAS")
    idioma = st.selectbox("Idioma / Language", ["Português", "English"])
    st.markdown(" ")
    st.markdown("""
    <div class="footer">
        <p style="margin: 0; font-weight: 600;">Criado por Elton Ostetti</p>
        <p style="margin: 0;">Laboratório de Biofármacos</p>
        <p style="margin: 0;">Instituto Butantan</p>
    </div>
    """, unsafe_allow_html=True)


# --- CONTEÚDO PRINCIPAL ---

# Cabeçalho Moderno (Corrigido para Digestão Enzimática)
st.markdown("# Simulador de Digestão Enzimática de DNA")
st.markdown("Configure suas amostras abaixo para visualizar o resultado da digestão in silico.")
st.markdown(" ") # Espaço extra

# Estrutura para guardar dados para o relatório
relatorio_dados = []
dados_para_plotar = []
labels_eixo_x = []
nomes_ladders = []

cols = st.columns(2)

# LOOP DOS POÇOS (Agora estilizados como Cards pelo CSS)
for i in range(num_pocos):
    col_atual = cols[i % 2]
    with col_atual:
        # O st.expander agora parece um "Card" branco flutuante
        with st.expander(f"🔹 Poço {i+1}", expanded=(i==0)):
            tipo = st.radio(f"Conteúdo {i+1}:", ["Amostra", "Ladder"], key=f"t_{i}", horizontal=True, label_visibility="collapsed")

            rotulo_padrao = str(i+1)

            if tipo == "Ladder":
                lad = st.selectbox("Selecione o Ladder:", list(LADDERS.keys()), key=f"l_{i}")
                ladder_data = [(tam, "Ladder", tam) for tam in LADDERS[lad]]
                dados_para_plotar.append(ladder_data)

                rotulo_custom = st.text_input("Rótulo no Gel:", value="M", key=f"lbl_{i}")
                labels_eixo_x.append(rotulo_custom)
                nomes_ladders.append(lad)

                # Dados para CSV
                relatorio_dados.append({
                    "Poço": i+1,
                    "Identificação": rotulo_custom,
                    "Tipo": "Ladder",
                    "Detalhes": lad,
                    "Bandas (pb)": "; ".join([str(t) for t in LADDERS[lad]])
                })
            else:
                nomes_ladders.append(None)
                tab_f, tab_t = st.tabs(["📂 Arquivo", "📝 Texto Manual"])
                seq, nome_arquivo = "", ""

                with tab_f:
                    up = st.file_uploader("Upload DNA/Fasta", type=['dna', 'fasta', 'txt', 'fa'], key=f"u_{i}")
                    if up:
                        nome_arquivo, seq = processar_upload(up)
                        if nome_arquivo == "Erro":
                            st.error(seq); seq = ""
                with tab_t:
                    txt = st.text_area("Colar Sequência", height=70, key=f"tx_{i}")
                    if txt and not seq:
                        nome_t, seq_t = processar_texto_manual(txt)
                        if nome_t != "Seq Manual": nome_arquivo = nome_t
                        seq = seq_t

                st.markdown("---")
                c1, c2 = st.columns([1, 2])
                circ = c1.checkbox("Circular?", True, key=f"c_{i}")
                enz = c2.multiselect("Enzimas de Restrição", TODAS_ENZIMAS, key=f"e_{i}")

                val_rotulo = nome_arquivo if nome_arquivo else str(i+1)
                rotulo_custom = st.text_input("Rótulo no Gel:", value=val_rotulo[:12], key=f"lbl_{i}")
                labels_eixo_x.append(rotulo_custom)

                if seq:
                    try:
                        res = calcular_digestao(seq, enz, circ)
                        dados_para_plotar.append(res)

                        fragmentos_str = "; ".join([str(int(b[0])) for b in res])
                        desc_enzimas = ", ".join(enz) if enz else ("Circular Uncut" if circ else "Linear Uncut")
                        relatorio_dados.append({
                            "Poço": i+1,
                            "Identificação": rotulo_custom,
                            "Tipo": "Amostra",
                            "Detalhes": desc_enzimas,
                            "Bandas (pb)": fragmentos_str
                        })

                    except Exception as e:
                        dados_para_plotar.append([])
                        st.error(f"Erro no cálculo: {e}")
                else:
                    dados_para_plotar.append([])
                    relatorio_dados.append({
                        "Poço": i+1,
                        "Identificação": rotulo_custom,
                        "Tipo": "Vazio",
                        "Detalhes": "-",
                        "Bandas (pb)": "-"
                    })

st.markdown(" ") # Espaço antes do gel
st.markdown("### Resultado da Eletroforese")

if any(dados_para_plotar):

    # Configuração de Cores do Gráfico
    if "Neon" in estilo_gel:
        bg_color = '#111827'; text_color = 'white'; color_sample = '#00ff41'; color_ladder = '#ff9900' # Neon otimizado
    elif "Profissional" in estilo_gel:
        bg_color = '#000000'; text_color = 'white'; color_sample = 'white'; color_ladder = 'white'
    else:
        bg_color = 'white'; text_color = 'black'; color_sample = 'black'; color_ladder = 'black'

    min_view = 50 + (100 * (agarose - 0.5))
    max_view = 25000 / (agarose * 0.8)

    fig = go.Figure()

    for i, lista_bandas in enumerate(dados_para_plotar):
        x_center = i + 1
        eh_ladder = (nomes_ladders[i] is not None)
        cor_atual = color_ladder if eh_ladder else color_sample

        if lista_bandas:
             massa_total = sum([b[2] for b in lista_bandas]) if not eh_ladder else 1

        for (tam_aparente, tipo_banda, tam_real) in lista_bandas:
            if tam_aparente < min_view or tam_aparente > max_view: continue

            width = 2; opacity = 0.8
            if eh_ladder:
                if tam_aparente in [3000, 1000, 500]: width = 7; opacity = 1.0
                elif tam_aparente >= 5000: width = 5; opacity = 0.9
                else: width = 3; opacity = 0.7
            else:
                if tipo_banda == "Supercoiled": fracao = 0.7
                elif tipo_banda == "Nicked (Relaxed)": fracao = 0.3
                else: fracao = tam_real / massa_total if massa_total > 0 else 0.5
                width = 3 + (8 * fracao)
                opacity = 0.6 + (0.4 * fracao)

            largura_banda = 0.28

            fig.add_trace(go.Scatter(
                x=[x_center - largura_banda, x_center + largura_banda],
                y=[tam_aparente, tam_aparente],
                mode='lines+markers',
                line=dict(color=cor_atual, width=width),
                marker=dict(color=cor_atual, size=width, symbol='circle'),
                opacity=opacity,
                showlegend=False,
                hoverinfo='text',
                hovertext=f"<b>~{int(tam_aparente)} pb</b><br>Real: {tam_real}<br>{labels_eixo_x[i]}"
            ))

            if eh_ladder:
                fig.add_trace(go.Scatter(
                    x=[x_center - 0.45], y=[tam_aparente], mode="text",
                    text=[str(tam_aparente)], textposition="middle left",
                    textfont=dict(color=text_color, size=10),
                    showlegend=False, hoverinfo='skip'
                ))

    max_range = max(num_pocos, 15) + 0.5

    fig.update_layout(
        plot_bgcolor=bg_color, paper_bgcolor=bg_color,
        height=700, margin=dict(t=40, b=40, l=40, r=40),
        xaxis=dict(
            tickmode='array', tickvals=list(range(1, num_pocos + 1)),
            ticktext=labels_eixo_x,
            tickfont=dict(color=text_color, size=14, family='Arial Black'),
            showgrid=False, zeroline=False, range=[0.2, max_range]
        ),
        yaxis=dict(
            type='log',
            range=[math.log10(min_view), math.log10(max_view)],
            showgrid=False, zeroline=False, showticklabels=False
        )
    )

    # Container visual para o gráfico ficar bonito no tema branco
    st.plotly_chart(fig, use_container_width=True)

    # --- BOTÃO DISCRETO DE EXPORTAÇÃO ---
    with st.expander("📥 Exportar Dados do Relatório"):
        df_resultados = pd.DataFrame(relatorio_dados)
        csv = df_resultados.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="Baixar Tabela (.csv)",
            data=csv,
            file_name='gel_resultado.csv',
            mime='text/csv',
        )

else:
    # Card de aviso vazio bonito
    st.info("👋 Para começar, adicione amostras nos cartões acima.")

# --- RODAPÉ ---
st.markdown("""
<div class="footer">
    <p><b>BioSpark</b> | Instituto Butantan</p>
</div>
""", unsafe_allow_html=True)
