import streamlit as st
import plotly.graph_objects as go
import pandas as pd
import math
import re
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

# --- 2. SISTEMA DE TRADUÇÃO ---
TEXTS = {
    "header_title": { "PT": "Simulador de Gel de Agarose", "EN": "Agarose Gel Simulator" },
    "header_sub": { "PT": "Ferramenta in silico para Digestão Enzimática e PCR.", "EN": "In silico tool for Enzymatic Digestion and PCR." },
    "sidebar_config": { "PT": "CONFIGURAÇÕES", "EN": "SETTINGS" },
    "sidebar_wells": { "PT": "Número de Poços", "EN": "Number of Wells" },
    "sidebar_agarose": { "PT": "Agarose (%)", "EN": "Agarose (%)" },
    "sidebar_visual": { "PT": "VISUALIZAÇÃO", "EN": "VISUALIZATION" },
    "sidebar_theme": { "PT": "Tema do Gel", "EN": "Gel Theme" },
    "guide_title": { "PT": "Guia Rápido", "EN": "Quick Guide" },
    "guide_content": {
        "PT": """
        **🧬 Funcionalidades:**
        * **Digestão:** Simula o corte com enzimas.
        * **PCR:** Simula amplificação (suporta overhangs).
        
        **📂 Arquivos:** .dna, .fasta, .txt
        
        **🛠️ Como Usar:**
        1. Ajuste **Poços** e **Agarose**.
        2. Escolha **Amostra**, **PCR** ou **Ladder**.
        3. Faça upload ou cole a sequência.
        """,
        "EN": """
        **🧬 Features:**
        * **Digestion:** Simulates enzyme cuts.
        * **PCR:** Simulates amplification (supports overhangs).
        
        **📂 Files:** .dna, .fasta, .txt
        
        **🛠️ How to Use:**
        1. Set **Wells** and **Agarose**.
        2. Select **Sample**, **PCR** or **Ladder**.
        3. Upload or paste sequence.
        """
    },
    "well_title": { "PT": "Poço", "EN": "Well" },
    "opt_sample": { "PT": "Digestão", "EN": "Digestion" },
    "opt_ladder": { "PT": "Ladder", "EN": "Ladder" },
    "opt_pcr": { "PT": "PCR", "EN": "PCR" },
    "sel_ladder": { "PT": "Selecione o Ladder:", "EN": "Select Ladder:" },
    "label_gel": { "PT": "Rótulo:", "EN": "Label:" },
    "tab_file": { "PT": "📂 Upload Arquivo", "EN": "📂 Upload File" },
    "tab_text": { "PT": "📝 Digitar/Colar", "EN": "📝 Type/Paste" },
    "upload_label": { "PT": "Upload DNA", "EN": "Upload DNA" },
    "paste_label": { "PT": "Sequência", "EN": "Sequence" },
    "check_circular": { "PT": "Circular?", "EN": "Circular?" },
    "sel_enzymes": { "PT": "Enzimas", "EN": "Enzymes" },
    "pcr_fwd": { "PT": "Primer Forward", "EN": "Forward Primer" },
    "pcr_rev": { "PT": "Primer Reverse", "EN": "Reverse Primer" },
    "result_title": { "PT": "Resultado da Eletroforese", "EN": "Electrophoresis Result" },
    "export_expander": { "PT": "Exportar Dados", "EN": "Export Data" },
    "btn_download": { "PT": "Baixar .csv", "EN": "Download .csv" },
    "empty_msg": { "PT": "Para começar, adicione amostras nos cartões acima.", "EN": "To start, add samples in the cards above." },
    "created_by": { "PT": "Desenvolvido por", "EN": "Developed by" },
    "lab_name": { "PT": "Laboratório de Biofármacos", "EN": "Biopharmaceuticals Laboratory" },
    "institute": { "PT": "Instituto Butantan", "EN": "Butantan Institute" },
    "pref_lang": { "PT": "Idioma / Language", "EN": "Language" },
    "report_bug": { "PT": "🐛 Reportar Problema", "EN": "🐛 Report Bug" },
    "warn_multiple": { "PT": "⚠️ Múltiplos sítios de ligação!", "EN": "⚠️ Multiple binding sites!" },
    "warn_no_product": { "PT": "Nenhum produto (Verifique orientação)", "EN": "No product (Check orientation)" },
}

# --- 3. ESTILO CSS (TURQUESA + COMPACTO) ---
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');

    .stApp {
        background: linear-gradient(180deg, #F0F9FF 0%, #FFFFFF 100%);
        font-family: 'Inter', sans-serif;
    }

    section[data-testid="stSidebar"] {
        background-color: #E0F7FA;
        border-right: 1px solid #B2EBF2;
    }

    div[data-baseweb="slider"] div[class*="StyledThumb"] {
        background-color: #0F766E !important;
        border-color: #0F766E !important;
    }
    div[data-baseweb="slider"] div[class*="StyledTrack"] > div {
        background-color: #0F766E !important;
    }
    div[data-baseweb="slider"] div[class*="StyledTrack"] {
        background-color: #B2EBF2 !important;
    }

    h1, h2, h3 {
        color: #0F172A !important;
        font-weight: 700 !important;
        letter-spacing: -0.02em;
    }
    
    .stExpander {
        background-color: #FFFFFF;
        border-radius: 6px !important;
        border: 1px solid #E2E8F0 !important;
        box-shadow: 0 1px 2px 0 rgba(0, 0, 0, 0.05) !important;
        margin-bottom: 0.5rem;
    }
    
    .stExpander:hover {
        border-color: #0F766E !important;
        box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1) !important;
    }

    div[data-testid="stExpanderDetails"] {
        padding: 0.5rem !important;
    }

    .stButton > button {
        background-color: #0F766E;
        color: white;
        border-radius: 6px;
        font-weight: 500;
        border: none;
    }
    .stButton > button:hover {
        background-color: #0d6e66;
        color: white;
    }
    
    span[data-baseweb="checkbox"] div {
        background-color: #0F766E !important;
    }
    
    button[data-baseweb="tab"] {
        font-size: 14px !important;
        font-weight: 500 !important;
    }

    /* RODAPÉ DA SIDEBAR COMPACTO */
    .sidebar-footer {
        margin-top: 20px;
        padding-top: 15px;
        border-top: 1px solid #B2EBF2;
        font-size: 11px;
        color: #111827; /* Preto suave para leitura */
        line-height: 1.6;
    }
    .sidebar-footer strong {
        color: #0F766E; /* Destaque sutil no nome */
    }
    .bug-report {
        font-size: 11px;
        color: #64748B;
        text-decoration: none;
        display: block;
        margin-top: 8px;
    }
    .bug-report:hover {
        color: #0F766E;
        text-decoration: underline;
    }
    
    .warning-text { color: #DC2626; font-weight: bold; font-size: 12px; margin: 2px 0; }
    .error-text { color: #B91C1C; font-size: 12px; margin: 2px 0; }
    
    /* Footer Principal (Página) */
    .footer {
        width: 100%;
        text-align: center;
        padding: 20px 0;
        font-size: 11px;
        color: #64748B;
        border-top: 1px solid #CBD5E1;
        margin-top: 40px;
        opacity: 0.8;
    }
</style>
""", unsafe_allow_html=True)

# --- 4. BACKEND ---

TODAS_ENZIMAS = sorted([str(e) for e in CommOnly])

LADDERS = {
    "1kb Plus DNA Ladder": [100, 200, 300, 400, 500, 650, 850, 1000, 1650, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 12000],
    "1kb DNA Ladder (Genérico)": [250, 500, 750, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 8000, 10000],
    "100bp DNA Ladder": [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1517, 2017],
    "High Mass": [1000, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 20000, 48500]
}

def clean_sequence(seq):
    if not seq: return ""
    return re.sub(r'[^a-zA-Z]', '', seq).upper()

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
        
        seq_final = clean_sequence(seq_limpa)
        if len(seq_final) > 0 and any(c not in "ATGCNRYKMSWBDHV" for c in seq_final[:100]): 
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
            return "Seq Manual", clean_sequence(texto)
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
        return [(tamanho_total, "Uncut", tamanho_total)]
        
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

def smart_pcr_search(template, fwd_raw, rev_raw, eh_circular):
    fwd = clean_sequence(fwd_raw)
    rev = clean_sequence(rev_raw)
    template = template.upper()
    diag = {'fwd_found': False, 'rev_found': False, 'products': 0}
    
    if len(fwd) < 10 or len(rev) < 10: return [], diag

    SEED = 12
    
    # 1. Busca FORWARD (Sempre 5'->3' Sense)
    fwd_seed = fwd[-SEED:]
    fwd_matches = [m.start() for m in re.finditer(fwd_seed, template)]
    if fwd_matches: diag['fwd_found'] = True
    
    # 2. Busca REVERSE (Inteligente)
    rev_matches = []
    
    # Hipótese A: Input 5'->3' (padrão)
    rev_seed_std = rev[-SEED:]
    rev_rc_std = str(Seq(rev_seed_std).reverse_complement())
    matches_std = [m.start() for m in re.finditer(rev_rc_std, template)]
    
    # Hipótese B: Input 3'->5' (visual)
    rev_seed_inv = rev[:SEED]
    rev_compl_inv = str(Seq(rev_seed_inv).complement())
    matches_inv = [m.start() for m in re.finditer(rev_compl_inv, template)]
    
    if len(matches_inv) > 0: rev_matches = matches_inv
    elif len(matches_std) > 0: rev_matches = matches_std
        
    if rev_matches: diag['rev_found'] = True
    
    produtos = []
    
    # 3. Calcular Produtos
    for f_pos in fwd_matches:
        f_end = f_pos + len(fwd_seed) 
        
        for r_pos in rev_matches:
            if r_pos > f_pos:
                dist = r_pos - f_end
                if dist >= 0:
                    produtos.append(len(fwd) + len(rev) + dist)
            elif eh_circular and r_pos < f_pos:
                dist = (len(template) - f_end) + r_pos
                produtos.append(len(fwd) + len(rev) + dist)
                
    diag['products'] = len(produtos)
    return [(p, "PCR", p) for p in sorted(produtos, reverse=True)], diag

# --- 5. INTERFACE DO USUÁRIO ---

if 'lang' not in st.session_state:
    st.session_state.lang = "PT"

with st.sidebar:
    st.markdown("""
    <div style="text-align: left; margin-bottom: 10px;">
        <h1 style="font-family: 'Inter', sans-serif; font-weight: 800; color: #0F766E; font-size: 26px; letter-spacing: -1px; margin:0;">
            BioSpark
        </h1>
        <p style="font-size: 10px; color: #0F766E; opacity: 0.7; margin:0; text-transform: uppercase; letter-spacing: 1px;">Studio</p>
    </div>
    """, unsafe_allow_html=True)
    st.markdown("---")

    lang = st.session_state.lang

    st.caption(TEXTS["sidebar_config"][lang])
    
    num_pocos = st.slider(TEXTS["sidebar_wells"][lang], 1, 15, 4)
    agarose = st.slider(TEXTS["sidebar_agarose"][lang], 0.5, 2.0, 1.0, 0.1)
    
    st.divider()
    
    st.caption(TEXTS["sidebar_visual"][lang])
    estilo_gel = st.selectbox(
        TEXTS["sidebar_theme"][lang], 
        ["Profissional (Dark P&B)", "Publicação (Light P&B)", "Neon (Verde/Laranja)"]
    )
    
    st.markdown("---")
    
    with st.expander(f"ℹ️ {TEXTS['guide_title'][lang]}"):
        st.markdown(TEXTS["guide_content"][lang])
    
    st.markdown("---")
    
    st.caption(TEXTS["pref_lang"][lang])
    idioma_selecionado = st.selectbox("Lang", ["Português", "English"], label_visibility="collapsed")
    
    novo_lang = "PT" if idioma_selecionado == "Português" else "EN"
    if novo_lang != st.session_state.lang:
        st.session_state.lang = novo_lang
        st.rerun()

    # RODAPÉ DA SIDEBAR COMPACTO E LIMPO
    st.markdown(f"""
    <div class="sidebar-footer">
        <strong>{TEXTS['created_by'][lang]} Elton Ostetti</strong><br>
        {TEXTS['lab_name'][lang]}<br>
        {TEXTS['institute'][lang]}<br>
        FAPESP • USP
        <a class="bug-report" href="mailto:e.ostetti.proppg@proppg.butantan.gov.br?subject=Bug%20Report%20BioSpark">
            {TEXTS['report_bug'][lang]}
        </a>
    </div>
    """, unsafe_allow_html=True)

# --- ÁREA PRINCIPAL ---

st.markdown(f"# {TEXTS['header_title'][lang]}")
st.markdown(TEXTS["header_sub"][lang])
st.markdown(" ")

relatorio_dados = []
dados_para_plotar = []
labels_eixo_x = []
nomes_ladders = [] 

cols = st.columns(4)

for i in range(num_pocos):
    col_atual = cols[i % 4]
    with col_atual:
        with st.expander(f"🔹 {TEXTS['well_title'][lang]} {i+1}", expanded=(i==0)):
            opcoes_tipo = [TEXTS['opt_sample'][lang], TEXTS['opt_pcr'][lang], TEXTS['opt_ladder'][lang]]
            tipo_display = st.radio("Tipo", options=opcoes_tipo, key=f"t_{i}", horizontal=True, label_visibility="collapsed")
            
            if tipo_display == TEXTS['opt_ladder'][lang]: tipo = "Ladder"
            elif tipo_display == TEXTS['opt_pcr'][lang]: tipo = "PCR"
            else: tipo = "Amostra"
            
            if tipo == "Ladder":
                lad = st.selectbox(TEXTS['sel_ladder'][lang], list(LADDERS.keys()), key=f"l_{i}")
                ladder_data = [(tam, "Ladder", tam) for tam in LADDERS[lad]]
                dados_para_plotar.append(ladder_data)
                
                rotulo_custom = st.text_input(TEXTS['label_gel'][lang], value="M", key=f"lbl_{i}")
                labels_eixo_x.append(rotulo_custom)
                nomes_ladders.append(lad)
                
                relatorio_dados.append({
                    "Poço": i+1, "Tipo": "Ladder", "Detalhes": lad, "Bandas (pb)": "; ".join([str(t) for t in LADDERS[lad]])
                })
            
            else: 
                nomes_ladders.append(None)
                tab_f, tab_t = st.tabs([TEXTS['tab_file'][lang], TEXTS['tab_text'][lang]])
                seq, nome_arquivo = "", ""
                
                with tab_f:
                    up = st.file_uploader(TEXTS['upload_label'][lang], type=['dna', 'fasta', 'txt', 'fa'], key=f"u_{i}", label_visibility="collapsed")
                    if up: nome_arquivo, seq = processar_upload(up)
                        
                with tab_t:
                    txt = st.text_area(TEXTS['paste_label'][lang], height=70, key=f"tx_{i}", label_visibility="collapsed", placeholder="ATGC...")
                    if txt and not seq: nome_arquivo, seq = processar_texto_manual(txt)
                
                st.markdown("---")
                val_rotulo = nome_arquivo if nome_arquivo else str(i+1)
                
                if tipo == "Amostra":
                    circ = st.checkbox(TEXTS['check_circular'][lang], True, key=f"c_{i}")
                    enz = st.multiselect(TEXTS['sel_enzymes'][lang], TODAS_ENZIMAS, key=f"e_{i}")
                    rotulo_custom = st.text_input(TEXTS['label_gel'][lang], value=val_rotulo[:10], key=f"lbl_{i}")
                    labels_eixo_x.append(rotulo_custom)

                    if seq:
                        res = calcular_digestao(seq, enz, circ)
                        dados_para_plotar.append(res)
                        fragmentos_str = "; ".join([str(int(b[0])) for b in res])
                        desc_enzimas = ", ".join(enz) if enz else "Uncut"
                        relatorio_dados.append({
                            "Poço": i+1, "Tipo": "Digestão", "Detalhes": desc_enzimas, "Bandas (pb)": fragmentos_str
                        })
                    else:
                        dados_para_plotar.append([])
                        relatorio_dados.append({"Poço": i+1, "Tipo": "Vazio", "Bandas (pb)": "-"})

                elif tipo == "PCR":
                    fwd = st.text_input(TEXTS['pcr_fwd'][lang], key=f"fwd_{i}", placeholder="5'->3'")
                    rev = st.text_input(TEXTS['pcr_rev'][lang], key=f"rev_{i}", placeholder="3'->5' (Visual) ou 5'->3'")
                    circ = st.checkbox(TEXTS['check_circular'][lang], False, key=f"cp_{i}")
                    
                    rotulo_custom = st.text_input(TEXTS['label_gel'][lang], value=f"PCR-{i+1}", key=f"lbl_{i}")
                    labels_eixo_x.append(rotulo_custom)
                    
                    if seq and fwd and rev:
                        res, diag = smart_pcr_search(seq, fwd, rev, circ)
                        dados_para_plotar.append(res)
                        
                        if not diag['fwd_found']: 
                            st.markdown(f"<p class='error-text'>{TEXTS['diag_fwd_fail'][lang]}</p>", unsafe_allow_html=True)
                        if not diag['rev_found']: 
                            st.markdown(f"<p class='error-text'>{TEXTS['diag_rev_fail'][lang]}</p>", unsafe_allow_html=True)
                        if diag['products'] > 1: 
                            st.markdown(f"<p class='warning-text'>{TEXTS['warn_multiple'][lang]}</p>", unsafe_allow_html=True)
                        elif diag['products'] == 0 and diag['fwd_found'] and diag['rev_found']: 
                            st.warning(TEXTS['warn_no_product'][lang])
                        
                        fragmentos_str = "; ".join([str(int(b[0])) for b in res])
                        relatorio_dados.append({
                            "Poço": i+1, "Tipo": "PCR", "Detalhes": f"Fwd:{fwd[:5]}.. Rev:{rev[:5]}..", "Bandas (pb)": fragmentos_str
                        })
                    else:
                        dados_para_plotar.append([])
                        relatorio_dados.append({"Poço": i+1, "Tipo": "Vazio", "Bandas (pb)": "-"})

st.markdown(" ") 
st.markdown(f"### {TEXTS['result_title'][lang]}")

if any(dados_para_plotar):
    
    if "Neon" in estilo_gel: bg, txt, c_samp, c_lad = '#111827', 'white', '#00ff41', '#ff9900'
    elif "Profissional" in estilo_gel: bg, txt, c_samp, c_lad = '#000000', 'white', 'white', 'white'
    else: bg, txt, c_samp, c_lad = 'white', 'black', 'black', 'black'

    # CORREÇÃO EIXO Y (100pb VISÍVEL)
    # Subtraímos um valor maior (min_view * 0.7) para garantir margem inferior
    min_view_base = 50 + (100 * (agarose - 0.5))
    min_view = min_view_base * 0.7 
    max_view = 25000 / (agarose * 0.8)
    
    # Eixo X Fixo (15 poços)
    max_range = max(num_pocos, 15) + 0.5

    fig = go.Figure()

    for i, bands in enumerate(dados_para_plotar):
        x = i + 1
        is_ladder = (relatorio_dados[i].get("Tipo") == "Ladder")
        color = c_lad if is_ladder else c_samp
        
        for (size, type, real) in bands:
            if size < (min_view * 0.9) or size > (max_view * 1.1): continue

            width = 2; opacity = 0.8
            if is_ladder:
                if size in [500, 1000, 3000]: width=5; opacity=1.0
            else:
                if type == "Supercoiled": width=4; opacity=0.7
                elif type == "PCR": width=3; opacity=0.9
            
            # Largura Fixa (0.25)
            fig.add_trace(go.Scatter(
                x=[x-0.25, x+0.25], y=[size, size], mode='lines',
                line=dict(color=color, width=width), opacity=opacity,
                hoverinfo='text', hovertext=f"{int(size)} pb"
            ))
            
            if is_ladder:
                fig.add_trace(go.Scatter(x=[x-0.35], y=[size], mode="text", text=[str(size)], textfont=dict(color=txt, size=9), showlegend=False))

    fig.update_layout(
        plot_bgcolor=bg, paper_bgcolor=bg, height=600,
        margin=dict(t=30, b=80, l=80, r=40),
        xaxis=dict(
            tickmode='array', tickvals=list(range(1, num_pocos + 1)),
            ticktext=labels_eixo_x,
            tickfont=dict(color=txt, size=14, family='Arial'),
            showgrid=False, zeroline=False, range=[0.5, max_range] 
        ),
        yaxis=dict(
            type='log',
            range=[math.log10(min_view), math.log10(max_view)],
            showgrid=False, zeroline=False, showticklabels=False
        )
    )
    
    st.plotly_chart(fig, use_container_width=True)
    
    with st.expander(f"📥 {TEXTS['export_expander'][lang]}"):
        df = pd.DataFrame(relatorio_dados)
        csv = df.to_csv(index=False).encode('utf-8')
        st.download_button(
            label=TEXTS['btn_download'][lang],
            data=csv,
            file_name='gel_result.csv',
            mime='text/csv',
        )

else:
    st.info(TEXTS['empty_msg'][lang])

# --- RODAPÉ PRINCIPAL ---
st.markdown("""
<div class="footer">
    <p><b>BioSpark</b></p>
</div>
""", unsafe_allow_html=True)
