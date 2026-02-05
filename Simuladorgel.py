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
    "sidebar_wells": { "PT": "Nº de Poços", "EN": "Well Count" },
    "sidebar_agarose": { "PT": "Agarose (%)", "EN": "Agarose (%)" },
    "sidebar_visual": { "PT": "VISUALIZAÇÃO", "EN": "VISUALIZATION" },
    "sidebar_theme": { "PT": "Tema do Gel", "EN": "Gel Theme" },
    "guide_title": { "PT": "Guia Rápido", "EN": "Quick Guide" },
    "guide_content": {
        "PT": """
        **Digestão:** Upload DNA + Enzimas.
        **PCR:** Cole os primers (Suporta overhangs).
        **Ladder:** Escolha o marcador.
        """,
        "EN": """
        **Digestion:** Upload DNA + Enzymes.
        **PCR:** Paste primers (Supports overhangs).
        **Ladder:** Select marker.
        """
    },
    "well_title": { "PT": "Poço", "EN": "Well" },
    "opt_sample": { "PT": "Digestão", "EN": "Digestion" },
    "opt_ladder": { "PT": "Ladder", "EN": "Ladder" },
    "opt_pcr": { "PT": "PCR", "EN": "PCR" },
    "sel_ladder": { "PT": "Selecione o Ladder:", "EN": "Select Ladder:" },
    "label_gel": { "PT": "Rótulo:", "EN": "Label:" },
    "tab_file": { "PT": "📂 Arquivo", "EN": "📂 File" },
    "tab_text": { "PT": "📝 Texto", "EN": "📝 Text" },
    "upload_label": { "PT": "Upload DNA", "EN": "Upload DNA" },
    "paste_label": { "PT": "Sequência", "EN": "Sequence" },
    "check_circular": { "PT": "Circular?", "EN": "Circular?" },
    "sel_enzymes": { "PT": "Enzimas", "EN": "Enzymes" },
    "pcr_fwd": { "PT": "Primer Forward", "EN": "Forward Primer" },
    "pcr_rev": { "PT": "Primer Reverse", "EN": "Reverse Primer" },
    "result_title": { "PT": "Resultado da Eletroforese", "EN": "Electrophoresis Result" },
    "export_expander": { "PT": "Exportar Dados", "EN": "Export Data" },
    "btn_download": { "PT": "Baixar .csv", "EN": "Download .csv" },
    "empty_msg": { "PT": "Adicione amostras para começar.", "EN": "Add samples to start." },
    "created_by": { "PT": "Dev.", "EN": "Dev." },
    "report_bug": { "PT": "✉️ Reportar Erro", "EN": "✉️ Report Bug" },
    "warn_multiple": { "PT": "⚠️ Múltiplos sítios de ligação!", "EN": "⚠️ Multiple binding sites!" },
    "warn_no_product": { "PT": "Nenhum produto.", "EN": "No product." },
    "diag_fwd_fail": { "PT": "❌ Fwd não anela", "EN": "❌ Fwd failed" },
    "diag_rev_fail": { "PT": "❌ Rev não anela", "EN": "❌ Rev failed" },
}

# --- 3. ESTILO CSS (CORRIGIDO E COMPACTO) ---
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
        padding-top: 1rem;
    }
    
    /* Compactar Sidebar */
    section[data-testid="stSidebar"] .block-container {
        padding-top: 1rem;
        padding-bottom: 1rem;
        gap: 0.5rem;
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
        font-size: 13px !important;
        padding: 5px 10px !important;
    }

    /* RODAPÉ LATERAL - COR ESCURA (#333) */
    .sidebar-footer {
        margin-top: 15px;
        padding-top: 15px;
        border-top: 1px solid #CBD5E1;
        font-size: 11px;
        color: #333333; 
        line-height: 1.4;
    }
    .sidebar-footer strong {
        color: #111827;
        font-weight: 600;
    }
    .bug-report {
        font-size: 10px;
        color: #4B5563;
        text-decoration: none;
        margin-left: 5px;
    }
    .bug-report:hover {
        color: #0F766E;
        text-decoration: underline;
    }
    
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
    
    .warning-text { color: #DC2626; font-weight: bold; font-size: 12px; margin: 2px 0; }
    .error-text { color: #B91C1C; font-size: 12px; margin: 2px 0; }
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
            except Exception as e: return "Erro", f"Erro .dna: {str(e)}"
        
        bytes_data = input_data.getvalue()
        try: conteudo = bytes_data.decode("utf-8")
        except: conteudo = bytes_data.decode("latin-1")
        
        if ">" in conteudo:
            try:
                iterator = SeqIO.parse(StringIO(conteudo), "fasta")
                record = next(iterator)
                return record.id, str(record.seq).upper()
            except: pass
            
        seq_limpa = "".join([l.strip() for l in conteudo.splitlines() if not l.startswith(">")])
        seq_final = clean_sequence(seq_limpa)
        if len(seq_final) > 0 and any(c not in "ATGCNRYKMSWBDHV" for c in seq_final[:100]):
             return "Erro", "Arquivo inválido."
        return nome_sugerido, seq_final
    except Exception as e: return "Erro", str(e)

def processar_texto_manual(texto):
    if ">" in texto:
        try:
            iterator = SeqIO.parse(StringIO(texto), "fasta")
            record = next(iterator)
            return record.id, str(record.seq).upper()
        except: pass
    return "Seq Manual", clean_sequence(texto)

def calcular_digestao(sequencia, enzimas, eh_circular):
    if not sequencia or sequencia.startswith("Erro"): return []
    seq_obj = Seq(sequencia)
    tamanho_total = len(seq_obj)
    
    if eh_circular and not enzimas:
        return [(tamanho_total * 1.4, "Nicked", tamanho_total), (tamanho_total * 0.7, "Supercoiled", tamanho_total)]
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
        if len(locais) == 1: fragmentos.append(tamanho_total)
        else:
            for i in range(len(locais)-1): fragmentos.append(locais[i+1] - locais[i])
            fragmentos.append((tamanho_total - locais[-1]) + locais[0])
            
    return [(frag, "Fragmento", frag) for frag in sorted(fragmentos, reverse=True)]

def smart_pcr_search(template, fwd_raw, rev_raw, eh_circular):
    fwd = clean_sequence(fwd_raw)
    rev = clean_sequence(rev_raw)
    template = template.upper()
    diag = {'fwd_found': False, 'rev_found': False, 'products': 0}
    
    if len(fwd) < 10 or len(rev) < 10: return [], diag

    SEED = 12
    fwd_seed = fwd[-SEED:]
    fwd_matches = [m.start() for m in re.finditer(fwd_seed, template)]
    if fwd_matches: diag['fwd_found'] = True
    
    rev_matches = []
    # Hipótese A: 5'->3' (Standard)
    rev_seed_std = rev[-SEED:]
    rev_rc_std = str(Seq(rev_seed_std).reverse_complement())
    matches_std = [m.start() for m in re.finditer(rev_rc_std, template)]
    
    # Hipótese B: 3'->5' (Visual)
    rev_seed_inv = rev[:SEED]
    rev_compl_inv = str(Seq(rev_seed_inv).complement())
    matches_inv = [m.start() for m in re.finditer(rev_compl_inv, template)]
    
    if len(matches_inv) > 0: rev_matches = matches_inv
    elif len(matches_std) > 0: rev_matches = matches_std
    if rev_matches: diag['rev_found'] = True
    
    produtos = []
    for f_pos in fwd_matches:
        f_end = f_pos + len(fwd_seed)
        for r_pos in rev_matches:
            if r_pos > f_pos:
                dist = r_pos - f_end
                if dist >= 0: produtos.append(len(fwd) + len(rev) + dist)
            elif eh_circular and r_pos < f_pos:
                dist = (len(template) - f_end) + r_pos
                produtos.append(len(fwd) + len(rev) + dist)
                
    diag['products'] = len(produtos)
    return [(p, "PCR", p) for p in sorted(produtos, reverse=True)], diag

# --- 5. INTERFACE ---

if 'lang' not in st.session_state: st.session_state.lang = "PT"

with st.sidebar:
    st.markdown("""<div style="text-align: left; margin-bottom: 10px;"><h1 style="color: #0F766E; margin:0; font-size:24px;">BioSpark</h1><p style="font-size: 10px; color: #0F766E; margin-top:-2px;">STUDIO</p></div>""", unsafe_allow_html=True)
    st.markdown("---")
    lang = st.session_state.lang
    
    st.caption(TEXTS["sidebar_config"][lang])
    num_pocos = st.slider(TEXTS["sidebar_wells"][lang], 1, 15, 4)
    agarose = st.slider(TEXTS["sidebar_agarose"][lang], 0.5, 2.0, 1.0, 0.1)
    
    st.markdown("---")
    st.caption(TEXTS["sidebar_visual"][lang])
    estilo_gel = st.selectbox(TEXTS["sidebar_theme"][lang], ["Profissional (Dark P&B)", "Publicação (Light P&B)", "Neon (Verde/Laranja)"])
    
    st.markdown("---")
    with st.expander(f"ℹ️ {TEXTS['guide_title'][lang]}"): st.markdown(TEXTS["guide_content"][lang])
    
    st.markdown("---")
    new_lang = st.selectbox("Lang", ["Português", "English"], index=0 if lang=="PT" else 1, label_visibility="collapsed")
    if (new_lang == "Português" and lang != "PT") or (new_lang == "English" and lang != "EN"):
        st.session_state.lang = "PT" if new_lang == "Português" else "EN"
        st.rerun()

    # RODAPÉ LATERAL COMPACTO E ESCURO
    st.markdown(f"""
    <div class="sidebar-footer">
        {TEXTS['created_by'][lang]} <strong>Elton Ostetti</strong>
        <a class="bug-report" href="mailto:e.ostetti.proppg@proppg.butantan.gov.br?subject=Bug%20Report%20BioSpark">{TEXTS['report_bug'][lang]}</a>
        <br>
        <span style="opacity:0.9;">FAPESP • USP • Instituto Butantan</span>
    </div>
    """, unsafe_allow_html=True)

# --- MAIN ---
st.markdown(f"# {TEXTS['header_title'][lang]}")
st.markdown(TEXTS["header_sub"][lang])
st.markdown(" ")

cols = st.columns(4)
relatorio = []
plot_data = []
labels_x = []

for i in range(num_pocos):
    with cols[i % 4]:
        with st.expander(f"🔹 {TEXTS['well_title'][lang]} {i+1}", expanded=(i==0)):
            tipo_disp = st.radio("Tipo", [TEXTS['opt_sample'][lang], TEXTS['opt_pcr'][lang], TEXTS['opt_ladder'][lang]], horizontal=True, label_visibility="collapsed", key=f"t{i}")
            
            if tipo_disp == TEXTS['opt_ladder'][lang]:
                lad = st.selectbox("Ladder", list(LADDERS.keys()), key=f"l{i}")
                plot_data.append([(t, "Ladder", t) for t in LADDERS[lad]])
                lbl = st.text_input(TEXTS['label_gel'][lang], "M", key=f"lb{i}")
                labels_x.append(lbl)
                relatorio.append({"Poço": i+1, "Tipo": "Ladder", "Detalhes": lad, "Bandas": str(LADDERS[lad])})
            
            else:
                tab1, tab2 = st.tabs([TEXTS['tab_file'][lang], TEXTS['tab_text'][lang]])
                seq, name = "", ""
                with tab1:
                    f = st.file_uploader("Upload", type=['dna','fasta','txt'], key=f"u{i}", label_visibility="collapsed")
                    if f: name, seq = processar_upload(f)
                with tab2:
                    t = st.text_area("Texto", height=70, key=f"tx{i}", label_visibility="collapsed")
                    if t and not seq: name, seq = processar_texto_manual(t)
                
                st.markdown("---")
                label = st.text_input(TEXTS['label_gel'][lang], (name[:10] if name else str(i+1)), key=f"lb{i}")
                labels_x.append(label)
                
                if tipo_disp == TEXTS['opt_sample'][lang]: # Digestão
                    c = st.checkbox(TEXTS['check_circular'][lang], True, key=f"c{i}")
                    e = st.multiselect(TEXTS['sel_enzymes'][lang], TODAS_ENZIMAS, key=f"e{i}")
                    if seq:
                        res = calcular_digestao(seq, e, c)
                        plot_data.append(res)
                        relatorio.append({"Poço": i+1, "Tipo": "Digestão", "Bandas": ";".join([str(int(x[0])) for x in res])})
                    else: 
                        plot_data.append([])
                        relatorio.append({"Poço": i+1, "Tipo": "Vazio"})
                        
                else: # PCR
                    fwd = st.text_input(TEXTS['pcr_fwd'][lang], key=f"fw{i}")
                    rev = st.text_input(TEXTS['pcr_rev'][lang], key=f"rv{i}")
                    c = st.checkbox(TEXTS['check_circular'][lang], False, key=f"cp{i}")
                    
                    if seq and fwd and rev:
                        res, diag = smart_pcr_search(seq, fwd, rev, c)
                        plot_data.append(res)
                        if not diag['fwd_found']: st.markdown(f"<p class='error-text'>{TEXTS['diag_fwd_fail'][lang]}</p>", unsafe_allow_html=True)
                        if not diag['rev_found']: st.markdown(f"<p class='error-text'>{TEXTS['diag_rev_fail'][lang]}</p>", unsafe_allow_html=True)
                        if diag['products'] > 1: st.markdown(f"<p class='warning-text'>{TEXTS['warn_multiple'][lang]}</p>", unsafe_allow_html=True)
                        elif diag['products'] == 0 and diag['fwd_found'] and diag['rev_found']: st.warning(TEXTS['warn_no_product'][lang])
                        relatorio.append({"Poço": i+1, "Tipo": "PCR", "Bandas": ";".join([str(int(x[0])) for x in res])})
                    else:
                        plot_data.append([])
                        relatorio.append({"Poço": i+1, "Tipo": "Vazio"})

# --- PLOTAGEM (LÓGICA RESTAURADA) ---
st.markdown(" ")
if any(plot_data):
    if "Neon" in estilo_gel: bg, txt, c_samp, c_lad = '#111827', 'white', '#00ff41', '#ff9900'
    elif "Profissional" in estilo_gel: bg, txt, c_samp, c_lad = '#000000', 'white', 'white', 'white'
    else: bg, txt, c_samp, c_lad = 'white', 'black', 'black', 'black'

    # CORREÇÃO EIXO Y PARA 100bp
    min_view = 20 # Trava em 20 para 100 aparecer bem acima
    max_view = 25000 / (agarose * 0.8)
    
    # LARGURA 15 POÇOS FIXA
    max_range = max(num_pocos, 15) + 0.5

    fig = go.Figure()
    for i, bands in enumerate(plot_data):
        x = i + 1
        is_ladder = (relatorio[i].get("Tipo") == "Ladder")
        color = c_lad if is_ladder else c_samp
        
        for (size, type, real) in bands:
            if size < (min_view * 0.9) or size > (max_view * 1.1): continue
            
            width = 2; opacity = 0.8
            if is_ladder:
                if size in [500, 1000, 3000]: width=5; opacity=1.0
            else:
                if type == "Supercoiled": width=4; opacity=0.7
                elif type == "PCR": width=3; opacity=0.9
            
            # VISUAL HALTERE (LINES+MARKERS)
            fig.add_trace(go.Scatter(
                x=[x-0.28, x+0.28], y=[size, size], mode='lines+markers',
                line=dict(color=color, width=width), 
                marker=dict(color=color, size=width, symbol='circle'), # AS BOLINHAS
                opacity=opacity,
                hoverinfo='text', hovertext=f"{int(size)} pb",
                showlegend=False
            ))
            
            if is_ladder:
                fig.add_trace(go.Scatter(x=[x-0.4], y=[size], mode="text", text=[str(size)], textfont=dict(color=txt, size=9), showlegend=False))

    fig.update_layout(
        plot_bgcolor=bg, paper_bgcolor=bg, height=600,
        margin=dict(t=30, b=30, l=50, r=30),
        xaxis=dict(tickvals=list(range(1, num_pocos+1)), ticktext=labels_x, showgrid=False, tickfont=dict(color=txt), range=[0.5, max_range]),
        yaxis=dict(type='log', range=[math.log10(min_view), math.log10(max_view)], showgrid=False, showticklabels=False)
    )
    st.plotly_chart(fig, use_container_width=True)
    
    df = pd.DataFrame(relatorio)
    csv = df.to_csv(index=False).encode('utf-8')
    with st.expander(f"📥 {TEXTS['export_expander'][lang]}"):
        st.download_button(TEXTS['btn_download'][lang], csv, "gel_results.csv", "text/csv")
else:
    st.info(TEXTS['empty_msg'][lang])

st.markdown("""<div class="footer"><p><b>BioSpark</b></p></div>""", unsafe_allow_html=True)
