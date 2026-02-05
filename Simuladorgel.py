import streamlit as st
import plotly.graph_objects as go
import pandas as pd
import math
import re
from Bio.Seq import Seq
from Bio.Restriction import RestrictionBatch, Analysis, CommOnly
from io import StringIO, BytesIO
from Bio import SeqIO

# --- 1. CONFIGURAÇÃO E ESTADO DA SESSÃO ---
st.set_page_config(page_title="BioSpark Studio", layout="wide", page_icon="🧬")

# Inicializa variáveis de estado (Memória) se não existirem
if 'lab_notebook' not in st.session_state:
    st.session_state.lab_notebook = ""

# Garante que os estados dos poços existam para receber o Batch Upload
for i in range(15):
    if f"seq_{i}" not in st.session_state: st.session_state[f"seq_{i}"] = ""
    if f"lbl_{i}" not in st.session_state: st.session_state[f"lbl_{i}"] = ""

# --- 2. FUNÇÕES AUXILIARES ---
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

# Função CALLBACK para o Batch Upload
def handle_batch_upload():
    uploaded_files = st.session_state.batch_uploader
    if uploaded_files:
        for idx, file in enumerate(uploaded_files):
            if idx >= 15: break # Limite de poços
            name, seq = processar_upload(file)
            if name != "Erro":
                st.session_state[f"seq_{idx}"] = seq
                st.session_state[f"lbl_{idx}"] = name[:10]

# --- 3. LÓGICA BIOLÓGICA (V1.0 ESTÁVEL) ---
def calcular_digestao(sequencia, enzimas, eh_circular):
    if not sequencia or sequencia.startswith("Erro"): return []
    seq_obj = Seq(sequencia)
    tamanho_total = len(seq_obj)
    
    if eh_circular and not enzimas:
        return [(tamanho_total * 1.4, "Nicked (Relaxed)", tamanho_total), (tamanho_total * 0.7, "Supercoiled", tamanho_total)]
    if not enzimas: return [(tamanho_total, "Linear", tamanho_total)]
    
    rb = RestrictionBatch(enzimas)
    analise = Analysis(rb, seq_obj, linear=not eh_circular)
    cortes = analise.full()
    locais = sorted(list(set([local for lista in cortes.values() for local in lista])))
    
    if not locais: return [(tamanho_total, "Uncut", tamanho_total)]
    fragmentos = []
    if not eh_circular:
        prev = 0
        for cut in locais:
            fragmentos.append(cut - prev); prev = cut
        fragmentos.append(tamanho_total - prev)
    else:
        if len(locais) == 1: fragmentos.append(tamanho_total)
        else:
            for i in range(len(locais)-1): fragmentos.append(locais[i+1] - locais[i])
            fragmentos.append((tamanho_total - locais[-1]) + locais[0])
    return [(frag, "Fragmento", frag) for frag in sorted(fragmentos, reverse=True)]

def calcular_pcr_biologico(sequencia, fwd_seq, rev_seq, eh_circular):
    if not sequencia or sequencia.startswith("Erro"): return [], False
    template = sequencia.upper()
    fwd = "".join(fwd_seq.split()).upper()
    rev = "".join(rev_seq.split()).upper()
    if len(fwd) < 10 or len(rev) < 10: return [], False

    # Lógica de Seed 15pb (Estável V1.0)
    SEED_SIZE = 15
    fwd_seed = fwd[-SEED_SIZE:] if len(fwd) > SEED_SIZE else fwd
    rev_seed = rev[-SEED_SIZE:] if len(rev) > SEED_SIZE else rev
    
    fwd_matches = [m.start() for m in re.finditer(fwd_seed, template)]
    rev_seed_rc = str(Seq(rev_seed).reverse_complement())
    rev_matches = [m.start() for m in re.finditer(rev_seed_rc, template)]
    
    produtos = []
    for f_pos in fwd_matches:
        f_3prime_end = f_pos + len(fwd_seed)
        for r_pos in rev_matches:
            if r_pos > f_pos:
                dist = r_pos - f_3prime_end
                if dist >= 0: produtos.append(len(fwd) + len(rev) + dist)
            elif eh_circular and r_pos < f_pos:
                dist = (len(template) - f_3prime_end) + r_pos
                produtos.append(len(fwd) + len(rev) + dist)
                
    return [(p, "PCR Product", p) for p in sorted(produtos, reverse=True)], len(produtos) > 1

# --- 4. TRADUÇÃO E DADOS ---
TODAS_ENZIMAS = sorted([str(e) for e in CommOnly])
LADDERS = {
    "1kb Plus DNA Ladder": [100, 200, 300, 400, 500, 650, 850, 1000, 1650, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 12000],
    "1kb DNA Ladder (Genérico)": [250, 500, 750, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 8000, 10000],
    "100bp DNA Ladder": [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1517, 2017],
    "High Mass": [1000, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 20000, 48500]
}
TEXTS = {
    "header_title": { "PT": "Simulador de Biologia Molecular", "EN": "Molecular Biology Simulator" },
    "header_sub": { "PT": "Digestão Enzimática e PCR In Silico.", "EN": "In Silico Enzymatic Digestion and PCR." },
    "sidebar_config": { "PT": "CONFIGURAÇÕES", "EN": "SETTINGS" },
    "sidebar_wells": { "PT": "Nº de Poços", "EN": "Well Count" },
    "sidebar_agarose": { "PT": "Agarose (%)", "EN": "Agarose (%)" },
    "batch_label": { "PT": "⚡ Upload Rápido em Lote (Arrastar múltiplos arquivos aqui preenche os poços)", "EN": "⚡ Smart Batch Upload" },
    "notebook_tab": { "PT": "📔 Caderno Lab", "EN": "📔 Lab Notebook" },
    "settings_tab": { "PT": "⚙️ Ajustes", "EN": "⚙️ Settings" },
    "well_title": { "PT": "Poço", "EN": "Well" },
    "opt_sample": { "PT": "Digestão", "EN": "Digestion" },
    "opt_ladder": { "PT": "Ladder", "EN": "Ladder" },
    "opt_pcr": { "PT": "PCR", "EN": "PCR" },
    "sel_ladder": { "PT": "Selecione o Ladder:", "EN": "Select Ladder:" },
    "label_gel": { "PT": "Rótulo:", "EN": "Label:" },
    "tab_file": { "PT": "📂 Upload", "EN": "📂 File" },
    "tab_text": { "PT": "📝 Digitar", "EN": "📝 Type" },
    "upload_label": { "PT": "Arquivo", "EN": "File" },
    "paste_label": { "PT": "Sequência", "EN": "Sequence" },
    "check_circular": { "PT": "Circular?", "EN": "Circular?" },
    "sel_enzymes": { "PT": "Enzimas", "EN": "Enzymes" },
    "pcr_fwd": { "PT": "Fwd (5'->3')", "EN": "Fwd (5'->3')" },
    "pcr_rev": { "PT": "Rev (5'->3')", "EN": "Rev (5'->3')" },
    "result_title": { "PT": "Resultado da Eletroforese", "EN": "Electrophoresis Result" },
    "export_expander": { "PT": "Exportar Dados", "EN": "Export Data" },
    "btn_download": { "PT": "Baixar .csv", "EN": "Download .csv" },
    "empty_msg": { "PT": "Adicione amostras.", "EN": "Add samples." },
    "created_by": { "PT": "Desenvolvido por", "EN": "Developed by" },
    "report_bug": { "PT": "✉️ Reportar Problema", "EN": "✉️ Report Bug" },
    "warn_multiple": { "PT": "⚠️ Múltiplos sítios!", "EN": "⚠️ Multiple sites!" },
    "warn_no_product": { "PT": "Sem produto.", "EN": "No product." },
    "ack_title": { "PT": "Apoio e Afiliação", "EN": "Support & Affiliation" }
}

# --- 5. ESTILO CSS (TURQUESA + MINIMALISTA) ---
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');
    .stApp { background: linear-gradient(180deg, #F0F9FF 0%, #FFFFFF 100%); font-family: 'Inter', sans-serif; }
    section[data-testid="stSidebar"] { background-color: #E0F7FA; border-right: 1px solid #B2EBF2; }
    div[data-baseweb="slider"] div[class*="StyledThumb"] { background-color: #0F766E !important; border-color: #0F766E !important; }
    div[data-baseweb="slider"] div[class*="StyledTrack"] > div { background-color: #0F766E !important; }
    div[data-baseweb="slider"] div[class*="StyledTrack"] { background-color: #B2EBF2 !important; }
    h1, h2, h3 { color: #0F172A !important; font-weight: 700 !important; letter-spacing: -0.02em; }
    .stExpander { background-color: #FFFFFF; border-radius: 6px !important; border: 1px solid #E2E8F0 !important; box-shadow: 0 1px 2px 0 rgba(0,0,0,0.05) !important; margin-bottom: 0.5rem; }
    .stExpander:hover { border-color: #0F766E !important; }
    .stButton > button { background-color: #0F766E; color: white; border-radius: 6px; font-weight: 500; border: none; }
    .stButton > button:hover { background-color: #0d6e66; color: white; }
    span[data-baseweb="checkbox"] div { background-color: #0F766E !important; }
    button[data-baseweb="tab"] { font-size: 13px !important; padding: 5px 10px !important; }
    .footer { width: 100%; text-align: center; padding: 20px 0; font-size: 11px; color: #64748B; border-top: 1px solid #CBD5E1; margin-top: 40px; opacity: 0.8; }
    .sidebar-footer { margin-top: 20px; padding-top: 15px; border-top: 1px solid #B2EBF2; font-size: 11px; color: #333333; line-height: 1.5; }
    .sidebar-footer strong { color: #0F766E; }
    .bug-report { font-size: 11px; color: #64748B; text-decoration: none; margin-top: 5px; display: inline-block; }
    .bug-report:hover { color: #0F766E; text-decoration: underline; }
    .warning-text { color: #DC2626; font-weight: bold; font-size: 12px; }
</style>
""", unsafe_allow_html=True)

# --- 6. INTERFACE PRINCIPAL ---
if 'lang' not in st.session_state: st.session_state.lang = "PT"
lang = st.session_state.lang

with st.sidebar:
    st.markdown("""<div style="text-align: left; margin-bottom: 20px;"><h1 style="color: #0F766E; margin:0; font-size:24px;">BioSpark</h1><p style="font-size: 10px; color: #0F766E; margin-top:-2px;">STUDIO v2.0</p></div>""", unsafe_allow_html=True)
    
    # Abas na Sidebar: Configuração vs Caderno
    tab_conf, tab_note = st.tabs([TEXTS["settings_tab"][lang], TEXTS["notebook_tab"][lang]])
    
    with tab_conf:
        num_pocos = st.slider(TEXTS["sidebar_wells"][lang], 1, 15, 4)
        agarose = st.slider(TEXTS["sidebar_agarose"][lang], 0.5, 2.0, 1.0, 0.1)
        estilo_gel = st.selectbox("Tema", ["Profissional (Dark P&B)", "Publicação (Light P&B)", "Neon (Verde/Laranja)"])
        st.markdown("---")
        idioma = st.selectbox("Lang", ["Português", "English"], label_visibility="collapsed")
        if (idioma == "Português" and lang != "PT") or (idioma == "English" and lang != "EN"):
            st.session_state.lang = "PT" if idioma == "Português" else "EN"
            st.rerun()

    with tab_note:
        st.caption("Suas anotações não somem ao recarregar.")
        st.session_state.lab_notebook = st.text_area("Notas", value=st.session_state.lab_notebook, height=300, label_visibility="collapsed", placeholder="Digite aqui...")

    st.markdown(f"""<div class="sidebar-footer"><strong>{TEXTS['created_by'][lang]} Elton Ostetti</strong><br><a class="bug-report" href="mailto:e.ostetti.proppg@proppg.butantan.gov.br?subject=Bug%20Report%20BioSpark">{TEXTS['report_bug'][lang]}</a><br><br><strong>{TEXTS['ack_title'][lang]}</strong><br>FAPESP • USP • Instituto Butantan</div>""", unsafe_allow_html=True)

st.markdown(f"# {TEXTS['header_title'][lang]}")
st.markdown(TEXTS["header_sub"][lang])

# ÁREA DE SMART BATCH UPLOAD (ARRASTAR E SOLTAR)
with st.expander(TEXTS['batch_label'][lang]):
    st.file_uploader("Arraste arquivos aqui para preencher sequencialmente", type=['dna', 'fasta', 'txt'], accept_multiple_files=True, key="batch_uploader", on_change=handle_batch_upload)

st.markdown(" ")

# --- LÓGICA DE POÇOS COM ESTADO ---
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
                
                rotulo_custom = st.text_input(TEXTS['label_gel'][lang], value="M", key=f"lbl_input_{i}")
                labels_eixo_x.append(rotulo_custom)
                nomes_ladders.append(lad)
                relatorio_dados.append({"Poço": i+1, "Tipo": "Ladder", "Detalhes": lad, "Bandas": str(LADDERS[lad])})
            
            else:
                nomes_ladders.append(None)
                tab_f, tab_t = st.tabs([TEXTS['tab_file'][lang], TEXTS['tab_text'][lang]])
                
                # Sincronização do Upload Individual com o Session State do Batch
                with tab_f:
                    up = st.file_uploader(TEXTS['upload_label'][lang], type=['dna', 'fasta', 'txt'], key=f"u_{i}", label_visibility="collapsed")
                    if up: 
                        n, s = processar_upload(up)
                        if n != "Erro":
                            st.session_state[f"seq_{i}"] = s
                            st.session_state[f"lbl_{i}"] = n[:10]

                with tab_t:
                    # O text_area lê do session_state (seq_i)
                    seq_input = st.text_area(TEXTS['paste_label'][lang], value=st.session_state.get(f"seq_{i}", ""), height=70, key=f"tx_area_{i}", label_visibility="collapsed", placeholder="ATGC...")
                    # Se o usuário digitar, atualiza o estado
                    if seq_input != st.session_state.get(f"seq_{i}", ""):
                        n_man, s_man = processar_texto_manual(seq_input)
                        st.session_state[f"seq_{i}"] = s_man

                st.markdown("---")
                
                # Rótulo sincronizado
                val_rotulo = st.session_state.get(f"lbl_{i}", str(i+1))
                rotulo_custom = st.text_input(TEXTS['label_gel'][lang], value=val_rotulo, key=f"lbl_final_{i}")
                labels_eixo_x.append(rotulo_custom)
                
                seq_final = st.session_state.get(f"seq_{i}", "")

                if tipo == "Amostra":
                    circ = st.checkbox(TEXTS['check_circular'][lang], True, key=f"c_{i}")
                    enz = st.multiselect(TEXTS['sel_enzymes'][lang], TODAS_ENZIMAS, key=f"e_{i}")
                    if seq_final:
                        try:
                            res = calcular_digestao(seq_final, enz, circ)
                            dados_para_plotar.append(res)
                            relatorio_dados.append({"Poço": i+1, "Tipo": "Digestão", "Bandas": str(res)})
                        except: dados_para_plotar.append([])
                    else:
                        dados_para_plotar.append([])
                        relatorio_dados.append({"Poço": i+1, "Tipo": "Vazio", "Bandas": "-"})

                elif tipo == "PCR":
                    fwd = st.text_input(TEXTS['pcr_fwd'][lang], key=f"fwd_{i}", placeholder="ATGC...")
                    rev = st.text_input(TEXTS['pcr_rev'][lang], key=f"rev_{i}", placeholder="ATGC...")
                    circ = st.checkbox(TEXTS['check_circular'][lang], False, key=f"cp_{i}")
                    
                    if seq_final and fwd and rev:
                        res, tem_inesp = calcular_pcr_biologico(seq_final, fwd, rev, circ)
                        dados_para_plotar.append(res)
                        if tem_inesp: st.markdown(f"<p class='warning-text'>{TEXTS['warn_multiple'][lang]}</p>", unsafe_allow_html=True)
                        if not res: st.warning(TEXTS['warn_no_product'][lang])
                        relatorio_dados.append({"Poço": i+1, "Tipo": "PCR", "Bandas": str(res)})
                    else:
                        dados_para_plotar.append([])
                        relatorio_dados.append({"Poço": i+1, "Tipo": "Vazio", "Bandas": "-"})

# --- GRÁFICO (PLOTLY) ---
st.markdown(" ") 
if any(dados_para_plotar):
    if "Neon" in estilo_gel: bg, txt, c_samp, c_lad = '#111827', 'white', '#00ff41', '#ff9900'
    elif "Profissional" in estilo_gel: bg, txt, c_samp, c_lad = '#000000', 'white', 'white', 'white'
    else: bg, txt, c_samp, c_lad = 'white', 'black', 'black', 'black'

    min_view = 25 
    max_view = 25000 / (agarose * 0.8)
    max_range = max(num_pocos, 15) + 0.5

    fig = go.Figure()
    for i, lista_bandas in enumerate(dados_para_plotar):
        x_center = i + 1
        eh_ladder = (nomes_ladders[i] is not None)
        cor_atual = color_ladder if eh_ladder else color_sample
        
        for (tam_aparente, tipo_banda, tam_real) in lista_bandas:
            if tam_aparente < (min_view * 0.9) or tam_aparente > (max_view * 1.1): continue
            width = 2; opacity = 0.8
            if eh_ladder:
                if tam_aparente in [500, 1000, 3000]: width = 7; opacity = 1.0
            else:
                if tipo_banda == "Supercoiled": width = 4; opacity = 0.7
                elif tipo_banda == "PCR Product": width = 3; opacity = 0.9
            
            # VISUAL PALITO/HALTERE RESTAURADO
            fig.add_trace(go.Scatter(
                x=[x_center - 0.28, x_center + 0.28], y=[tam_aparente, tam_aparente],
                mode='lines+markers', # As bolinhas estão aqui!
                line=dict(color=cor_atual, width=width),
                marker=dict(color=cor_atual, size=width, symbol='circle'),
                opacity=opacity, showlegend=False,
                hoverinfo='text', hovertext=f"<b>~{int(tam_aparente)} pb</b>"
            ))
            
            if eh_ladder:
                fig.add_trace(go.Scatter(x=[x_center - 0.4], y=[tam_aparente], mode="text", text=[str(tam_aparente)], textfont=dict(color=txt, size=9), showlegend=False))

    fig.update_layout(
        plot_bgcolor=bg, paper_bgcolor=bg, height=600,
        margin=dict(t=40, b=40, l=40, r=40),
        xaxis=dict(tickvals=list(range(1, num_pocos+1)), ticktext=labels_eixo_x, showgrid=False, tickfont=dict(color=txt), range=[0.5, max_range]),
        yaxis=dict(type='log', range=[math.log10(min_view), math.log10(max_view)], showgrid=False, showticklabels=False)
    )
    st.plotly_chart(fig, use_container_width=True)
    
    with st.expander(f"📥 {TEXTS['export_expander'][lang]}"):
        df = pd.DataFrame(relatorio_dados)
        csv = df.to_csv(index=False).encode('utf-8')
        st.download_button(TEXTS['btn_download'][lang], csv, 'gel_result.csv', 'text/csv')
else:
    st.info(TEXTS['empty_msg'][lang])

st.markdown("""<div class="footer"><p><b>BioSpark</b></p></div>""", unsafe_allow_html=True)
