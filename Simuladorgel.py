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

# --- 2. SISTEMA DE TRADUÇÃO (DICIONÁRIO) ---
TEXTS = {
    "header_title": {
        "PT": "Simulador de Gel de Agarose",
        "EN": "Agarose Gel Simulator"
    },
    "header_sub": {
        "PT": "Ferramenta in silico para Digestão Enzimática e PCR.",
        "EN": "In silico tool for Enzymatic Digestion and PCR."
    },
    "sidebar_config": {
        "PT": "CONFIGURAÇÕES",
        "EN": "SETTINGS"
    },
    "sidebar_wells": {
        "PT": "Número de Poços",
        "EN": "Number of Wells"
    },
    "sidebar_agarose": {
        "PT": "Agarose (%)",
        "EN": "Agarose (%)"
    },
    "sidebar_visual": {
        "PT": "VISUALIZAÇÃO",
        "EN": "VISUALIZATION"
    },
    "sidebar_theme": {
        "PT": "Tema do Gel",
        "EN": "Gel Theme"
    },
    "guide_title": {
        "PT": "Guia Rápido",
        "EN": "Quick Guide"
    },
    "guide_content": {
        "PT": """
        **Modos de Uso:**
        * **Digestão:** Upload do DNA + Enzimas.
        * **PCR:** Upload do Molde + Primers.
          - Forward: Insira 5'->3'.
          - Reverse: Insira 3'->5' (Exatamente como lido na fita antisense).
          - Suporta Overhangs (detecta pela extremidade 3').
        * **Ladder:** Marcadores de peso molecular.
        """,
        "EN": """
        **Modes:**
        * **Digestion:** DNA Upload + Enzymes.
        * **PCR:** Template Upload + Primers.
          - Forward: Input 5'->3'.
          - Reverse: Input 3'->5' (As read on antisense strand).
          - Supports Overhangs (3' seed detection).
        * **Ladder:** Molecular weight markers.
        """
    },
    "well_title": {
        "PT": "Poço",
        "EN": "Well"
    },
    "opt_sample": {
        "PT": "Digestão",
        "EN": "Digestion"
    },
    "opt_ladder": {
        "PT": "Ladder",
        "EN": "Ladder"
    },
    "opt_pcr": {
        "PT": "PCR",
        "EN": "PCR"
    },
    "sel_ladder": {
        "PT": "Selecione o Ladder:",
        "EN": "Select Ladder:"
    },
    "label_gel": {
        "PT": "Rótulo:",
        "EN": "Label:"
    },
    "tab_file": {
        "PT": "Upload Arquivo",
        "EN": "Upload File"
    },
    "tab_text": {
        "PT": "Digitar/Colar",
        "EN": "Type/Paste"
    },
    "upload_label": {
        "PT": "Arraste seu arquivo aqui",
        "EN": "Drag your file here"
    },
    "paste_label": {
        "PT": "Cole a sequência (FASTA ou pura)",
        "EN": "Paste sequence (FASTA or raw)"
    },
    "check_circular": {
        "PT": "Circular?",
        "EN": "Circular?"
    },
    "sel_enzymes": {
        "PT": "Enzimas",
        "EN": "Enzymes"
    },
    "pcr_fwd": {
        "PT": "Primer Forward (5'->3')",
        "EN": "Forward Primer (5'->3')"
    },
    "pcr_rev": {
        "PT": "Primer Reverse (3'->5')",
        "EN": "Reverse Primer (3'->5')"
    },
    "result_title": {
        "PT": "Resultado da Eletroforese",
        "EN": "Electrophoresis Result"
    },
    "export_expander": {
        "PT": "Exportar Dados",
        "EN": "Export Data"
    },
    "btn_download": {
        "PT": "Baixar .csv",
        "EN": "Download .csv"
    },
    "empty_msg": {
        "PT": "Para começar, adicione amostras nos cartões acima.",
        "EN": "To start, add samples in the cards above."
    },
    "created_by": {
        "PT": "Criado por",
        "EN": "Created by"
    },
    "lab_name": {
        "PT": "Laboratório de Biofármacos",
        "EN": "Biopharmaceuticals Laboratory"
    },
    "institute": {
        "PT": "Instituto Butantan",
        "EN": "Butantan Institute"
    },
    "pref_lang": {
        "PT": "Idioma / Language",
        "EN": "Language"
    },
    "report_bug": {
        "PT": "🐛 Reportar Problema",
        "EN": "🐛 Report Bug"
    },
    "warn_multiple": {
        "PT": "⚠️ INESPECIFICIDADE: Múltiplos sítios de ligação!",
        "EN": "⚠️ NON-SPECIFIC: Multiple binding sites!"
    },
    "warn_no_product": {
        "PT": "Nenhum produto (Verifique se o 3' anela corretamente)",
        "EN": "No product (Check if 3' anneals correctly)"
    }
}

# --- 3. ESTILO CSS (TURQUESA + MINIMALISTA) ---
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

    /* SLIDER PERSONALIZADO */
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
    
    /* ABAS VISÍVEIS */
    button[data-baseweb="tab"] {
        font-size: 13px !important;
        padding: 10px !important;
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
    
    .bug-report {
        font-size: 11px;
        color: #64748B;
        text-decoration: none;
        margin-top: 5px;
        display: inline-block;
    }
    .bug-report:hover {
        color: #0F766E;
        text-decoration: underline;
    }
    
    .warning-text {
        color: #DC2626;
        font-weight: bold;
        font-size: 12px;
        margin-top: 5px;
    }
</style>
""", unsafe_allow_html=True)

# --- 4. BACKEND (LÓGICA BIOLÓGICA & PCR PRO) ---

TODAS_ENZIMAS = sorted([str(e) for e in CommOnly])

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

def calcular_pcr_biologico(sequencia, fwd_seq, rev_seq, eh_circular):
    """
    Simula PCR com lógica de SEED 3' e Overhangs.
    Fwd: 5'->3'. Seed = final da string.
    Rev: 3'->5'. Seed = INÍCIO da string (pois 3' está na esquerda).
    """
    if not sequencia or sequencia.startswith("Erro"): return [], False
    template = sequencia.upper()
    fwd = "".join(fwd_seq.split()).upper()
    rev = "".join(rev_seq.split()).upper()
    
    if len(fwd) < 10 or len(rev) < 10: return [], False 

    SEED_SIZE = 15 # Tamanho da semente para anelamento
    
    # --- 1. Definir a Semente (Parte que anela no 3') ---
    # FWD (5'->3'): A parte que polimeriza é o final da string (3').
    fwd_seed = fwd[-SEED_SIZE:] if len(fwd) > SEED_SIZE else fwd
    
    # REV (3'->5'): O usuário digitou ao contrário. O 3' é o começo da string (índice 0).
    # Então a semente são os primeiros caracteres.
    rev_seed_3to5 = rev[:SEED_SIZE] if len(rev) > SEED_SIZE else rev
    
    # Para buscar no template (que é 5'->3'), precisamos inverter a semente do reverso
    # Ex: User 3'-ATGC-5' -> Busca 5'-CGTA-3' no template? Não.
    # Se o Primer é 3'-ATGC-5', ele pareia com 5'-TACG-3' no template.
    # Então buscamos a sequência exata COMPLEMENTAR no template?
    # O usuário pediu "Não inverter". Assumimos que ele digitou a sequência que PAREIA.
    # Se ele digitou 3'->5', a sequência correspondente no template (5'->3') é simplesmente a string invertida.
    # Ex: Primer 3'-A T G C-5' pareia com Template 5'-T A C G-3'.
    # Se o usuário digitou A T G C... ele quer que pareie com T A C G?
    # Vamos assumir a interpretação mais direta de "3'->5' visual":
    # O usuário vê T-A-C-G na fita de baixo (lendo esq-dir). Ele digita T-A-C-G.
    # Isso corresponde a A-T-G-C na fita de cima.
    # Então buscamos a string invertida no template 5'->3'.
    rev_seed_search = rev_seed_3to5[::-1]
    
    # --- 2. Encontrar Sítios ---
    fwd_matches = [m.start() for m in re.finditer(fwd_seed, template)]
    rev_matches = [m.start() for m in re.finditer(rev_seed_search, template)]
    
    produtos = []
    
    # --- 3. Parear Sítios ---
    for f_pos in fwd_matches:
        # Posição onde termina o anelamento do FWD (3' efetivo no template)
        f_3prime_end = f_pos + len(fwd_seed)
        
        for r_pos in rev_matches:
            # r_pos é o início do match no template.
            # Como a semente foi invertida, r_pos corresponde ao 5' do sítio no template
            # que pareia com o 3' do primer reverso.
            
            # Cálculo de distância interna (DNA a ser copiado entre os primers)
            
            # CASO LINEAR (Forward antes do Reverse)
            if r_pos > f_pos:
                distancia_interna = r_pos - f_3prime_end
                
                if distancia_interna >= 0:
                    # Tamanho Final = Tamanho Primer Fwd Completo + Tamanho Primer Rev Completo + Miolo
                    tamanho_total = len(fwd) + len(rev) + distancia_interna
                    produtos.append(tamanho_total)
            
            # CASO CIRCULAR (Reverse aparece "antes" no array linear, cruzando a origem)
            elif eh_circular and r_pos < f_pos:
                dist_fim = len(template) - f_3prime_end
                dist_inicio = r_pos
                distancia_interna = dist_fim + dist_inicio
                
                tamanho_total = len(fwd) + len(rev) + distancia_interna
                produtos.append(tamanho_total)
                
    tem_inespecificidade = len(produtos) > 1
    
    if not produtos: return [], False
    
    # Retorna lista de tuplas para plotagem
    return [(p, "PCR Product", p) for p in sorted(produtos, reverse=True)], tem_inespecificidade

# --- 5. INTERFACE DO USUÁRIO ---

if 'lang' not in st.session_state:
    st.session_state.lang = "PT"

with st.sidebar:
    st.markdown("""
    <div style="text-align: left; margin-bottom: 20px;">
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

    # Rodapé com Email de Reporte
    st.markdown(f"""
    <div style="font-size: 11px; color: #334155; line-height: 1.4; margin-top: 15px;">
        <strong>{TEXTS['created_by'][lang]} Elton Ostetti</strong><br>
        {TEXTS['lab_name'][lang]}<br>
        {TEXTS['institute'][lang]}<br><br>
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
                    "Poço": i+1,
                    "Identificação": rotulo_custom,
                    "Tipo": "Ladder",
                    "Detalhes": lad,
                    "Bandas (pb)": "; ".join([str(t) for t in LADDERS[lad]])
                })
            
            else:
                nomes_ladders.append(None)
                tab_f, tab_t = st.tabs([f"📂 {TEXTS['tab_file'][lang]}", f"📝 {TEXTS['tab_text'][lang]}"])
                seq, nome_arquivo = "", ""
                
                with tab_f:
                    up = st.file_uploader(TEXTS['upload_label'][lang], type=['dna', 'fasta', 'txt', 'fa'], key=f"u_{i}", label_visibility="collapsed")
                    if up: 
                        nome_arquivo, seq = processar_upload(up)
                        if nome_arquivo == "Erro": 
                            st.error(seq); seq = ""
                with tab_t:
                    txt = st.text_area(TEXTS['paste_label'][lang], height=70, key=f"tx_{i}", label_visibility="collapsed", placeholder="ATGC...")
                    if txt and not seq: 
                        nome_t, seq_t = processar_texto_manual(txt)
                        if nome_t != "Seq Manual": nome_arquivo = nome_t
                        seq = seq_t
                
                st.markdown("---")
                val_rotulo = nome_arquivo if nome_arquivo else str(i+1)
                
                if tipo == "Amostra":
                    circ = st.checkbox(TEXTS['check_circular'][lang], True, key=f"c_{i}")
                    enz = st.multiselect(TEXTS['sel_enzymes'][lang], TODAS_ENZIMAS, key=f"e_{i}")
                    rotulo_custom = st.text_input(TEXTS['label_gel'][lang], value=val_rotulo[:10], key=f"lbl_{i}")
                    labels_eixo_x.append(rotulo_custom)

                    if seq:
                        try:
                            res = calcular_digestao(seq, enz, circ)
                            dados_para_plotar.append(res)
                            
                            fragmentos_str = "; ".join([str(int(b[0])) for b in res])
                            desc_enzimas = ", ".join(enz) if enz else "Uncut"
                            relatorio_dados.append({
                                "Poço": i+1,
                                "Identificação": rotulo_custom,
                                "Tipo": "Digestão",
                                "Detalhes": desc_enzimas,
                                "Bandas (pb)": fragmentos_str
                            })
                        except Exception as e:
                            dados_para_plotar.append([])
                            st.error("Error")
                    else:
                        dados_para_plotar.append([])
                        relatorio_dados.append({"Poço": i+1, "Tipo": "Vazio", "Bandas (pb)": "-"})

                elif tipo == "PCR":
                    fwd = st.text_input(TEXTS['pcr_fwd'][lang], key=f"fwd_{i}", placeholder="ATGC... (5'->3')")
                    # REVERSE: 3'->5' (Sem template, apenas ordem)
                    rev = st.text_input(TEXTS['pcr_rev'][lang], key=f"rev_{i}", placeholder="ATGC... (3'->5')")
                    circ = st.checkbox(TEXTS['check_circular'][lang], False, key=f"cp_{i}")
                    
                    rotulo_custom = st.text_input(TEXTS['label_gel'][lang], value=f"PCR-{i+1}", key=f"lbl_{i}")
                    labels_eixo_x.append(rotulo_custom)
                    
                    if seq and fwd and rev:
                        try:
                            res, tem_inespecificidade = calcular_pcr_biologico(seq, fwd, rev, circ)
                            dados_para_plotar.append(res)
                            
                            if tem_inespecificidade:
                                st.markdown(f"<p class='warning-text'>{TEXTS['warn_multiple'][lang]}</p>", unsafe_allow_html=True)
                            
                            if not res:
                                st.warning(TEXTS['warn_no_product'][lang])
                            
                            fragmentos_str = "; ".join([str(int(b[0])) for b in res])
                            relatorio_dados.append({
                                "Poço": i+1,
                                "Identificação": rotulo_custom,
                                "Tipo": "PCR",
                                "Detalhes": f"Fwd(3'):..{fwd[-5:]} / Rev(3'):..{rev[:5]}",
                                "Bandas (pb)": fragmentos_str if res else "Nenhum"
                            })
                        except Exception as e:
                            dados_para_plotar.append([])
                            st.error(f"Error: {e}")
                    else:
                        dados_para_plotar.append([])
                        relatorio_dados.append({"Poço": i+1, "Tipo": "Vazio", "Bandas (pb)": "-"})

st.markdown(" ") 
st.markdown(f"### {TEXTS['result_title'][lang]}")

if any(dados_para_plotar):
    
    if "Neon" in estilo_gel:
        bg_color = '#111827'; text_color = 'white'; color_sample = '#00ff41'; color_ladder = '#ff9900'
    elif "Profissional" in estilo_gel: 
        bg_color = '#000000'; text_color = 'white'; color_sample = 'white'; color_ladder = 'white'
    else: 
        bg_color = 'white'; text_color = 'black'; color_sample = 'black'; color_ladder = 'black'

    # Margem de segurança para bandas pequenas
    min_view_calc = 50 + (100 * (agarose - 0.5))
    min_view = min_view_calc * 0.8
    max_view = 25000 / (agarose * 0.8)

    fig = go.Figure()

    for i, lista_bandas in enumerate(dados_para_plotar):
        x_center = i + 1
        eh_ladder = (nomes_ladders[i] is not None)
        cor_atual = color_ladder if eh_ladder else color_sample

        if lista_bandas:
             massa_total = sum([b[2] for b in lista_bandas]) if not eh_ladder else 1
        
        for (tam_aparente, tipo_banda, tam_real) in lista_bandas:
            # Filtro visual mais permissivo
            if tam_aparente < (min_view * 0.9) or tam_aparente > (max_view * 1.1): continue

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
        height=700,
        # Margens ajustadas para não cortar o eixo Y
        margin=dict(t=40, b=80, l=80, r=40),
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
    
    st.plotly_chart(fig, use_container_width=True)
    
    with st.expander(f"📥 {TEXTS['export_expander'][lang]}"):
        df_resultados = pd.DataFrame(relatorio_dados)
        csv = df_resultados.to_csv(index=False).encode('utf-8')
        st.download_button(
            label=TEXTS['btn_download'][lang],
            data=csv,
            file_name='gel_result.csv',
            mime='text/csv',
        )

else:
    st.info(TEXTS['empty_msg'][lang])

# --- RODAPÉ ---
st.markdown("""
<div class="footer">
    <p><b>BioSpark</b></p>
</div>
""", unsafe_allow_html=True)
