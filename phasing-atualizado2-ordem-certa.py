import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

arquivo = "USE-ESSE-PARA-GRAFICOS.xlsx"
df = pd.read_excel(arquivo, sheet_name="Comparacao")

df = df.copy()
df["PRU"] = pd.to_numeric(df["PRU"], errors="coerce")
df = df.dropna(subset=["PRU", "diplotipo_minha_norm"])


# ------------------------------------------------------------------
# PARSE BÁSICO DO HAPLOTOPO
# ------------------------------------------------------------------
def parse_star_and_suffix(hap_str):
    """
    Recebe algo como '1CG', '2CG', '1TG', '17CG'
    Retorna (star_num, suffix), ex.: ('1', 'CG') ou ('1', 'TG')
    """
    if pd.isna(hap_str):
        return (None, None)
    hap_str = str(hap_str).strip().upper()
    m = re.match(r"^(\d+)([A-Z]+)$", hap_str)
    if not m:
        return (None, None)
    star_num, suffix = m.group(1), m.group(2)
    return star_num, suffix


# ------------------------------------------------------------------
# GRUPO A: SOMENTE ESTRELAS (*1/*2, etc.)
# ------------------------------------------------------------------
def diplotipo_estrelas_somente(diplotipo_norm):
    """
    Converte '1CG|2CG' -> '*1/*2'
    (ignora completamente TG/CG/TA etc, usa só o número)
    """
    if pd.isna(diplotipo_norm):
        return None
    parts = str(diplotipo_norm).upper().split("|")
    if len(parts) != 2:
        return None
    star1, _ = parse_star_and_suffix(parts[0])
    star2, _ = parse_star_and_suffix(parts[1])
    if star1 is None or star2 is None:
        return None
    estrelas = sorted([f"*{star1}", f"*{star2}"])
    return "/".join(estrelas)


# ------------------------------------------------------------------
# GRUPO B: CYP2C19 + TG  (COM TEXTO "CG ou TA" QUANDO HOUVER TG)
# ------------------------------------------------------------------
def diplotipo_com_TG(diplotipo_norm):
    """
    Converte, por exemplo:
        '1CG|1CG'  -> '*1/*1'
        '1CG|1TG'  -> '*1(CG ou TA)/*1TG'
        '2CG|17CG' -> '*2/*17'
        '2TG|17CG' -> '*2TG/*17(CG ou TA)'
    Regra:
      - Se NÃO houver nenhum TG: apenas *n/*m
      - Se houver TG em 1 alelo: o alelo NÃO-TG recebe '(CG ou TA)'
      - Se houver TG nos 2 alelos: *nTG/*mTG (sem '(CG ou TA)')
    """
    if pd.isna(diplotipo_norm):
        return None
    parts = str(diplotipo_norm).upper().split("|")
    if len(parts) != 2:
        return None

    parsed = []
    for hap in parts:
        star, suffix = parse_star_and_suffix(hap)
        if star is None:
            return None
        is_TG = (suffix == "TG")
        parsed.append((star, is_TG))

    (star1, tg1), (star2, tg2) = parsed
    any_TG = tg1 or tg2

    def fmt(star, is_TG, any_TG):
        base = f"*{star}"
        if is_TG:
            return base + "TG"
        else:
            # Só adiciona "(CG ou TA)" se existir TG no outro alelo
            if any_TG:
                return base + "(CG ou TA)"
            else:
                return base

    label1 = fmt(star1, tg1, any_TG)
    label2 = fmt(star2, tg2, any_TG)

    # Ordena apenas dentro do diplótipo, para ter representação consistente
    labels = sorted([label1, label2])
    return "/".join(labels)


df["grupo_estrelas"] = df["diplotipo_minha_norm"].apply(diplotipo_estrelas_somente)
df["grupo_TG"]       = df["diplotipo_minha_norm"].apply(diplotipo_com_TG)

df1 = df.dropna(subset=["grupo_estrelas"]).copy()
df2 = df.dropna(subset=["grupo_TG"]).copy()


# ------------------------------------------------------------------
# FUNÇÕES PARA ORDENAR PELO FENÓTIPO (LENTO -> NORMAL -> RÁPIDO)
# ------------------------------------------------------------------

# Mapeamento de atividade aproximado:
#   *2, *3 = perda de função (lento)
#   *1     = normal
#   *17    = ganho de função (rápido)
activity_map = {
    "2": 0.0,
    "3": 0.0,
    "1": 1.0,
    "17": 2.0,
}
lof_set = {"2", "3"}  # alelos de perda de função
tg_bonus = 2.0        # "peso extra" para alelos com TG (mais rápidos)


def sort_key_estrelas(label):
    """
    label tipo '*1/*2'
    Ordena por:
      1) soma da atividade (menor = mais lento, vai primeiro)
      2) dentro da mesma soma, diplótipo com mais LOF (mais *2/*3) vem mais à esquerda
      3) por fim, ordena pelos números das estrelas
    Isso gera, por exemplo:
      *1/*2  (mais lento)
      *2/*17
      *1/*1
      *1/*17
      *17/*17 (mais rápido)
    e se existir *3, ele entra junto com os lentos.
    """
    parts = label.replace("*", "").split("/")
    stars = [p for p in parts if p]

    activities = [activity_map.get(s, 1.0) for s in stars]
    total_activity = sum(activities)
    lof_count = sum(1 for s in stars if s in lof_set)
    stars_sorted = sorted(int(s) for s in stars)

    # menor atividade primeiro; mais LOF primeiro; depois número da estrela
    return (total_activity, -lof_count, stars_sorted)


def sort_key_TG(label):
    """
    label tipo:
        '*1/*1'
        '*1(CG ou TA)/*1TG'
        '*2/*17'
        '*2(CG ou TA)/*17TG'
        '*1TG/*17TG'
    Regras:
      1) primeiro diplótipos sem TG
         depois com 1 TG
         depois com 2 TG
      2) dentro de cada grupo, ordenar pela atividade (incluindo 'bônus' de TG)
      3) em empates, mais LOF à esquerda
    """
    # Remove o "(CG ou TA)" para facilitar o parse
    clean = label.replace("(CG ou TA)", "")
    alleles = clean.split("/")

    stars = []
    tg_flags = []

    for a in alleles:
        a = a.strip()
        m = re.match(r"^\*(\d+)(TG)?$", a)
        if not m:
            continue
        star = m.group(1)
        is_tg = m.group(2) == "TG"
        stars.append(star)
        tg_flags.append(is_tg)

    n_TG = sum(tg_flags)

    activities = []
    lof_count = 0
    for star, is_tg in zip(stars, tg_flags):
        base = activity_map.get(star, 1.0)
        if is_tg:
            base += tg_bonus
        activities.append(base)
        if star in lof_set:
            lof_count += 1

    total_activity = sum(activities)
    stars_sorted = sorted(int(s) for s in stars)

    # n_TG determina o "bloco" (0 TG -> 1 TG -> 2 TG)
    return (n_TG, total_activity, -lof_count, stars_sorted)


# ------------------------------------------------------------------
# FUNÇÃO GERAL DO GRÁFICO
# ------------------------------------------------------------------
def grafico_pru_por_grupo(df_in, coluna_grupo, titulo, sort_key_func=None):
    """
    df_in: dataframe com colunas [PRU, coluna_grupo]
    coluna_grupo: 'grupo_estrelas' ou 'grupo_TG'
    titulo: título do gráfico
    sort_key_func: função que recebe o rótulo do grupo e retorna a chave de ordenação
    """
    grupos = (
        df_in.groupby(coluna_grupo)["PRU"]
        .apply(list)
        .to_dict()
    )

    # Define ordem: se não tiver função, usa ordem alfabética
    if sort_key_func is None:
        grupos_ordenados = dict(sorted(grupos.items(), key=lambda x: x[0]))
    else:
        grupos_ordenados = dict(
            sorted(grupos.items(), key=lambda x: sort_key_func(x[0]))
        )

    group_labels = []
    data_values = []

    # Só entram diplótipos que REALMENTE existem (N>0)
    for geno, valores in grupos_ordenados.items():
        n = len(valores)
        group_labels.append(f"(N={n}) {geno}")
        data_values.append(valores)

    # Monta vetor "long" apenas para o jitter
    categories = []
    for label, valores in zip(group_labels, data_values):
        categories.extend([label] * len(valores))

    data_combined = np.concatenate([np.array(v) for v in data_values])
    category_feature = np.array(categories)

    df_plot = pd.DataFrame({
        "variavel_numerica": data_combined,
        "variavel_categorica": category_feature
    })

    fig, ax = plt.subplots(figsize=(12, 6))

    boxplot_data = []
    positions = np.arange(1, len(group_labels) + 1)

    for i, label in enumerate(group_labels):
        y = df_plot[df_plot["variavel_categorica"] == label]["variavel_numerica"]
        x = np.random.normal(i + 1, 0.04, size=len(y))
        ax.scatter(x, y, alpha=0.8)
        boxplot_data.append(y.values)

    for i, data in enumerate(boxplot_data):
        if data.size > 0:
            ax.boxplot(
                data,
                positions=[positions[i]],
                widths=0.6,
                patch_artist=True,
                boxprops=dict(facecolor=(1, 1, 1, 0), edgecolor='black'),
                medianprops=dict(color='darkblue', linewidth=2)
            )

    # Linha de corte PRU = 208
    ax.axhline(y=208, linestyle='--', linewidth=1.5)

    ax.set_ylim(0, 400)

    ax.set_title(titulo)
    ax.set_ylabel("Unidades de reação P2Y12 (PRU)")

    ax.set_xticks(positions)
    ax.set_xticklabels(group_labels, rotation=45, ha='right')

    ax.set_aspect(aspect=0.02)
    fig.subplots_adjust(bottom=0.25)

    plt.show()


# ------------------------------------------------------------------
# CHAMADAS FINAIS (A) e (B)
# ------------------------------------------------------------------

# (A) CYP2C19 – somente estrelas, ordem por fenótipo
grafico_pru_por_grupo(
    df1,
    "grupo_estrelas",
    titulo="CYP2C19",
    sort_key_func=sort_key_estrelas
)

# (B) CYP2C19 + TG – mantém blocos 0 TG, 1 TG, 2 TG e aplica '(CG ou TA)'
grafico_pru_por_grupo(
    df2,
    "grupo_TG",
    titulo="CYP2C19 + CYP2C:TG",
    sort_key_func=sort_key_TG
)
