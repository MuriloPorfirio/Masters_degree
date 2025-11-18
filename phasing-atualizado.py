import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

arquivo = "USE-ESSE-PARA-GRAFICOS.xlsx"
df = pd.read_excel(arquivo, sheet_name="Comparacao")

df = df.copy()
df["PRU"] = pd.to_numeric(df["PRU"], errors="coerce")
df = df.dropna(subset=["PRU", "diplotipo_minha_norm"])


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

def diplotipo_com_TG(diplotipo_norm):
    """
    Converte '1CG|1TG' -> '*1/*1TG'
            '2CG|17CG' -> '*2/*17'
            '1TG|17CG' -> '*1TG/*17'
    Considera TG se o sufixo for 'TG'.
    """
    if pd.isna(diplotipo_norm):
        return None
    parts = str(diplotipo_norm).upper().split("|")
    if len(parts) != 2:
        return None

    def rotulo_hap(hap_str):
        star, suffix = parse_star_and_suffix(hap_str)
        if star is None:
            return None
        base = f"*{star}"
        if suffix == "TG":
            return base + "TG"
        else:
            return base

    label1 = rotulo_hap(parts[0])
    label2 = rotulo_hap(parts[1])
    if label1 is None or label2 is None:
        return None
    labels = sorted([label1, label2])
    return "/".join(labels)

df["grupo_estrelas"] = df["diplotipo_minha_norm"].apply(diplotipo_estrelas_somente)
df["grupo_TG"]       = df["diplotipo_minha_norm"].apply(diplotipo_com_TG)

df1 = df.dropna(subset=["grupo_estrelas"]).copy()
df2 = df.dropna(subset=["grupo_TG"]).copy()

def grafico_pru_por_grupo(df_in, coluna_grupo, titulo):
    """
    df_in: dataframe com colunas [PRU, coluna_grupo]
    coluna_grupo: 'grupo_estrelas' ou 'grupo_TG'
    titulo: título do gráfico
    """
    grupos = (
        df_in.groupby(coluna_grupo)["PRU"]
        .apply(list)
        .to_dict()
    )

    grupos_ordenados = dict(sorted(grupos.items(), key=lambda x: x[0]))

    group_labels = []
    data_values = []

    for geno, valores in grupos_ordenados.items():
        n = len(valores)
        group_labels.append(f"(N={n}) {geno}")
        data_values.append(valores)

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

    ax.axhline(y=208, linestyle='--', linewidth=1.5)

    # limites do eixo Y
    ax.set_ylim(0, 400)

    ax.set_title(titulo)
    ax.set_ylabel("Unidades de reação P2Y12 (PRU)")

    ax.set_xticks(positions)
    ax.set_xticklabels(group_labels, rotation=45, ha='right')

    ax.set_aspect(aspect=0.02)
    fig.subplots_adjust(bottom=0.25)

    plt.show()


grafico_pru_por_grupo(df1, "grupo_estrelas", titulo="CYP2C19")

grafico_pru_por_grupo(df2, "grupo_TG", titulo="CYP2C19 + CYP2C:TG")
