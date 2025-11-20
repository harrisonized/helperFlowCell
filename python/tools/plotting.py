import os
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from typing import List, Dict, Any
from .list_tools import items_in_a_not_b
from .math import (unpaired_t_test, fishers_lsd, tukey_multiple_comparisons,
                   bonferroni_multiple_comparisons, apply_unpaired_t_test,
                   generate_lognormal_data)


def save_fig(fig, dirpath: str = 'figures', filename: str = 'fig',
             height: int = 500, width: int = 800, scale: int = 3,
             save_html: bool = False):
    """Save a Plotly figure as PNG and optionally HTML."""
    png_dir = os.path.join(dirpath, 'png')
    os.makedirs(png_dir, exist_ok=True)

    fig.write_image(os.path.join(png_dir, f'{filename}.png'),
                    height=height, width=width, scale=scale)

    if save_html:
        html_dir = os.path.join(dirpath, 'html')
        os.makedirs(html_dir, exist_ok=True)
        fig.write_html(os.path.join(html_dir, f'{filename}.html'))


def plot_scatter(df: pd.DataFrame, x: str, y: str, group_by: str = None,
                 size: str = None, xlabel: str = None, ylabel: str = None,
                 title: str = None, xmin: float = None, xmax: float = None,
                 ymin: float = None, ymax: float = None,
                 yaxis_type: str = 'linear', color: str = '#1f77b4',
                 color_discrete_map: dict = None, hover_data: list = None,
                 mode: str = 'markers') -> go.Figure:
    """Create a scatter plot using Plotly."""
    if xmin is None:
        xmin = df[x].min()
    if xmax is None:
        xmax = df[x].max()
    xrange = [xmin - (xmax - xmin) * 0.1, xmax + (xmax - xmin) * 0.1]

    if ymin is None:
        ymin = df[y].min()
    if ymax is None:
        ymax = df[y].max()

    if yaxis_type == 'linear':
        yrange = [ymin - (ymax - ymin) * 0.1, ymax + (ymax - ymin) * 0.1]
    else:
        yrange = [np.log10(ymin) if ymin > 0 else None,
                  np.log10(ymax) if ymax > 0 else None]

    fig = go.Figure()

    if group_by is None:
        fig.add_trace(go.Scatter(
            x=df[x], y=df[y], mode=mode,
            marker=dict(size=5, color=color)
        ))
    else:
        for group in df[group_by].unique():
            mask = df[group_by] == group
            group_color = color_discrete_map.get(group) if color_discrete_map else None

            hovertext = []
            for idx in df[mask].index:
                text = '<br>'.join([f"{col}={df.loc[idx, col]}"
                                    for col in ([x, y] + (hover_data or []))])
                hovertext.append(text)

            fig.add_trace(go.Scatter(
                x=df.loc[mask, x], y=df.loc[mask, y],
                mode=mode, name=group,
                marker=dict(size=5, color=group_color),
                hovertext=hovertext
            ))

    fig.update_layout(
        title=dict(text=title, x=0.5),
        xaxis=dict(title=xlabel, range=xrange),
        yaxis=dict(title=ylabel, range=yrange, type=yaxis_type,
                   showgrid=True, gridcolor='#E4EAF2', zeroline=False),
        plot_bgcolor='rgba(0,0,0,0)',
        showlegend=group_by is not None,
        hovermode='closest'
    )

    return fig


def plot_violin(df: pd.DataFrame, x: str, y: str, group_by: str = None,
                xlabel: str = None, ylabel: str = None, title: str = None,
                ymin: float = None, ymax: float = None,
                legend_order: list = None, xaxis_angle: int = 45,
                yaxis_type: str = 'linear', color: str = '#1f77b4',
                color_discrete_map: dict = None, hover_data: list = None,
                sort: bool = True, descending: bool = True,
                violinmode: str = 'group') -> go.Figure:
    """Create a violin plot using Plotly."""
    if ymin is None:
        ymin = df[y].min()
    if ymax is None:
        ymax = df[y].max()

    if yaxis_type == 'linear':
        yrange = [ymin - (ymax - ymin) * 0.1, ymax + (ymax - ymin) * 0.1]
    else:
        yrange = [np.log10(ymin) if ymin > 0 else None,
                  np.log10(ymax) if ymax > 0 else None]

    if sort:
        order = df.groupby(x)[y].median().sort_values(ascending=not descending).index

    fig = go.Figure()

    if group_by is None:
        if sort:
            df = df.set_index(x).loc[order].reset_index()

        fig.add_trace(go.Violin(
            x=df[x], y=df[y],
            box_visible=True, meanline_visible=True,
            points='all', jitter=0.2, pointpos=-1,
            marker=dict(size=5, color=color)
        ))
    else:
        if legend_order is None:
            legend_order = sorted(df[group_by].unique())

        for group in legend_order:
            mask = df[group_by] == group
            group_df = df[mask]

            if sort:
                group_df = group_df.set_index(x).reindex(order).reset_index()

            group_color = color_discrete_map.get(group) if color_discrete_map else None

            fig.add_trace(go.Violin(
                x=group_df[x], y=group_df[y],
                name=group, legendgroup=group, scalegroup=group,
                box_visible=True, meanline_visible=True,
                points='all', jitter=0.2, pointpos=-1,
                marker=dict(size=5, color=group_color)
            ))

    fig.update_layout(
        title=dict(text=title, x=0.5),
        xaxis=dict(title=xlabel, tickangle=xaxis_angle),
        yaxis=dict(title=ylabel, range=yrange, type=yaxis_type,
                   showgrid=True, gridcolor='#E4EAF2', zeroline=False),
        violinmode=violinmode,
        plot_bgcolor='rgba(0,0,0,0)',
        showlegend=group_by is not None,
        hovermode='closest'
    )

    return fig


def get_significance_code(p: float, digits: int = 3) -> str:
    """Convert p-value to significance symbol."""
    if pd.isna(p):
        return 'NA'

    p = abs(p)
    if p > 0.20:
        return 'ns'
    elif p > 0.05:
        return f'{p:.{digits}f}'
    elif p > 0.01:
        return '*'
    elif p > 0.001:
        return '**'
    elif p > 0.0001:
        return '***'
    else:
        return '****'


def generate_base_level(n_groups: int) -> list:
    """Generate starting positions for significance brackets."""
    if n_groups == 1:
        return [0]
    elif n_groups == 2:
        return [1]

    mid = n_groups // 2
    if n_groups % 2 == 0:
        nlevels = list(range(1, mid + 1)) + list(range(mid - 1, 0, -1))
    else:
        nlevels = list(range(1, mid + 1)) + list(range(mid, 0, -1))

    res = []
    value = 1
    nreps = n_groups - 1

    for step in nlevels:
        res.extend([value] * nreps)
        value += step
        nreps -= 1

    return res


def plot_multiple_comparisons(df: pd.DataFrame, x: str, y: str,
                              xlabel: str = None, ylabel: str = None,
                              title: str = None, ymin: float = 0,
                              xaxis_angle: int = 60, test: str = 't_test',
                              show_numbers: bool = False, digits: int = 4,
                              show_brackets: bool = True,
                              custom_group_order: list = None):
    """Create boxplot with significance brackets using matplotlib."""
    import matplotlib.pyplot as plt
    import seaborn as sns
    from itertools import combinations

    df = df[~df[y].isna()].copy()

    if custom_group_order:
        group_names = [g for g in custom_group_order if g in df[x].unique()]
        df = df[df[x].isin(group_names)]
        df[x] = pd.Categorical(df[x], categories=group_names, ordered=True)
    else:
        group_names = sorted(df[x].unique())

    n_groups = len(group_names)

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6))

    # Create boxplot
    colors = plt.cm.Dark2(np.linspace(0, 1, min(8, n_groups)))
    if n_groups > 8:
        colors = plt.cm.tab20(np.linspace(0, 1, n_groups))

    bp = ax.boxplot([df[df[x] == g][y].dropna() for g in group_names],
                    positions=range(1, n_groups + 1),
                    patch_artist=True, widths=0.6)

    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    # Add jittered points
    for i, group in enumerate(group_names, 1):
        group_data = df[df[x] == group][y].dropna()
        jitter = np.random.uniform(-0.15, 0.15, len(group_data))
        ax.scatter(i + jitter, group_data, alpha=0.6, s=20, color='black')

    # Add mean line
    for i, group in enumerate(group_names, 1):
        mean_val = df[df[x] == group][y].mean()
        ax.hlines(mean_val, i - 0.3, i + 0.3, colors='black',
                  linestyles='dashed', linewidth=1)

    if show_brackets and n_groups > 1:
        # Compute p-values
        if test == 't_test':
            pval_result = apply_unpaired_t_test(
                df, index_cols=items_in_a_not_b(list(df.columns), [x, y]),
                group_name=x, metric=y, custom_group_order=custom_group_order
            )
            pvals = {}
            if pval_result is not None:
                for col in pval_result.columns:
                    if '-' in col:
                        pvals[col] = pval_result[col].iloc[0] if len(pval_result) > 0 else np.nan
        elif test == 'fishers_lsd':
            pvals = fishers_lsd(df, group=x, metric=y, custom_group_order=custom_group_order)
        elif test == 'tukey':
            pvals = tukey_multiple_comparisons(df, group=x, metric=y, custom_group_order=custom_group_order)
        elif test == 'bonferroni':
            pvals = bonferroni_multiple_comparisons(df, group=x, metric=y, custom_group_order=custom_group_order)
        else:
            raise ValueError("Choose from test='t_test', 'fishers_lsd', 'tukey', or 'bonferroni'")

        # Calculate bracket positions
        pairs = list(combinations(range(n_groups), 2))
        bracket_params = pd.DataFrame(pairs, columns=['left', 'right'])
        bracket_params['dist'] = bracket_params['right'] - bracket_params['left']
        bracket_params = bracket_params.sort_values('dist')
        bracket_params['base_level'] = generate_base_level(n_groups)
        bracket_params['level'] = (bracket_params['base_level'] +
                                   (bracket_params['left']) % bracket_params['dist'])

        # Compute starting height
        y_data = df[y].dropna()
        h_low = y_data.max() * 1.1
        space = h_low / 10

        # Draw brackets
        for _, row in bracket_params.iterrows():
            left = int(row['left']) + 1
            right = int(row['right']) + 1
            level = row['level']
            col_name = f"{group_names[int(row['right'])]}-{group_names[int(row['left'])]}"

            if show_numbers:
                pval_text = f"{pvals.get(col_name, np.nan):.{digits}f}"
            else:
                pval_text = get_significance_code(pvals.get(col_name, np.nan))

            y_pos = h_low + (level - 1) * space

            # Draw bracket
            ax.plot([left + 0.1, left + 0.1, right - 0.1, right - 0.1],
                    [y_pos - space * 0.1, y_pos, y_pos, y_pos - space * 0.1],
                    'k-', linewidth=0.8)

            # Add text
            ax.text((left + right) / 2, y_pos + space * 0.1, pval_text,
                    ha='center', va='bottom', fontsize=8 if n_groups <= 6 else 6)

    ax.set_xticks(range(1, n_groups + 1))
    ax.set_xticklabels(group_names, rotation=xaxis_angle, ha='right')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_ylim(bottom=ymin)

    plt.tight_layout()
    return fig


def plot_modal_histograms(df: pd.DataFrame, group: str = 'group',
                          value: str = 'value', xlabel: str = None,
                          ylabel: str = None, title: str = None,
                          x_ticks: list = None, colors: list = None,
                          nbins: int = 100, spar: float = 0.33,
                          show_bins: bool = False, max_scale: float = 8289.7211,
                          pos_decades: float = 4.5, lin_width: float = 2.25,
                          extra_neg: float = 0, width_basis: float = 0.25):
    """Plot modal histograms mimicking FlowJo style."""
    from scipy.interpolate import UnivariateSpline
    from scipy.ndimage import gaussian_filter1d

    if x_ticks is None:
        x_ticks = [10, 250, 500, 1000, 2000, 4000, 8000]

    # Simple logicle-like transform (simplified)
    def logicle_transform(x, w=lin_width, t=max_scale, m=pos_decades):
        """Simplified logicle transform."""
        x = np.asarray(x)
        result = np.zeros_like(x, dtype=float)
        mask = x > 0
        result[mask] = np.log10(x[mask] + 1) / m
        return result

    fig, ax = plt.subplots(figsize=(10, 4))

    groups = df[group].unique()
    n_groups = len(groups)

    if colors is None or len(colors) < n_groups:
        colors = plt.cm.tab10(np.linspace(0, 1, n_groups))

    for i, grp in enumerate(groups):
        grp_data = df[df[group] == grp][value].dropna()

        if len(grp_data) == 0:
            continue

        # Transform data
        transformed = logicle_transform(grp_data)

        # Compute histogram
        hist, bins = np.histogram(transformed, bins=nbins, density=True)
        bin_centers = (bins[:-1] + bins[1:]) / 2

        # Normalize to mode
        hist = hist / hist.max() if hist.max() > 0 else hist

        # Smooth
        hist_smooth = gaussian_filter1d(hist, sigma=2)

        # Plot
        color = colors[i] if isinstance(colors[i], str) else colors[i]
        ax.fill_between(bin_centers, hist_smooth, alpha=0.3, color=color, label=grp)
        ax.plot(bin_centers, hist_smooth, color=color, linewidth=1.2)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_ylim(0, 1)
    ax.legend(loc='upper right')

    plt.tight_layout()
    return fig
