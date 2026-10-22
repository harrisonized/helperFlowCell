min_cells <- 3  # do not set this less than 3

# exclude from ribbon plots or polar area charts
composite_gates <- c(
    "Single Cells", "Live Cells", "CD45+",
    "Non-DCs", "DCs",
    "CD115+", "Monocytes",
    "T cells", "MCs",
     "CD45-hi"
)

cell_type_to_color <- list(
    'Ly6C-hi Monocytes'  = '#2ca02c',  # green
    'Ly6C-hi Mono'       = '#2ca02c',  # green
    'Ly6C-int Monocytes' = '#3b82f6',  # blue 500
    'Ly6C-int Mono'      = '#3b82f6',  # blue 500
    'Ly6C-lo Monocytes'  = '#ef4444',  # red 500
    'Ly6C-lo Mono'       = '#ef4444',  # red 500
    'Neutrophils'        = '#1f77b4',  # blue
    'Neu'                = '#1f77b4',  # blue
    'MoMacs'             = '#ff7f0e',  # orange
    'Macrophages'        = '#d62728',  # red
    'AM'                 = '#d62728',  # red
    'IM'                 = '#f59e0b',  # amber 500
    'Eosinophils'        = "#10b981",  # emerald 500
    'DCs'                = '#17becf',  # cyan
    'DC1'                = '#17becf',  # cyan
    'DC2'                = '#6366f1',  # indigo 500
    'pDCs'               = '#a855f7',  # purple 500
    'NK Cells'           = '#bcbd22',  # olive
    'B Cells'            = '#8c564b',  # brown
    'B cells'            = '#8c564b',  # brown
    'CD4+ T cells'       = '#e377c2',  # pink
    'CD8+ T cells'       = '#9467bd',  # purple
    'Other'              = '#7f7f7f'   # gray
)

custom_group_order <- c(
    'F, WT', 'F, het', 'F, homo',
    'M, WT', 'M, hemi'
)


# ----------------------------------------------------------------------
# Examples

# custom_group_order <- c(
#     'F, DMSO, WT', 'F, DMSO, het', 'F, DMSO, homo',
#     'F, R848, WT', 'F, R848, homo',
#     'M, DMSO, WT', 'M, DMSO, hemi',
#     'M, R848, WT', 'M, R848, hemi'
# )


# custom_group_order <- c(
#     'F, NT, WT', 'F, NT, het',
#     'F, PBS, WT', 'F, PBS, het',
#     'F, TGA, WT', 'F, TGA, het',
# )


# custom_group_order <- c(
#     'F, WT', 'F, het', 'F, homo',
#     'M, WT', 'M, hemi'
# )
