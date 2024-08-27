# Lists for flow analysis

## Objects
## id_cols
## numerical_cols
## ignored_cell_types
## cell_type_replacements
## flowjo_metadata_cols


id_cols <- c('filepath', 'filename', 'fcs_name')

numerical_cols <- c('Count', 'Cells', 'Cells/Single Cells',
    'Cells/Single Cells/Single Cells', 'Cells/Single Cells/Single Cells/Live Cells'
)

ignored_cell_types <- c(
    'Large Cells', 'Small Cells', 'Myeloid Cells', 'Non-Dendritic Cells', 'Non-Macrophages',
    'CD3-CD19-', 'Ly6C-Ly6G-', 'F4_80-CD11b-', 'Ly6G-', 'NK1.1-'
)

cell_type_replacements <- c('B cells'='B Cells',
    'NK cells'='NK Cells',
    'Ly6C-lo Monocytes'='Ly6C-int Monocytes'
)

#' to be deprecated
#' 
flowjo_metadata_cols <- c(
    'organ',
    'mouse_id',
    'treatment_group',
    'strain'
)
