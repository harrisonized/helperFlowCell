id_cols = ['filepath', 'filename', 'fcs_name']

initial_gates = [
    'Ungated', 'Cells', 'Cells/Single Cells',
    'Cells/Single Cells/Single Cells',
    'Cells/Single Cells/Single Cells/Live Cells'
]

cell_type_spell_check = {
    'B cells': 'B Cells',
    'NK cells': 'NK Cells'
}

cell_type_ignore = [
    'Large Cells', 'Small Cells',
    'CD3-CD19-', 'Ly6C-Ly6G-', 'F4_80-CD11b-', 'Ly6G-', 'NK1.1-',
    'Myeloid Cells', 'Non-Dendritic Cells', 'Non-Macrophages',
    'mNeonGreen-',
    'CD19-CD3-', 'CD19-CD90-', 'CD11b+', 'CD117-',
    'CD4+ T cells', 'CD8+ T cells'
]

mouse_db_ignore = [
    'system_id', 'rack', 'position', 'alive',
    'pcr_confirmation', 'color', 'ignore'
]
