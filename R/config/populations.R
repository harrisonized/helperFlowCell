# Flow populations

## Objects
## bm_populations
## pb_populations
## pc_populations
## spleen_populations
## populations_for_organ


bm_populations <- c(
    'Cells/Single Cells/Single Cells/Live Cells',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/B cells',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/T cells',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/T cells/CD4+ T cells',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/T cells/CD8+ T cells',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/Myeloid Cells/Ly6G-/NK cells',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/Myeloid Cells/Neutrophils',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/Myeloid Cells/Ly6G-/NK1.1-/Non-Dendritic Cells/Ly6C-hi Monocytes',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/Myeloid Cells/Ly6G-/NK1.1-/Non-Dendritic Cells/Ly6C-int Monocytes',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/Myeloid Cells/Ly6G-/NK1.1-/Non-Dendritic Cells/Ly6C-/Macrophages',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/Myeloid Cells/Ly6G-/NK1.1-/Dendritic Cells'
)

pb_populations <- c(
    'Cells/Single Cells/Single Cells/Live Cells',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/B cells',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/T cells',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/T cells/CD4+ T cells',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/T cells/CD8+ T cells',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/Myeloid Cells/Ly6G-/NK cells',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/Myeloid Cells/Neutrophils',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/Myeloid Cells/Ly6G-/NK1.1-/Non-Dendritic Cells/Monocytes',  # rename this
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/Myeloid Cells/Ly6G-/NK1.1-/Non-Dendritic Cells/Ly6C-lo',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/Myeloid Cells/Ly6G-/NK1.1-/Non-Dendritic Cells/Ly6C-/Macrophages',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/Myeloid Cells/Ly6G-/NK1.1-/Dendritic Cells'
)

pc_populations <- c(
    'Cells/Single Cells/Single Cells/Live Cells',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/B cells',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/T cells',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/T cells/CD4+ T cells',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/T cells/CD8+ T cells',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/Myeloid Cells/Ly6G-/NK cells',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/Myeloid Cells/Neutrophils',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/CD3-CD19-/Ly6G-Ly6C-/NK1.1-/F4_80-hi Macrophages',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/CD3-CD19-/Ly6G-Ly6C-/NK1.1-/F4_80-int Macrophages'
)

spleen_populations <- c(
    'Cells/Single Cells/Single Cells/Live Cells',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/B cells',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/T cells',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/T cells/CD4+ T cells',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/T cells/CD8+ T cells',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/Myeloid Cells/Ly6G-/NK cells',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/Myeloid Cells/Neutrophils',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/Myeloid Cells/Ly6G-/NK1.1-/Non-Dendritic Cells/Ly6C+ Monocytes',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/Myeloid Cells/Ly6G-/NK1.1-/Non-Dendritic Cells/Ly6C-lo',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/Myeloid Cells/Ly6G-/NK1.1-/Non-Dendritic Cells/Ly6C-/Macrophages',
    'Cells/Single Cells/Single Cells/Live Cells/CD45+/Myeloid Cells/Ly6G-/NK1.1-/Dendritic Cells'
)  # same as pb


populations_for_organ <- list2env(list(
    'bm'= bm_populations,
    'pb'= pb_populations,
    'pc'= pc_populations,
    'spleen'= spleen_populations
))
