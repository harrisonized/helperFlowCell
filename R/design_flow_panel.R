## Assist in designing flow panels

wd = dirname(this.path::here())  # wd = '~/github/R/harrisonRTools'
library('readxl')
library('tidyr')
suppressMessages(library('dplyr'))
library('optparse')
library('logr')
source(file.path(wd, 'R', 'functions', 'file_io.R'))  # read_excel_or_csv
source(file.path(wd, 'R', 'functions', 'df_tools.R'))  # rename_columns
source(file.path(wd, 'R', 'functions', 'text_tools.R'))  # title_to_snake_case
source(file.path(wd, 'R', 'functions', 'list_tools.R'))  # find_first_match_index, multiple_replacement


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-file"), default='data/flow/panel.txt',
                metavar='data/flow/panel.txt', type="character",
                help="specify the antibodies to use in your flow panel"),

    make_option(c("-a", "--antibody-inventory"), default='data/flow/antibody_inventory.xlsx',
                metavar='data/flow/antibody_inventory.xlsx', type="character",
                help="antibody inventory table"),

    make_option(c("-c", "--instrument-config"), default='data/flow/instrument_config.csv',
                metavar='data/flow/instrument_config.csv', type="character",
                help="instrument configuration file"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt[['troubleshooting']]
# troubleshooting = TRUE

if (!troubleshooting) {
    # note: except for the last file, all files are output in the troubleshooting folder
    data_dir = file.path(wd, dirname(opt[['input-file']]), 'troubleshooting')
    if (!dir.exists(file.path(data_dir, 'lists'))) {
        dir.create(file.path(data_dir, 'lists'), recursive=TRUE)
    }
}

# Start Log
start_time = Sys.time()
log <- log_open(paste0("design_flow_panel-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Configs

col_to_new_col = c(
    'id'='channel_id',
    'fluorochrome_detected'='fluorophore',
    'fluorochrome'='fluorophore'
)


antibody_replacements <- c(

    # special characters
    'α' = 'a',
    'β' = 'b',
    'γ' = 'g',
    '™' = '',
    '([A-Za-z0-9],)\\s*([A-Za-z0-9])' = '\\1 \\2',  # comma-separated list

    # animals
    '[Gg][Oo][Aa][Tt]' = 'goat',
    '[Mm][Oo][Uu][Ss][Ee]' = 'mouse',
    '[Rr][Aa][Tt]' = 'rat',
    '[Dd][Oo][Nn][Kk][Ee][Yy]' = 'donkey',
    '[Bb][Oo][Vv][Ii][Nn][Ee]' = 'bovine',
    '[Ss][Hh][Ee][Ee][Pp]' = 'sheep',

    # general replacements
    '[Bb][Cc][Ll1]-{0,1}' = 'Bcl-',
    '.*^[Cc][Dd]' = 'CD',
    'FoxP3' = 'Foxp3',
    'gammaDelta' = 'gd',
    'INFg' = 'IFNg',
    'Immunoglobulin ' = 'Ig',
    '.*^[Ii][Ll]' = 'IL',
    '(IL)([0-9]+)' = '\\1-\\2',
    '.*^Ly-' = 'Ly',
    'MHC {0,1}(I{1,2})' = 'MHC\\1',
    '.*^[Nn][Oo]tch' = 'Notch',
    '.*^[Oo]nly' = 'Only',
    '[Rr][Oo][Rr][gγy][yt]' = 'RORgt',
    '.*^[Tt][Cc][Rr]' = 'TCR',
    '[Tt][Cc][Rr][Bb-]\\w*' = 'TCRb',
    '.*^[Tt][Dd][Tt]' = 'TdT',
    '[Tt][Gg][Ff][Bb-]\\w*' = 'TGFb',
    'Unlabel{1,2}ed' = 'Unlabeled',
    'Va(lpha|) {0,1}([0-9]+)' = 'Va\\2',
    'Vbeta' = 'Vb',
    
    # custom exact replacements
    '^\\(CXCR4\\)$' = 'CXCR4',
    ' Fixable Viability Kit' = '',
    'NK cell Pan' = 'CD49b',
    'TCRb TCRcb' = 'TCRb, TCRcb',
    'Vb8.1 Vb8.2' = 'Vb8.1, 2',
    ' {0,1}\\(Tonegawa nomenclat\\)' = ''
)


fluorophore_replacements <- c(

    # special characters
    '®' = '',
    '/' = '-',
    '^ ([A-Za-z]+)' = '\\1',

    # general replacements
    '(Alexa) {0,1}(Fluor|) {0,1}' = 'AF',
    'A(F|) {0,1}([0-9]+)' = 'AF\\2',
    '[Aa][Pp][Cc]' = 'APC',
    'APC {0,1}-{0,1}([A-Za-z]+)' = 'APC-\\1',
    '[Ff]ire {0,1}([0-9]+)' = 'Fire \\1',
    '^APC-[Ff]ire$' = 'APC-Fire 750',
    '^Atto-' = 'ATTO ',
    '[Bb][Ii][Oo](tin|)(,|)( |)([A-Za-z]+|)' = 'Biotin',
    '(BU[Vv]) {0,1}([0-9]+)' = 'BUV\\2',
    '(B[Vv]) {0,1}([0-9]+)' = 'BV\\2',
    '([Dd][Ll]|Dy[Ll]ight) {0,1}-{0,1}([0-9]+)' = 'DL\\2',
    'e[Ff]((l|)uor|) {0,1}([0-9]+)' = 'eFluor \\3',
    'eVolve {0,1}([0-9]+)' = 'eVolve \\1',
    'Fluos' = 'Annexin-V-FLUOS',
    'FITC-AF488' = 'FITC',  # was previously 'FITC/AF488'
    'Indo {0,1}1' = 'Indo-1',
    'Maybe ' = '',
    '(Pac Blue|PB)' = 'Pacific Blue',
    '[Pp][Ee] {0,1}-{0,1}([A-QS-Za-qs-z])' = 'PE-\\1',
    '^PE ' = '^PE-',
    '^[Pp][Ee]$' = 'PE',
    '^PE-Dazzle$' = 'PE-Dazzle 594',
    'PerCP {0,1}-{0,1}([A-Za-z]+)' = 'PerCP-\\1',
    '(Zenon {0,1}|)(pHrodo) (iFL|) {0,1}([A-Za-z]+)' = '\\2 \\4',
    '(Ultra-LEAF |)([Pp]urified|[Pp]ure|[Uu]nlabeled)' = 'Purified',
    '(Q[D|d])(ot|) ([0-9]+)' = 'QD\\3',
    '^RPE$' = 'R-PE',
    '^RPM$' = 'R-PE',
    'red' = 'Red',
    'Tx{0,1}Re{0,1}d{0,1}' = 'Texas Red',
    'Vio([0-9]+)' = 'Vio \\1',
    'YG([0-9]+)' = 'YG \\1',

    # custom exact replacements
    '^eFluor 506 and eFluor 780 \\(APC-Cy7\\)$' = 'APC-Cy7',
    '^eFluor 660 \\(APC\\)$' = 'APC',
    '^PE-Vio {0,1}([0-9]+) \\(txRed\\)$' = 'PE-Vio \\1',

    # possible typos
    '^Fire 750$' = 'APC-Cy7',
    '^BUV427$' = 'BV421'
)

fixable_dyes <- c("DAPI", "Zombie UV", "Zombie Aqua")


# ----------------------------------------------------------------------
# Instrument Config

instr_cfg <- read_excel_or_csv(file.path(wd, opt[['instrument-config']]))
colnames(instr_cfg) <- unlist(lapply(colnames(instr_cfg), title_to_snake_case))  # format columns
instr_cfg <- rename_columns(instr_cfg, col_to_new_col)

# standardize names
instr_cfg[['fluorophore']] = multiple_replacement(
    instr_cfg[['fluorophore']], fluorophore_replacements, func='gsub'
)

# order list
instr_cfg <- instr_cfg[
    wrapr::orderv(instr_cfg[, c('excitation', 'bandpass_filter')], decreasing=TRUE), 
]

available_fluorophores <- separate_rows(instr_cfg, 'fluorophore', sep=', ')[['fluorophore']]

# save
if (!troubleshooting) {
    filepath = file.path(data_dir, '_instrument_config.csv')
    write.table(instr_cfg, file = filepath, row.names = FALSE, sep=',')

    filepath = file.path(data_dir, 'lists', 'available_fluorophores.txt')
    write.table(sort(unique(available_fluorophores)), filepath,
                row.names = FALSE, col.names = FALSE, quote = FALSE)
}

instr_cfg_expanded <- separate_rows(instr_cfg, 'fluorophore', sep=', ')


# ----------------------------------------------------------------------
# Antibody Inventory

ab_inv <- read_excel_or_csv(file.path(wd, opt[['antibody-inventory']]))
ab_inv <- ab_inv[, 1:(find_first_match_index('\\.{3}\\d{2}', colnames(ab_inv))-1)]  # filter extra columns
colnames(ab_inv) <- unlist(lapply(colnames(ab_inv), title_to_snake_case))  # column names
colnames(ab_inv) <- unlist(lapply(colnames(ab_inv), function(x) gsub('[.]', '', x)))  # column nmaes

# remove 'c2 a0', the "no-break space"
# see: https://stackoverflow.com/questions/68982813/r-string-encoding-best-practice
for (col in colnames(ab_inv)) {
    ab_inv[(!is.na(ab_inv[[col]]) & (ab_inv[[col]] == enc2utf8("\u00a0"))), col] <- NA
    ab_inv[[col]] = unlist(
        lapply(ab_inv[[col]], function(x) gsub(paste0('.*^', enc2utf8("\u00a0"), '+'), '', x))
    )    
}

# standardize names
for (col in c('antibody', 'alternative_name')) {
    ab_inv[[col]] = multiple_replacement(ab_inv[[col]], antibody_replacements, func='gsub')
}
ab_inv[['fluorophore']] = multiple_replacement(ab_inv[['fluorophore']], fluorophore_replacements, func='gsub')

# fill in missing fixable_dyes
mask <- (ab_inv[['antibody']] %in% fixable_dyes) & (is.na(ab_inv[['fluorophore']]))
ab_inv[mask, 'fluorophore'] <- ab_inv[mask, 'antibody']

unavailable_fluorophores = items_in_a_not_b(
    sort(unique( ab_inv[['fluorophore']] )),
    sort(unique( available_fluorophores ))
)

# save
if (!troubleshooting) {

    filepath = file.path(data_dir, '_antibody_inventory.csv')
    write.table(ab_inv, file = filepath, row.names = FALSE, sep=',')

    filepath = file.path(data_dir, 'lists', 'antibodies.txt')
    write.table(sort(unique(ab_inv[['antibody']])), filepath,
                row.names = FALSE, col.names = FALSE, quote = FALSE)

    filepath = file.path(data_dir, 'lists', 'alternative_names.txt')
    write.table(sort(unique(ab_inv[['alternative_name']])), filepath,
                row.names = FALSE, col.names = FALSE, quote = FALSE)

    filepath = file.path(data_dir, 'lists', 'all_fluorophores.txt')
    write.table(sort(unique(ab_inv[['fluorophore']])), filepath,
                row.names = FALSE, col.names = FALSE, quote = FALSE)

    filepath = file.path(data_dir, 'lists', 'unavailable_fluorophores.txt')
    write.table(sort(unique(unavailable_fluorophores)), filepath,
                row.names = FALSE, col.names = FALSE, quote = FALSE)
}


# ----------------------------------------------------------------------
# Build matrix


# Read data
panel <- read_csv_from_text(file.path(wd, opt[['input-file']]), encoding='UTF-8')

# merge channel_ids
panel = merge(
    ab_inv[(ab_inv[['antibody']] %in% panel[[1]]),  # selected antibodies only
           c('antibody', 'company', 'catalog_no', 'fluorophore')],
    instr_cfg_expanded[, c('fluorophore', 'channel_id')],  # channels
    by='fluorophore', suffixes=c('', '_'),
    all.x=TRUE, all.y=FALSE,
    na_matches = 'never'
)
panel <- fillna(panel, c('channel_id'), val='other')  # Biotin


# compute metadata for the panel
ab_counts <- panel %>%
    group_by(antibody, channel_id) %>%
    summarise(
        # num_antibodies=n(),  # don't necessarily need this
        num_fluorophores=n_distinct(fluorophore),
        fluorophores = toString(unique(fluorophore)),
        most_common_fluorophore=names(table(fluorophore))[which.max(table(fluorophore))],
        .groups='keep'
    )
instr_metadata <- ab_counts %>%
    group_by(channel_id) %>%
    summarise(
        fluorophores = toString(unique(most_common_fluorophore)),
        representative_fluorophore=names(
            table(most_common_fluorophore)[which.max(table(most_common_fluorophore))]
        )
    )
instr_metadata <- merge(
    instr_metadata,
    instr_cfg[, c('channel_id', 'laser', 'pmt')],
    by='channel_id', suffixes=c('', '_'),
    all.x=TRUE, all.y=FALSE,
    na_matches = 'never',
    sort=FALSE
)


# ----------------------------------------------------------------------
# Build matrix

# pivot into table where rows are antibodies and columns are non-overlapping colors
df <- pivot_wider(
    ab_counts,
    id_cols = 'antibody',
    names_from = 'channel_id',
    values_from = 'num_fluorophores',
    values_fill = 0
) %>% ungroup()


# Note: it is important that rowSums are computed before colSums

# Sort rows
num_channels_per_ab <- df %>%
    select_if(is.numeric) %>%
    mutate_at(
        colnames(select_if(df, is.numeric)),
        function(x) ifelse(x >= 1, 1, 0)
    ) %>% rowSums()
df[['num_channels_per_ab']] <- num_channels_per_ab
df <- df[order(df[['num_channels_per_ab']]), ]

# Sort columns
num_abs_per_channel <- df %>%
    select_if(is.numeric) %>%
    mutate_at(
        colnames(select_if(df, is.numeric)),
        function(x) ifelse(x >= 1, 1, 0)
    ) %>% colSums()
df <- append_dataframe(df, dataframe_row_from_named_list(num_abs_per_channel), reset_index=FALSE)
df <- df[ , order(unlist(df['count', ]), decreasing=TRUE)]  # sort columns

# Move id columns to the front
id_cols <- c('antibody', 'num_channels_per_ab', 'other')
channel_cols <- items_in_a_not_b(colnames(df), id_cols)
df <- df[ , c('antibody', 'num_channels_per_ab', channel_cols, 'other')]


# ----------------------------------------------------------------------
# Add additional data

# add comma-separated fluorophore list to each row
# for each antibody, aggregate fluorophores into comma-separated list
abs_vs_fluorophores <- ab_counts %>%
    group_by(antibody) %>%
    summarise(fluorophores = toString(unique(most_common_fluorophore)))
df = merge(
    df,
    abs_vs_fluorophores,
    by='antibody', suffixes=c('', '_'),
    all.x=TRUE, all.y=FALSE,
    na_matches = 'never',
    sort=FALSE
)

# add instrument metadata and available fluorophores to each column
df <- append_dataframe(
    df,
    stranspose(instr_metadata[
        items_in_a_not_b(colnames(instr_metadata), 'representative_fluorophore')
    ], 'channel_id'),
    infront=TRUE, reset_index=FALSE
)

# Optional: Rename column names from channel_id (1, 2, 3, etc.) to representative_fluorophore
fluorophores_for_channels = instr_metadata[['representative_fluorophore']]
names(fluorophores_for_channels) = instr_metadata[['channel_id']]
df <- rename_columns(df, fluorophores_for_channels[channel_cols])

# last row is the count
row.names(df)[nrow(df)] <- 'count'
df <- reset_index(df)  # export index information without the frameshift


# save
if (!troubleshooting) {
    filepath = file.path(dirname(data_dir), 'antibodies_vs_channels.csv')
    write.table(df, file = filepath, row.names = FALSE, sep=',')
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
