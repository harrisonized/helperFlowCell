# Regexes to standardize antibody and fluorophore names

## Objects
## col_to_new_col
## antibody_replacements
## fluorophore_replacements


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
    '7-AAD' = '7AAD',
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
    '[Ss][Cc][Aa](-|)\\d+' = 'Sca1',
    '.*^[Tt][Cc][Rr]' = 'TCR',
    '[Tt][Cc][Rr][Bb-]\\w*' = 'TCRb',
    '.*^[Tt][Dd][Tt]' = 'TdT',
    '[Tt][Gg][Ff][Bb-]\\w*' = 'TGFb',
    'Unlabel{1,2}ed' = 'Unlabeled',
    'Va(lpha|) {0,1}([0-9]+)' = 'Va\\2',
    'Vbeta' = 'Vb',
    'XCR-1' = 'XCR1',
    
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
    '7-AAD' = '7AAD',
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
    'e[Ff]((l|)uor|) {0,1}([0-9]+)' = 'eF\\3',
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
    '^eF(luor| |)506 and eF(luor| |)780 \\(APC-Cy7\\)$' = 'APC-Cy7',
    '^eF(luor| |)660 \\(APC\\)$' = 'APC',
    '^PE-Vio {0,1}([0-9]+) \\(txRed\\)$' = 'PE-Vio \\1',

    # possible typos
    '^Fire 750$' = 'APC-Cy7',
    '^BUV427$' = 'BV421'
)
