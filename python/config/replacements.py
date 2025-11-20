instr_cfg_colreps = {
    'id': 'channel_id',
    'fluorochrome_detected': 'fluorophore',
    'fluorochrome': 'fluorophore'
}

ab_inv_colreps = {
    'fluorophor': 'fluorophore',
    'catalog_nr': 'catalog_no'
}

antibody_replacements = {
    # special characters
    'α': 'a',
    'β': 'b',
    'γ': 'g',
    '™': '',
    r'([A-Za-z0-9],)\s*([A-Za-z0-9])': r'\1 \2',

    # animals
    r'[Gg][Oo][Aa][Tt]': 'goat',
    r'[Mm][Oo][Uu][Ss][Ee]': 'mouse',
    r'[Rr][Aa][Tt]': 'rat',
    r'[Dd][Oo][Nn][Kk][Ee][Yy]': 'donkey',
    r'[Bb][Oo][Vv][Ii][Nn][Ee]': 'bovine',
    r'[Ss][Hh][Ee][Ee][Pp]': 'sheep',

    # general replacements
    '7-AAD': '7AAD',
    r'[Bb][Cc][Ll1]-{0,1}': 'Bcl-',
    r'.*^[Cc][Dd]': 'CD',
    'FoxP3': 'Foxp3',
    'gammaDelta': 'gd',
    'INFg': 'IFNg',
    'Immunoglobulin ': 'Ig',
    r'.*^[Ii][Ll]': 'IL',
    r'(IL)([0-9]+)': r'\1-\2',
    r'.*^Ly-': 'Ly',
    r'MHC {0,1}(I{1,2})': r'MHC\1',
    r'.*^[Nn][Oo]tch': 'Notch',
    r'.*^[Oo]nly': 'Only',
    r'[Rr][Oo][Rr][gγy][yt]': 'RORgt',
    r'[Ss][Cc][Aa](-|)\d+': 'Sca1',
    r'.*^[Tt][Cc][Rr]': 'TCR',
    r'[Tt][Cc][Rr][Bb-]\w*': 'TCRb',
    r'.*^[Tt][Dd][Tt]': 'TdT',
    r'[Tt][Gg][Ff][Bb-]\w*': 'TGFb',
    r'Unlabel{1,2}ed': 'Unlabeled',
    r'Va(lpha|) {0,1}([0-9]+)': r'Va\2',
    'Vbeta': 'Vb',
    'XCR-1': 'XCR1',

    # custom exact replacements
    r'^\(CXCR4\)$': 'CXCR4',
    ' Fixable Viability Kit': '',
    'NK cell Pan': 'CD49b',
    'TCRb TCRcb': 'TCRb, TCRcb',
    'Vb8.1 Vb8.2': 'Vb8.1, 2',
    r' {0,1}\(Tonegawa nomenclat\)': ''
}

fluorophore_replacements = {
    # special characters
    '®': '',
    '/': '-',
    r'^ ([A-Za-z]+)': r'\1',

    # general replacements
    '7-AAD': '7AAD',
    r'(Alexa) {0,1}(Fluor|) {0,1}': 'AF',
    r'A(F|) {0,1}([0-9]+)': r'AF\2',
    r'\(allophycocyanin\)( |)': '',
    r'[Aa][Pp][Cc]': 'APC',
    r'APC( |-|)([CcYyVioFfire0-9]+)': r'APC-\2',
    r'^APC-[Ff]ire( |-|)([0-9]+)': r'APC-Fire \2',
    r'^APC-[Ff]ire$': 'APC-Fire 750',
    r'^(A|a)(T|t)(T|t)(O|o)( |-|)([0-9]+)': r'ATTO \6',
    r'[Bb][Ii][Oo](tin|)(,|)( |)([A-Za-z]+|)': 'Biotin',
    r'BD Horizon( |)[Bb][Uu][Vv]': 'BUV',
    r'(BU[Vv]) {0,1}([0-9]+)': r'BUV\2',
    r'(B[Vv]) {0,1}([0-9]+)': r'BV\2',
    r'Brilliant Violet( |)': 'BV',
    r'([Dd][Ll]|Dy[Ll]ight) {0,1}-{0,1}([0-9]+)': r'DL\2',
    r'e[Ff]((l|)uor|) {0,1}([0-9]+)': r'eF\3',
    r'eVolve {0,1}([0-9]+)': r'eVolve \1',
    'Fluos': 'Annexin-V-FLUOS',
    'FITC-AF488': 'FITC',
    r'Fluorescein \(FITC\)': 'FITC',
    r'Indo {0,1}1': 'Indo-1',
    'Maybe ': '',
    r'(Pac Blue|PB)': 'Pacific Blue',
    r'^[Pp][Ee]$': 'PE',
    r'^PE ': 'PE-',
    r'[Pp][Ee] {0,1}-{0,1}([A-QS-Za-qs-z]+)': r'PE-\1',
    r'PE(-|/)Dazzle( |)([A-Za-z0-9 ]+|)': r'PE-Dazzle \3',
    r'^PE-Dazzle( |)$': 'PE-Dazzle 594',
    r'^PE-[Ff]ire( |)([0-9]+)$': r'PE-Fire \2',
    r'^RPE$': 'PE',
    r'^RPM$': 'PE',
    r'PE-\(R-PE, R-phycoerythrin\)': 'PE',
    r'PerCP {0,1}-{0,1}([CcYyeF0-9]+)': r'PerCP-\1',
    r'Super Bright( |)': 'SB',
    r'(Zenon {0,1}|)(pHrodo) (iFL|) {0,1}([A-Za-z]+)': r'\2 \4',
    r'(Ultra-LEAF |)([Pp]urified|[Pp]ure|[Uu]nlabeled)': 'Purified',
    r'(Q[D|d])(ot|) ([0-9]+)': r'QD\3',
    'red': 'Red',
    r'Spark NIR( +|)([0-9]+)': r'SNIR\2',
    r'Tx{0,1}Re{0,1}d{0,1}': 'Texas Red',
    r'Vio([0-9]+)': r'Vio \1',
    r'YG([0-9]+)': r'YG \1',
    r'^Zombie Aq$': 'Zombie Aqua',

    # custom exact replacements
    r'^eF(luor| |)506 and eF(luor| |)780 \(APC-Cy7\)$': 'APC-Cy7',
    r'^eF(luor| |)660 \(APC\)$': 'APC',
    r'^PE-Vio {0,1}([0-9]+) \(txRed\)$': r'PE-Vio \1',

    # possible typos
    r'^Fire 750$': 'APC-Cy7',
    r'^BUV427$': 'BV421'
}
