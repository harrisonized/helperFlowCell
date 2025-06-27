## helperFlowCell

This is a general purpose library to assist with handling flow cytometry data.

1. `build_design_matrix.R` takes a text input of markers (`panel.csv`) and compares it to reference files (`'antibody_inventory.xlsx'` and `instrument_config.xlsx'`) and outputs a succinct Excel spreadsheet in which each row is an antibody and each column is a distinct channel on your flow instrument. Positive cells are shaded gray for easy identification.
2. `plot_spectra.R` takes a text input of markers (`panel.csv`) and plots the fluorophores against channels available on your instrument (`instrument_config.xlsx'`) to show how well the fluorophores align with the detectors.
3. `reshape_compensation.R` takes the compensation exported from Flowjo (eg. `Acquisition-defined.csv`) and generates `compensation_matrix.txt` so that it can be imported into FACSDiva. It also allows you to rearrange the channels by providing a list in `fluorochrome_order.txt'`.
4. `analyze_populations.R` is the main analysis script. It takes raw exported counts from Flowjo's table editor (`Table.csv`) and user-defined `metadata.csv`, merges and reshapes them into a single table, then for each cell type in each organ, calculates two-sided unpaired T test against all groups. If provided, `mice.csv` can also be merged in, enabling the calculation of a `weeks_old` parameter.
5. `quantify_mfi.R` and `analyze_fluorescence_gates.R` are variations of `analyze_populations` that are still under development.


## Installation

Install the following packages in R:

```R
install.packages('import')
install.packages("optparse")
install.packages("logr")
install.packages("zeallot")  # %<-% operator
install.packages('magrittr')  # %>% operator
install.packages("progress")

# file IO
install.packages("openxlsx")
install.packages("readxl")
install.packages("XML")
install.packages("jsonlite")

# data
install.packages('plyr')
install.packages('dplyr')
install.packages('reshape2')
install.packages('tidyr')
install.packages('wrapr')
install.packages("stringr")
install.packages("stringi")

# plotting
install.packages('ggplot2')
install.packages('superb')
install.packages('cowplot')
install.packages('plotly')
install.packages('htmlwidgets')
```

## Data Requirements

Set up the following files and place them in the `ref` directory.

1. `antibody_inventory.xlsx`: This file should have at minimum the antibody and fluorophore columns. You should probably also include other identifying information, such as the company and catalog_no, but for the purposes of this repo, that's optional.

2. `instrument_config.csv`: This file should have the following columns: id, pmt, laser, fluorochromes, where each row is a separate laser/detector combination and the fluorochromes column should contain a comma-separated list of all the fluorophores that laser/detector combo should be able to detect. For the purposes of this repo, it's okay to leave fluorophores out, but do not duplicate a fluorophore, even if it is detected across multiple channels.

3. `spectra.csv`: Download the the fluorescence intensity data from [FPBase](https://www.fpbase.org/spectra/) by selecting your desired fluorophores, then clicking the "Share" button on the bottom right of the lower toolbar. An example file was included, since this is publicly available data.

 
## Getting Started
 
1. All scripts are meant to be run from the command line. For example:

    ```bash
    Rscript R/design_flow_panel.R
    ```

## Contributing

Please create a Github issue or reach out to me directly.

## Copyright

This code is copyright by Harrison Wang in 2023.
