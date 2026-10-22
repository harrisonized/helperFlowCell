## helperFlowCell

The purpose of this package is to help with analyses being done by FlowJo 10.

**Upstream** of your flow cytometry experiment, `build_design_matrix.R` produces a succinct Excel spreadsheet to assist with selecting your panel based on your instrument configuration and antibody inventory. `plot_spectra.R` plots the fluorophores against channels available on your instrument separated by laser. Finally, `reshape_compensation.R` allows you to copy the perfect compensation from the previous experiment back into the FACSDiva for your next experiment.

**Downstream** of your analysis, `analyze_counts.R` takes exported counts data along with some metadata you define and generates presentation-ready box plots and ribbon diagrams along with statistics. MFI data can be analzyed using `analyze_mfi.R` and `analyze_tvd.R`.


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
install.packages('pandoc')  # save static plotly images
install.packages('ggprism')
install.packages('svglite')
update.packages('systemfonts')
```

There may be extra libraries I did not account for here. If so, simply install them. In python, create a conda environment ('r-reticulate') and install plotly and kaleido:

```bash
conda create -n r-reticulate python=3.9
conda activate r-reticulate
pip install plotly==5.3.1
pip3 install -U kaleido==0.2.1  # important! higher versions result in some plots being blank at export
pip install numpy==1.24.2
```

If you have a Linux machine, additionally you need to install pandoc via your terminal:

```bash
sudo apt install pandoc
```

To help R find your conda environment, create the `.Renviron` hidden file on in your home directory and include the following three lines:

```bash
LD_LIBRARY_PATH="/path/to/miniconda3/envs/r-reticulate/lib/libpython3.9.dylib"
RETICULATE_MINICONDA_PATH="/path/to/miniconda3"
RETICULATE_PYTHON_ENV="/path/to/miniconda3/envs/r-reticulate/"
```

Note that if you have a higher version of plotly and kaleido in your conda environment, you will see the following entirely unhelpful error. See this [Stack Overflow link](https://stackoverflow.com/questions/71118601/saving-a-plotly-image-not-working-with-kaleido-even-though-it-is-installed) for more details.

```R
Error in save_image(): ! {reticulate} wasn't able to find a Python environment. ℹ If you have an existing Python installation, use reticulate::use_python() to inform {reticulate} of it. ℹ To have {reticulate} install Python for you, reticulate::install_python().
Backtrace:
    ▆ 1. └─save_fig(...)
      2. ├─base::suppressWarnings(...)
      3. │ └─base::withCallingHandlers(...)=
      4. └─plotly::save_image(...)
      5. └─plotly::kaleido()
      6. └─rlang::abort(...)
      Error in save_image(fig, file = file.path(file.path(dirpath, "png"), paste0(filename, : Calls: save_fig
```

Lastly, put this in your bashrc:

```bash
# mac
RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools"

# linux
RSTUDIO_PANDOC="/usr/lib/rstudio/resources/app/bin/quarto/bin/tools/x86_64"
```

## Reference Files

To use `build_design_matrix.R` and `plot_spectra.R`, set up the following files and place them in the `ref` directory.

1. `antibody_inventory.xlsx`: This file should have at minimum the antibody and fluorophore columns. You should probably also include other identifying information, such as the company and catalog_no, but for the purposes of this repo, that's optional.

2. `instrument_config.csv`: This file should have the following columns: id, pmt, laser, fluorochromes, where each row is a separate laser/detector combination and the fluorochromes column should contain a comma-separated list of all the fluorophores that laser/detector combo should be able to detect. For the purposes of this repo, it's okay to leave fluorophores out, but do not duplicate a fluorophore, even if it is detected across multiple channels.

3. `spectra.csv`: An example file was included, since this is publicly available data. If there is a fluorophore you did not find, download it from [FPBase](https://www.fpbase.org/spectra/) by selecting your desired fluorophores, then clicking the "Share" button on the bottom right of the lower toolbar, then copy it into the `spectra.csv` file.


## Getting Started

### Set Up Files for Analysis

To get the counts, open Flowjo's Table Editor and drag in all the gates you want for the organ that you want. Change the metric to Count instead of "Freq. of Parent". The fastest way to do this is to open the workspace file using Sublime or other text editor, then go to the TableEditor section. Under your table of interest, replace `statistic="fj.stat.freqofparent"` with `statistic="Count"` for all rows.

To get the mfi and sdev, do the same thing as above in Sublime, but instead, use regex to search and replace this block: `name=.*?rename=""`. For mfi, replace with `name="Geometric Mean : GFP"  identifier="Geometric Mean"  statistic="Geometric Mean"  parameterName="Comp-GFP-A :: GFP"  rename=""`. For sdev, replace with `name="Robust SD : GFP"  identifier="Robust SD"  statistic="Robust SD"  parameterName="Comp-GFP-A :: GFP"  rename=""`. Make sure you name your file so you know what kind of mean you use. For example, I name my mfi files `gmfi.csv` and `rsdev.csv` to denote geometric mean and robust sd, respectively.

Finally, set up the metadata file with the following column names:

```
fcs_name
mouse_id
sex	
treatment
zygosity
organ
is_excluded
```

Additionlly, if you have raw counts, you can include the following columns to get statistics on absolute counts:

```
total_cells_per_organ
total_viable_cells
viability
```
Structure your directory as follows:

```R
helperFlowCell
└─ data
    ├─ flow-counts
    │   └─ exp1-counts.csv  # required
    ├─ flow-mfi
    │   └─ exp1-mfi.csv   # only for analyze_mfi and analyze_tvd
    ├─ flow-sdev
    │   └─ exp1-sdev.csv  # only for analyze_tvd
    ├─ flow-metadata
    │   └─ exp1-metadata.csv  # required
    └─ mouse
        └─ exp1-mouse-db.csv  # optional
```

Note that in each directory, you can place multiple csv files and the scripts will append all the csvs in the directory.


### How to Run

After setting up the files to be analyzed, run your desired script from the command line. For example:

```bash
cd path/to/helperFlowCell
Rscript R/analyze_counts.R
```

## Contributing

Please create a Github issue or reach out to me directly.

## Copyright

This code is copyright by Harrison Wang in 2023.
