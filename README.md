## helperFlowCell

This tool helps with designing flow panels if you have a large collection of antibodies. Currently, there are two scripts:

1. `design_flow_panel.R` takes your `panel.csv` along with the reference files and outputs a succinct Excel spreadsheet in which each row is an antibody and each column is a distinct channel on your flow instrument. Positive cells are shaded gray for easy identification.
2. `plot_spectra.R` automatically plots the fluorophores listed in your `panel.csv` file.


## Installation

Install the following packages in R:

```R
install.packages('dplyr')
install.packages('wrapr')
install.packages("readxl")
install.packages("openxlsx")
install.packages('cowplot')
install.packages("optparse")
```

## Data Requirements

Set up the following files and place them in the `ref` directory.

1. `antibody_inventory.xlsx`: This file should have at minimum the antibody and fluorophore columns. You should probably also include other identifying information, such as the company and catalog_no, but for the purposes of this repo, that's optional.

2. `instrument_config.csv`: This file should have the following columns: id, pmt, laser, fluorochromes, where each row is a separate laser/detector combination and the fluorochromes column should contain a comma-separated list of all the fluorophores that laser/detector combo should be able to detect. For the purposes of this repo, it's okay to leave fluorophores out, but do not duplicate a fluorophore, even if it is detected across multiple channels.

3. `spectra.csv`: Download the the fluorescence intensity data from [FPBase](https://www.fpbase.org/spectra/) by selecting your desired fluorophores, then clicking the "Share" button on the bottom right of the lower toolbar. An example file was included, since this is publicly available data.

In the future, if I can find publicly available example data for antibody\_inventory.csv and instrument\_config.csv, I will include those here.

 
## Getting Started
 
1. Create a file named `panel.csv` that contains two columns: `antibody` and a `fluorophore`. The first column is important for `design_flow_panel.R`, and you can leave the second column blank at first. Run the panel. 

	```bash
	Rscript R/design_flow_panel.R
	```
	
2. After you have decided on which fluorophore combinations to use, enter that information into the `fluorophore` column of `panel.csv`, then run `plot_spectra.R`. This script is based on this [tutorial](https://bradyajohnston.github.io/posts/2022-09-03-plotting-fluorescence/) by Brady Johnston.

	```bash
	Rscript R/plot_spectra.R
	```

## Contributing

Please create a Github issue or reach out to me directly.

## Copyright

This code is copyright by Harrison Wang in 2023.