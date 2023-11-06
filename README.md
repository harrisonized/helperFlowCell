## helperFlowCell

This tool helpss with designing flow panels if you have a large collection of antibodies. It outputs a succinct table in which each row is an antibody and each column is a distinct channel on your flow instrument.

## Installation

Install the following packages in R:

```R
install.packages('dplyr')
install.packages('wrapr')
install.packages("readxl")
install.packages("openxlsx")
install.packages("optparse")
```

## Data Requirements

Set up the following two files and place it in the `ref` directory.

1. **antibody\_inventory.xlsx**: This file should have at minimum the antibody and fluorophore columns. You should probably also include other identifying information, such as the company and catalog_no, but for the purposes of this repo, that's optional.

2. **instrument\_config.csv**: This file should have the following columns: id, pmt, laser, fluorochromes, where each row is a separate laser/detector combination and the fluorochromes column should contain a comma-separated list of all the fluorophores that laser/detector combo should be able to detect. For the purposes of this repo, it's okay to leave fluorophores out, but do not duplicate a fluorophore, even if it is detected across multiple channels.

In the future, I will provide example files.
 
## Getting Started
 
1. In a text file named panel.txt, add a list of antibodies you'd like to build your panel for. Each line/row should be a unique antibody (eg. CD3 one one line, CD4 on the next, etc.)

2. Run the panel. You can use various options, like -i, -a, and -c to point to filenames different from the default option.

	```bash
	Rscript R/design_flow_panel.R
	```

## Contributing

Please create a Github issue or reach out to me directly.

## Copyright

This code is copyright by Harrison Wang in 2023.