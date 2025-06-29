# Biome4all

<img src="www/hutch.png" alt="Biome4all Logo" width="150" height="150" />

**Biome4all** is a Shiny application designed to help visualize and analyze microbiome data.

---

## Features

* **Data Upload**

  * Upload a metadata CSV (first column: sample IDs)
  * Upload a taxon abundance CSV (first column: `tax_name`)

* **Alpha Diversity**

  * Compute Observed, Shannon, Simpson, and Evenness indices
  * Boxplots grouped by any metadata field
  * Downloadable alpha diversity table

* **Stacked Bar Plot**

  * Display relative abundances per sample
  * Facet by any metadata column (e.g., Treatment, Timepoint)
  * Dynamic Viridis palette based on number of top taxa

* **Reads per Sample**

  * Bar chart of total reads by sample

* **Heatmap**

  * Interactive heatmap of raw abundance matrix (Viridis colors)

* **PCoA (Bray–Curtis)**

  * Distance calculation on transposed abundance (samples × taxa)
  * 2D ordination with points, centroids, and 95% ellipses
  * Color-coded by any metadata field

---

## Installation

1. Clone or download this repository:

   ```bash
   git clone https://github.com/yourusername/biome4all.git
   cd biome4all
   ```

2. Ensure you have R (≥ 4.0) installed.

3. Install required R packages (run in R console):

   ```r
   install.packages(c(
     "shiny", "shinythemes", "ggplot2", "plotly", "dplyr",
     "tibble", "readr", "tidyr", "vegan", "pheatmap", "viridis", "ape"
   ))
   ```

4. Place the app’s logo (`hutch.png`) in the `www/` folder (it’s already referred to in the README).

---

## Usage

1. In the project directory, launch R or RStudio.

2. Run the app:

   ```r
   library(shiny)
   runApp("app.R")
   ```

3. In your browser, click **Upload files** and select:

   * **Metadata CSV** (first column = sample IDs, additional columns = metadata fields)
   * **Abundance CSV** (first column = `tax_name`, subsequent columns = counts per sample)

4. Navigate through tabs:

   * **Data**: Preview raw tables
   * **Alpha Diversity**: Choose `Group by` field and `Diversity index` → view boxplot → download table
   * **Stacked Bar**: Adjust number of taxa and facet by metadata → view interactive bar plot
   * **Reads per Sample**: View total reads bar chart
   * **Heatmap**: View raw abundance heatmap
   * **PCoA (Bray–Curtis)**: Choose `Color/Group by` field → view ordination plot

---

## File Structure

```
biome4all/
├── app.R
├── README.md
└── www/
    └── hutch.png       # Logo displayed in UI and README
```

* **app.R**: Main Shiny application script
* **README.md**: This documentation file
* **www/hutch.png**: Logo displayed at the top of the app and in README

---

## Contributing

Contributions, bug reports, and pull requests are welcome:

1. Fork this repository.
2. Create a new branch: `git checkout -b feature-name`.
3. Make your changes and commit: `git commit -m "Description of change"`.
4. Push to your fork: `git push origin feature-name`.
5. Open a Pull Request on GitHub.

Please ensure your code follows consistent style and includes minimal but clear comments.

---

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

[![DOI](https://zenodo.org/badge/994345371.svg)](https://doi.org/10.5281/zenodo.15570995)

---

## Contact / Assistance

If you need help, have questions, or want to report an issue, please contact:

> **Théo Ghelfenstein-Ferreira**
> Email: [theo.ferreira@aphp.fr](mailto:theo.ferreira@aphp.fr)

Feel free to reach out for assistance with installation, data formatting, or troubleshooting.

