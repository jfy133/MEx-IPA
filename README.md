# MEx-IPA

**M**(alt)**Ex**(tract)-**I**(nteractive)**P**(lotting)A(pp)

Interactive viewer of MALT-Extract results. 

[MaltExtract](https://github.com/rhuebler/MaltExtract) is a part of the [HOPS pipeline](https://github.com/rhuebler/HOPS) ([Huebler et al. 2019 bioRxiv](https://doi.org/10.1101/534198)).
It extracts various metrics of from alignments as stored
in MEGAN6 ([Huson et al. 2016 PloS Comp. Bio.](https://doi.org/10.1371/journal.pcbi.1004957)) RMA6 files.
It is particularly designed for identifiying authentic ancient DNA of 
a given organism.

This R shiny app is a modified implementation of the metric visualisation 
script of the HOPS pipeline, with additional interactive plot functionality.

It will display the following distributions of the alignments of a given sample
and taxonomic node:
  * General statistics table of the node (e.g. mean alignents to a node)
  * C to T miscorporation(a.k.a damage plots)
  * Read length
  * Edit distance
  * Percent identity
  * Percentage of reference covered (breadth coverage)
  * Fold coverage (depth coverage)

It also provides two comparison panes - for a given characteristic, a) for a
given node, the plots for every sample with data on that node, and b) for a 
given sample, plots for every node. This allows for fast detection of potential
samples or taxa of interest from large datasets.  

![](Example display of MEx-IPA)

## Preparation
This shiny app has been tested on `R` version 3.6, and requires the following 
packages:

 * shinyWidgets (tested v0.4.8)
 * shinycustomloader (tested v0.9.0)
 * patchwork (tested v0.0.1)
 * DT (tested v0.6)
 * plotly (tested v4.9.0)
 * data.table (tested v1.12.2)
 * tidyverse (tested v1.2.1)
 * shiny (tested v1.3.2)

To install these packages in R:

```r
install.packages(c("shiny", "tidyverse", "data.table", "plotly", "DT", 
	"patchwork", "shinycustomloader", "shinyWidgets"))
```

To download the app, change to a directory where you wish to install the 
app, and run the following in your terminal

```bash
git clone https://github.com/jfy133/MEx-IPA.git
```

This will download all the app files from this github repository to a
directory named MEX-IPA.

## Running

Once installed there are two methods of running the app. 

1) Run the following command from your terminal

```bash
R -e "shiny::runApp('/<path>/<to>/MEx-IPA')"
```

![How to run in a Terminal](assets/images/02-rstudio_instructions.png)

and then go to the IP address given once loaded in your internet browser.

2 ) Alternatively, you can use Rstudio to load either the `server.R` or `ui.R` 
file, and press the 'Run App' button in the top right hand corner of the 
code pane.

![How to run in Rstudio](assets/images/03-rstudio_instructions.png)

> NB: For both methods, you must run on the app on an a machine that is able to 
> access your MaltExtract results directory via a directory path e.g. as listed 
> in from your file explorer. It does not take remote URLs.

3) Type or copy and paste the path to your MaltExtract results directory. Once
the RunSummary.txt file is found in this directory, the 'Run Visualisation' 
button will appear.

4) Customise the display based on the options on the left. Remove from
file name allows you to remove a given string (e.g. file suffix) from the
sample names. Interactive plots allow you to hover over lines or bars to get
specific alignment counts for that position or bin.

> NB: If you have a very large MaltExtract results run, plots may take a 
> long time to display - particularly on the multi-plot panels.

5) You can save a PDF report of the non-interactive single sample displays with 
the 'download PDF' report button. The PDF can be modified by loading into 
vector-based image software such as [inkscape](https://inkscape.org/). 

> NB: If loading into inkscape, load with 'internal import' so text is 'stored
> as text'. Text boxes can be made editable by selecting the box and going 
> Text > Remove Manual Kerns

6) You can download the alignment stats table in a range of tabular formats 
using the 'Download Table' button above the table.

## FAQs

### I just see a white screen where the plots should be!

This can be because it takes a long time to render the plots - particularly
when you are displaying many plots at once or if loading the data
from a remote server into your own computer. Be patient.

### What is the difference between the 'everything at 0' plots and the 'No Input data found' messages

This is something to do with the way MaltExtract saves results. The both 
essentially mean the same thing: there are no alignments to that specific node.

### Why is the interactive plots a bit ugly with titles and axis labels?
The `plotly` package is unfortunately not fully compatible with all `ggplot` 
functions. I tried different ways of hacking around it but failed. However,
this does not remove from the utility of being able to get fine resolution
statistics of each characteristic. All report downloads are 

## References

Huebler, R. et al. HOPS: Automated detection and authentication of pathogen DNA in archaeological remains. bioRxiv 534198 (2019). doi:10.1101/53419

Huson, D. H. et al. MEGAN Community Edition - Interactive Exploration and Analysis of Large-Scale Microbiome Sequencing Data. PLoS Comput. Biol. 12, e1004957 (2016).

## Change Log

**0.4**
  * New rewrite with dynamically generate plots and an option for interactive plots

<0.4 - Experimental alpha versions