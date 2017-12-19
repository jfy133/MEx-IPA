# MALT-Extract_iViewer

Interactive viewer of MALT-Extract results. Currently allows you to see damage patterns, edit distance and median fragment length of a taxon compared all average lengths of every taxon in a sample.

Requires the R packages tidyverse and gridExtra.

To run on your own data, change the directory in the script (where indicated) before running in either Rstudio or 

```
R -e "shiny::runApp('/<path>/<to>/MALT-Extract_iViewer')"
```

and then go to the IP address given once loaded in your internet browers.
