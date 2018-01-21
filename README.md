# MALT-Extract_iViewer

Interactive viewer of MALT-Extract results. Currently allows you to see damage patterns, edit distance and median fragment length of a taxon compared all average lengths of every taxon in a sample.

Requires the R packages shiny, tidyverse, and gridExtra.

Input is a zipped MALT-Extract Beta v1.2 results directory. The largest file size is 50MB, but can be increased in the 'options' on line 10 of the `app.R` file. 

To run,

```
R -e "shiny::runApp('/<path>/<to>/MALT-Extract_iViewer')"
```

and then go to the IP address given once loaded in your internet browser.
