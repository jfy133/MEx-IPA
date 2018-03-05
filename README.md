# MALT-Extract_iViewer

Interactive viewer of MALT-Extract results. Currently allows you to see damage patterns, edit distance and a read length distribution.

Requires the R packages shiny, tidyverse, and gridExtra.

Input is a zipped MALT-Extract Beta v1.2 results directory. The largest file size is 50MB, but can be increased in the 'options' on line 10 of the `app.R` file. 

To run, either:

Load the app in Rstudio and press 'Run' in the top right of the code pane

OR

run the following from your terminal:

```
R -e "shiny::runApp('/<path>/<to>/MALT-Extract_iViewer')"
```

and then go to the IP address given once loaded in your internet browser.
