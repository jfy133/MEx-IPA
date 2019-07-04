app <- ShinyDriver$new("../")
app$snapshotInit("LoadData_SingleSampleDropdown_ChangeSampleDropdown_ChangeTaxonDropDown_ChangeFilterDropdown")

app$setInputs(select_dir = "/home/fellows/Documents/Scripts/shiny_web_apps/MEx-IPA/dev/test_data/output_dairymicrobes_archive")
app$setInputs(submit = "click")
app$snapshot()
app$snapshot()
app$setInputs(selected_file = "CTW002.A0101.SG1.1_S0_L001_R1_001.fastq.combined.fq.prefixed.extractunmapped.bam.rma6")
# Input 'oute0a0bc6ec9d4a5c4_rows_current' was set, but doesn't have an input binding.
# Input 'oute0a0bc6ec9d4a5c4_rows_all' was set, but doesn't have an input binding.
app$snapshot()
app$setInputs(selected_file = "CTW001.A0101.SG1.1_S0_L001_R1_001.fastq.combined.fq.prefixed.extractunmapped.bam.rma6")
# Input 'out9648179001b0268f_rows_current' was set, but doesn't have an input binding.
# Input 'out9648179001b0268f_rows_all' was set, but doesn't have an input binding.
app$snapshot()
app$setInputs(selected_node = "Enterococcus_casseliflavus")
# Input 'out36f0de81c3090c58_rows_current' was set, but doesn't have an input binding.
# Input 'out36f0de81c3090c58_rows_all' was set, but doesn't have an input binding.
app$snapshot()
app$setInputs(selected_node = "Enterococcus_casseliflavus_EC20")
app$snapshot()
app$setInputs(selected_node = "Debaryomyces_hansenii")
# Input 'out16fae399487179a1_rows_current' was set, but doesn't have an input binding.
# Input 'out16fae399487179a1_rows_all' was set, but doesn't have an input binding.
app$setInputs(`selected_node-selectized` = "rumin")
app$setInputs(selected_node = "Lactobacillus_ruminis_ATCC_27782")
# Input 'out254d33ff112daa22_rows_current' was set, but doesn't have an input binding.
# Input 'out254d33ff112daa22_rows_all' was set, but doesn't have an input binding.
app$snapshot()
app$setInputs(selected_filter = "default")
# Input 'out1e5c3c40ec0187b6_rows_current' was set, but doesn't have an input binding.
# Input 'out1e5c3c40ec0187b6_rows_all' was set, but doesn't have an input binding.
app$snapshot()
app$setInputs(selected_filter = "ancient")
app$snapshot()
