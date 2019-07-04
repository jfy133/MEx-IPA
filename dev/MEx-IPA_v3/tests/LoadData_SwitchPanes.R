app <- ShinyDriver$new("../")
app$snapshotInit("LoadData_SwitchPanes")

app$setInputs(select_dir = "/home/fellows/Documents/Scripts/shiny_web_apps/MEx-IPA/dev/test_data/output_dairymicrobes_archive")
app$snapshot()
app$setInputs(submit = "click")
app$snapshot()
app$snapshot()
app$snapshot()
