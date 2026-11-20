scripts <- sort(list.files(pattern = "^0[1-5]_.*[.]R$"))
stopifnot(length(scripts) == 5L)
for (script in scripts) {
  message("Running ", script)
  source(script, local = new.env(parent = globalenv()))
}
