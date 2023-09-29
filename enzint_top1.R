datasets = list.files(path = 'datasets', pattern = 'PXD')
for (dataset in datasets) {
  print(dataset)
  for (i in 0:9) {
    file = paste('datasets/', dataset, '/crux-output/open_5', i, '.make-pin.pin', sep = '')
    df = read_delim(file)
    rank = substring(df$SpecId, nchar(df$SpecId), nchar(df$SpecId))
    df = df[rank == '1',]
    df$enzInt = 0
    write_delim(df, paste('datasets/', dataset, '/crux-output/open_1_', i, '.make-pin.pin', sep = ''), delim = '\t', quote = "none")
  }
}

datasets = list.files(path = 'datasets', pattern = 'PXD')
for (dataset in datasets) {
  print(dataset)
  for (i in 0:9) {
    file = paste('datasets/', dataset, '/crux-output/narrow_5', i, '.make-pin.pin', sep = '')
    df = read_delim(file)
    rank = substring(df$SpecId, nchar(df$SpecId), nchar(df$SpecId))
    df = df[rank == '1',]
    df$enzInt = 0
    write_delim(df, paste('datasets/', dataset, '/crux-output/narrow_1_', i, '.make-pin.pin', sep = ''), delim = '\t', quote = "none")
  }
}
