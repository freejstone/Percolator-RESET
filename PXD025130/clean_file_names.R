library(tidyverse)

narrow_file = read_delim('conga.narrow.filtered.txt')
open_file = read_delim('conga.open.filtered.txt')

narrow_file = narrow_file %>% select(-drop_scan, -database)
open_file = open_file %>% select(-drop_scan, -database)
names(narrow_file) = gsub('_', ' ', names(narrow_file))
names(open_file) = gsub('_', ' ', names(open_file))

write_delim(narrow_file, 'narrow.filtered.txt', delim = '\t', na = '', quote = 'none')

write_delim(open_file, 'open.filtered.txt', delim = '\t', na = '', quote = 'none')
