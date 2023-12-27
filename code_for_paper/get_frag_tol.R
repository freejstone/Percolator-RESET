#get fragment tolerance from search param files and update open_fragger
library(xfun)

ids = list.files()

for (id in ids) {
  get_param_file = readLines(paste(id, 'crux-output', 'narrow_50.tide-search.params.txt', sep = '/'))
  fragment_tol = get_param_file[grepl('mz-bin-width=', get_param_file)][1]
  fragment_tol = sub('.*=', '', fragment_tol)
  gsub_file(paste(id, 'closed_fragger.params', sep = '/'), 'fragment_mass_tolerance = 20', paste('fragment_mass_tolerance =', fragment_tol))
  gsub_file(paste(id, 'open_fragger.params', sep = '/'), 'fragment_mass_tolerance = 20', paste('fragment_mass_tolerance =', fragment_tol))
  
  gsub_file(paste(id, 'closed_fragger.params', sep = '/'), 'fragment_mass_units = 1', 'fragment_mass_units = 0')
  gsub_file(paste(id, 'open_fragger.params', sep = '/'), 'fragment_mass_units = 1', 'fragment_mass_units = 0')
  
}
