# GA

## Memory usage

`simulate` formerly stored the full diagnostic structure in `resp.diag`, which could be large. It now retains only the `cav_frac_t` time series, reducing memory needs of the response object. For a representative run with 1000 time steps and 3 stories, this change lowers memory from about 168 kB to 8 kB (~95% reduction).

