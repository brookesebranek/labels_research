peaks_and_window_table<-data.table::foverlaps(peaks_table, window_table, nomatch=NULL)
track_labels_and_window_table<-data.table::foverlaps(track_labels_table, window_table, nomatch = NULL)
coverage_and_window_table<-data.table::foverlaps(coverage_table, window_table, nomatch=NULL)
errors_and_window_table<-data.table::foverlaps(errors_table, window_table, nomatch=NULL)
