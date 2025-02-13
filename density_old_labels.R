library(data.table)
labelcounts <- fread ("wc -l Downloads/data/*/samples/*/*/problems/*/labels.bed",col.names = c("number_of_labels" ,"file"))
labelcounts[order(number_of_labels)]
int.pattern <- list("[0-9]*", as.integer)
start.end.pattern <- list(":", start=int.pattern, "-", end=int.pattern)
label_counts_with_positions = nc::capture_first_df(labelcounts, file=start.end.pattern, nomatch.error = FALSE)
label_counts_with_positions[, bases := end-start]
 str(label_counts_with_positions)
 label_counts_with_positions[, density := number_of_labels/bases]
 str(label_counts_with_positions)
 label_counts_with_positions[order(density)]
 label_counts_with_positions[order(number_of_labels)]
 
