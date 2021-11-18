#creating individual ggplots first for simple visualization
library(ggplot2)
#ggplot for the coverage table
ggplot2::ggplot()+geom_path(aes(x=chromStart, y=matches), data=coverage_and_window_table)
#ggplot for the peak detection
ggplot2::ggplot()+geom_segment(aes(x=chromStart, y=mean, xend=chromEnd, yend=mean), color="orange", size=3, data= peaks_and_window_table)
#ggplot for the label detection
ggplot2::ggplot()+geom_rect(aes(xmin=chromStart, ymin=-Inf, xmax=chromEnd, ymax=Inf, fill=annotation), alpha= 0.25, data=track_labels_and_window_table)
#ggplot for the errors
ggplot2::ggplot()+geom_rect(aes(xmin=chromStart, ymin=-Inf, xmax=chromEnd, ymax=Inf,  linetype=status), fill="transparent", color="black", size=1, data=errors_and_window_table)+
  scale_linetype_manual(
    "error type",
    values=c(
      correct=0,
      "false negative"=3,
      "false positive"=1))
#to add them all to be one concise image, we use + 
ggplot2::ggplot()+geom_path(
  aes(x=chromStart, 
      y=matches), 
    color="gray50", data=coverage_and_window_table)+
  geom_segment(
    aes(x=chromStart, 
        y=mean, 
        xend=chromEnd, 
        yend=mean), 
    color="orange", size=3, data= peaks_and_window_table)+
  geom_rect(aes(
    xmin=chromStart, 
    ymin=-Inf, 
    xmax=chromEnd, 
    ymax=Inf, 
    fill=annotation), alpha= 0.25, data=track_labels_and_window_table)+
  geom_rect(aes(
    xmin=chromStart, 
    ymin=-Inf, 
    xmax=chromEnd, 
    ymax=Inf,  
    linetype=status), fill="transparent", color="black", size=1, data=errors_and_window_table)+
  scale_linetype_manual(
    "error type",
    values=c(
      correct=0,
      "false negative"=3,
      "false positive"=1))

