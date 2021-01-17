# load the data ---
source(file = "load_data/load_data_final.R")

# count ---
med_count <- apply(count_mRNA, 2, median)
IQR_count <- apply(count_mRNA, 2, IQR)

data_count <- data.frame(med = med_count, IQR = IQR_count)

ggplot_count <- ggplot(data = data_count, aes(x = med, y = IQR)) +
  geom_point(color = "#3c3c3c") + xlab("Median") + ylab("IQR") + 
  geom_smooth(method = "loess") + ggtitle("Count data") +
  theme_Publication() + 
  theme(axis.text.x = element_text(size = 12))
# ggplot_count

# log2-CPM ---
logCPM_data <- log.cpm(count_mRNA)

med_logCPM <- apply(logCPM_data, 2, median)
IQR_logCPM <- apply(logCPM_data, 2, IQR)

plot(med_logCPM, IQR_logCPM)
data_logCPM <- data.frame(med = med_logCPM, IQR = IQR_logCPM)

ggplot_logCPM <- ggplot(data = data_logCPM, aes(x = med, y = IQR)) +
  geom_point(color = "#3c3c3c") + xlab("Median") + ylab("IQR") + 
  geom_smooth(method = "loess") + ggtitle("Log2-CPM data") +
  theme_Publication()+ 
  theme(axis.text.x = element_text(size = 12))

# ggplot_logCPM

# CPM ---
med_CPM <- apply(2^logCPM_data, 2, median)
IQR_CPM <- apply(2^logCPM_data, 2, IQR)

plot(med_CPM, IQR_CPM)
data_CPM <- data.frame(med = med_CPM, IQR = IQR_CPM)

ggplot_CPM <- ggplot(data = data_CPM, aes(x = med, y = IQR)) +
  geom_point(color = "#3c3c3c") + xlab("Median") + ylab("IQR") + 
  geom_smooth(method = "loess") + ggtitle("CPM data") +
  theme_Publication()+ 
  theme(axis.text.x = element_text(size = 12))

# ggplot_CPM

# VST ---
vst_data <- t(vst(t(round(count_mRNA))))

med_vst <- apply(vst_data, 2, median)
IQR_vst <- apply(vst_data, 2, IQR)

data_vst <- data.frame(med = med_vst, IQR = IQR_vst)

ggplot_vst <- ggplot(data = data_vst, aes(x = med, y = IQR)) +
  geom_point(color = "#3c3c3c") + xlab("Median") + ylab("IQR") + 
  geom_smooth(method = "loess") + ggtitle("VST data") +
  theme_Publication()+ 
  theme(axis.text.x = element_text(size = 12))

# ggplot_vst


# arrange the plot ---
final_plot <- ggarrange(ggplot_count, ggplot_CPM, ggplot_logCPM, ggplot_vst, 
          labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)
# final_plot

ggsave(final_plot, filename = "pdf/mean_variance.pdf")
print("Figure saved in 'pdf/mean_variance.pdf' file")
