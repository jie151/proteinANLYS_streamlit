library(venneuler)
library(eulerr)

my_data <- readLines("./venn_data.txt", n=2)

colnames <- unlist(strsplit(my_data[1],split = ',;'))
data <- unlist(strsplit(my_data[2],split = ',;'))
data <- as.numeric(data)

venn_data <- setNames(data, colnames)
print(venn_data)

vd <- euler(venn_data)
png("./file/image/plot1_10_2.png")
plot(vd,
     fills = list(fill = c("#fbb4ae", "#b3cde3", "#ccebc5"), alpha = 0.6),
     labels = list(col = "red", font = 8),
     edges = list(col = "black", lex = 2),
     quantities = TRUE)

dev.off()
