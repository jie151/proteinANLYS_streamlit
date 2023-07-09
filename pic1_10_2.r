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
     fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
     alpha = 0.45,
     labels = list(col = "red", font = 8),
     edges = list(col = "black", lex = 2),
     quantities = TRUE,
     legend = TRUE)

dev.off()