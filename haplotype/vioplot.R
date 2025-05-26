setwd("./")
install.packages("ggpubr")
library(ggpubr)
data <- read.table("data.txt", sep = "\t" ,header = T)


#PH
p1 <- ggviolin(data, x = "Hap", y = "PH", fill = "Hap",ylab = "Plant height (cm)",
         add = "boxplot",add.params = list(fill="white")) + 
  stat_compare_means()+
  theme_bw()
