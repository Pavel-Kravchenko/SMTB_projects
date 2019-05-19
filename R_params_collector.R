library(readr)
library(data.table)
library(magrittr)
library(reshape2)
library(gplots)


#args <- commandArgs(trailingOnly = T)

#first_dir <- toString(args[1]) # first dir
dir_name <- "/home/pavel/Desktop/Work/Project/Simulation/result_big/"
print(dir_name)

#output_dir <- toString(args[3])
output_dir <- "/home/pavel/Desktop/Work/Project/Simulation/result_big/output"

"%+%" <- function(...){
  paste0(...)
}



if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Dir already exists!")
}

dir_name_files <- dir(dir_name)
print("dir_name_files")
print(dir_name_files)
dir_name_files <- grep(".+.csv", dir_name_files, value=TRUE)
dir_name_files


Total <- data.frame(generation=NaN, Ne=NaN, recombination=NaN, mutation=NaN, Kendall=NaN, len_before_5=NaN, len_before_10=NaN, len_before_15=NaN, len_before_20=NaN)
Total <- na.omit(Total)
Total



detach_package <- function(pkg, character.only = FALSE)
{
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}

detach_package(dplyr)
detach(package:dplyr)
#detach("package:dplyr", unload=TRUE)
#name_list = list()
for (file_name in dir_name_files){
  print(file_name)
  
  
  LD_txt <- read_csv(dir_name %+% file_name)  
  if (nrow(LD_txt) == 0){
    next
  }
  
  LD_txt <- LD_txt[order(LD_txt$Len, decreasing = F),] 
  x <- LD_txt$Len
  y <- LD_txt$D
  
  
  a <- LD_txt$Len
  b <- LD_txt$`R^2`
  
  #wind <- as.integer(round(length(x)/50, 0))
  wind <- length(x)*0.1
  
  if (wind == 0) {wind <- 2}
  
  par(mfrow=c(2,2))
  plot(x, y, col=grey(.7),  main = "D plot", xlab ="Distance", ylab = "D", ylim=c(-0.25, 0.25))
  
  f <- rep(1/wind, wind)
  y_lag <- filter(y, f, sides=1)
  lines(x, y_lag, col="red")
  
  plot(a, b, col=grey(.7),  main = "R^2 plot", xlab ="Distance", ylab = "R^2", ylim=c(0, 1))
  
  q <- rep(1/wind, wind)
  b_lag <- filter(b, q, sides=1)
  lines(a, b_lag, col="blue")
  
  hist(y, breaks = "Sturges", main = "D histogram", xlab ="LD", ylab = "Frequency", pch=16, col = "orange")
  hist(x, breaks = "Sturges", main = "Len histogram", xlab ="Len", ylab = "Frequency", pch=16, col = "green") 
  
  # #-------------------------------------------------------------------------
  # if (length(df1$x) > 50000) { df1 <- df1[1:50000,]}
  # #if (length(df1$x) > 50000) {df1 <- df1[sample(nrow(df1), 50000), ]}
  
  #cor_test$estimate[[1]]
  cor_test = tryCatch({
    cor.test(a, b, method = "kendall")
  }, error = function(e) {
    return(NaN)
  })
  
  cor_test_estimate = tryCatch({
    cor_test$estimate[[1]]
  }, error = function(e) {
    return(NaN)
  })
  

  file_name_split <- strsplit(file_name, split='_', fixed=TRUE)
  print(file_name_split)
  
  calculate_first <- function(a, b_lag, y_value){
    tmp <- which((b_lag <= y_value) == TRUE)
    if (length(tmp) == 0){
      tmp <- a[length(a)]
      return(NaN)
    }
    else{
      tmp <- a[min(tmp)]
      return(tmp)
    }
    
  }
  

  curr <- data.frame(generation=file_name_split[[1]][3], Ne=file_name_split[[1]][4], recombination=file_name_split[[1]][5], mutation=file_name_split[[1]][6], Kendall=cor_test_estimate, len_before_5=calculate_first(a, b_lag,0.05), len_before_10=calculate_first(a, b_lag,0.1), len_before_15=calculate_first(a, b_lag,0.15), len_before_20=calculate_first(a, b_lag,0.2))
  print(curr)
  Total <- rbind(Total, curr)
  
}

#Total <- na.omit(Total)
#Total



Ne_val = 4000
mutation_val = 1
library(dplyr)
data <- Total %>% filter(Ne == Ne_val, mutation == mutation_val) %>% 
  select(c(generation, Kendall, recombination))
length(unique(data["generation"][[1]]))
length(unique(data["recombination"][[1]]))

#################################################
m = matrix(c(0)*length(unique(data["generation"][[1]]))*length(unique(data["recombination"][[1]])), 
           nrow=length(unique(data["generation"][[1]])), ncol=length(unique(data["recombination"][[1]])))
m
colnames(m) <- as.character(sort(as.numeric(as.character(unique(data["recombination"][[1]])))))
rownames(m) <- as.character(sort(as.numeric(as.character(unique(data["generation"][[1]])))))

m

for (i in 1:nrow(data)){
  #print(i)
  m[as.character(as.numeric(as.character(data$generation[[i]]))), as.character(as.numeric(as.character(data$recombination[[i]])))] = data$Kendall[[i]]
}
m

heatmap.2(m,
          col = colorpanel(100,"red","yellow","green"),
          margins = c(5, 10),
          trace = "none", 
          lhei = c(2, 8),
          scale = c("none"),
          symbreaks = min(m, na.rm=TRUE),
          na.color="blue",
          dendrogram = "none", 
          Colv = FALSE,
          cellnote=round(m,2),
          notecex=1.0,
          Rowv = F,
          notecol="black",
          main = "R^2 heatmap Kendall. Ne: " %+% Ne_val, 
          xlab = '% of recombination', 
          ylab = 'generations')
################################################
data <- Total %>% filter(Ne == Ne_val, mutation == mutation_val) %>% 
  select(c(generation, len_before_5, recombination))
length(unique(data["generation"][[1]]))
length(unique(data["recombination"][[1]]))

m = matrix(c(0)*length(unique(data["generation"][[1]]))*length(unique(data["recombination"][[1]])), 
           nrow=length(unique(data["generation"][[1]])), ncol=length(unique(data["recombination"][[1]])))
m
colnames(m) <- as.character(sort(as.numeric(as.character(unique(data["recombination"][[1]])))))
rownames(m) <- as.character(sort(as.numeric(as.character(unique(data["generation"][[1]])))))

m

for (i in 1:nrow(data)){
  #print(i)
  m[as.character(as.numeric(as.character(data$generation[[i]]))), as.character(as.numeric(as.character(data$recombination[[i]])))] = data$len_before_5[[i]]
}
m




heatmap.2(m,
          col = colorpanel(100,"red","yellow","green"),
          margins = c(5, 10),
          trace = "none", 
          lhei = c(2, 8),
          scale = c("none"),
          symbreaks = min(m, na.rm=TRUE),
          na.color="blue",
          dendrogram = "none", 
          Colv = FALSE,
          cellnote=round(m,2),
          notecex=1.0,
          Rowv = F,
          notecol="black",
          main = "R^2 heatmap len_before_5. Ne: " %+% Ne_val, 
          xlab = '% of recombination', 
          ylab = 'generations')

#################################################
data <- Total %>% filter(Ne == Ne_val, mutation == mutation_val) %>% 
  select(c(generation, len_before_10, recombination))
length(unique(data["generation"][[1]]))
length(unique(data["recombination"][[1]]))

m = matrix(c(0)*length(unique(data["generation"][[1]]))*length(unique(data["recombination"][[1]])), 
           nrow=length(unique(data["generation"][[1]])), ncol=length(unique(data["recombination"][[1]])))
m
colnames(m) <- as.character(sort(as.numeric(as.character(unique(data["recombination"][[1]])))))
rownames(m) <- as.character(sort(as.numeric(as.character(unique(data["generation"][[1]])))))

m

for (i in 1:nrow(data)){
  #print(i)
  m[as.character(as.numeric(as.character(data$generation[[i]]))), as.character(as.numeric(as.character(data$recombination[[i]])))] = data$len_before_10[[i]]
}
m




heatmap.2(m,
          col = colorpanel(100,"red","yellow","green"),
          margins = c(5, 10),
          trace = "none", 
          lhei = c(2, 8),
          scale = c("none"),
          symbreaks = min(m, na.rm=TRUE),
          na.color="blue",
          dendrogram = "none", 
          Colv = FALSE,
          cellnote=round(m,2),
          notecex=1.0,
          Rowv = F,
          notecol="black",
          main = "R^2 heatmap len_before_10. Ne: " %+% Ne_val, 
          xlab = '% of recombination', 
          ylab = 'generations')


#################################################
data <- Total %>% filter(Ne == Ne_val, mutation == mutation_val) %>% 
  select(c(generation, len_before_15, recombination))
length(unique(data["generation"][[1]]))
length(unique(data["recombination"][[1]]))

m = matrix(c(0)*length(unique(data["generation"][[1]]))*length(unique(data["recombination"][[1]])), 
           nrow=length(unique(data["generation"][[1]])), ncol=length(unique(data["recombination"][[1]])))
m
colnames(m) <- as.character(sort(as.numeric(as.character(unique(data["recombination"][[1]])))))
rownames(m) <- as.character(sort(as.numeric(as.character(unique(data["generation"][[1]])))))

m

for (i in 1:nrow(data)){
  #print(i)
  m[as.character(as.numeric(as.character(data$generation[[i]]))), as.character(as.numeric(as.character(data$recombination[[i]])))] = data$len_before_15[[i]]
}
m




heatmap.2(m,
          col = colorpanel(100,"red","yellow","green"),
          margins = c(5, 10),
          trace = "none", 
          lhei = c(2, 8),
          scale = c("none"),
          symbreaks = min(m, na.rm=TRUE),
          na.color="blue",
          dendrogram = "none", 
          Colv = FALSE,
          cellnote=round(m,2),
          notecex=1.0,
          Rowv = F,
          notecol="black",
          main = "R^2 heatmap len_before_15. Ne: " %+% Ne_val, 
          xlab = '% of recombination', 
          ylab = 'generations')

####################################################
data <- Total %>% filter(Ne == Ne_val, mutation == mutation_val) %>% 
  select(c(generation, len_before_20, recombination))
length(unique(data["generation"][[1]]))
length(unique(data["recombination"][[1]]))

m = matrix(c(0)*length(unique(data["generation"][[1]]))*length(unique(data["recombination"][[1]])), 
           nrow=length(unique(data["generation"][[1]])), ncol=length(unique(data["recombination"][[1]])))
m
colnames(m) <- as.character(sort(as.numeric(as.character(unique(data["recombination"][[1]])))))
rownames(m) <- as.character(sort(as.numeric(as.character(unique(data["generation"][[1]])))))

m

for (i in 1:nrow(data)){
  #print(i)
  m[as.character(as.numeric(as.character(data$generation[[i]]))), as.character(as.numeric(as.character(data$recombination[[i]])))] = data$len_before_20[[i]]
}
m




heatmap.2(m,
          col = colorpanel(100,"red","yellow","green"),
          margins = c(5, 10),
          trace = "none", 
          lhei = c(2, 8),
          scale = c("none"),
          symbreaks = min(m, na.rm=TRUE),
          na.color="blue",
          dendrogram = "none", 
          Colv = FALSE,
          cellnote=round(m,2),
          notecex=1.0,
          Rowv = F,
          notecol="black",
          main = "R^2 heatmap len_before_20. Ne: " %+% Ne_val, 
          xlab = '% of recombination', 
          ylab = 'generations')
dev.off()

