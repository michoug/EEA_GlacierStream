suppressMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
fileList <- args[2]
fileout <- args[3]

dat <- read_tsv(filename, comment = "##", show_col_types = FALSE)

names(dat)[1] <- "query"

ecNumb <- c("3.2.1.20", "3.2.1.21", "3.2.1.14", "3.4.11.1", "3.1.3.1")
koNumb <- "ko:K03553"

dat_f <- dat %>%
  filter(EC %in% ecNumb | KEGG_ko %in% koNumb) %>%
  select(c(query, EC, KEGG_ko))


write_tsv(as.data.frame(dat_f$query), fileList, col_names = F)

write_tsv(dat_f, fileout)