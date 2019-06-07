args = commandArgs(trailingOnly = TRUE)

dir = toString(args[1])
file_in = toString(args[2])
file_out = toString(args[3])

setwd(dir)

saveRDS(read.table(file_in), file_out)

