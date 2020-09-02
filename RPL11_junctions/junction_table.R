library(data.table)
library(dplyr)


# go to directory with all the tidy bed files
setwd('tidy/')
# get paths to the files with RPL11.txt at the end of the file name
files <- list.files(pattern ="RPL11.txt$")

## make columns with start & end position, strand, and total count
all_info = data.table()
for(file in files){
  name = basename(file) # gets the file name stored as "name"
  name = gsub("_RPL11.txt", "", name) # removes _RPL11.txt extension from the name
  
  input = read.table(file) # reads the file as "input"
  input_dt = as.data.table(input) # changes it to a data.table
  names(input_dt) <- c("Start", "Stop", "Count", "Strand") # gives names to the columns
  
  all_info <- rbind(all_info, input_dt) # appends the input file to the end of all_info, making one big data.table
  
}
# get the total count for each junction and store it as the total_count column in junction_table
junction_table = all_info[,.(total_count=sum(Count)), by=.(Start, Stop, Strand)]
# sort the table by start position
junction_table = junction_table[order(Start),]


# add the labels corresponding to the bar plot
label_table = data.table(labels = list('skip', 'b', 'b2', 'a', 'c', 'd', 'normal', '?', '93neg', '93pos', '93pos2', 'long', 'paper'))
junction_table = cbind(junction_table, label_table)

# put data for individual samples as a column
for(file in files){
  name = basename(file) # stores the file name as "name"
  name = gsub("_RPL11.txt", "", name) # removes _RPL11.txt from the name
  
  input = read.table(file) # reads in the file at "input"
  input_dt = as.data.table(input) # makes the "input" as data.tale
  names(input_dt) <- c("Start", "Stop", "Count", "Strand") # gives the input column names
  
  junction_table = left_join(junction_table, input_dt, 
                             by = c("Start", "Stop", "Strand")) # takes the counts in each junction in the input file, adds it onto junction_table on the right
  
  names(junction_table)[ncol(junction_table)] <- name # makes the column name the name of the sample
  
}

# make the positions relative to the normal start & end positions
junction_table = mutate(junction_table, Start = Start - 23694791, Stop = Stop - 23695797)

# fill NAs with 0s
junction_table[is.na(junction_table)] <- 0

# makes the junction_table in writeable form
df <- apply(junction_table,2,as.character)

# write the final table into a file
write.table(df, 
            "all_junction_table.txt", 
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE) 
