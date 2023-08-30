rm(list = ls())
load("~/Desktop/verisense_count_steps/output/output_input/meta/basic/meta_4800014_90001_0_0.cwa.RData")
timestamp_per_5_secs <- M[["metashort"]][["timestamp"]]
steps_per_5_secs <- M[["metashort"]][["step_count"]]
num_minutes = length(steps_per_5_secs) %/% 12
timestamp_per_mins = rep('', num_minutes)
steps_per_mins = rep(0, num_minutes)
for(i in 1:num_minutes){
  timestamp_per_mins[i] = timestamp_per_5_secs[(i -1) * 12 + 1]
  for(j in 1:12){
    steps_per_mins[i] = steps_per_mins[i] + steps_per_5_secs[(i - 1) * 12 + j]
  }
}

file_name = "~/Desktop/verisense_count_steps/steps_per_min.csv"
my_data <- data.frame(timestamp_per_mins, steps_per_mins)
write.csv(my_data, file = file_name, row.names = FALSE)