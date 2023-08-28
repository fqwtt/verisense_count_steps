verisense_count_steps <- function(input_data=runif(500,min=-1.5,max=1.5), coeffs=c(0,0,0)) {
  # by Matthew R Patterson, mpatterson@shimmersensing.com
  ## Find peaks of RMS acceleration signal according to Gu et al, 2017 method
  # This method is based off finding peaks in the summed and squared acceleration signal
  # and then using multiple thresholds to determine if each peak is a step or an artefact.
  # An additional magnitude threshold was added to the algorithm to prevent false positives 
  # in free living data. 
  #
  # returns sample location of each step
  
  # 这种方法是基于在加速度信号的求和与平方后找到峰值，
  # 然后使用多个阈值来判断每个峰值是步伐还是干扰。
  # 为了防止在自由生活数据中出现假阳性，
  # 在算法中添加了一个额外的幅度阈值。
  # 返回每个步伐的采样位置
  
  # 采样率，这里是15HZ，可以手动设置。表示加速度计每秒采集15个数据点
  fs = 15 # temporary for now, this is manually set
  
  # 计算input_data的欧几里得范数，第一列的平方+第二列的平方+第三列的平方得到
  acc <- sqrt(input_data[,1]^2 + input_data[,2]^2 + input_data[,3]^2)
  
  # sd()函数用于计算一组数值数据的标准差
  # 这里如果标准差小于0.025，则认为没有步数
  if (sd(acc) < 0.025) {
    # acceleration too low, no steps
    
    # 计算加速度数据acc的持续时间（以秒为单位）
    num_seconds = round(length(acc) / fs)
    
    # 创建一个长度为num_seconds的向量steps_per_sec，并将其所有元素初始化为0。这表示在整个时间段内，每秒的步数都为0
    steps_per_sec = rep(0,num_seconds)
  } else {
    # Search for steps
    # Thresholds
    # 设置一些阈值，coeffs从1到8分别对应myscript.R里的parameters= c(3, 5, 15, -0.5, 3, 4, 0.001, 1.2)里的八个数字
    k <- coeffs[[1]] # 3
    period_min <- coeffs[[2]] # 5
    period_max <- coeffs[[3]] # 15
    sim_thres <- coeffs[[4]] # -0.5 # similarity threshold 相似性阈值
    cont_win_size <- coeffs[[5]] # 3  # continuity window size 连续窗口大小
    cont_thres <- coeffs[[6]] #4  # continuity threshold 连续性阈值
    var_thres <- coeffs[[7]] # 0.001  # variance threshold 方差阈值
    mag_thres <- coeffs[[8]] # 1.2  # 幅度阈值
    
    # find the peak rms value is every range of k
    half_k <- round(k/2) # 计算k的一半并四舍五入到最接近的整数
    segments <- floor(length(acc) / k) # 计算加速度数据acc的长度除以k的结果向下取整，将结果赋值给segments。表示将加速度数据分成segments个部分，以便在每个部分中查找峰值。
    peak_info <- matrix(NA,nrow=segments,ncol=5) # 创建一个矩阵peak_info，其行数等于segments，列数为5。矩阵的所有元素初始化为NA。这个矩阵可能用于存储在每个数据段中找到的峰值信息。
    # peak_info[,1] - peak location
    # peak_info[,2] - acc magnitude
    # peak_info[,3] - periodicity (samples)
    # peak_info[,4] - similarity
    # peak_info[,5] - continuity
    
    # for each segment find the peak location
    # 在加速度计数据（acc）的每个分段中查找峰值。
    # 它使用for循环遍历每个分段，并在每个分段中找到最大值（即峰值）
    for (i in 1:segments) {
      start_idx <- (i-1) * k + 1 # 计算当前分段的起始索引
      end_idx <- start_idx + (k-1) # 计算当前分段的结束索引
      tmp_loc_a <- which.max(acc[start_idx:end_idx]) # 在当前分段中找到最大值的位置
      tmp_loc_b <- (i-1) * k + tmp_loc_a # 将最大值的位置转换为原始加速度数据acc中的索引
      
      # only save if this is a peak value in range of -k/2:+K/2
      # 计算以tmp_loc_b为中心、长度为k的窗口的起始和结束索引。
      # 这里添加了边界检查，以确保索引在acc的范围内
      start_idx_ctr <- tmp_loc_b - half_k
      if (start_idx_ctr < 1) {
        start_idx_ctr <- 1
      }
      end_idx_ctr <- tmp_loc_b + half_k
      if (end_idx_ctr > length(acc)) {
        end_idx_ctr <- length(acc)
      }
      
      # 在以tmp_loc_b为中心的窗口中找到最大值的位置
      check_loc <- which.max(acc[start_idx_ctr:end_idx_ctr])
      
      # 检查最大值是否位于窗口的中心位置。如果是，则将其视为峰值，并将其索引和值存储在peak_info矩阵中
      if (check_loc == (half_k + 1)) {
        peak_info[i,1] <- tmp_loc_b
        peak_info[i,2] <- max(acc[start_idx:end_idx])
      }
    }
    
    # 检查peak_info矩阵的第一列是否包含NA值，仅保留值不是NA的行
    peak_info <- peak_info[is.na(peak_info[,1])!=TRUE,] # get rid of na rows
    
    # filter peak_info[,2] based on mag_thres
    # 根据mag_thres（幅度阈值）过滤peak_info矩阵中的峰值。
    # 即：仅保留那些第二列值（即峰值幅度）大于mag_thres的行。
    peak_info <- peak_info[peak_info[,2] > mag_thres,]
    
    # 88行，93行对peak_info做了两次过滤，要检查过滤后的peak_info矩阵中是否至少有两个步伐
    #（每个步伐由5个元素表示，因此至少需要10个元素）
    if (length(peak_info) > 10) {  # there must be at least two steps
      num_peaks <- length(peak_info[,1]) # 计算peak_info矩阵中的峰值数量
      
      no_steps = FALSE # 初始化一个名为no_steps的变量，其值为FALSE，用于表示是否检测到步伐。
      
      # 检查是否检测到至少两个峰值
      if (num_peaks > 2) {
        # Calculate Features (periodicity, similarity, continuity)
        
        # 计算相邻峰值之间的时间间隔，并将结果存储在peak_info矩阵的第三列中
        peak_info[1:(num_peaks-1),3] <- diff(peak_info[,1]) # calculate periodicity
        
        # 根据period_min过滤峰值，仅保留那些周期性大于period_min的峰值
        peak_info <- peak_info[peak_info[,3] > period_min,] # filter peaks based on period_min
        
        # 根据period_max过滤峰值，仅保留那些周期性小于period_max的峰值
        peak_info <- peak_info[peak_info[,3] < period_max,]   # filter peaks based on period_max 
      } else {
        no_steps = TRUE # 如果检测到的峰值数量小于等于2，则将no_steps变量设置为TRUE，表示未检测到步伐
      }
    } else {
      no_steps = TRUE # # 如果peak_info中的元素数量小于等于10，则将no_steps变量设置为TRUE，表示未检测到步伐
    }
    
    # 检查peak_info矩阵是否为空，即是否没有检测到峰值
    # 检查peak_info矩阵中的所有元素是否为NA，即是否所有峰值都被过滤掉了
    # 检查no_steps变量是否为TRUE，即是否未检测到步伐
    # 如果满足上述任一条件，则执行和29、32行一样的操作
    if ( length(peak_info)==0 || length(peak_info) == sum(is.na(peak_info)) || no_steps == TRUE) {
      # no steps found
      num_seconds = round(length(acc) / fs)
      steps_per_sec = rep(0,num_seconds)
    } else {
      # calculate similarity
      # 计算peak_info矩阵中的峰值数量
      num_peaks <- length(peak_info[,1]) 
      
      # 计算相邻峰值幅度之间的差异（二阶差分），并将结果的绝对值取负数。将结果存储在peak_info矩阵的第四列中，作为峰值之间的相似性度量
      peak_info[1:(num_peaks-2),4] <- -abs(diff(peak_info[,2],2)) # calculate similarity
      
      # 根据sim_thres过滤峰值，仅保留那些相似性大于sim_thres的峰值
      peak_info <- peak_info[peak_info[,4] > sim_thres,]  # filter based on sim_thres
      
      # 上一行代码可能会导致peak_info矩阵的第一列出现NA值。这行代码检查第一列中的每个元素是否为NA，并仅保留不包含NA值的行。
      peak_info <- peak_info[is.na(peak_info[,1])!=TRUE,] # previous statement can result in an NA in col-1
      
      
      # calculate continuity
      # 计算峰值之间的连续性，并根据一系列阈值和参数对峰值进行过滤
      if (length(peak_info[,3]) > 5) {
        end_for <- length(peak_info[,3])-1
        for (i in cont_thres:end_for) {
          # for each bw peak period calculate acc var
          v_count <- 0 # count how many windows were over the variance threshold
          for (x in 1:cont_thres) {
            if (var(acc[peak_info[i-x+1,1]:peak_info[i-x+2,1]]) > var_thres) {
              v_count = v_count + 1
            }
          }
          if (v_count >= cont_win_size) {
            peak_info[i,5] <- 1 # set continuity to 1, otherwise, 0
          } else {
            peak_info[i,5] <- 0
          }
        }
      } 
      peak_info <- peak_info[peak_info[,5]==1,1] # continuity test - only keep locations after this
      peak_info <- peak_info[is.na(peak_info)!=TRUE] # previous statement can result in an NA in col-1
      
      if (length(peak_info)==0) {
        # no steps found
        num_seconds = round(length(acc) / fs)
        steps_per_sec = rep(0,num_seconds)
      } else {
      
        # debug plot
        # is_plot = F
        # if (is_plot) {
        #   library(ggplot2)
        #   library(plotly)
        #   acc.df <- data.frame(acc=acc, det_step=integer(length(acc)))
        #   acc.df$det_step[peak_info] <- 1  # to plot annotations, prepare a 0/1 column on dataframe
        #   acc.df$idx <- as.numeric(row.names(acc.df))
        #   pl <- ggplot(data=acc.df,aes(x=idx,y=acc)) 
        #   pl2 <- pl + geom_line()
        #   pl3 <- pl2 + geom_point(data=subset(acc.df,det_step==1),aes(x=idx,y=acc),color='red',size=1,alpha=0.7)
        #   pl4 <- ggplotly(pl3)
        #   print(pl4)  
        # }
        
        # for GGIR, output the number of steps in 1 second chunks
        start_idx_vec <- seq(from=1,to=length(acc),by=fs)
        steps_per_sec <- table(factor(findInterval(peak_info, start_idx_vec), levels = seq_along(start_idx_vec)))
        steps_per_sec <- as.numeric(steps_per_sec)
      }
    }
  }
        
  return(steps_per_sec)

  # # 计算每分钟的步数
  # num_minutes <- 24 * 60
  # steps_per_minute <- numeric(num_minutes)
  # steps_per_sec_df <- data.frame(steps_per_sec)
  # write.csv(steps_per_sec_df, file = "~/Desktop/steps_per_sec_output.csv", row.names = FALSE, col.names = FALSE)
  # for (i in 1:num_minutes) {
  #   start_idx <- (i - 1) * 60  + 1
  #   end_idx <- i * 60
  #   if (end_idx > length(steps_per_sec)) {
  #     end_idx <- length(steps_per_sec)
  #   }
  #   steps_per_minute[i] <- sum(steps_per_sec[start_idx:end_idx])
  # }
  # 
  # # 返回一个包含24 * 60列的数据框
  # return(data.frame(steps_per_minute))
}
