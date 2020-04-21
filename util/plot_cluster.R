# get cluster solution on specific time series
timeSeriSplitPlot <- function(k) {
  split_point <- fisher$run(k)$cluster %>% 
    extract(-length(.))
  
  plot(fisher$seq_data, type = 'l', main = sprintf('k = %s', k))
  split_point %>% 
    sapply(function(x) {
      x[1]
    }) %>% {
      abline(v = ., col = 'red', lwd = 2)
    }
}

