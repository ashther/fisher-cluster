
# TODO find Inflection point of convex curve, confirm best cluster
FisherClust <- R6::R6Class(
  classname = 'FisherClust',
  
  public = list(
    data = NULL, 
    is_matrix = FALSE, 
    dmtx = NULL, 
    loss_mtx = NULL, 
    loss_idx = NULL, 
    
    # accept numeric vector, matrix and sparse matrix
    initialize = function(data, verbose = TRUE) {
      self$data <- data
      
      if (any(class(data) %in% c('dgCMatrix', 'dgTMatrix'))) {
        len <- nrow(data) - 1
        self$is_matrix <- TRUE
      } else if (any(class(data) %in% c('numeric', 'integer'))) {
        len <- length(data) - 1
      } else {
        rlang::abort(sprintf(paste0(
          'invalid data type, accept numeric vector, matrix or sparse matrix, ', 
          'not %s'
        ), class(data)[1]))
      }
      
      self$dmtx <- matrix(0, len, len, dimnames = list(2:(len + 1), seq_len(len)))
      
      len <- len - 1
      self$loss_mtx <- matrix(
        0, len, len, 
        dimnames = list(rownames(self$dmtx)[-1], colnames(self$dmtx)[-1])
      )
      
      self$loss_idx <- self$loss_mtx
    }, 
    
    build = function() {
      private$dmtx_create()
      private$loss_create()
    }, 
    
    run = function(k = 3) {
      private$get_cluster(k)
    }, 
    
    plot = function(k1, k2 = NULL) {
      if (is.null(k2)) {
        k2 <- ceiling(nrow(self$loss_mtx) / 2)
      } else {
        k2 <- min(k2, nrow(self$loss_mtx))
      }
      
      p1 <- purrr::map_dbl(seq_len(k2), ~ private$get_cluster(.x)$loss) %>% 
        dplyr::tibble(id = seq_along(.), x = .) %>% 
        ggplot2::ggplot(ggplot2::aes(id, x)) + 
        ggplot2::geom_point(size = 3, alpha = 0.5)
      
      p2 <- tibble::enframe(private$get_cluster(k1)$cluster) %>% 
        tidyr::unnest(value) %>% 
        dplyr::mutate_all(as.integer) %>% 
        ggplot2::ggplot(ggplot2::aes(name, value)) + 
        ggplot2::geom_point(size = 3, alpha = 0.5)
      
      gridExtra::grid.arrange(p1, p2, ncol = 2)
    }
  ), 
  
  private = list(
    dij_caculate = function(i, j) {
      avg <- mean(self$seq_data[i:j], na.rm = TRUE)
      round(sum((self$seq_data[i:j] - avg) ^ 2), 3)
    }, 
    
    dmtx_create = function() {
      if (self$is_matrix) {
        self$dmtx <- dMtxCreateCpp(as.matrix(self$data))
        dimnames(self$dmtx) <- list(seq_len(nrow(self$dmtx)) + 1, seq_len(ncol(self$dmtx)))
      } else {
        for (k in seq_len(nrow(self$dmtx))) {
          for (n in k:nrow(self$dmtx)) {
            self$dmtx[n, k] <- private$dij_caculate(k, n + 1)
          }
        }
      }
    }, 
    
    loss_create = function() {
      for (k in colnames(self$loss_mtx)) {
        for (n in rownames(self$loss_mtx)[which(colnames(self$loss_mtx) == k):nrow(self$loss_mtx)]) {
          min_loss <- minLossSplitCpp(self$dmtx, self$loss_mtx, as.integer(n), as.integer(k))
          self$loss_mtx[n, k] <- min_loss$value
          self$loss_idx[n, k] <- min_loss$idx
        }
      }
    }, 
    
    get_cluster = function(k) {
      name_mtx <- c('1', '2', rownames(self$loss_idx))
      n <- name_mtx[length(name_mtx)]
      cluster <- list()
      
      if (k <= 1) {
        return(list(cluster = list('1' = name_mtx), loss = self$dmtx[n, 1]))
      }
      
      # if (k >= as.integer(n)) {
      #   return(list(cluster = setNames(as.list(name_mtx), name_mtx), loss = 0))
      # }
      
      loss <- ifelse(
        k >= as.integer(n), 0, self$loss_mtx[n, as.character(k)]
      )
      while (TRUE) {
        if (k >= as.integer(n)) {
          for (i in rev(seq_len(as.integer(n)))) {
            cluster[[as.character(i)]] <- as.character(i)
          }
          break
        }
        
        if (!n %in% rownames(self$loss_idx) | k == 1) {
          cluster[['1']] <- name_mtx[1:which(name_mtx == n)]
          break
        }
        
        c_start <- self$loss_idx[n, as.character(k)]
        cluster[[as.character(k)]] <- name_mtx[
          which(name_mtx == as.character(c_start)):which(name_mtx == n)
        ]
        k <- k - 1
        n <- name_mtx[which(name_mtx == as.character(c_start)) - 1]
      }
      return(list(cluster = cluster, loss = loss))
    }
  )
)

# seq_data <- c(9.3, 1.8, 1.9, 1.7, 1.5, 1.3, 1.4, 2.0, 1.9, 2.3, 2.1)
# fisher <- FisherClust$new(seq_data)
# private <- fisher$.__enclos_env__$private
# self <- fisher$.__enclos_env__$self
# fisher$build()
# fisher$run(3)
