
library(dplyr)
library(purrr)
library(text2vec)
library(jiebaR)

Rcpp::sourceCpp('src/minLossSplitCpp.cpp')
Rcpp::sourceCpp('src/dMtxCreate.cpp')
source('src/fisher_cluster.R')

seg <- worker(stop_word = 'dict/stop_words.utf8', write = 'NOFILE')
n_split <- 100

# readRDS('data/script.rds') %>% 
#   filter(title == '赛末点') %>%
#   pull(content) %>%
#   str_split('[:punct:](?=[\\w]+(（.*）)*：)') %>%
#   unlist() %>%
#   # paste0(collapse = '\n\n') %>%
#   writeLines('data/match_point.txt')
script <- readLines('data/match_point.txt') %>% 
  stringr::str_squish() %>% 
  grep('^$', ., invert = TRUE, value = TRUE) %>% 
  .[-1]
script <- tibble(id = seq_along(script), content = script) %>% 
  mutate(vid = id %/% floor(length(id) / n_split) + 1) %>% 
  group_nest(vid) %>% 
  mutate(words = map(data, ~ segment(.x$content, seg)))

it <- itoken(script$words)
v <- create_vocabulary(it)
  # prune_vocabulary(
  #   term_count_min = 2, 
  #   doc_count_min = 2
  # )
vectorizer <- vocab_vectorizer(v)
dtm <- create_dtm(it, vectorizer, type = 'dgCMatrix')

fc <- FisherClust$new(dtm)
fc$build()
clt <- fc$run(20)

result <- tibble::enframe(clt$cluster, 'cluster', 'vid') %>% 
  tidyr::unnest(vid) %>% 
  mutate(vid = as.integer(vid)) %>% 
  inner_join(script) %>% 
  select(cluster, vid, data) %>% 
  tidyr::unnest(data)

