library(purrr)
library(dplyr)

# generate data
data <- rbeta(n = 10, shape1 = 80, shape2 = 80)

prob_k1 <- rbeta(n = 10, shape1 = 80, shape2 = 10)
prob_k2 <- 1-prob_k1


# perform operations on prob_k and data in a data.frame
operations_df <- tibble(components = c('1', '2'),
                        probability = list(prob_k1, prob_k2)) %>%
  
  # sum over list column
  mutate(n = map_dbl(probability, sum)) %>%
  
  # mean for each row, using list column and a single 1-element vector
  mutate(mu = map2_dbl(probability, n, ~ (1/.y) * sum(data * .x))) 

operations_df
  
# this doesn't work
# variance for each row, using list column, and two 1-element vectors
operations_df %>%
  mutate(var = pmap_dbl(probability, n, mu, ~ (1/(..2-1)) * sum(..1 * data^2) - ..3^2))

# this does work
(1/(operations_df$n[1]-1)) * sum(operations_df$probability[[1]] * data^2) - operations_df$mu[1]^2
(1/(operations_df$n[2]-1)) * sum(operations_df$probability[[2]] * data^2) - operations_df$mu[2]^2

# breaking it up into two map2 calls works:
operations_df %>%
  mutate(var = map2_dbl(n, probability, ~ (1/(.x-1)) * sum(.y * data^2))) %>%
  mutate(var = map2_dbl(var, mu, ~ .x - .y^2))
