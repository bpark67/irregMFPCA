#' Helper function to format simulated data into a long format
#'
#' @param dat original data
#'
#' @return long formatted data
#' @noRd
#'
long_format = function(dat){
  # Grab number of cases; last column is the time index
  # Grab length of time index as well; by default should be 51
  N = ncol(dat) - 1; l = length(t)

  # Pull out first column, which correponds to first case
  long_data = dat[, c(1, N + 1)] %>%
    na.omit() %>%
    mutate(id = "1") %>%
    `colnames<-`(c("value", "t", "id")) %>%
    select("id", "t", "value")

  # Iterate over remaining columns
  for(i in 2:N){
    long_data = long_data %>%
      rbind(dat[, c(i, N + 1)] %>%
              na.omit() %>%
              mutate(id = i) %>%
              `colnames<-`(c("value", "t", "id")) %>%
              select("id", "t", "value")
      )
  }

  long_data$id = as.numeric(long_data$id)
  long_data = long_data %>% arrange(id)
  return(long_data)
}

#' # Function that formats data into a tibble
#'
#' @param dat original data
#' @param tindex time vector
#'
#' @return a tibble formatted data
#' @export
#' @import tidyverse
#' @import brolgar
#'
#'
tibble_format = function(dat, tindex = seq(0, 1, length.out = 51)){
  dat = dat %>% t() %>% as.data.frame()
  # Grab number of cases, length of time index, and number of components
  N = ncol(dat) - 1; l = length(tindex); c = nrow(dat)/l

  # Split data into each component's dataframe
  chunks = split(1:nrow(dat), ceiling(1:nrow(dat)/l))
  list_dat = lapply(chunks, function(c) dat[c,])

  # Format into long format (Similar to reading "Walli" data)
  list_long_dat = lapply(list_dat, long_format)
  long_data = list_long_dat %>%
    reduce(left_join, by = c("id", "t")) %>%
    `colnames<-`(c("id", "t", paste0("component", 1:c)))

  df = as_tsibble(long_data,
                  key = "id",
                  index = "t",
                  regular = F)

  df = df %>%
    add_n_obs()

  return(df)
}
