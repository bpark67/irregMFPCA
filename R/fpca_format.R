#' Helper function to fill vector from a tibble
#'
#' @param v vector
#' @param target target to fill
#' @param source source to fill from
#'
#' @return a filled vector
#' @export
#'
vector_filler = function(v, target, source){
  N = length(v)
  uid = unique(source$id)
  for(k in 1:N){
    v[k] = source %>%
      filter(id == uid[k]) %>%
      select(target)  %>%
      pull() %>%
      list()
  }
  return(v)
}


#' Format functional data ready for FPCA analysis
#'
#' @param data original data
#'
#' @return an fpca formatted dataset
#' @export
#'
fpca_format = function(data){
  # Grab the number of components
  C = ncol(data)-3

  # Convert to FPCA friendly format
  df_tib = data %>% as_tibble()

  # Define variables
  uid = unique(df_tib$id)

  N = length(uid)
  Time = rep(0,N)

  for(c in 1:C){
    assign(paste0("Component", c), rep(0, N))
  }

  # Fill list with vectors of observation for each person
  Time = vector_filler(Time, target = "t", source = df_tib)

  for(c in 1:C){
    assign(paste0("Component", c),
           vector_filler(get(paste0("Component", c)),
                         target = paste0("component", c),
                         source = df_tib))
  }

  # Fill dataframe with each vector
  return_tibble = tibble(uid, Time)

  for(c in 1:C){
    return_tibble = return_tibble %>%
      add_column(v = get(paste0("Component", c)), .name_repair = "universal") %>%
      suppressWarnings()
  }

  colnames(return_tibble) = c("uid", "Time", paste0("Component", 1:C))
  # return(tibble(uid, Component1, Component2, Time))
  return(return_tibble)
}
