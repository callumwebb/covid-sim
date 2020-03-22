
# Script for generating infection data, and animating results

library(tidyverse)
library(gganimate)
library(kernlab)
library(EpiDynamics)


#' Calculate the number of individuals in each group at every
#' time step
#'
#' @param sir_results results of SEIR function
#' @param pop_size population size
#'
#' @return data frame of SEIR values for each time step
calc_state_frames <- function(seir_results, pop_size) {
  numbers <- seir_results %>% 
    mutate(
      S = round(S * pop_size),
      E = round(E * pop_size),
      I = round(I * pop_size),
      R = pop_size - S - I - E,
      needs_hosp = round(pop_size * needs_hosp),
      treated = round(pop_size * treated),
      untreated = round(pop_size * untreated)
    ) 
  numbers
}


#' Initialise visualisation matrices.
#' 
#' We will use three binary matrices to track the state of the pandemic.
#' S - is susceptible to infection
#' E - is exposed
#' I - is currently infectious
#' R - is "recovered" (includes those that die)
#'
#' @param p_size population size, must equal nrow * ncol
#' @param nrow number of rows
#' @param ncol number of columns
#' @param e0 number of exposed
#'
#' @return a list of three matrices
init_matrices <- function(p_size, nrow, ncol, e0 = floor(p_size * 0.1)) {
  if (p_size != nrow * ncol) stop("population must be the product of nrow and ncol")
  
  S <- matrix(1L, nrow = nrow, ncol = ncol)
  E <- matrix(0L, nrow = nrow, ncol = ncol)
  I <- matrix(0L, nrow = nrow, ncol = ncol)
  R <- matrix(0L, nrow = nrow, ncol = ncol)
  
  initial_exposed <- sample(1:p_size, e0)
  E[initial_exposed] <- 1L 
  S[initial_exposed] <- 0L
  
  list(s = S, e = E, i = I, r = R)
}


#' Visualise a binary matrix
#'
#' @param mat matrix to visualise
#'
#' @return a ggplot2 plot
plot_matrix <- function(mat) {
  p_data <- expand.grid(1:nrow(mat), 1:ncol(mat)) %>% 
    cbind(as.vector(mat)) %>% 
    set_names(c("y", "x", "state"))
  p <- ggplot(p_data, aes(x = x, y = y, fill = state)) +
    geom_raster() +
    theme_void() +
    coord_fixed()
  p
}


#' Calculate likelihood matrix for susceptible -> exposed transition
#' 
#' Uses kernel density estimate based on infected locations
#'
#' @param S suseptible matrix
#' @param I infectious matrix
#' @param bandwidth bandwidth of kernel
#' @param min_prob lower bound of probability for a susceptible individual,
#'   to allow a non-zero chance of a new spatial cluster forming
#'
#' @return a matrix of probability weights, for the purpose of random sampling
calc_exposure_likelihood <- function(S, I, bandwidth = 1, min_prob = 0.2) {
  H <- diag(c(bandwidth, bandwidth))
  S_locs <- which(S == 1L, arr.ind = TRUE)
  
  # If no infectious currently, all are equally likely
  if (sum(I) == 0) {
    result <- matrix(0, nrow = nrow(S), ncol = ncol(S))
    result[S_locs] <- 1
    result
  } else {
    infectious <- which(I == 1L, arr.ind = TRUE)
    kde <- ks::kde(
      x = infectious, 
      H = H,
      eval.points = S_locs
    ) 
    
    # normalise output to range [min_prob, 1], assume fixed lower bound of 0
    # from kde output
    max_val <- max(kde$estimate)
    min_val <- min(kde$estimate)
    vals <- ((kde$estimate - min_val) / (max_val - min_val)) * 
      (1 - min_prob) + min_prob
    
    result <- matrix(0, nrow = nrow(S), ncol = ncol(S))
    result[S_locs] <- vals
    result
  }
}


update_matrices <- function(S, E, I, R, state_frame, bandwidth = 0.1, min_prob = 0.0) {
  browser()
  # S -> E -> I -> R
  # But more than one transition can occur in a day, so the transitions need to be
  # processed in order, with calculations including upstream state changes
  
  # S -> E === S_E
  exposure_likelihood <- calc_exposure_likelihood(S, I, bandwidth, min_prob)
  S_E <- sum(S) - state_frame$S
  new_E <- sample(1:length(S), size = S_E, prob = exposure_likelihood)
  S[new_E] <- 0L
  E[new_E] <- 1L
  
  # E -> I
  dR <- state_frame$R - sum(R)
  dI <- state_frame$I - sum(I)
  E_I <- dI + dR
  new_I <- sample(1:length(E), size = E_I, prob = E)
  E[new_I] <- 0L
  I[new_I] <- 1L
  
  # I -> R
  new_R <- sample(1:length(I), size = dR, prob = I)
  I[new_R] <- 0L
  R[new_R] <- 1L
  
  list(s = S, e = E, i = I, r = R)
}


gen_visual_matrix <- function(S, E, I, R) {
  result <- matrix(NA_character_, nrow = nrow(S), ncol = ncol(S))
  result[which(S == 1)] <- "S"
  result[which(E == 1)] <- "E"
  result[which(I == 1)] <- "I"
  result[which(R == 1)] <- "R"
  result
}




