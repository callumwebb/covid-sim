
# Script for generating infection data, and animating results

library(tidyverse)
library(gganimate)
library(kernlab)
library(EpiDynamics)


#' Calculate the number of individuals in each group at every
#' time step
#'
#' @param sir_results results of SIR function
#' @param pop_size population size
#'
#' @return data frame of SIR values for each time step
calc_state_frames <- function(sir_results, pop_size) {
  numbers <- sir_results %>% 
    mutate(
      S = round(S * pop_size),
      I = round(I * pop_size),
      R = pop_size - S - I
    )
  numbers
}


#' Initialise visualisation matrices.
#' 
#' We will use three binary matrices to track the state of the pandemic.
#' S - is susceptible to infection
#' I - is currently infected
#' R - is "recovered" (includes those that die)
#'
#' @param p_size population size, must equal nrow * ncol
#' @param nrow number of rows
#' @param ncol number of columns
#' @param i0 number of infected at time zero
#'
#' @return a list of three matrices
init_matrices <- function(p_size, nrow, ncol, i0 = floor(p_size * 0.1)) {
  if (p_size != nrow * ncol) stop("population must be the product of nrow and ncol")
  
  S <- matrix(1L, nrow = nrow, ncol = ncol)
  I <- matrix(0L, nrow = nrow, ncol = ncol)
  R <- matrix(0L, nrow = nrow, ncol = ncol)
  
  initial_infected <- sample(1:p_size, i0)
  I[initial_infected] <- 1L 
  S[initial_infected] <- 0L
  
  list(s = S, i = I, r = R)
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


#' Calculate likelihood matrix for susceptible -> infected transition
#' 
#' Uses kernel density estimate based on infected locations
#'
#' @param S suseptible matrix
#' @param I infected matrix
#' @param bandwidth bandwidth of kernel
#' @param min_prob lower bound of probability for a susceptible individual,
#'   to allow a non-zero chance of a new spatial cluster forming
#'
#' @return a matrix of probability weights, for the purpose of random sampling
calc_infection_likelihood <- function(S, I, bandwidth = 1, min_prob = 0.2) {
  H <- diag(c(bandwidth, bandwidth))
  S_locs <- which(S == 1L, arr.ind = TRUE)
  
  infected <- which(I == 1L, arr.ind = TRUE)
  kde <- ks::kde(
    x = infected, 
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

