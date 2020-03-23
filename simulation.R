
# Functions for generating infection data, and animating results

library(tidyverse)
library(gganimate)
library(kernlab)
library(EpiDynamics)
library(grid)
library(gridExtra)


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
      E = round(E * pop_size),
      I = round(I * pop_size),
      R = round(R * pop_size),
      S = pop_size - E - I - R,
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


update_matrices <- function(S, E, I, R, MT, state_frame, bandwidth = 0.1, min_prob = 0.0) {
  # S -> E -> I -> R
  # But more than one transition can occur in a day, so the transitions need to be
  # processed in order, with calculations including upstream state changes
  
  # S -> E === S_E
  exposure_likelihood <- calc_exposure_likelihood(S, I, bandwidth, min_prob)
  S_E <- sum(S) - state_frame$S
  if (S_E > 0) {
    new_E <- sample(1:length(S), size = S_E, prob = exposure_likelihood)
    S[new_E] <- 0L
    E[new_E] <- 1L  
  }

  # E -> I
  dR <- state_frame$R - sum(R)
  dI <- state_frame$I - sum(I)
  E_I <- dI + dR
  if (E_I > 0) {
    new_I <- sample(1:length(E), size = E_I, prob = E)
    E[new_I] <- 0L
    I[new_I] <- 1L  
  }
  
  # Check those making the E -> I transition to see how 
  # many that would require hospitalisation can get it,
  # based on the current number of I cases requiring 
  # hospitalisation
  if (E_I > 0) {
    # This is a rough way of doing this, but it is for illustrative purposes
    new_cases_req_treat <- round(E_I * hospitalisation_rate)
    beds_max <- round(beds_per_1000pop / 1E3 * 
      (sum(S) + sum(E) + sum(I) + sum(R)))
    beds_occupied <- round(min(state_frame$I * hospitalisation_rate, beds_max))
    beds_available <- max(0, beds_max - beds_occupied)
    if (new_cases_req_treat > beds_available) {
      missed_treatment <- new_cases_req_treat - beds_available
      MT[sample(1:length(I), size = missed_treatment, prob = I)] <- 1L
    } 
  }
  
  # I -> R
  if (dR > 0) {
    new_R <- sample(1:length(I), size = dR, prob = I)
    I[new_R] <- 0L
    R[new_R] <- 1L 
  }
  
  list(s = S, e = E, i = I, r = R, mt = MT)
}


#' Incrementally generate ggplot2 plots for each time step (state frame)
#'
#' @param state_frames 
#' @param matrices 
#' @param bandwidth 
#' @param min_prob 
#' @param title 
#' @param pos 
#'
#' @return a list of ggplot2 plots
gen_sim_plots <- function(state_frames, matrices, bandwidth, min_prob, title, pos = "l") {
  # cat("\nrendering frame: 1\n")
  plots <- vector(mode = "list", length = nrow(state_frames))
  if (pos == "l") {
    margin_l <- 1
    margin_r <- 0.5
  } else {
    margin_l <- 0.5
    margin_r <- 1
  }
  plots[[1]] <- plot_frame(matrices, title, margin_l = margin_l, margin_r = margin_r)
  
  # Initialise missed treatment matrix
  matrices[["mt"]] <- matrix(0L, nrow = nrow(matrices$s), ncol = ncol(matrices$s))
  for (i in seq_len(nrow(state_frames) - 1)) {
    # cat(paste0("rendering frame:", i + 1, "\n"))
    matrices <- update_matrices(
      matrices$s,
      matrices$e,
      matrices$i,
      matrices$r,
      matrices$mt,
      state_frames[i + 1, ],
      bandwidth = bandwidth,
      min_prob = min_prob
    )

    p <- plot_frame(matrices, title, margin_l = margin_l, margin_r = margin_r) 
    plots[[i + 1]] <- p
  }
  plots
}


#' Turn a single state of S, E, I, R and MT matrices into a visualisation
#'
#' @param matrices 
#' @param title 
#' @param margin_l 
#' @param margin_r 
#'
#' @return a ggplot2 plot
plot_frame <- function(matrices, title, margin_l, margin_r) {
  state_mat <- gen_visual_matrix(matrices$s, matrices$e, matrices$i, matrices$r, matrices$mt) 
  
  p_data <- expand.grid(1:nrow(state_mat), 1:ncol(state_mat)) %>% 
    cbind(as.vector(state_mat)) %>% 
    set_names(c("y", "x", "state"))

  p <- ggplot(p_data, aes(x = x, y = y, fill = state)) +
    geom_raster() +
    coord_fixed() +
    scale_x_continuous(expand = expand_scale(0, 0)) +
    scale_y_continuous(expand = expand_scale(0, 0)) +
    scale_fill_manual(
      values = c(
        "S" = "#EEEEEE", 
        "E" = "#FFEE93", 
        "I" = "#FF5b5b",
        "R" = "#454545",
        "MT" = "#FF3030"
      )
    ) +
    ggtitle(title) +
    theme_void(base_size = 20) +
    theme(
      legend.position = "none",
      plot.title = element_text(
        hjust = 0,
        margin = margin(5, 0, 5, 0, unit = "pt")
      ),
      plot.margin = margin(5, margin_r, 5, margin_l, unit = "pt"),
    )
  p
}


#' Generate a single matrix from S, E, I, R, MT matrices that represents
#' colouring for each cell
#'
#' @param S 
#' @param E 
#' @param I 
#' @param R 
#' @param MT 
#'
#' @return a character matrix
gen_visual_matrix <- function(S, E, I, R, MT) {
  result <- matrix(NA_character_, nrow = nrow(S), ncol = ncol(S))
  result[which(S == 1)] <- "S"
  result[which(E == 1)] <- "E"
  result[which(I == 1)] <- "I"
  result[which(R == 1)] <- "R"
  result[which(MT == 1)] <- "MT"
  result
}


#' Save a list of png frames
#'
#' @param plots_A 
#' @param plots_B 
#' @param out_path 
#' @param break_at 
#'
#' @return
save_frames <- function(plots_A, plots_B, out_path = "./png/", break_at = 200) {
  for (i in seq_along(plots_A)) {
    if (i >= break_at) break
    cat(paste0("Saving frame ", i, "\n"))
    
    lay <- rbind(
      c(1, 1),
      c(2, 3),
      c(4, 4),
      c(5, 5)
    )
    title_grob <- textGrob(
      label = "#FlattenTheCurve",  hjust = 0, x = 0.02,
      gp = gpar(col = "black", fontsize = 32, fontface = "bold")
    )
    day_grob <- textGrob(
      label = paste0("Day: ", i, "\n"),
      gp = gpar(col = "black", fontsize = 22, fontfamily = "Helvetica Neue", fontface = "bold"),
      x = 0.02, hjust = 0, y = 0
    )
    blurb_grob <- textGrob(
      "People remaining red didn't get a hospital bed when they needed it",
      gp = gpar(col = "black", fontsize = 20, fontfamily = "Helvetica Neue", fontface = "italic"),
      hjust = 0, 
      x = 0.02, vjust = 0
    )
    gs <- list(title_grob, plots_A[[i]], plots_B[[i]], day_grob, blurb_grob)
    
    png(
      filename = paste0(out_path, formatC(i, width = 4, format = "d", flag = "0"), ".png"), 
      width = 700, height = 500, units = "px"
    )
    grid.arrange(grobs = gs, layout_matrix = lay, heights = c(0.6, 4, 0.5, 0.5))
    dev.off()
    # ggsave(
    #   paste0(out_path, formatC(i, width = 4, format = "d", flag = "0"), ".png"), 
    #   arrangeGrob(plots_A[[i]], plots_B[[i]], width = 4000, units = "px")
    # )
  }
}


gen_line_plots <- function(seir_baseline, seir_distancing, out_path = "./png/lines/") {
  plots <- vector(mode = "list", length = nrow(seir_baseline))
  for (i in seq_len(nrow(seir_baseline))) {
    plots[[i]] <- plot_lines_frame(seir_baseline, seir_distancing, i)
  }
  plots
}


plot_lines_frame <- function(seir_A, seir_B, current_t) {
  df_A <- seir_A %>% 
    pivot_longer(
      S:R,
      names_to = "variable",
      values_to = "value"
    ) %>% 
    mutate(scenario = "bau")
  df_B <- seir_B %>% 
    pivot_longer(
      S:R,
      names_to = "variable",
      values_to = "value"
    ) %>% 
    mutate(scenario = "isolate")
  df_lines <- bind_rows(df_A, df_B) %>% 
    filter(time <= current_t)
  
  df_points <- data.frame(
    t = current_t,
    variable = "I",
    scenario = c("bau", "isolate"),
    value = c(seir_A[current_t, "I"], seir_B[current_t, "I"])
  )
  
  df_infectious_line <- df_lines %>% 
    filter(variable == "I")
  
  # browser()
  
  p <- ggplot(
    df_infectious_line, 
    aes(x = time, y = value, colour = as.factor(scenario), group = as.factor(scenario))
  ) +
    geom_line() +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      name = "infectious",
      limits = c(0, 0.1)
    ) +
    scale_x_continuous(
      name = "day",
      limits = c(min(seir_A$time), max(seir_B$time)),
    ) +
    theme_minimal(base_size = 20) +
    theme(legend.position = "none")
  
  p
}


