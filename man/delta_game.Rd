% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/game_func.R
\name{delta_game}
\alias{delta_game}
\title{Simulates a 2-arm bandit game with binary or gaussian payoffs}
\usage{
delta_game(
  type,
  seed = NULL,
  trial_num,
  pay1,
  pay2,
  pay1_sd,
  pay2_sd,
  ev0 = 0.5,
  softmax_type = "norm",
  alphas_setup = "single",
  alpha,
  alphaM = NULL,
  alpha_obs = NULL,
  alpha_obsM = NULL,
  temp,
  p_observe = 0,
  p_observeA = 0.5,
  decay_obs = NULL,
  decay_bool = F,
  decay_function = "",
  growth_function = "hyper",
  inf_bonus = 0,
  decay_inf = 0,
  decay_temp = 0,
  growth_temp = 0
)
}
\arguments{
\item{type}{'binom': binary outcome; 'norm': gaussian outcome}

\item{seed}{seed}

\item{trial_num}{number of trials within a game}

\item{pay1}{payoff mean for arm 1}

\item{pay2}{payoff mean for arm 2}

\item{pay1_sd}{payoff sd for arm 1}

\item{pay2_sd}{payoff sd for arm 2}

\item{ev0}{numeric vector (2) specifying initial value belief (at the start of the game)}

\item{softmax_type}{'norm': with inverse temperature par (higher->more precision) 'inv': with temperature par (higher->less precision)}

\item{alphas_setup}{string vector specifying the number of learning rates. Options: 'single', 'obs', 'val' or 'valobs'}

\item{alpha}{learning rate for choice trials}

\item{alpha_obs}{learning rate for observational trials}

\item{temp}{temperature parameter (either inverse or not, dependent on the softmax_type argument)}

\item{p_observe}{probability of observing any arm}

\item{p_observeA}{probability of observing arm 1}

\item{decay_obs}{decay of p_observe parameter}

\item{decay_function}{'hyper': hyperbolic; 'exp': exponential}

\item{growth_function}{'hyper': hyperbolic}

\item{inf_bonus}{information bonus (not fully implemented)}

\item{decay_inf}{information bonus decay (not fully implemented)}

\item{decay_temp}{decay of the temp softmax parameter}

\item{growth_temp}{growth of the temp softmax parameter}
}
\value{
tibble containing a simulated game with choices, outcomes and expected values
}
\description{
Simulates a 2-arm bandit game with binary or gaussian payoffs
}
