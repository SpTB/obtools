
#' Simulates a 2-arm bandit game with binary or gaussian payoffs
#'
#' @param type 'binom': binary outcome; 'norm': gaussian outcome
#' @param seed
#' @param trial_num number of trials within a game
#' @param pay1 payoff mean for arm 1
#' @param pay2 payoff mean for arm 2
#' @param pay1_sd payoff sd for arm 1
#' @param pay2_sd payoff sd for arm 2
#' @param ev0 numeric vector (2) specifying initial value belief (at the start of the game)
#' @param softmax_type 'norm': with inverse temperature par (higher->more precision) 'inv': with temperature par (higher->less precision)
#' @param alpha learning rate for choice trials
#' @param alpha_obs learning rate for observational trials
#' @param temp temperature parameter (either inverse or not, dependent on the softmax_type argument)
#' @param p_observe probability of observing any arm
#' @param p_observeA probability of observing arm 1
#' @param decay decay of p_observe parameter
#' @param decay_function 'hyper': hyperbolic; 'exp': exponential
#' @param growth_function 'hyper': hyperbolic
#' @param inf_bonus information bonus (not fully implemented)
#' @param decay_inf information bonus decay (not fully implemented)
#' @param decay_temp decay of the temp softmax parameter
#' @param growth_temp growth of the temp softmax parameter
#'
#' @return
#' @export
#'
#' @examples
delta_game <- function(type, seed = NULL, trial_num, pay1, pay2, pay1_sd, pay2_sd, ev0 = 0.5, softmax_type = 'norm', alpha, alpha_obs=NULL, temp, p_observe=0, p_observeA=0.5, decay=NULL, decay_function='', growth_function='hyper',inf_bonus=0, decay_inf=0, decay_temp=0, growth_temp=0) {
#decay_funcion: currently can take arguments: 'exp' (exponential), 'hyper' (hyperbolic), 'qhyper' (quasi-hyperbolic)

  if(!is.null(seed)) set.seed(seed)
  if (is.null(alpha_obs)) alpha_obs = alpha #if not specified, set alpha for observation trials (alpha_obs) to the same value as alpha for choices (alpha)

   choice = reject = observed_outcomes = observation  =rep(-1,trial_num)
   p_observe = rep(p_observe, trial_num)
   inf_bonus = rep(inf_bonus, trial_num)
   temp      = rep(temp, trial_num)
   if (!is.null(decay)) {

     if (growth_function!='hyper') stop ("Only hyperbolic growth is currently implemented")

     for (t in 1:max(trial_num)) {
       if (decay_function=='') {
         stop('Please provide decay_function argument')
       }
       else if (decay_function=='exp') {
         p_observe[t] = wztools::decay_exp(A0=p_observe[1], k=decay,      t=t-1)
         inf_bonus[t] = wztools::decay_exp(A0=inf_bonus[1], k=decay_inf,  t=t-1)
         temp     [t] = wztools::decay_exp(A0=temp     [1], k=decay_temp, t=t-1)
       }
       else if (decay_function=='hyper') {
         p_observe[t] = wztools::decay_hyper (A0=p_observe[1], k=decay,      t=t-1)
         inf_bonus[t] = wztools::decay_hyper (A0=inf_bonus[1], k=decay_inf,  t=t-1)
         if (softmax_type=='norm') {
           temp[t] = wztools::growth_hyper(A0=temp[1], k=growth_temp, t=t-1)
         }
         else if (softmax_type=='inv') {
           temp[t] = wztools::decay_hyper(A0=temp [1], k=decay_temp, t=t-1)
         }
       }
       else if (decay_function=='qhyper') {
        stop('Not implemented yet!')
       }
       else {
         stop('decay_function argument needs to be either "exp", "hyper" or "qhyper"')
       }
     }
   }

  if (type=='binom') {
    outcomes = cbind(rbinom(trial_num, 1, pay1),
                     rbinom(trial_num, 1, pay2))

  } else if (type=='norm') {
    outcomes = cbind(rnorm(trial_num, pay1, pay1_sd),
                     rnorm(trial_num, pay2, pay2_sd))
  }

  ev = ev_after = matrix(ev0, trial_num, 2)


  for (t in 1:trial_num) {

    #update
    if(t>1) {
      if (choice[t-1]>-1) {
        ev[t,choice[t-1]] = ev[t-1,choice[t-1]] + pe * alpha
        ev[t,reject[t-1]] = ev[t-1,reject[t-1]]
      } else {
        ev[t,observation[t-1]] = ev[t-1,observation[t-1]] + pe * alpha_obs
        ev[t,abs(observation[t-1]-3)] = ev[t-1,abs(observation[t-1]-3)]
      }

    }
    #choose
    did_choose = rbinom(1, 1, (1-p_observe[t]))
    if (did_choose) {
      evt = ev[t,]
      # #!!!information bonus stuff not implemented
      # freqa = sum(c(observation, choice)==1) #how many times up to this point option A was chosen
      # freqb = sum(c(observation, choice)==2)
      # if (is.na(freqa)) freqa = 0
      # if (is.na(freqb)) freqb = 0
      # inf_diff = freqa-freqb
      #TODO: potentially, inf bonus can be scaled by the difference in inf
      #A = -(inf_bonus[t] * inf_diff) #bonus goes to the less chosen option
      A = 0
      evt[1] = evt[1] + A #adding the inf bonus (minus values substract bonus from evt[1], which is equivalent to adding the bonus for evt[2])
      if (softmax_type=='norm') prob1 = wztools::softmax(evt, temp[t]) else prob1 = wztools::softmax_inverse(evt, temp[t])# probability of choosing bandit 1
      choice[t] = rbinom(1,1, (1-prob1)) +1 # 1-> left; 2->right
      if (choice[t] == 1) reject[t] = 2
      if (choice[t] == 2) reject[t] = 1
      #update (delta)
      observed_outcomes[t] = outcomes[t, choice[t]]
      pe = observed_outcomes[t] - ev[t,choice[t]]
    } else {
      choice[t] = reject[t] = -1
      observation[t] = rbinom(1,1, p_observeA) + 1
      observed_outcomes[t] = outcomes[t, observation[t]]
      pe = observed_outcomes[t] - ev[t,observation[t]]
    }

    #update 2 (after)
    if (choice[t]>-1) {
      ev_after[t,choice[t]] = ev_after[t,choice[t]] + pe * alpha
      ev_after[t,reject[t]] = ev_after[t,reject[t]]
    } else {
      ev_after[t,observation[t]] = ev_after[t,observation[t]] + pe * alpha_obs
      ev_after[t,abs(observation[t]-3)] = ev_after[t,abs(observation[t]-3)]
    }

  } # t loop

  out = tibble::tibble(
    ev1 = ev[,1],
    ev2 = ev[,2],
    ev1_after = ev_after[,1],
    ev2_after = ev_after[,2],
    choice = choice,
    observation = observation,
    outcome = observed_outcomes,
    trial = 1:trial_num)
}


kalman_game<-function(trial_num, pay1, pay2,
                      choice_rule='softmax', cond_type,
                      temp, ev0, ev_var0, sigN, sigC, binom_outcomes, sd_gauss){
  #!currently unnecessary (could be useful if rewards distributed normally)
  #!currently doesn't support observed trials

  p1 = sample(pay1,1) # sample mean payment for bandit 1
  p2 =sample(pay2[pay2!=p1],1) # sample mean for bandit 2

  choice = reject = rep(-1,trial_num)
  if (binom_outcomes) {
    outcomes = cbind(rbinom(trial_num, 1, p1),
                     rbinom(trial_num, 1, p2)
                     )
  } else {
    outcomes = cbind(round(rnorm(trial_num, p1, sd_gauss)),
                     round(rnorm(trial_num, p2, sd_gauss))
    )
  }

  observed_outcomes = rep(-1, trial_num)
  ev = matrix(ev0, trial_num, 2)
  ev_var = matrix(ev_var0, trial_num, 2)

  # update (kalman)
  for (t in 1:trial_num) {

    #update beliefs
    if (t>1) {
      ev[t,] = updated_preds[[1]]
      ev_var[t,] = updated_preds[[2]]
    }
    #choose
    prob1 = wztools::softmax(ev[t,], temp) # probability of choosing bandit 1
    choice[t] = rbinom(1,1, (1-prob1)) +1 # 1-> left; 2->right


    if (choice[t] == 1) reject[t] = 2
    if (choice[t] == 2) reject[t] = 1
    observed_outcomes[t] = outcomes[t, choice[t]]

    #fit kalman
    updated_preds = wztools::kalman(vals = ev[t,], vars = ev_var[t,],
                           num_chosen = choice[t],
                           reward = observed_outcomes[t],
                           sigma_noise = sigN,
                           sigma_change = sigC)
  } # t-loop

  out = tibble::tibble(
    ev1 = ev[,1],
    ev2 = ev[,2],
    ev_var1 = ev_var[,1],
    ev_var2 = ev_var[,2],
    choice = choice,
    outcome = observed_outcomes,
    trial = 1:trial_num)

}

beta_softmax_opt_game <- function(trial_num, pay1, pay2, temp) {
  #obsolete: just a special case of beta_softmax_game
  choice = reject = rep(-1,trial_num)
  outcomes = cbind(rbinom(trial_num, 1, pay1),
                   rbinom(trial_num, 1, pay2))
  observed_outcomes = rep(-1, trial_num)

  #initial beliefs
  alphas = betas = matrix(1.0, trial_num, 2)
  ev             = matrix(0.5, trial_num, 2)

  for (t in 1:trial_num) {

    #update
    if(t>1) {
      alphas[t,choice[t-1]] = alphas[t-1,choice[t-1]] + outcomes[t-1, choice[t-1]]
      alphas[t,reject[t-1]] = alphas[t-1,reject[t-1]]
      betas [t,choice[t-1]] = betas [t-1,choice[t-1]] + (1-abs(outcomes[t-1, choice[t-1]]))
      betas [t,reject[t-1]] = betas [t-1,reject[t-1]]

      ev[t,choice[t-1]] = (alphas[t,choice[t-1]]/(alphas[t,choice[t-1]] + betas[t,choice[t-1]]))
      ev[t,reject[t-1]] = (alphas[t,reject[t-1]]/(alphas[t,reject[t-1]] + betas[t,reject[t-1]]))
    }
    #choose
    prob1 = wztools::softmax(ev[t,], temp) # probability of choosing bandit 1
    choice[t] = rbinom(1,1, (1-prob1)) +1 # 1-> left; 2->right
    if (choice[t] == 1) reject[t] = 2
    if (choice[t] == 2) reject[t] = 1

    observed_outcomes[t] = outcomes[t, choice[t]]


  } # t loop

  out = tibble::tibble(
    ev1 = ev[,1],
    ev2 = ev[,2],
    choice = choice,
    outcome = observed_outcomes,
    alpha1  = alphas[,1],
    beta1  = betas[,1],
    alpha2  = alphas[,2],
    beta2  = betas[,2],
    trial = 1:trial_num)
}

#' Title
#'
#' @param trial_num number of trials per game
#' @param pay1 payoff probability of the better arm
#' @param pay2 payoff probability of the worse arm
#' @param temp inverse temperature of the sofmax choice function
#' @param game_type str: 'blocked' or 'er'. 'blocked' means mixed game choice sequence is blocked: 10 observe + 10 choose. 'er' (event-related) means mixed game choice sequence is randomized
#' @param cond_type str: 'mixed' or 'free' 'free':
#' @param aW 2-num vector: weights for the alpha parameter (positive update)  (1st element for free-choice cond; 2nd <- for observe cond)
#' @param bW 2-num vector: weight vector for the beta parameter (negative update)
#'
#' @return generates a single game using beta distribution to represent values
#' @export
#'
#' @examples
beta_softmax_game <- function(trial_num, pay1, pay2, temp, game_type='blocked',
                              cond_type='mixed', aW=c(1,1), bW=c(1,1),
                              gamma=1, udb=0, usb=0, kf=1) {

  choice = reject = rep(-1,trial_num)
  outcomes = cbind(rbinom(trial_num, 1, pay1),
                   rbinom(trial_num, 1, pay2))
  observed_outcomes = rep(-1, trial_num)

  #ratings
  rating = rep(-1, trial_num)
  shape1 = shape2 = rep(1, trial_num)
  ratEV = rep(.5, trial_num)
  # rated = sample(c(rep(-1, round(trial_num/3*2)), rep(1, round(trial_num/6)), rep(2, round(trial_num/6)))) #due to rounding the vector will probably be too short
  # if(length(rating)>length(rated)) {
  #   rated = c(rep(-1, length(rating)-length(rated)), rated) #fill it with -1s at the start
  # } else if (length(rating)<length(rated)) {
  #   rated = rated[1:length(rating)]
  # }
  ##
  rated = sample(rep(c(0,0,0,1,1), (length(rating)/5) ))
  rated[1] = 0
  #rated[length(rated)]=0 #no first or last rating

  #initial beliefs
  alphas = betas = matrix(1.0, trial_num, 2)
  ev             = matrix(0.5, trial_num, 2)

  #set choice sequence
  if (cond_type=='free') {
    observedTrials = rep(0,trial_num) #no observed trials in the free condition
  } else if (game_type=='blocked') {
    observedTrials = rep(c(1,0),each=trial_num/2) #if mixed and blocked, first half observed
  } else if (game_type=='er') {
    observedTrials = sample(rep(c(1,0),each=trial_num/2)) # if mixed and event-related, randomize
  }

  free = as.numeric(!observedTrials)

  for (t in 1:trial_num) {

    #update
    if(t>1) {

      if (observedTrials[t-1]==1) { #observing weights
        aW_trial = aW[2]
        bW_trial  =bW[2]
      } else { #choosing weights
        aW_trial = aW[1]
        bW_trial  =bW[1]
      }

      alphas[t,choice[t-1]] = alphas[t-1,choice[t-1]] + (outcomes[t-1, choice[t-1]]*aW_trial)
      alphas[t,reject[t-1]] = alphas[t-1,reject[t-1]]
      betas [t,choice[t-1]] = betas [t-1,choice[t-1]] + ((1-abs(outcomes[t-1, choice[t-1]])) * bW_trial)
      betas [t,reject[t-1]] = betas [t-1,reject[t-1]]

      #explore bonus TODO: no negatives?
      ud = udb * ((alphas[t,1] + betas[t,1]) - (alphas[t,2] + betas[t,2])) /10
      us = usb * (alphas[t,2] + betas[t,2] + alphas[t,1] + betas[t,1]) /10

      ev[t,choice[t-1]] = (alphas[t,choice[t-1]]/(alphas[t,choice[t-1]] + betas[t,choice[t-1]]))
      ev[t,reject[t-1]] = (alphas[t,reject[t-1]]/(alphas[t,reject[t-1]] + betas[t,reject[t-1]]))

      #update difference distribution pars
      shape1[t] = alphas[t,1] + betas[t,2]
      shape2[t] = alphas[t,2] + betas[t,1]
      ratEV[t] =  shape1[t] /  (shape1[t] + shape2[t])
      #rate (needs to happen AFTER update!)
      if (rated[t-1]>0)  {
        rating[t-1] = rbeta(n=1, shape1=shape1[t],shape2=shape2[t])
      }

    }
    #choose
    if (observedTrials[t]==1) { # observed (random) case
      choice[t] = rbinom(n=1,size=1, prob=.5)+1 # 1-> left; 2->right;
    } else { #softmax choice
      #add uncert bonus to the EVB vector (bad arm; if +: increase bad arm, if -: decrease)
      if (t==1) evb = ev[t,] else evb = c(ev[t,1] ,ev[t,2] +us + ud)
      prob1 = wztools::softmax(evb, temp) # probability of choosing bandit 1
      choice[t] = rbinom(1,1, (1-prob1)) +1 # 1-> left; 2->right
    }
    if (choice[t] == 1) reject[t] = 2
    if (choice[t] == 2) reject[t] = 1

    observed_outcomes[t] = outcomes[t, choice[t]]

  } # t loop

  out = tibble::tibble(
    free = !observedTrials,
    payB = rep(pay2,trial_num),
    ev1 = ev[,1],
    ev2 = ev[,2],
    choice = choice,
    outcome = observed_outcomes,
    rated = rated,
    rating = rating,
    alpha1  = alphas[,1],
    beta1  = betas[,1],
    alpha2  = alphas[,2],
    beta2  = betas[,2],
    alpha3 = shape1,
    beta3 = shape2,
    ratEV = ratEV,
    trial = 1:trial_num)
}

beta_game <- function(trial_num, pay1, pay2,
                      game_type='blocked', choice_rule,  cond_type='mixed',
                      aW=c(1,1), bW=c(1,1), gamma=1, temp, udb=0, usb=0, kf=1, ...) {

  choice = reject = rep(-1,trial_num)
  outcomes = cbind(rbinom(trial_num, 1, pay1),
                   rbinom(trial_num, 1, pay2))
  observed_outcomes = rep(-1, trial_num)

  #ratings
  rating = rep(-1, trial_num)
  shape1 = shape2 = rep(1, trial_num)
  ratEV = rep(.5, trial_num)
  # rated = sample(c(rep(-1, round(trial_num/3*2)), rep(1, round(trial_num/6)), rep(2, round(trial_num/6)))) #due to rounding the vector will probably be too short
  # if(length(rating)>length(rated)) {
  #   rated = c(rep(-1, length(rating)-length(rated)), rated) #fill it with -1s at the start
  # } else if (length(rating)<length(rated)) {
  #   rated = rated[1:length(rating)]
  # }
  ##
  rated = sample(rep(c(0,0,0,1,1), (length(rating)/5) ))
  rated[1] = 0

  #initial beliefs
  alphas = betas = matrix(1.0, trial_num, 2)
  ev             = matrix(0.5, trial_num, 2)

  #set choice sequence
  if (cond_type=='free') {
    observedTrials = rep(0,trial_num) #no observed trials in the free condition
  } else if (game_type=='blocked') {
    observedTrials = rep(c(1,0),each=trial_num/2) #if mixed and blocked, first half observed
  } else if (game_type=='er') {
    observedTrials = sample(rep(c(1,0),each=trial_num/2)) # if mixed and event-related, randomize
  }

  free = as.numeric(!observedTrials)

  for (t in 1:trial_num) {

    #update
    if(t>1) {

      if (observedTrials[t-1]==1) { #observing weights
        aW_trial = aW[2]
        bW_trial  =bW[2]
      } else { #choosing weights
        aW_trial = aW[1]
        bW_trial  =bW[1]
      }

      alphas[t,choice[t-1]] = alphas[t-1,choice[t-1]] + (outcomes[t-1, choice[t-1]]*aW_trial)
      alphas[t,reject[t-1]] = alphas[t-1,reject[t-1]]
      betas [t,choice[t-1]] = betas [t-1,choice[t-1]] + ((1-abs(outcomes[t-1, choice[t-1]])) * bW_trial)
      betas [t,reject[t-1]] = betas [t-1,reject[t-1]]

      #explore bonus TODO: no negatives?
      ud = udb * ((alphas[t,1] + betas[t,1]) - (alphas[t,2] + betas[t,2])) /10
      us = usb * (alphas[t,2] + betas[t,2] + alphas[t,1] + betas[t,1]) /10

      ev[t,choice[t-1]] = (alphas[t,choice[t-1]]/(alphas[t,choice[t-1]] + betas[t,choice[t-1]]))
      ev[t,reject[t-1]] = (alphas[t,reject[t-1]]/(alphas[t,reject[t-1]] + betas[t,reject[t-1]]))

      #update difference distribution pars
      shape1[t] = alphas[t,1] + betas[t,2]
      shape2[t] = alphas[t,2] + betas[t,1]
      ratEV[t] =  shape1[t] /  (shape1[t] + shape2[t])
      #rate (needs to happen AFTER update!)
      if (rated[t-1]>0)  {
        rating[t-1] = rbeta(n=1, shape1=shape1[t],shape2=shape2[t])
      }

    }
    #choose
    if (observedTrials[t]==1) { # observed (random) case
      choice[t] = rbinom(n=1,size=1, prob=.5)+1 # 1-> left; 2->right;
    } else {
      if (choice_rule=='softmax') {
        #add uncert bonus to the EVB vector (bad arm; if +: increase bad arm, if -: decrease)
        if (t==1) evb = ev[t,] else evb = c(ev[t,1] ,ev[t,2] +us + ud)
        probA = wztools::softmax(evb, temp) # probability of choosing bandit 1
      } else if (choice_rule=='thompson') {
        probA = wztools::thompson(alphas[t,], betas[t,])
      }

      choice[t] = rbinom(1,1, (1-probA)) +1 # 1-> left; 2->right
    }
    if (choice[t] == 1) reject[t] = 2
    if (choice[t] == 2) reject[t] = 1

    observed_outcomes[t] = outcomes[t, choice[t]]

  } # t loop

  out = tibble::tibble(
    free = !observedTrials,
    payB = rep(pay2,trial_num),
    ev1 = ev[,1],
    ev2 = ev[,2],
    choseA = 1-(choice-1),
    choice = choice,
    outcome = observed_outcomes,
    rated = rated,
    rating = rating,
    alpha1  = alphas[,1],
    beta1  = betas[,1],
    alpha2  = alphas[,2],
    beta2  = betas[,2],
    alpha3 = shape1,
    beta3 = shape2,
    ratEV = ratEV,
    trial = 1:trial_num)
}

