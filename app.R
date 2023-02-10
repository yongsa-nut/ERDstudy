#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.

library(shiny)
library(ggplot2)
library(cowplot)

## --- Backend Simulation --- ##

## Functions for simulation of the model ##

softmax <- function(par){
  par1 <- sort(par, decreasing = TRUE)
  Lk <- par1[1]
  for (k in 1:(length(par)-1)) {
    Lk <- max(par1[k+1], Lk) + log1p(exp(-abs(par1[k+1] - Lk))) 
  }
  val <- exp(par - Lk)
  return(val)
}

normalized <- function(x){
  if(sum(x)==0) return(x)
  return(x/sum(x))
}

cal_p_o <- function(prior_state, prior_obs, new_obs, time){
  ##Cal calculate P(O)  = Sum_{s}P(O|S)P(S)
  p_o_s <- obs_diff_num(prior_obs[[time]], new_obs)
  return(sum(prior_state*p_o_s))
}

gen_new_state <- function(prior_state, prior_obs, new_obs, time){
  ## P(S|O) = norm(P(O|S)*P(S))
  p_o_s <- obs_diff_num(prior_obs[[time]], new_obs)
  p_s   <- normalized(prior_state*p_o_s)
  ##if new_obs == time, then that means you get it on time so you know the outcome.
  if(new_obs == time){
    p_s[1:length(p_s)] <- 0
    p_s[time] <- 1
  }
  return(p_s)
}

obs_diff_num <- function(diff_dist, new_obs){
  ##Generate ERD from diff
  post_obs <- rep(0,10) #gen up to 10 times
  for(i in c(1:length(diff_dist))){
    ind <- new_obs + ind_to_diff[[i]]
    ind <- max(min(ind,10),1) ##anything above or below max/min are set to max/min
    post_obs[[ind]] <- post_obs[[ind]] + diff_dist[[i]]
  }
  return(post_obs)
}

reward_bel <- function(bel, switch_costs){
  ## The mapping is not necessary one to one.
  return(sum(bel[7:length(bel)])*switch_costs[runway])
}

zero_out_state <- function(state, time){
  if(time <= 0){
    return(normalized(state))
  }
  state[1:time] <- 0
  return(normalized(state))
}


Q_wait <- function(prior_state, prior_obs, switch_costs, time=1, depth=0){
  prior_state <- zero_out_state(prior_state, time-1) ## clean up in the case of prior state not clean before
  ##3  terminal conditions:
  if(depth == 0 || time == runway || (1 %in% prior_state)){
    return(reward_bel(prior_state, switch_costs))
  }
  q_temp <- 0
  
  ##Loop through possible observations upto 8
  ## Note this is quite slow and a lot of redundant 
  p_o_list   <- rep(0,8-time)
  q_max_list <- rep(0,8-time)
  ##zero out this time step since it's impossible in the future
  prior_state <- zero_out_state(prior_state, time = time)
  for(o in c((time+1):8)){
    ## Cal P(o) (Don't need to update prior_state to zero out time+1 because time+1 is still possible next)
    p_o_list[o - time] <- cal_p_o(prior_state, prior_obs, o, time+1)
    ## Step to new state (gen_post)
    new_state <- gen_new_state(prior_state, prior_obs, o, time+1)
    ## Cal q_wait of the new state
    q_wait <- Q_wait(new_state, prior_obs, switch_costs, time+1, depth-1)
    ## Cal q_switch of the new state
    q_switch <- switch_costs[time+1]
    ## pick the max and multiply by P(o) then sum for Q
    q_max_list[o - time] <- max(q_wait,q_switch)
  }
  q_temp <- sum(normalized(p_o_list)*q_max_list)
  
  return(q_temp)
}
##Q_wait(prior_state_norm, prior_obs_norm, sub_switch_costs, time = 5, depth = 1)

prob_decision <- function(prior_state, prior_obs, switch_costs, time = 1, depth=0, new_obs){
  ## Gen Probability for each action, wait = 1, switch = 2
  q_switch <- switch_costs[time]
  new_state <- gen_new_state(prior_state, prior_obs, new_obs, time)
  q_wait   <- Q_wait(new_state, prior_obs, switch_costs, time, depth)
  return(softmax(c(q_wait,q_switch)))
}

sim_one_seq <- function(prior_state, prior_obs, switch_costs, obs, depth){
  probs <- list()
  states <- list(prior_state)
  ##Only up to the last obs as it's either runway or delivered
  for(i in c(1:(length(obs)-1))){
    probs[[i]] <- prob_decision(prior_state, prior_obs, switch_costs, time=i, depth=depth, new_obs = obs[[i]])
    prior_state <- gen_new_state(prior_state, prior_obs, obs[[i]], time = i)
    states[[i+1]] <- prior_state
  }
  return(list(probs = probs,states = states))
}

update_obs <- function(prior, obs, alpha=1){
  ## Update prior obs_t from obs (diff between deliver day and msg) 
  ## obs in order of trial (length(obs) == current_trial)
  if(length(obs)==0) return(prior)
  
  post <- prior
  
  for(i in c(1:length(obs))){
    ind <- which(ind_to_diff == obs[[i]])
    post[ind] <- post[ind] + 1*alpha
    alpha <- alpha*alpha
  }
  return(post)
}

get_week_seq <- function(seqs, week){
  ## Gen week from a full trial x week
  ## in case you do not receive msgs every week
  ## Don't update when you receive the product because you don't receive msg.
  ## Essentially, you ignore the last one since it's either runway or received
  temp <- c()
  for(i in c(1:length(seqs))){
    if(week < length(seqs[[i]])){
      temp <- c(temp,seqs[[i]][[week]])
    }
  }
  return(temp)
  
  temp <- sapply(c(1:length(seqs)), function(x) {seqs[[x]][week]} ) 
  return(temp[!is.na(temp)])
}

update_obs_all_time <- function(prior, obs_list, alpha=1){
  ##Update all obs_t (t = week) given a full list
  for(i in c(1:runway)){
    obs_week <- get_week_seq(obs_list,i) 
    prior[[i]] <- update_obs(prior[[i]], obs_week, alpha)
  }
  return(prior)
}

update_state <- function(prior, outcome, beta=1){
  ## Update State given the outcome (delivered date)
  if(length(outcome)==0) return(prior)
  
  for(i in outcome){
    prior[i] <- prior[i] + 1*beta
    beta <- beta*beta
  }
  return(prior)
}

sim_n_seq <- function(prior_state, prior_obs, switch_cost, obs_list, outcome, num_trial=1, depth=0, alpha=1, beta=1){
  ##Simulate a whole experiment 
  ## return a list of probs of decision to wait, list of prob state, list of prob obs from trial 1 to num_trial
  ##general priors vs specific priors
  prior_state_temp <- normalized(prior_state)
  prior_obs_temp   <- lapply(c(1:runway), function(x) normalized(prior_obs[[x]]))
  results_probs  <- list()
  results_states <- list()
  prob_state_list <- list(prior_state_temp)
  prob_obs_list   <- list(prior_obs_temp)
  
  diff_obs_list <- diff_cal(obs_list, outcome)
  
  for(i in c(1:min(num_trial, length(obs_list)))){
    result_temp <- sim_one_seq(prior_state_temp, prior_obs_temp, switch_cost, obs_list[[i]], depth)
    results_probs[[i]]  <- result_temp$probs
    results_states[[i]] <- result_temp$states
    
    ##Update/Learning
    prior_state_temp <- normalized(update_state(prior_state, outcome[1:i], beta))
    prior_obs_temp   <- lapply(c(1:runway), function(x) normalized(update_obs_all_time(prior_obs, diff_obs_list[1:i], alpha)[[x]]) )
    
    prob_state_list[[i+1]] <- prior_state_temp
    prob_obs_list[[i+1]]   <- prior_obs_temp
  }
  return(list(prob_wait = results_probs, states = results_states, prob_state = prob_state_list, prob_obs = prob_obs_list))
}

## Visualization ## 

### Things to visualize: 
#### Within trial: P(s|o) and P(wait)
##### -> Test within trial 
#### After trial: updated P(o|s) and P(s)
##### -> Test the update

show_bel <- function(bel, type){
  ## Plot a distribution 
  bel <- normalized(bel)
  if(type=="Obs"){
    dat_temp <- data.frame(diff = seq(-3,3), prob = bel)
    g <- ggplot(dat_temp, aes(x=diff, y=prob)) + geom_line() +
      scale_x_continuous(breaks = seq(-3,3)) + ylim(0,1)
  }
  if(type=="State"){
    dat_temp <- data.frame(week = seq(1,10), prob = bel)
    g <- ggplot(dat_temp, aes(x=week, y = prob)) + geom_line() + 
      scale_x_continuous(breaks = seq(1,10)) + ylim(0,1)
  }
  return(g)
}

plot_bel_trial <- function(bel_state, bel_obs, title="Belief Distribution"){
  ##state + 6 weeks (runway)
  g_state  <- show_bel(bel_state,"State")  + ggtitle("P(S)")
  g_obs_w1 <- show_bel(bel_obs[[1]],"Obs") + ggtitle("P(O|S): week 1")
  g_obs_w2 <- show_bel(bel_obs[[2]],"Obs") + ggtitle("P(O|S): week 2")
  g_obs_w3 <- show_bel(bel_obs[[3]],"Obs") + ggtitle("P(O|S): week 3")
  g_obs_w4 <- show_bel(bel_obs[[4]],"Obs") + ggtitle("P(O|S): week 4")
  g_obs_w5 <- show_bel(bel_obs[[5]],"Obs") + ggtitle("P(O|S): week 5")
  ## Make no decision at week 6
  #g_obs_w6 <- show_bel(bel_obs[[6]],"Obs") + ggtitle("P(O|S): week 6") 
  g <- plot_grid(g_state, g_obs_w1, g_obs_w2, g_obs_w3, g_obs_w4, g_obs_w5, nrow=1)
  title <- ggdraw() + draw_label(title, fontface = 'bold', size = 14)
  g <- plot_grid(title, g, ncol=1, rel_heights = c(0.1,1))
  
  return(g)
}

plot_population_decision <- function(prob_waits, n_pop){
  wait_pop <- c()
  max_pop <- n_pop
  for(i in 1:length(prob_waits)){
    wait <- floor(n_pop*prob_waits[[i]][1])
    n_pop <- wait
    wait_pop <- c(wait_pop, wait)
  }
  wait_pop <- c(wait_pop, n_pop)
  dat_temp <- data.frame( week = c(1:(length(prob_waits)+1)), N_wait = wait_pop)
  g <- ggplot(dat_temp, aes(x = week, y=N_wait)) + geom_bar(stat = "identity") + 
    ylim(0,max_pop) + ggtitle("Number of people who decide to wait (N=100)")
  return(g)
}

plot_prob_wait <- function(prob_waits){
  dat_temp <- data.frame( week = c(1:length(prob_waits)), 
                          prob_of_wait = sapply(c(1:length(prob_waits)), function(x) prob_waits[[x]][[1]] ) )
  g <- ggplot(dat_temp ,aes(x = week, y = prob_of_wait)) + geom_line() + ylim(0,1) + 
    scale_x_continuous(breaks = c(1:length(prob_waits))) + ggtitle("Probability of wait")
  return(g)
}

plot_decision <- function(prob_waits, n_pop, title="Trial 0"){
  g_prob <- plot_prob_wait(prob_waits)
  g_popu <- plot_population_decision(prob_waits, n_pop)
  g <- plot_grid(g_prob, g_popu, nrow = 1) 
  title <- ggdraw() + draw_label(title, fontface = 'bold', size = 20)
  g <- plot_grid(title, g, ncol=1, rel_heights = c(0.1,1))
  return(g)
}


plot_states <- function(prob_states, title="P(S)/Subjective ERD Trial #0"){
  
  g_list <- list()
  for(i in c(1:length(prob_states))){
    g_list[[i]] <-  show_bel(prob_states[[i]],"State")+ 
      ggtitle(paste0("P(S|O) Week: ",i))
  }
  
  g <- plot_grid(plotlist = g_list, nrow=1)
  title <- ggdraw() + draw_label(title, fontface = 'bold', size = 14)
  g <- plot_grid(title, g, ncol=1, rel_heights = c(0.1,1))
  
  return(g)
}

## Predefined Parameters ## 
diff_cal <- function(msgs, delivered){
  out_msgs <- msgs
  for(t in c(1:7)){
    seq_len <- length(out_msgs[[t]])
    for(i in c(1:seq_len)){
      out_msgs[[t]][[i]] <- delivered[t] - msgs[[t]][[i]]
    }
  }
  return(out_msgs)
}

delivered_week   <- c(6,5,6,7,8,7,6)
delivered_week_rv <- c(7,8,7,6,5,6,6)

pushback_seq <- list(c(5,5,5,6,6,6), c(3,4,4,5,5),
                     c(4,4,4,5,6,6), c(5,5,5,6,7,7),
                     c(6,6,6,6,7,8), c(6,6,6,6,6,7),
                     c(5,5,5,5,6,6))
PB_LowC_seq <- pushback_seq

random_seq <- list(c(5,5,7,7,6,6), c(7,7,6,6,5),
                   c(5,5,5,6,6,6), c(7,7,6,6,7,7),
                   c(6,6,7,7,7,8), c(4,4,6,6,7,7),
                   c(5,5,5,7,7,6))

acc_seq  <- list(c(5,5,6,6,6,6), c(4,4,5,5,5),
                 c(6,6,6,6,6,6), c(6,6,6,7,7,7),
                 c(7,7,8,8,8,8), c(7,7,7,7,7,7),
                 c(5,5,5,6,6,6))

PB_RV_seq <- list(c(5,5,5,6,7,7), c(6,6,6,6,7,8),
                  c(6,6,6,6,6,7), c(5,5,5,6,6,6),
                  c(3,4,4,5,5), c(4,4,4,5,6,6),
                  c(5,5,5,5,6,6))

diff_push_back_seq <- diff_cal(pushback_seq, delivered_week)
diff_PB_LowC_seq <- diff_cal(PB_LowC_seq, delivered_week)
diff_random_seq <- diff_cal(random_seq, delivered_week)
diff_acc_seq <- diff_cal(acc_seq, delivered_week)
diff_PB_RV_seq <- diff_cal(PB_RV_seq, delivered_week_rv)

## Parameter Declalation ## 

runway <- 6
switch_costs <- -1*c(37500, 40000, 45000, 55000, 70000, 100000)
switch_costs_l <- -1*c(20000,24000, 30000, 38000, 50000, 100000)

sub_switch_costs   <- -1*c(3.75, 4, 4.5, 5.5, 7, 10)
sUb_switch_costs_l <- -1*c(2, 2.4, 3, 3.8, 5, 10)

ind_to_diff   <-  c(-3,-2,-1,0,1,2,3)
prior_obs     <- c(0.01,1/9, 2/9, 1, 2/9, 1/9, 0.01) 
prior_obs_time <- list(prior_obs,prior_obs,prior_obs,prior_obs,prior_obs,prior_obs)

##state is from t = 1 to 10?
prior_state   <- c(1,1,1,1,1,1,1,1,1,1)
prior_state_norm <- normalized(prior_state)
prior_obs_norm   <- lapply(c(1:6), function(x) normalized(prior_obs_time[[x]]))

## Other tunable parameters
alpha = 1
beta  = 1
depth = 0

## Variables for GUI

ws_type <- c("Push Back", "Random", "Acc", "Push Back Low Cost", "Push Back Reserved")

## --- FrontEnd --- ##
parsePriorState <- function(stateString){
  prior_state <- as.numeric(unlist(strsplit(stateString,",")))
  if( length(prior_state) != 10){
    return(FALSE)
  }
  return(normalized(prior_state))
}

parsePriorObs <- function(obsString){
  prior_obs <- as.numeric(unlist(strsplit(obsString,",")))
  if(length(prior_obs) != 7){
    return(FALSE)
  }
  ##repeat prior_obs for 6 time steps
  prior_obs <- normalized(prior_obs)
  prior_obs_time <- list(prior_obs, prior_obs, prior_obs,
                         prior_obs, prior_obs, prior_obs)
  return(prior_obs_time)
}

parseSWcost <- function(sw_cost){
  sw_cost <- as.numeric(unlist(strsplit(sw_cost,",")))
  if(length(sw_cost) != 6){
    return(FALSE)
  }
  return(-1*sw_cost)
}

erdPrint <- function(erd_seqs){
  str <- ""
  for(i in 1:length(erd_seqs)){
    str <- paste0(str,"Trial #",i,": ",toString(erd_seqs[[i]]), "\n")
  }
  return(str)
}

# Define UI for application that draws a histogram
ui <- fluidPage(
    titlePanel("ERD Simulation"),
    fluidRow(
      column(4, 
             h4("Wholesaler: ERD Sequences (#Trial = 7, Runway = 6)"),
             selectInput("wsType","WS type",ws_type),
             div(style="width:300px;",verbatimTextOutput("erdSequences"))
             ),
      column(4, 
             h4("Agent's main parameters"),
             ## Prior State 
             textInput("priorStateInput","Prior State (10 values from 1 to 10)",value="1,1,1,1,1,1,1,1,1,1"),
             ## Prior Observation
             textInput("priorObsInput","Prior Obs (7 values from -3 to 3)", value="0.01, 0.11, 0.22, 1, 0.22, 0.11, 0.01"),
             ## Switching cost 
             textInput("swcostInput","Switching Cost (Positive from 1 to 6/runway)", value="3.75, 4, 4.5, 5.5, 7, 10")
             ),
      column(4,
             h4("Hyperparameters"),
             numericInput("alphaInput", "P(S|O) Learnign Rate", value = 1),
             numericInput("betaInput", "P(S) Learning Rate", value = 1),
             numericInput("depthInput", "Depth", value = 0),
             actionButton("runButton","Run"),
             )
    ),
    hr(),
    plotOutput('plotBel0', height="250px"),
    plotOutput('plotDecision1'),
    plotOutput('plotERD1', height="250px"),
    plotOutput('plotBel1', height="250px"),
    
    plotOutput('plotDecision2'),
    plotOutput('plotERD2', height="250px"),
    plotOutput('plotBel2', height="250px"),
    
    plotOutput('plotDecision3'),
    plotOutput('plotERD3', height="250px"),
    plotOutput('plotBel3', height="250px"),
    
    plotOutput('plotDecision4'),
    plotOutput('plotERD4', height="250px"),
    plotOutput('plotBel4', height="250px"),
    
    plotOutput('plotDecision5'),
    plotOutput('plotERD5', height="250px"),
    plotOutput('plotBel5', height="250px"),
    
    plotOutput('plotDecision6'),
    plotOutput('plotERD6', height="250px"),
    plotOutput('plotBel6', height="250px"),
    
    plotOutput('plotDecision7'),
    plotOutput('plotERD7', height="250px"),
    plotOutput('plotBel7', height="250px")
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  simResults <-  eventReactive(input$runButton,{
    ##Setting up parameters 
    ### Wholesaler ERD sequence
    ws <- input$wsType
    obs_seq <- switch(
      ws, 
      "Push Back" = pushback_seq,
      "Random"    = random_seq,
      "Acc"       = acc_seq,
      "Push Back Low Cost" = PB_LowC_seq,
      "Push Back Reserved" = PB_RV_seq,
    )
    
    ### Agent's parameters 
    prior_state    <- parsePriorState(input$priorStateInput)
    prior_obs_time <- parsePriorObs(input$priorObsInput)
    sw_cost        <- parseSWcost(input$swcostInput)
    
    ### Hyperparameters
    beta  <- input$betaInput
    alpha <- input$alphaInput
    depth <- input$depthInput
    
    validate(
      need(prior_state,"Prior State's length is not 10,"),
      need(prior_obs_time,"Prior obs's length is not 7."),
      need(sw_cost,"Switch Cost's length is not 6.")
    )
    
    sim_n_seq(prior_state, prior_obs_time, sw_cost, 
              obs_seq, delivered_week, num_trial=7, 
              depth = depth, alpha = alpha, beta = beta)
  })
  output$plotBel0 <- renderPlot({
    result <- simResults()
    plot_bel_trial(result$prob_state[[1]], result$prob_obs[[1]],
                   "Prior Belief Distributions (P(S) and P(O|S))")
  })
  
  output$plotDecision1 <- renderPlot({
    result <- simResults()
    plot_decision(result$prob_wait[[1]], 100, paste0("Trial 1: (", toString(pushback_seq[[1]]),")"))
  })
  output$plotERD1 <- renderPlot({
    result <- simResults()
    plot_states(result$states[[1]], "P(S)/Subjective ERD Trial #1")
  })
  output$plotBel1 <- renderPlot({
    result <- simResults()
    plot_bel_trial(result$prob_state[[2]], result$prob_obs[[2]], "Belief Distribution Trial #1")
  })
  
  output$plotDecision2 <- renderPlot({
    result <- simResults()
    plot_decision(result$prob_wait[[2]], 100, paste0("Trial 2: (", toString(pushback_seq[[2]]),")"))
  })
  output$plotERD2 <- renderPlot({
    result <- simResults()
    plot_states(result$states[[2]], "P(S)/Subjective ERD Trial #2")
  })
  output$plotBel2 <- renderPlot({
    result <- simResults()
    plot_bel_trial(result$prob_state[[3]], result$prob_obs[[3]], "Belief Distribution Trial #2")
  })
  
  output$plotDecision3 <- renderPlot({
    result <- simResults()
    plot_decision(result$prob_wait[[3]], 100, paste0("Trial 3: (", toString(pushback_seq[[3]]),")"))
  })
  output$plotERD3 <- renderPlot({
    result <- simResults()
    plot_states(result$states[[3]], "P(S)/Subjective ERD Trial #3")
  })
  output$plotBel3 <- renderPlot({
    result <- simResults()
    plot_bel_trial(result$prob_state[[4]], result$prob_obs[[4]], "Belief Distribution Trial #3")
  })
  
  output$plotDecision4 <- renderPlot({
    result <- simResults()
    plot_decision(result$prob_wait[[4]], 100, paste0("Trial 4: (", toString(pushback_seq[[4]]),")"))
  })
  output$plotERD4 <- renderPlot({
    result <- simResults()
    plot_states(result$states[[4]], "P(S)/Subjective ERD Trial #4")
  })
  output$plotBel4 <- renderPlot({
    result <- simResults()
    plot_bel_trial(result$prob_state[[5]], result$prob_obs[[5]], "Belief Distribution Trial #4")
  })
  
  output$plotDecision5 <- renderPlot({
    result <- simResults()
    plot_decision(result$prob_wait[[5]], 100, paste0("Trial 5: (", toString(pushback_seq[[5]]),")"))
  })
  output$plotERD5 <- renderPlot({
    result <- simResults()
    plot_states(result$states[[5]], "P(S)/Subjective ERD Trial #5")
  })
  output$plotBel5 <- renderPlot({
    result <- simResults()
    plot_bel_trial(result$prob_state[[6]], result$prob_obs[[6]], "Belief Distribution Trial #5")
  })
  
  output$plotDecision6 <- renderPlot({
    result <- simResults()
    plot_decision(result$prob_wait[[6]], 100, paste0("Trial 6: (", toString(pushback_seq[[6]]),")"))
  })
  output$plotERD6 <- renderPlot({
    result <- simResults()
    plot_states(result$states[[6]], "P(S)/Subjective ERD Trial #6")
  })
  output$plotBel6 <- renderPlot({
    result <- simResults()
    plot_bel_trial(result$prob_state[[7]], result$prob_obs[[7]], "Belief Distribution Trial #6")
  })
  
  output$plotDecision7 <- renderPlot({
    result <- simResults()
    plot_decision(result$prob_wait[[7]], 100, paste0("Trial 7: (", toString(pushback_seq[[7]]),")"))
  })
  output$plotERD7 <- renderPlot({
    result <- simResults()
    plot_states(result$states[[7]], "P(S)/Subjective ERD Trial #7")
  })
  output$plotBel7 <- renderPlot({
    result <- simResults()
    plot_bel_trial(result$prob_state[[8]], result$prob_obs[[8]], "Belief Distribution Trial #7")
  })
  
  
  
  output$erdSequences <- renderText({
    ws <- input$wsType
    obs_seq <- switch(
      ws, 
      "Push Back" = pushback_seq,
      "Random"    = random_seq,
      "Acc"       = acc_seq,
      "Push Back Low Cost" = PB_LowC_seq,
      "Push Back Reserved" = PB_RV_seq,
    )
    erdPrint(obs_seq)
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
