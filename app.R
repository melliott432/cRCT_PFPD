library(shiny)
library(ggplot2)
require(geeCRT)
require(dplyr)
require(logitnorm)
require(mvtnorm)
require(R2jags)
require(rjags)
require(viridis)

# ── Placeholder functions – replace these with your own ──────────────────────

simulate_data <- function(run_name, parms){
  
  # Required packages
  require(geeCRT)
  require(dplyr)
  require(logitnorm)
  require(mvtnorm)
  
  # Define parameters
  n_village = parms$n_village # number of clusters
  n_ind = parms$n_ind # common cluster size
  n_arm = parms$n_arm # number of arms
  pii = parms$pii # follow-up rate
  tau = parms$tau # attrition ICC
  delta = parms$delta # treatment effect
  sd_mult = parms$sd_mult # common standard deviation
  p1 = parms$p1 # baseline or control outcome proportion
  n_ind_total = n_village*n_ind # total number of individuals in the trial
  beta1 = logit(delta + p1) # logit-transformed treatment effect 
  beta0 = logit(p1) # logit-transformed baseline
  
  # Set seed for reproducibility
  seed = as.numeric(pii*100 + tau*100 + delta*100) + n_village + n_ind
  set.seed(seed)
  
  ##############################################################################
  ## ATTRITION
  
  # Beta distribution for cluster-level proportions
  a = (pii*(1-tau))/tau
  b = (1/tau) - 1 - a
  pr_cluster = rbeta(n = n_village, shape1 = a, shape2 = b) # Cluster-level follow-up rate
  village_flwup_rate = rep(pr_cluster, each = n_ind)
  
  # Binomial distribution for individual attrition indicators
  followup_final_indic <- c() # R = 1 is observed, R = 0 is missing
  for(i in 1:n_village){
    Rvec <- rbinom(n = n_ind, size = 1, prob = pr_cluster[i])
    followup_final_indic <- c(followup_final_indic, Rvec)
  }
  
  ##############################################################################
  ## ADHERENCE
  
  # Create std devs for random effects
  village_sd = 1*sd_mult # village-level RE
  #clinic_sd = 1*sd_mult # clinic-level RE
  #household_sd = 1*sd_mult # household-level RE
  individual_sd = 1*sd_mult # individual-level RE
  
  # Create ID vectors
  ind_id = c(1:n_ind_total)
  arm_id = c(rep(c(1:n_arm), each = n_ind_total/n_arm))
  village_id = c(rep(c(1:n_village), each = n_ind))
  
  # Intervention (arm) effect
  arm_eff = c(rep(beta0, length = n_ind_total/2), rep(beta1, length = n_ind_total/2))
  
  # Village (cluster) effect
  village_eff = rnorm(n_village, mean = 0, sd = village_sd)
  village_eff = rep(village_eff, each = n_ind)
  
  # Individual effect
  individual_eff = rnorm(n_ind_total, mean = 0, sd = individual_sd)
  
  # Sample baseline adherence indicators
  adherence_base_logit = village_eff + individual_eff
  adherence_base_prob = invlogit(adherence_base_logit)
  adherence_base_indic = rbinom(n_ind_total, 1, prob = adherence_base_prob)
  
  # Sample follow-up adherence indicators
  adherence_final_logit = village_eff + individual_eff + arm_eff
  adherence_final_prob = invlogit(adherence_final_logit)
  adherence_final_indic = rbinom(n_ind_total, 1, prob = adherence_final_prob)
  
  # Set missing observations to NA
  adherence_final_indic[which(followup_final_indic == 0)] = NA
  
  ##############################################################################
  
  # Collate data set
  sim_trial_df <- data.frame(ind_id, village_id, arm_id, 
                             individual_eff, village_eff, arm_eff, 
                             adherence_base_logit, adherence_base_prob, adherence_base_indic, 
                             adherence_final_logit, adherence_final_prob, adherence_final_indic,
                             village_flwup_rate, followup_final_indic)
  
  return(sim_trial_df)
  
}



run_model <- function(df) {
  
  ##############################################################################
  #### Format data
  
  not_missing = which(is.na(df$adherence_final_indic) == FALSE) # index which observations are not missing
  n_village = length(unique(df$village_id)) # number of villages
  
  # select non-missing outcome observations
  village <- df$village_id[not_missing] # village IDs
  Y <- df$adherence_final_indic[not_missing] # adherence outcome
  X <- df$arm_id[not_missing] - 1 # arm assignment
  
  # calculate cluster follow-up proportions 
  Ri <- df |>
    group_by(village_id) |> # group by village
    summarise(Rprop = sum(followup_final_indic)/length(followup_final_indic)) # calculate proportion
  Ri <- Ri$Rprop
  Ri <- ifelse(Ri>0.999, 0.999, Ri) # set bounds so JAGS doesn't break
  Ri <- ifelse(Ri<0.001, 0.001, Ri)
  R <- df$followup_final_indic # pull individual follow up indicators
  
  # when clusters are unequal sizes
  # create index for individual follow-up variable
  index <- c()
  index[1] <- 1
  for(i in 1:n_village){index[i] <- min(which(df$village_id == i))}
  index[n_village+1] <- (nrow(df)+1)
  
  jagsData <- list("Y" = Y, # outcome
                   "X" = X, # arm
                   "R" = R, # follow-up
                   "index" = index, # cluster sizes
                   "village" = village, # village id
                   "J" = length(Y), # number of outcome observations
                   "N" = length(unique(df$village_id)), # number of villages
                   "pisq" = pi^2) # 3.14^2 (because JAGS doesn't have this)
  
  ############################################################################
  ### Define model for estimation
  
  model_string <- # model as string for JAGS to read
    
    "model {
        
        ### Likelihood ###
        
        for(j in 1:J){ # Individual level
          Y[j] ~ dbern(ilogit(logitp[j]))
          logitp[j] = beta0 + beta1*X[j] + eta[village[j]] 
        }
      
        for(n in 1:N){ # Cluster level
            for(m in index[n]:(index[n+1]-1)){
              R[m] ~ dbern(r[n])
            }
          r[n] ~ dbeta(a, b) T(0.001,0.999)
          eta[n] ~  dnorm(0, prec_eta)
        }
      
        ### Priors ###
        
        beta0 ~ dunif(-10.0, 10.0)
        beta1 ~ dunif(-10.0, 10.0)
        a ~ dunif(0, 100.0)
        b ~ dunif(0, 100.0)
        prec_eta ~ dunif(0, 1000)
        
        p1 = exp(beta0)/(exp(beta0)+1)
        p2 = exp(beta0+beta1)/(exp(beta0+beta1)+1)
        delta = p2 - p1
        tau = 1/(a+b+1)
        pii = a/(a+b)
        rho_y = (1/prec_eta) / ((1/prec_eta) + (pisq/3))
        
      }"
  
  # Compile and initialize JAGS model object
  model <- jags.model(file = textConnection(model_string), data = jagsData, n.chains = 4)
  
  # Run warm up iterations
  update(model, n.iter = 4000)
  
  # Sample from posterior
  posterior_sample <- coda.samples(model, n.iter = 5000,
                                   variable.names = c("pii", "tau", "delta", "a", "b", "rho_y", "beta0", "beta1"))
  
  # Extract posterior samples
  posterior <- as.data.frame(rbind(posterior_sample[[1]], posterior_sample[[2]], posterior_sample[[3]], posterior_sample[[4]]))
  
  return(posterior)
  
}

CalcPower <- function(alpha, n, m, n_arms, pii, tau, rho, p1, p2){
  
  z.a2 = qnorm(alpha/2)
  
  k = n/n_arms
  
  G <- 1 + (((m*pii)-1)*rho)+(1-pii)*(1+(m-1)*tau)*rho
  
  z.b <- sqrt(k*m*pii*(1/G)*((p1-p2)^2)*(((p2*(1-p2))+(p1*(1-p1)))^(-1))) + z.a2 
  
  power <- pnorm(z.b)
  
  return(power)
}



run_power <- function(parms, estimates){
  
  B = 1000
  
  power.dist <- data.frame( "alpha" = rep(parms$alpha,B),
                            "n" = rep(parms$n_village,B),
                            "m" = rep(parms$n_ind,B),
                            "n_arms" = rep(parms$n_arm,B),
                            "pii" = sample(estimates$pii, B, replace = T),
                            "tau" = sample(estimates$tau, B, replace = T),
                            "rho" = sample(estimates$rho_y, B, replace = T),
                            "delta" = sample(estimates$delta, B, replace = T),
                            "beta0" = sample(estimates$beta0, B, replace = T))
  
  power.dist <- power.dist |>
    mutate(p1 = exp(beta0)/(exp(beta0)+1),
           p2 = p1 + delta,
           power = CalcPower(alpha = alpha, n = n, m = m, n_arms = n_arms, pii = pii, tau = tau, rho = rho, p1 = p1, p2 = p2))
  
  return(power.dist)
}



# ── UI ───────────────────────────────────────────────────────────────────────

ui <- navbarPage(
  
  title = "Power Calculator",
  
  tags$head(tags$style(HTML("
    body { font-family: sans-serif; background-color: #f5f5f5; }

    .section-box {
      background: #e8e8e8;
      border: 1px solid #ccc;
      border-radius: 4px;
      padding: 14px 16px;
      margin-bottom: 14px;
      box-sizing: border-box;
      width: 100%;
      overflow: hidden;
    }
    .section-box h4 { margin-top: 0; font-size: 15px; font-weight: normal; color: #33498A; }
    .section-box p  { font-size: 12px; color: #555; margin: 2px 0 10px 0; }
    .subsection-label {
      font-weight: bold;
      font-size: 13px;
      margin: 10px 0 4px 0;
    }

    /* inline numeric inputs */
   .inline-inputs {
      display: flex;
      align-items: center;
      gap: 6px;
      margin-bottom: 6px;
      width: 100%;
      box-sizing: border-box;
    }
    .inline-inputs label {
      flex: 1;
      font-size: 13px;
      white-space: normal;
      word-wrap: break-word;
      min-width: 0;
    }
    .inline-inputs .form-group {
      margin: 0;
      flex-shrink: 0;
    }
    .inline-inputs .form-group input[type=text],
    .inline-inputs .form-group input[type=number] {
      width: 70px !important;
      min-width: 0;
      padding: 2px 4px;
      font-size: 13px;
      box-sizing: border-box;
    }

    /* Run button */
    #run_btn {
  width: 100%;
  background-color: #33498A;
  color: white;
  font-size: 14px;   /* ← reduce this, was 18px */
  font-weight: bold;
  border: none;
  border-radius: 4px;
  padding: 14px;
  cursor: pointer;
  letter-spacing: 1px;
}
  "))),

  # ── Tab 1: Landing page ─────────────────────────────────────────────────
  tabPanel("Home",
           fluidPage(
             div(style = "max-width: 800px; margin: 40px auto; line-height: 1.6;",
                 h2("Welcome to the Push-forward Power Calculator"),
                 p("This tool implements a Bayesian simulation method to estimate statistical power for cluster-randomized trials
          under various sample size scenarios, accounting for attrition."),
                 h3("How to use this app"),
                 h4("From a previous study or historical data, you will need to enter the following information:"),
                 p("1. A range for the total number of clusters in your study. The app will create a sequence from min to max by 2's."),
                 p("2. A range for the cluster size (number of individuals in each cluster). The app will create a sequence from min to max by 5's."),
                 p("3. The outcome in the control arm as a proportion (i.e. adherence rate = 0.60)"),
                 p("4. The raw effect size for the binary outcome (i.e. if you expect the adherence rate to increase from 0.6 to 0.7, you would enter 0.1 here)"),
                 p("5. The standard deviation for cluster-level random effects"),
                 p("6. The follow-up rate"),
                 p("7. The intraclass correlation coefficient for missingness (i.e. how correlated are two individuals in the same cluster, in terms of attrition)"),
                 h3("Getting started"),
                 p("Head to the ", strong("Calculator"), " tab to begin.")
             )
           )
  ),
  
  # ── Tab 2: Your existing app ─────────────────────────────────────────────
  tabPanel("Calculator",
           fluidRow(

    # ── Left sidebar ─────────────────────────────────────────────────────────
    column(4,

      # Section 1 — Sample size
      div(class = "section-box",
        h4(strong("1) Enter sample size information")),

        div(class = "inline-inputs",
          tags$label("Min number of clusters"),
          numericInput("n_clusters_lo", NULL, value = "", width = "70px")
        ),
        div(class = "inline-inputs",
          tags$label("Max number of clusters"),
          numericInput("n_clusters_hi", NULL, value = "", width = "70px")
        ),
        div(class = "inline-inputs",
            tags$label("Min cluster size"),
            numericInput("cluster_size_lo", NULL, value = "", width = "70px")
        ),
        div(class = "inline-inputs",
            tags$label("Max cluster size"),
            numericInput("cluster_size_hi", NULL, value = "", width = "70px")
        )
      ),

      # Section 2 — Previous study info
      div(class = "section-box",
        h4(strong("2) Enter information from previous study")),

        div(class = "subsection-label", tags$u("Outcome")),
        div(class = "inline-inputs",
          tags$label("Proportion in control arm"),
          numericInput("p_control", NULL, value = NA, min = 0, max = 1, step = 0.01, width = "80px")
        ),
        div(class = "inline-inputs",
          tags$label("Effect size"),
          numericInput("delta", NULL, value = NA, min = 0, max = 1, step = 0.01, width = "80px")
        ),
        div(class = "inline-inputs",
            tags$label("Standard deviation of cluster-level random effects"),
            numericInput("sd_village", NULL, value = NA, min = 0, step = 0.01, width = "80px")
        ),

        div(class = "subsection-label", tags$u("Attrition")),
        div(class = "inline-inputs",
          tags$label("Follow-up rate"),
          numericInput("followup_rate", NULL, value = NA, min = 0, max = 1, step = 0.01, width = "80px")
        ),
        div(class = "inline-inputs",
          tags$label("Intraclass correlation"),
          numericInput("icc_attrition", NULL, value = NA, min = 0, max = 1, step = 0.01, width = "80px")
        ),
      ),

      # Run button
      actionButton("run_btn", "Calculate power")
    ),

    # ── Main panel — plot ─────────────────────────────────────────────────────
    column(8,
      plotOutput("power_plot", height = "520px")
    )
  )
)
)

# ── Server ───────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  output$cluster_size_warning <- renderUI({
    val <- as.integer(input$cluster_size_lo)
    if (!is.na(val) && val %% 2 != 0) {
      tags$small("⚠️ Min cluster size must be an even number",
                 style = "color: red; font-size: 11px; display: block;")
    }
  })
  
  output$cluster_size_warning2 <- renderUI({
    val <- as.integer(input$cluster_size_lo)
    if (!is.na(val) && val %% 2 != 0) {
      tags$small("⚠️ Max cluster size must be an even number",
                 style = "color: red; font-size: 11px; display: block;")
    }
  })
  
  output$run_time_warning <- renderUI({
    val <- nrow(ssGrid)
    if (val > 5) {
      tags$small("This may take a few minutes...",
                 style = "color: red; font-size: 11px; display: block;")
    }
  })
  
  # Reactive: only fires when Run is clicked
  power_values <- eventReactive(input$run_btn, {

    n_clusters_min   <- input$n_clusters_lo
    n_clusters_max  <- input$n_clusters_hi 
    cluster_size_min <- input$cluster_size_lo
    cluster_size_max <- input$cluster_size_hi

    # Create grid of sample size combinations to run
    n_clusters = seq(n_clusters_min, n_clusters_max, 2)
    cluster_size = seq(cluster_size_min, cluster_size_max, 5)
    ssGrid <- expand.grid(n_clusters = n_clusters, cluster_size = cluster_size)
    
    results_list <- list()
    
    withProgress(message = "Running, this may take a few minutes...",
                 value = 0, {
                   
    for(i in 1:nrow(ssGrid)){
                     
    incProgress(1/nrow(ssGrid), detail = paste("Scenario", i, "of", nrow(ssGrid)))                
    # Run the pipeline
      
    parms <- list("n_village" = ssGrid$n_clusters[i], # number of clusters
                  "n_ind" = ssGrid$cluster_size[i], # common cluster size
                  "n_arm" = 2, # number of arms,
                  "alpha" = 0.05,
                  "pii" = input$followup_rate, # follow-up rate
                  "tau" = input$icc_attrition, # attrition ICC
                  "delta" = input$delta, # treatment effect
                  "sd_village" = input$sd_village, # common standard deviation
                  "p1" = input$p_control # prop in control arm
                  )
    
    df    <- simulate_data(parms)
    
    estimates  <- run_model(df)
    
    power  <- run_power(parms, estimates)
    
    results_list[[i]] <- data.frame(
      n     = ssGrid$n_clusters[i],
      m     = ssGrid$cluster_size[i],
      power = power
    )
    
    }
                 })
    power.df <- do.call(rbind, results_list)
    
  })

  # Plot
  output$power_plot <- renderPlot({

    req(power_values())
    power.df <- power_values()

    ggplot(power.df, aes(x = power, colour = as.factor(m)
                         ,linetype = as.factor(m), fill = as.factor(m)
    )) +
      facet_wrap(~ n, labeller = labeller(
        n = ~ paste("n =", .),
        .multi_line = FALSE
      )) +
      geom_density(mapping = aes(y = after_stat(scaled)), lwd = 0.5, alpha = 0.2,
                   key_glyph = "path") +
      labs(title = "Push-forward Power Distributions",
           x = "Power",
           y = "Density",
           color = "Cluster Size",
           fill = "Cluster Size",
           linetype = "Cluster Size") +
      xlim(0, 1) +
      scale_color_manual(values = viridis(7)) +
      scale_linetype_manual(values = c("twodash",  "dashed", "longdash", "solid")) +
      scale_fill_manual(values = viridis(6)) +
      guides(color = guide_legend(override.aes = list(fill = NA, alpha = NA))) +
      theme(legend.position = "right", axis.text.x = element_text(angle = 90, hjust=1))
    
  })
}

# ── Launch ───────────────────────────────────────────────────────────────────

shinyApp(ui = ui, server = server)
