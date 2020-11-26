
set.seed(1)
seed <- 1L
n_threads <- 1L

# 
# This section: set up data and run model as needed
#

# Set up model (see mcstate 'sir_models' vignette and mcstate's scripts/generate_vignette_data.R)

# We use the SIR model in the dust package
gen_sir <- dust::dust_example("sir")
dt <- 0.25
inv_dt <- 1/dt
end_step <- 100
sir <- gen_sir$new(data = list(dt = dt, S_ini = 1000, I_ini = 10, beta = 0.2, gamma = 0.1),
                   step = 0,
                   n_particles = 10L,
                   n_threads = n_threads,
                   seed = seed)
steps <- seq(0, end_step/dt, by = 1)
dust_simulate <- dust::dust_iterate(sir, steps)


# Set up data and particle filter as in mcstate vignette
# Read in real data which was fitted to
incidence <- read.table("sir_incidence_data.csv", header = TRUE, sep = ",")
true_history <- readRDS("sir_true_history.rds")

sir_data <- mcstate::particle_filter_data(data = incidence,
                                          time = "day",
                                          rate = 1/dt)

case_compare <- function(state, prev_state, observed, pars = NULL) {
  exp_noise <- 1e6
  incidence_modelled <-
    prev_state[1, , drop = TRUE] - state[1, , drop = TRUE]
  incidence_observed <- observed$cases
  lambda <- incidence_modelled +
    rexp(n = length(incidence_modelled), rate = exp_noise)
  dpois(x = incidence_observed, lambda = lambda, log = TRUE)
}

# Single run of the particle filter
n_particles <- 100
filter <- mcstate::particle_filter$new(data = sir_data,
                                       model = gen_sir,
                                       n_particles = n_particles,
                                       compare = case_compare,
                                       seed = 1L)
filter$run(save_history = TRUE, pars = list(dt = dt))

# MCMC run
beta <- mcstate::pmcmc_parameter("beta", 0.2, min = 0)
gamma <- mcstate::pmcmc_parameter("gamma", 0.1, min = 0)

# Tuned to get good mixing (see vignette for details)
proposal_matrix <- matrix(c(0.0002203173, 0.0001507378, 0.0001507378, 0.0001265065),
                          nrow = 2, ncol = 2, byrow = TRUE)
mcmc_pars <- mcstate::pmcmc_parameters$new(list(beta = beta, gamma = gamma), proposal_matrix)

n_steps <- 2000
n_chains <- 1
pmcmc_run <-
  mcstate::pmcmc(
    mcmc_pars,
    filter,
    n_steps,
    save_state = TRUE,
    save_trajectories = TRUE,
    progress = TRUE,
    n_chains = n_chains
  )
processed_chains <- mcstate::pmcmc_thin(pmcmc_run, burnin = 500, thin = 2)

n_steps <- 5000
# Multiple chains
pmcmc_1 <-
  mcstate::pmcmc(
    mcmc_pars,
    filter,
    n_steps,
    save_state = TRUE,
    save_trajectories = TRUE,
    progress = TRUE,
    n_chains = n_chains
  )
pmcmc_2 <-
  mcstate::pmcmc(
    mcmc_pars,
    filter,
    n_steps,
    save_state = TRUE,
    save_trajectories = TRUE,
    progress = TRUE,
    n_chains = n_chains
  )
pmcmc_3 <-
  mcstate::pmcmc(
    mcmc_pars,
    filter,
    n_steps,
    save_state = TRUE,
    save_trajectories = TRUE,
    progress = TRUE,
    n_chains = n_chains
  )
pmcmc_4 <-
  mcstate::pmcmc(
    mcmc_pars,
    filter,
    n_steps,
    save_state = TRUE,
    save_trajectories = TRUE,
    progress = TRUE,
    n_chains = n_chains
  )

# Forecast from MCMC run
forecast <- predict(processed_chains,
                    steps = seq(400, 800, 4),
                    prepend_trajectories = TRUE,
                    seed = processed_chains$predict$seed)
mini_forecast <- forecast$state[,sample.int(ncol(forecast$state), size = 10),]

# 
# This section: make plots from data above
#

# Generic function to plot particle trajectories from the SIR model
cols <- c(S = "#999966", I = "#8c8cd9", R = "#cc0044")
plot_trajectories <- function(history,
                              times,
                              alpha = 1) {
  par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
  col_alpha <- c(S = rgb(0x99/256, 0x99/256, 0x66/256, alpha = alpha), 
                 I = rgb(0x8c/256, 0x8c/256, 0xd9/256, alpha = alpha), 
                 R = rgb(0xcc/256, 0x00/256, 0x44/256, alpha = alpha))
  matplot(times, t(history[1, , -1]), type = "l",
          xlab = "Time", ylab = "Number of individuals",
          col = col_alpha[["S"]], lty = 1, ylim = range(history))
  matlines(times, t(history[2, , -1]), col = col_alpha[["I"]], lty = 1)
  matlines(times, t(history[3, , -1]), col = col_alpha[["R"]], lty = 1)
  legend("left", lwd = 1, col = cols, legend = names(cols), bty = "n")
}


# Figure Panel A: simulation of model
pdf("figure4a.pdf", width = 7, height = 5)
plot_trajectories(dust_simulate, steps[-1] * dt)
dev.off()

# Figure 4 Panel B: particle filter

pdf("figure4b.pdf", width = 7, height = 5)
plot_trajectories(filter$history(), seq_len(100), alpha = 0.05)
# Add observations
matpoints(seq_len(100), t(true_history[1:3, , 2:101]), pch = 21,
          bg = cols, col = "#000000", lwd = 1)
dev.off()

# Figure 4 Panel C: incidence curves

# Extract the cumulative incidence partition and convert to point incidence
incidence_modelled <- apply(filter$history()[4,,], 1, diff)
avg <- apply(incidence_modelled, 1, mean)

# Plot particles in grey, observations in black, particle avg as points
pdf("figure4c.pdf", width = 7, height = 5)
matplot(seq_len(100), incidence_modelled, type = "l",
        xlab = "Time", ylab = "Cases",
        col = rgb(0x99/256, 0x99/256, 0x99/256, alpha = 0.05), lty = 2,
        ylim = range(incidence_modelled))
matlines(seq_len(100), incidence$cases, lty = 1,
          col = "#000000", lwd = 1)
matpoints(seq_len(100), avg, pch = 20,
         col = "#000000")
legend("topright", lwd = 1, lty = c(2,1), col = c("#999999", "#000000"),
       legend = c("Modelled", "Observed"), bty = "n")
dev.off()

# Figure 4 Panel D: forecasts

pdf("figure4d.pdf", width = 7, height = 5)
plot_trajectories(mini_forecast[,,1:101], seq_len(100))
# Add nowcast observations
matpoints(seq_len(50), t(true_history[1:3, , 2:51]), pch = 21,
          bg = cols, col = "#000000", lwd = 1)
# Add forecast observations
matpoints(seq(51, 100), t(true_history[1:3, , 52:101]), pch = 18, cex = 1,
          lwd = 1, col = "#000000")
dev.off()


# Figure 5
cols <- rev(viridis::viridis(4))

# Combine chains in columns for plotting
beta <- cbind(pmcmc_1$pars[,'beta'], pmcmc_2$pars[,'beta'], pmcmc_3$pars[,'beta'], pmcmc_4$pars[,'beta'])
gamma <- cbind(pmcmc_1$pars[,'gamma'], pmcmc_2$pars[,'gamma'], pmcmc_3$pars[,'gamma'], pmcmc_4$pars[,'gamma'])

pdf("figure5.pdf", width = 14, height = 8)
par(mar = c(3, 3, 2, 1), 
    mgp = c(2, 0.5, 0), oma = c(1, 1, 1, 1),
    fig=c(0,1,0.67,1))
matplot(beta, type = "l", lty = 1, 
        xlab = "Iteration", 
        ylab = "beta", col = cols)

par(fig=c(0,1,0.33,0.67), new = TRUE)
matplot(gamma, type = "l", lty = 1, 
        xlab = "Iteration", 
        ylab = "gamma", col = cols)

par(fig=c(0.4,0.9,0,0.33), new = TRUE)
hist(c(beta)/c(gamma), main = "", breaks = 14, xlab = "R0")

par(fig=c(0,0.3,0,0.33), xpd=TRUE, new = TRUE)
plot.new()
legend("left", cex = 0.75, fill = cols, bty = "n", legend = paste("chain", seq_len(4)))
dev.off()
