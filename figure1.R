library(backports)

sirs_gen <- odin.dust::odin_dust("sirs.R")
sirs_times <- rep(NA_real_, 16)
for (i in seq(1, 16)) {
  sirs <- sirs_gen$new(data=list(I_ini=1), n_particles = 1000L, n_threads = as.integer(i), step = 0)
  sirs_times[i] <- bench::system_time(sirs$run(500000))['real']
}

volatility_gen <- dust::dust_example("volatility")
volatility_times <- rep(NA_real_, 16)
for (i in seq(1, 16)) {
  volatility <- volatility_gen$new(data=list(), n_particles = 1000L, n_threads = as.integer(i), step = 0)
  volatility_times[i] <- bench::system_time(volatility$run(500000))['real']
}

volatility_small_times <- rep(NA_real_, 16)
for (i in seq(1, 16)) {
  volatility <- volatility_gen$new(data=list(), n_particles = 100L, n_threads = as.integer(i), step = 0)
  volatility_small_times[i] <- bench::system_time(volatility$run(500000))['real']
}

sirs_gen <- dust::dust_example("sir")
sir_times <- rep(NA_real_, 16)
for (i in seq(1, 16)) {
  sir <- sir_gen$new(data=list(), n_particles = 1000L, n_threads = as.integer(i), step = 0)
  sir_times[i] <- bench::system_time(sir$run(1000))['real']
}

# sircovid vignette setup
# Read the data from the CSV
death_data <- read.table("ons_deaths.csv", 
                         sep=",", 
                         header = TRUE, 
                         row.names = NULL, 
                         stringsAsFactors = FALSE)
serology <- read.table("ons_serology.csv", 
                       sep=",", 
                       header = TRUE, 
                       row.names = NULL, 
                       stringsAsFactors = FALSE)
# Combine with the death data
i <- match(death_data$date, serology$date)
serology <- serology[i, ]
rownames(serology) <- NULL
data <- cbind(death_data, serology[setdiff(names(serology), "date")])

# Clean data
data[,'date'] <- as.Date(data[,'date'], format = "%d/%m/%y")

# Cut out outside date range
data <- data[data$date >= as.Date("2020-03-14"), ]

# Prepare for sircovid. Here only analyse up until restrictions ease
first_wave_data <- data[data$date < as.Date("2020-05-11"), ]

colnames(first_wave_data) <- c("date", "deaths", "npos_15_64", "ntot_15_64")
missing_cols <- c("icu", "hosp", "pillar2_cases", "pillar2_pos", "general", "deaths_comm", "deaths_hosp", "admitted", "new", "new_admitted", 'pillar2_tot', 'pillar2_over25_pos', 'pillar2_over25_tot', 'pillar2_over25_cases', 'react_pos', 'react_tot')
na_col <- as.data.frame(matrix(NA_integer_, 
                               nrow = nrow(first_wave_data),
                               ncol = length(missing_cols)), 
                        row.names = NULL)
colnames(na_col) <- missing_cols
first_wave_data <- cbind(first_wave_data, na_col)

steps_per_day <- 4L
pf_data <- sircovid::sircovid_data(first_wave_data, 
                                   start_date = "2020-01-01", 
                                   dt = 1 / steps_per_day)

n_chains <- 1L
n_steps <- 1e2
n_particles <- 96L

sircovid_times <- rep(NA_real_, 16)
for (i in seq(1, 16)) {
  sircovid_pf <- sircovid::carehomes_particle_filter(pf_data,
                                                     n_particles = n_particles,
                                                     n_threads = i,
                                                     seed = 1L)
  sircovid_times[i] <- bench::system_time(mcstate::pmcmc(mcmc_params,
                                                sircovid_pf,
                                                n_steps = n_steps,
                                                n_chains = n_chains,
                                                save_trajectories = FALSE,
                                                progress = TRUE))['real']
}

# Pre-computed times using the above
sirs_times <- c(83.935750246048,
          46.0407934188843,
          29.9706220626831,
          22.4376611709595,
          18.6930551528931,
          15.6459732055664,
          13.8609700202942,
          11.7859258651733,
          11.0854496955872,
          9.31847095489502,
          8.87862968444824,
          8.00534129142761,
          7.32798647880554,
          7.63635516166687,
          7.24844360351562,
          6.45346426963806)

volatility_times <- c(25.4822378158569,
                13.1624171733856,
                8.69704103469849,
                6.59937334060669,
                5.30535507202148,
                4.55759644508362,
                3.85768342018127,
                3.37826156616211,
                2.98937320709229,
                2.70262622833252,
                2.47994637489319,
                2.4098162651062,
                2.15113806724548,
                2.0420835018158,
                1.93813991546631,
                1.80801820755005)

volatility_small_times <- c(2.52029323577881,
                      1.3471953868866,
                      0.918437719345093,
                      0.654603481292725,
                      0.59818172454834,
                      0.486013650894165,
                      0.40122652053833,
                      0.404340982437134,
                      0.350444316864014,
                      0.293945074081421,
                      0.291177272796631,
                      0.272449254989624,
                      0.254233598709106,
                      0.245222568511963,
                      0.259147644042969,
                      0.228818655014038)

sir_times <- c(0.12324047088623,
         0.0791833400726318,
         0.0565016269683838,
         0.0436475276947021,
         0.0356009006500244,
         0.0268158912658691,
         0.0276618003845215,
         0.0255947113037109,
         0.0245544910430908,
         0.0219776630401611,
         0.0205354690551758,
         0.019913911819458,
         0.0193817615509033,
         0.0198855400085449,
         0.0186755657196045,
         0.0179836750030518)

sircovid_times <- c(199.46611,
              105.35823,
              69.58677,
              53.92774,
              46.00638,
              39.41239,
              34.47344,
              29.70441,
              28.17281,
              26.23955,
              24.16216,
              21.34760,
              21.62700,
              19.90814,
              19.61710,
              17.48254)
 

# Set y data
n_exp <- 4
times = matrix(c(sirs_times, sir_times, sircovid_times, volatility_times), 
                 nrow = n_exp, byrow = TRUE)
for (i in 1:n_exp) {
  times[i,] <- max(times[i,])/times[i,]
}

# Set x data
cores <- seq(1,16)

# Set line colours
cols <- viridis::viridis(n_exp)
cols <- c("#E69F00", "#E69F00", "#56B4E9", "#009E73")
lty = c(1,3,1,1,2)

# Plot lines on the same figure
pdf("figure1.pdf", width = 7, height = 5)
par(mar = c(4.1, 5.1, 0.5, 0.5))
matplot(cores, t(times), type = "l",
         xlab = "Cores", ylab = "Speedup",
        col = cols, lty = lty, lwd = 2, ylim = range(times),
        log = 'xy')
matlines(cores,seq(1,16), col = "#999999", lty = 2)
legend("topleft", lwd = 1, col = c(cols, "#999999"), lty = lty,
       legend = c("SIRS", "SIR (short)", "SIRCOVID", "Volatility","Maximum"), bty = "n")
dev.off()
