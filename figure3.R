set.seed(1)

# Set parameters
alpha <- 0.91
sigma <- 1
n_particles <- 10L
n_threads <- 1L
seed <- 1L
end_step <- 200
steps <- seq(0, end_step, by = 1)

# Load model from dust package
volatility <- dust::dust_example("volatility")
mod <- volatility$new(list(alpha = alpha, sigma = sigma), step = 0,
                      n_particles = n_particles, n_threads = n_threads, seed = seed)
mod$set_state(matrix(rnorm(n_particles, 0, 1), 1, n_particles))
res <- dust::dust_iterate(mod, steps)

# Plot all particle trajectories in gray at 80% opacity
pdf("figure3.pdf", width = 7, height = 5)
par(mar = c(4.1, 5.1, 0.5, 0.5))
matplot(steps, t(res[1, ,]), type = "l", xlab = "t", ylab = "x",
        col = rgb(red = 0.8, green = 0.8, blue = 0.8, alpha = 0.8),
        lty = 1, ylim = range(res))
# Overplot the second trajectory in black
matlines(steps, res[1,2,], col = "#000000", lty = 1)
# Add a loess regression line in red
matlines(steps, predict(loess(apply(res[1,,], MARGIN = 2, mean) ~ steps)), col = "#ff0000", lty = 1)
dev.off()

