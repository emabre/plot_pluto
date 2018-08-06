data <- read.table("energy_cons.dat")
# data <- c(1, 2, 3)

# set device driver for interactive plotting, I hope!
X11()
# plot default
plot(data)
