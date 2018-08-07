#!/usr/local/bin/Rint

# See: https://thesquareplanet.com/blog/interactive-r-scripts/
library(grDevices)
library(utils)
library(stats)
library(graphics)
library(datasets)
library(methods)
library(base)

data <- read.table("/home/konrad/simulazioni/sims_pluto/disch_outcap/out/energy_cons.dat")
# data <- c(1, 2, 3)

# set device driver for interactive plotting, I hope!
X11()

# plot
par(mfrow=c(3,1)) # all plots on one page

plot(data$t, data$Etot, type="o")
plot(data$t, data$E_tc_in, type="o")
plot(data$t, data$current, type="o")