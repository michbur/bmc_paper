#setwd("/home/tux/Work/paper/chipPCR/bmc_article/figures/")
setwd("/tmp")
require(chipPCR)
require(qpcR)

# Figures for chipPCR, The R Journal

###########################################################################
# Figure 2
#################
# Load the shiny package (chipPCR should already be loaded)
# Run from a R console following commands
require(shiny)

#Invoke the shiny AmpSim app in the default browser
runApp(paste(find.package("chipPCR")[1],"/AmpSim.gui", sep = ""))

###########################################################################

###########################################################################
# Figure 1
#################
pdf(file = "problems.pdf", width = 12, height = 8)
# Use AmpSim to generate an amplification curve with 40 cycles
# and a different Cq.
res.pos <- AmpSim(cyc = 1:40, noise = TRUE, b.eff = -12, nnl = 0.02)
res.pos[5, 2] <- res.pos[5, 2] * 6

res.low <- AmpSim(cyc = 1:40, noise = TRUE, b.eff = -20, bl = 0.5, 
		  ampl = 0.58, Cq = 33)
# Add missing value to res.low at cycle 31
res.low[31, 2] <- NA

res.neg <- AmpSim(cyc = 1:40, b.eff = -0.1, bl = 0.05, ampl = 0.4, Cq = 1, 
		  noise = FALSE, nnl = 0.5)
		      
res.pos.CPP <- cbind(1:40, CPP(res.pos[, 1], res.pos[, 2], 
		     bg.outliers = TRUE, smoother = TRUE, method = "smooth", 
		      method.norm = "minmax", method.reg = "lmrob")$y)
		      
res.low.NA <- cbind(1:40, CPP(res.low[, 1], res.low[, 2], smoother = TRUE, 
		    method = "smooth", bg.outliers = TRUE, method.norm = "minmax", 
		    method.reg = "lmrob")$y)
		      
res.neg.exc <- cbind(1:40, amptester(res.neg[, 2]))

par(mfrow = c(1,2), las = 0, bty = "n", cex.axis = 1.5, cex.lab = 1.5, 
    font = 2, cex.main = 1.8, oma = c(1,1,1,1))
plot(NA, NA, xlim = c(1,40), ylim = c(0, max(res.pos[, 2])), xlab = "Cycle", 
     ylab = "Raw fluorescence")
mtext("A", cex = 2, side = 3, adj = 0, font = 2)

lines(res.pos, lwd = 2)
lines(res.low, col = 2, lwd = 2)
arrows(38, min(res.low[, 2], na.rm = TRUE), 38, max(res.low[, 2], 
       na.rm = TRUE), code=3, lwd=3, angle=90, col="grey")
text(38, max(res.low[, 2], na.rm = TRUE) * 0.7,"SNR", cex=1.2)

arrows(29,0.42,31,0.51, lwd=2)
text(29,0.38, "NA", cex=1.2)

points(res.pos[5, 1], res.pos[5, 2], pch=21, cex=4, lwd=5, col="orange")
text(res.pos[5, 1], res.pos[5, 2] * 1.2, "Outlier", cex=1.2)

lines(res.neg, col = 4, lwd = 2)
text(20, mean(res.neg[, 2]) *0.9, "No amplification", cex=1.2, col = 
"blue")


plot(NA, NA, xlim = c(1,40), ylim = c(0, max(res.pos[, 2])), xlab = "Cycle", 
     ylab = "Pre-processed fluorescence")
abline(h = 0.03, lty = 2, lwd = 2)
mtext("B", cex = 2, side = 3, adj = 0, font = 2)

lines(res.pos.CPP, lwd = 2)
lines(res.low.NA, col = 2, lwd = 2)
lines(res.neg.exc, col = 4, lwd = 2)

legend(1, 1, c("Positive (outlier removed)", "Positive (scaled)", 
       "Negative", "Threshold line nof Cq"), 
       col = c("black", "red", "blue", "black"), lty = c(1,1,1,2), 
       lwd = 2, bty = "n")
lines(c(15.1,15.1), c(-1,0.03), lwd = 2, col = "black")
text(14, 0.06, "Cq")
lines(c(28.5,28.5), c(-1,0.03), lwd = 2, col = "red")
text(27, 0.06, "Cq", col = "red")


dev.off()

# Figure 3
#################
pdf(file = "AmpSim.pdf", width = 12, height = 6)
# Draw an empty plot for 40 cycles with user defined parameters.
par(las = 0, bty = "n", cex.axis = 1.2, cex.lab = 1.2, font = 2, 
    cex.main = 1.1, oma = c(1,1,1,1))
plot(NA, NA, xlim = c(1,40), ylim = c(0,1.1), xlab = "Cycle", ylab = "RFU")
colors <- rainbow(8)
# Create eight amplification curves. The approximate Cqs are synthesized 
# as temporary Cqs by adding a random value to a starting Cq of 25. Note: 
# ``noise'' is set TRUE with a level of nnl = 0.03. This adds some scatter 
# to the amplification curves.

sapply(1L:8, function(i) {
  Cq.tmp <- 25 + rnorm(1) * 5
  
  tmp <- AmpSim(1:40, Cq = Cq.tmp, noise = TRUE, nnl = 0.03)
  lines(tmp, col = colors[i], lwd = 2)
  
  # Add the approximate Cq values to the plot
  text(3, 1 - i / 10, paste("Cq ", round(Cq.tmp, 2)), col = colors[i])
})

###########################################################################
# Figure 7
#################
# Load MBmca package (v. 0.0.3-3 or later)
require(MBmca)

# Create an graphic device for two empty plots.
par(mfrow = c(1,2))
plot(NA, NA, xlim = c(1,45), ylim = c(0.01,1.1), xlab = "Cycles", 
     ylab = "Fluorescence", main = "")
mtext("A", cex = 2, side = 3, adj = 0, font = 2)

# Create a sequence of "targeted" Cq values (Cq.t) between 15 and 34 cycles.

Cq.t <- rep(seq(15, 34, 3.5), 3)

# In-silico experiment set up: Define the levels for the decadic dilutions
# with concentrations from 100 to 0.001 (six steps) as three replicates.

dilution <-rep(sapply((2:-4), function(i) {10^i}), 3)

# Create an empty matrix for the results of the concentration
# dependent Cq values.

ma.out <- matrix(data = NA, nrow = 45, ncol = length(Cq.t))

# Use AmpSim to simulate amplification curves at different concentrations. 
# The simulation is performed with the addition of some noise. This generates 
# unique (non-reproducebale) amplification curves, even under idential paramter 
# settings.

Cq.out <- vector()
# Simulate a qPCR reaction with AmpSim for 45 cycles and some noise.

for (i in 1L:18) {
      ma.out[1:45, i] <- AmpSim(cyc = c(1:45), b.eff = -50, bl = 0.001, 
				ampl = 1, Cq = Cq.t[i], noise = TRUE, 
				nnl = 0.02)[, 2]
      lines(1:45, ma.out[, i])
      tmpP <- mcaSmoother(1:45, ma.out[, i])
# Calculate the pseudo Second Derivative Maximum (SDM) (Cq) using 
# the diffQ2 function from the MBmca package.
      Cq.tmp <- diffQ2(tmpP, inder = TRUE, warn = FALSE)$xTm1.2.D2[1]
      abline(v = Cq.tmp)
      Cq.out <- c(Cq.out, Cq.tmp)
}

# Assign the calculated Cqs to the corresponding concentrations.
tmp <- data.frame(dilution[1:6], Cq.out[1:6], Cq.out[7:12],  Cq.out[13:18])
		  
# Determine the amplification efficiency by using the effcalc function.
effcalc(tmp[, 1], tmp[, 2:4], CI = TRUE)
mtext("B", cex = 2, side = 3, adj = 0, font = 2) 

dev.off()

###########################################################################
# Figure 4
#################
pdf(file = "SDM.pdf", width = 12, height = 8)
# Use AmpSim to generate an amplification curve with 40 cycles
# and an approximate Cq of 20 and assign it to the object isPCR.
# isPCR is an object of the class "data.frame".
isPCR <- AmpSim(cyc = 1:40, Cq = 20)

# Invoke the inder function for the object isPCR to interpolate 
# the derivatives of the simulated data as object res. The Nip 
# parameter was set to 5. This leads to smoother curves. res is
# an object of the class "der".
res <- inder(isPCR, Nip = 5)

# Plot the the object res and add descriptions to the elements.

par(las = 0, bty = "n", cex.axis = 1.5, cex.lab = 1.5, 
    font = 2, cex.main = 1.8, oma = c(1,1,1,1))

plot(isPCR, xlab = "Cycle", ylab = "RFU", ylim = c(-0.15,1),
     main = "", type = "b", pch = 20, lwd = 2)
colors <- rainbow(4)
# Add graphical elements for the dervatives and the calculated
# Cq values FDM, SDM, SDm and SDC.

  lines(res[, "x"], res[, "d1y"], col = "blue", lwd = 2)
  lines(res[, "x"], res[, "d2y"], col = "red", lwd = 2)
  
# Fetch the Cq values from res with the summary function
  summ <- summary(res, print = FALSE)
  
  abline(v = summ, col = colors, lwd = 2)
  text(15, 0.3, paste("FDM ~ ", round(summ["FDM"], 2)), 
       cex = 1.5, col = colors[1])
  text(15, 0.2, paste("SDM ~ ", round(summ["SDM"], 2)), 
       cex = 1.5, col = colors[2])
  text(15, - 0.1, paste("SDm ~ ", round(summ["SDm"], 2)), 
       cex = 1.5, col = colors[3])
  text(15, 0.7, paste("SDC ~ ", round(summ["SDC"], 2)), 
       cex = 1.5, col = colors[4])
       
  legend(1.1, 0.9, c("raw", "first derivative", "second derivative"), 
         col = c(1,4,2), lty = c(2,1,1), cex = 1.2, bty = "n")

# Summary of the object res.
summ
#     FDM      SDM      SDm      SDC 
#19.81407 19.03015 20.98995 19.98604
dev.off()

###########################################################################
pdf(file = "fixNA.pdf", width = 12, height = 12)

# Simulation of an ideal amplification curve with 40 cycles
# The other paramter of the AmpSim function are identical to
# the default.

res <- AmpSim(cyc = 1:40)

# Introduce a missing value (cycle 18) in the transition between 
# the background and the exponential phase.
res.NA <- res
res.NA[18, 2] <- NA

# Helper function to highlight the position of the missing value.
abliner <- function(x1 = 17.5, x2 = 18.5, y1 = 0.09, y2 = 0.14) {
		abline(v = c(x1, x2), col = "red")
		abline(h = c(y1, y2), col = "red")
	    }

par(las = 0, mfrow = c(2,2), bty = "n", cex.axis = 1.5, cex.lab = 1.5, 
    font = 2, cex.main = 1.8, oma = c(1,1,1,1))
plot(res, xlab = "Cycles", ylab = "RFI", type = "b", pch = 20, 
     main = "Without missing value")
abliner()
mtext("A", cex = 2, side = 3, adj = 0, font = 2)
res.NA.linear <- fixNA(res.NA[, 1], res.NA[, 2], spline = FALSE, 
		       verbose = FALSE)

plot(res.NA, xlab = "Cycles", ylab = "RFI", type = "b", pch = 20, 
     main = "With missing\n value during transition")
abliner()
mtext("B", cex = 2, side = 3, adj = 0, font = 2)

res.NA.spline <- fixNA(res.NA[, 1], res.NA[, 2], spline = TRUE, 
		       verbose = FALSE)
	       
plot(res.NA.linear, xlab = "Cycles", ylab = "RFI", type = "b", 
     pch = 20, main = "Linear imputed\n NA value during transition")
abliner()
mtext("C", cex = 2, side = 3, adj = 0, font = 2)
     
plot(res.NA.spline, xlab = "Cycles", ylab = "RFI", type = "b", 
      pch = 20, main = "Spline imputed\n NA value during transition")
abliner()
mtext("D", cex = 2, side = 3, adj = 0, font = 2)
par(mfrow = c(1,1))

dev.off()

###########################################

pdf(file = "bgmax.pdf", width = 12, height = 8)
par(las = 0, mfrow = c(2,1), bty = "n", cex.axis = 1.5, cex.lab = 1.5, 
    font = 2, cex.main = 1.8, oma = c(1,1,1,1))

res <- AmpSim(cyc = 1:40, Cq = 25)
plot(res, xlim = c(1,40), ylim = c(-0.1,1), xlab = "Cycles", 
     ylab = "RFI", 
     main = "Estimation of the Background Range\n in Absence of Noise", 
     type = "b", pch = 20)
background <- bg.max(res[, 1], res[, 2])
mtext("A", cex = 2, side = 3, adj = 0, font = 2)

points(background@.Data[, 3], col = "red", type = "b", pch = 20)
points(background@.Data[, 4], col = "blue", type = "b", pch = 20)
abline(v = background@bg.start)
  text(background@bg.start, 0.2, "Background start", pos = 4)
abline(v = background@bg.stop, col = "blue")
  text(background@bg.stop, 0.25, "Background stop", pos = 4, 
	col = "blue")
abline(v = background@amp.stop, col = "green")
  text(background@amp.stop, 0.3, "Plateau transition", pos = 4, 
	col = "green")
legend(4, 1, c("Raw data", "First derivative", "Second derivative"), 
       pch = rep(20,3), col = c(1,2,4), bty = "n")

res <- AmpSim(cyc = 1:40, Cq = 25, noise = TRUE)
plot(res, xlim = c(1,40), ylim = c(-0.1,1), xlab = "Cycles", 
     ylab = "RFI", 
     main = "Estimation of the Background Range\n in Presence of Noise", 
     type = "b", pch = 20)
mtext("B", cex = 2, side = 3, adj = 0, font = 2)
background <- bg.max(res[, 1], res[, 2])

points(background@.Data[, 3], col = "red", type = "b", pch = 20)
points(background@.Data[, 4], col = "blue", type = "b", pch = 20)
abline(v = background@bg.start)
text(background@bg.start, 0.2, "Background start", pos = 4)
abline(v = background@bg.stop, col = "blue")
text(background@bg.stop, 0.25, "Background stop", pos = 4, col = "blue")
abline(v = background@amp.stop, col = "green")
text(background@amp.stop, 0.3, "Plateau transition", pos = 4, col = 
"green")
legend(4, 1, c("Raw data", "First derivative", "Second derivative"), 
       pch = rep(20,3), col = c(1,2,4), bty = "n")
par(mfrow = c(1,1))
    
dev.off()

###########################################
pdf(file = "inder.pdf", width = 12, height = 6)


# First example
# Plot all data from C127EGHP and calculate the SDM (Second Derivative 
# Maximum) values with the diffQ2() function (Note: the inder parameter
# is set as TRUE)
# first plot the samples detected with EvaGreen and next the samples 
# detected with the Hydrolysis probe

require(MBmca)

pointer <- function (x, pos = 1, w = 5, stat = TRUE){
  xx <- pos + rep(seq(-0.1, 0.1, length.out = w), ceiling(length(x)/w))
  yy <- sort(x)
  points(xx[1:length(yy)], yy, pch = 19)
  
  if (stat == TRUE)
    x.median <- median(x, na.rm = T)
    x.mad <- mad(x, na.rm = T) * 2
    param <- c(length= 0, code = 3, pch = 15, cex = 1.5)
    arrows(xx[1] * 0.98, x.median, tail(xx, 1) * 1.02, 
	    x.median, param, lwd = 3)
    arrows(xx[1] * 1.01, x.median + x.mad, tail(xx, 1) * 0.99, 
	    x.median + x.mad, param, lwd = 2, lty = 2)
    arrows(xx[1] * 1.01, x.median - x.mad, tail(xx, 1) * 0.99, 
	    x.median - x.mad, param, lwd = 2, lty = 2)
}

amp.liner <- function(range, input, colors = "black") {
  sapply(range, function(i) {
	 lines(input[, 2], input[, i], col = colors, pch = 19)
	 tmpP <- mcaSmoother(input[, 2], input[, i])
	 SDM <- diffQ2(tmpP, inder = TRUE, warn = FALSE)[["xTm1.2.D2"]][1]
	 abline(v = SDM)
	 SDM
       }
  )
}


par(mfrow = c(1,3), las = 0, bty = "n", cex.axis = 1.5, 
    cex.lab = 1.5, font = 2, cex.main = 1.8, oma = c(1,1,1,1))
plot(NA, NA, xlim = c(1,40), ylim = c(0,10), xlab = "Cycle", 
      ylab = "Fluorescence", main = "EvaGreen")
mtext("A", cex = 2, side = 3, adj = 0, font = 2)

EG <- amp.liner(range = 3L:34, input = C127EGHP)
  
plot(NA, NA, xlim = c(1,40), ylim = c(0,10), xlab = "Cycle", 
      ylab = "Fluorescence", main = "Hydrolysis probe")
mtext("B", cex = 2, side = 3, adj = 0, font = 2)

HP <- amp.liner(range = 35L:66, input = C127EGHP)

plot(NA, NA, xlim = c(0.8,2.2), ylim = c(13,14), xaxt = "n", 
     xlab = "", ylab = "Cq (SDM, diffQ2)")
text(c(1.05,2), c(13.05,13.05), c("EG", "HP"), cex = 2)
mtext("C", cex = 2, side = 3, adj = 0, font = 2)
pointer(EG, pos = 1, w = 8)
pointer(HP, pos = 2, w = 8)

par(mfrow = c(1,1))
    
#dev.off()

###########################################

#pdf(file = "inder_fit.pdf", width =12, height = 6)

fit.amp <- function(cyc, fluo, plot = FALSE) {

	ampl <- quantile(fluo, 0.999)
	bl <- quantile(fluo, 0.001)
	Cq <- round(mean(cyc))
	b.eff <- 1
	
 	fit <- nls(fluo ~ bl + ampl / (1 + exp(- (cyc - Cq) / b.eff)), 
 		     start = list(Cq = Cq, b.eff = b.eff, ampl = ampl, 
 		     bl = bl)
 		     )
	
	res.pred <- data.frame(cyc, predict(fit))
	res <- inder(res.pred[, 1], res.pred[, 2])
	if (plot) {
	    lines(res[, 1], res[, 4])
	}
	summary(res)[2] # SDM
}

tmp <- C126EG595

out <- unlist(lapply(2L:ncol(tmp), function(x) fit.amp(tmp[, 1], tmp[, x])))

layout(matrix(c(1,2,1,3), 2, 2, byrow = TRUE))

plot(NA, NA, xlim = c(1,40), ylim = c(min(tmp[, 2L:97]), 
     max(tmp[, 2L:97])), xlab = "Cycle", ylab = "Raw fluorescence")
mtext("A", cex = 2, side = 3, adj = 0, font = 2)
lapply(2L:97, function(x) {
    lines(tmp[, 1], tmp[, x], col = ifelse(out[x - 1] < 15.5, 
	  "red", "black"), lwd = 2)
  }
)
abline(v = out)

plot(NA, NA, xlab = "Cycle", ylab = "RFU''(Cycle)", main = "", xlim = c(0,40), 
     ylim = c(-850, 850))
lapply(2L:97, function(x) {
    fit.amp(tmp[, 1], tmp[, x], plot = TRUE)
  }
)
mtext("B", cex = 2, side = 3, adj = 0, font = 2)

hist(out, xlab = "Cq (SDM)", main = "", 
     breaks = seq(14.8, 15.8, 0.05), col = rainbow(96))
abline(v = 15.5, lty = 2)
mtext("C", cex = 2, side = 3, adj = 0, font = 2)
dev.off()

###########################################
#
# Identical to code example in chipPCR
# manual v 0.0.7.
#
###########################################
pdf(file = "amptester.pdf", width = 16, height = 10)

# Arrange graphs in orthogonal matrix and set parameter for the plot.
par(las = 0, bty = "n", cex.axis = 1.5, cex.lab = 1.5, 
    font = 2, cex.main = 1.4, oma = c(1,1,1.5,1))
    
# Use of amptester for time-dependent measurements. Amplification curves 
# from the capillaryPCR data set were processed in a loop. The results of 
# amptester are added to the raw data. 

colors <- rainbow(8)
plot(NA, NA, xlim = c(0,80), ylim = c(0,1300), xlab = "Time [min]", 
     ylab = "Voltage (micro V)", main = "ccPCR")
mtext("B", cex = 2, side = 3, adj = 0, font = 2)
sapply(c(1,3,5,7), function(i) {
    xy.tmp <- cbind(capillaryPCR[1L:750, i], capillaryPCR[1:750, i + 1])
    res.ampt <- amptester(xy.tmp[, 2])
    res.bg <- summary(bg.max(xy.tmp[, 1], xy.tmp[, 2]))
    res.ampt <- ifelse(res.ampt@decisions[2] == TRUE && 
		       res.ampt@decisions[4] == TRUE, "positve", "negative"
		       )

    lines(xy.tmp[, 1], xy.tmp[, 2], type = "b", pch = 20, col = colors[i])
    text(75, max(na.omit(xy.tmp[, 2])), res.ampt, cex = 1.3, col = colors[i])
    
    arrows(res.bg["bg.stop"], i * 10, res.bg["amp.stop"], i * 10, 
	   col = colors[i], code = 3, lwd = 2, length = 0.1)
  }
)

###########################################

# Load the drc package for the five-parameter log-logistic fit.
require(drc)
# Load the testdat dataset from the qpcR package without 
# loading the entire package.
data(testdat, package = "qpcR")

# Arrange graphs in an orthogonal matrix and set the plot parmeters.
par(mfrow = c(4,6))

# Apply the the amptester function to all amplification curve data
# and write the results to the main of the plots.

sapply(2L:ncol(testdat), function(x) {
	res.ampt <- amptester(testdat[, x])
	# Use bg.max to indicate the region of the estimated background range 
	# and the plateau phase.
	res.bg <- summary(bg.max(testdat[, 1], testdat[, x]))
	
	# Make a logical connection by two tests (shap.noisy, lrt.test and 
	# tht.dec) of amptester to decide if an amplification reaction is 
	# positive or negative.
	decision <- ifelse(!res.ampt@decisions[1] &&
			    res.ampt@decisions[2] && 
			    res.ampt@decisions[4], 
			   "positive", "negative")

	plot(testdat[, 1], testdat[, x], xlab = "Cycle", ylab = "RFU", 
	     pch = 19, main = paste0("Decision: ", decision)
	)
	# Use the drm function with a five-parameter log-logistic fit model.
	lines(predict(drm(testdat[, x] ~ testdat[, 1], fct = LL.5())), col = 2)
	
	# Add vertical lines for the estimated start (green) and end (blue) 
	# of the exponential region in case of positive amplification 
	# reactions.
	if (decision == "positive") abline(v = c(res.bg["bg.stop"], 
					     res.bg["amp.stop"]), col = c(3,4))
	# Add sample label.
	text(10, max(testdat[, x]) * 0.95, colnames(testdat)[x], cex = 1.3, 
	     col = 4)
	}
)

#########################################################
data(testdat, package = "qpcR")
x <- testdat[1:45, 1]
y.pos <- testdat[1:45, 2]
y.neg <- testdat[1:45, 4]
nh <- trunc(length(x) * 0.20)
nt <- trunc(length(x) * 0.15)

y.pos.head <- head(y.pos, n = nh)
y.neg.head <- head(y.neg, n = nh)
y.pos.tail <- tail(y.pos, n = nt)
y.neg.tail <- tail(y.neg, n = nt)

lb.pos <- median(y.pos.head) + 2 * mad(y.pos.head)
ub.pos <- median(y.pos.tail) - 2 * mad(y.pos.tail)

lb.neg <- median(y.neg.head) + 2 * mad(y.neg.head)
ub.neg <- median(y.neg.tail) - 2 * mad(y.neg.tail)

res.shapiro.pos <- shapiro.test(y.pos)
res.shapiro.neg <- shapiro.test(y.neg)

res.wt.pos <- wilcox.test(head(y.pos, n = nh), tail(y.pos, n = nt), alternative = "less")
res.wt.neg <- wilcox.test(head(y.neg, n = nh), tail(y.neg, n = nt), alternative = "less")

###
RGt <- function(y) {
   ws <- ceiling((15 * length(y)) / 100)
    if (ws < 5) 
      ws <- 5
    if (ws > 15) 
      ws <- 15
    y.tmp <- na.omit(y[-c(1:5)])
    x <- 1:length(y.tmp)
    suppressWarnings(
      res.reg <- sapply(1L:(length(y.tmp)), function (i)  {
        round(summary(lm(y.tmp[i:c(i + ws)] ~ x[i:c(i + ws)]))[["r.squared"]], 4)
      }
      )
    )
    
    # Binarize R^2 values. Everything larger than 0.8 is positve
    res.LRt <- res.reg
    # Define the limits for the R^2 test
    res.LRt[res.LRt < 0.8] <- 0
    res.LRt[res.LRt >= 0.8] <- 1
    # Seek for a sequence of at least six positve values (R^2 >= 0.8)
    # The first five measurepoitns of the amplification curve are skipped
    # beacuse most technologies and probetechnologies tend to overshot
    # in the start (background) region.
    res.out <- sapply(5L:(length(res.LRt) - 6), function(i) {
	ifelse(sum(res.LRt[i:(i + 4)]) == 5, TRUE, FALSE)
      }
    )
    out <- cbind(1L:(length(y.tmp)), res.reg)
    #res.out
}
###

#par(mfrow = c(2,2))
layout(matrix(c(1,1,2,2,3,4,5,6), 2, 4, byrow = TRUE), respect = TRUE)
plot(x, y.pos, xlim = c(1, 50), ylim = c(-0.1, 15.5), xlab = "Cycle", ylab = "RFU", main = "Positive amplification", type = "b", pch = 19)

mtext("A", cex = 2, side = 3, adj = 0, font = 2)
abline(v = c(nh, length(x) - nt), lty = 3)

abline(h = lb.pos, lty = 2, col = "red")
text(4, 1.5, "Noise\nmedian + 2 * mad", col = "red", cex = 1.3)

abline(h = ub.pos, lty = 2, col = "green")
text(48, 13, "Signal\nmedian - 2 * mad", col = "green", cex = 1.3)

arrows(4.5, 12.5, 42.5, 12.5, length = 0.1, angle = 90, code = 3)
text(25, 14.5, paste("W = ", res.wt.pos$statistic, "\np-value = ", res.wt.pos$p.value))
text(5,5, paste("Fold change: \n", round(ub.pos/lb.pos, 2)))

plot(x, y.neg, xlim = c(1, 50), ylim = c(-0.03, 0.08), xlab = "Cycle", ylab = "RFU", main = "Negative amplification", type = "b", pch = 19)
mtext("B", cex = 2, side = 3, adj = 0, font = 2)
abline(v = c(nh, length(x) - nt), lty = 3)

abline(h = lb.neg, lty = 2, col = "red")
text(4, 0.06, "Noise\nmedian + 2 * mad", col = "red", cex = 1.3)

abline(h = ub.neg, lty = 2, col = "green")
text(48, -0.02, "Signal\nmedian - 2 * mad", col = "green", cex = 1.3)

arrows(4.5, 0.04, 42.5, 0.04, length = 0.1, angle = 90, code = 3)
text(25, 0.06, paste("W = ", res.wt.neg$statistic, "\np-value = ", res.wt.neg$p.value))
text(5, 0.075, paste("Fold change: \n", round(ub.neg/lb.neg, 2)))

qqnorm(y.pos, pch = 19, main = paste("W = ", res.shapiro.pos$statistic, "\np-value = ", res.shapiro.pos$p.value))
mtext("C", cex = 2, side = 3, adj = 0, font = 2)
qqline(y.pos, col = "orange", lwd = 2)

plot(RGt(y.pos), xlab = "Cycle", ylab = "R^2", main = "LRt", pch = 19, type = "b")
abline(h = 0.8, col = "black", lty = 2)

qqnorm(y.neg, pch = 19, main = paste("W = ", res.shapiro.neg$statistic, "\np-value = ", res.shapiro.neg$p.value))
mtext("D", cex = 2, side = 3, adj = 0, font = 2)
qqline(y.neg, col = "orange", lwd = 2)

plot(RGt(y.neg), xlab = "Cycle", ylab = "R^2", main = "LRt", pch = 19, type = "b")
abline(h = 0.8, col = "black", lty = 2)

dev.off()
###########################################

pdf(file = "MFIaggr.pdf", width = 12, height = 6)
# 
par(las = 0, bty = "n", cex.axis = 1.2, cex.lab = 1.2, 
    font = 2, cex.main = 1.2, oma = c(1,1,1,1))

plot(MFIaggr(VIMCFX96_60[, 1], VIMCFX96_60[, 2:ncol(VIMCFX96_60)], 
     llul = c(1,10)), CV = FALSE)

plot(MFIaggr(VIMCFX96_60[, 1], VIMCFX96_60[, 2:ncol(VIMCFX96_60)], 
     llul = c(1,40)), CV = FALSE)

dev.off()

###################################
pdf(file = "smoother.pdf", width = 11, height = 8)

# Simulate an amplification curve as object tmp and plot the data.

tmp <- AmpSim(cyc = c(1:35), bl = 0)

par(las = 0, bty = "n", cex.axis = 1.5, cex.lab = 1.5, 
    font = 2, cex.main = 1.8, oma = c(1,1,1,1), fig = c(0,1,0.55,1))
plot(tmp, type = "b", col = 1, pch = 20, xlab = "", ylab = "RFU", 
      main = "Raw data")
mtext("A", cex = 2, side = 3, adj = 0, font = 2)

# Apply all (parameter method = "all") smoothers/filter with the default 
# setting to the amplification curve of the object tmp. Smoothers / Filters:
# Savitzky-Golay smoothing filter
# locally-weighted polynomial regression
# moving average, windowsize 3
# cubic spline smooth
# standard cubic spline smooth
# Friedman's SuperSmoother
# weighted Whittaker smoothing with first order finite difference penalty
# weighted Whittaker smoothing with a second order finite difference penalty
res <- smoother(tmp[, 1], tmp[, 2], method = "all", CPP = FALSE)

# Calculate the difference between the ideal curve (tmp) and the smoothed curves
# (res) and asssign the results to the object res.out
res.out <- cbind(cycle = tmp[, 1], tmp[, 2] - res)

# Plot the smoothed curves
par(fig = c(0,1,0,0.65), new = TRUE)
plot(NA, NA, type = "b", col = 2, pch = 15, xlim = c(1,35), ylim = c(-0.1,0.1),
     xlab = "Cycle", ylab = "delta refMFI (raw - smoothed)",
     main = "Smoothed / Filtered data")
     
mtext("B", cex = 2, side = 3, adj = 0, font = 2) 
legend(1.5, 0.1, ncol = 2, colnames(res.out[, 2:9]), pch = 15:22, 
       lwd = 2, col = c(2:9))

# Plot the results.
sapply(2:9, function(i) {
      lines(res.out[, 1], res.out[, i], type = "b", col = i, pch = i + 13)
     }
)

dev.off()

############################################
pdf(file = "fixNA_CPP.pdf", width = 11, height = 8)
require(MBmca)
# Draw an empty plot for 45 cycles with user defined parameters.
layout(matrix(c(1,2,1,3), 2, 2, byrow = TRUE))
par(las = 0, bty = "n", cex.axis = 1.5, cex.lab = 1.5, font = 2, 
    cex.main = 1.1, oma = c(1,1,1,1))
    
plot(NA, NA, xlim = c(1,45), ylim = c(0, 0.7), xlab = "Cycle", 
     ylab = "refMFI", main = "Raw data")
mtext("A", cex = 2, side = 3, adj = 0, font = 2) 

# Take the first 45 cycles of the C54 data set and assign it to the
# object dat.
dat <- C54[1:45, ]

# Plot the amplification curve data.
apply(C54[, c(2:4)], 2, function(y) lines(C54[, 1], y))

# Draw an empty plot for the output of the preprocessed data.

plot(NA, NA, xlim = c(1,45), ylim = c(0, 0.55), xlab = "Cycle", 
     ylab = "refMFI", main = "Pre-processed data")
mtext("B", cex = 2, side = 3, adj = 0, font = 2) 

# Use CPP to preprocesses the data. CPP is used with the default parameters.
# This means: the slope of the background range is corrected by a linear 
# robust regression (MM-estimator), outliers are not removed, the curve is
# Smoothed be a Savitzky-Golay smoothing filter and NAs are removed by the 
# spline method. The results of CPP are assigned to the objects D1, D2, and D3.

D1 <- cbind(1:45, CPP(dat[, 1], dat[, 2])$y.norm)
D2 <- cbind(1:45, CPP(dat[, 1], dat[, 3])$y.norm)
D3 <- cbind(1:45, CPP(dat[, 1], dat[, 4])$y.norm)

# Plot the output of the preprocessed curves.
lines(D1, col = 1)
lines(D2, col = 2)
lines(D3, col = 3)

# Define the dilutions
dilution <- c(1E0, 1E-3, 1E-6)

# The function diffQ2 with the parameter "inder = TRUE" is used to calculate 
# the Cq value by means SDM. The results are assigned to the corresponding 
# concentration in the object res.

res <- data.frame(dilution, 
		  rbind(diffQ2(D1, inder = TRUE, warn = FALSE)$xTm1.2.D2[1], 
		        diffQ2(D2, inder = TRUE, warn = FALSE)$xTm1.2.D2[1], 
			diffQ2(D3, inder = TRUE, warn = FALSE)$xTm1.2.D2[1]
			)
		 )
# Finally the object res is used in effcalc to show the results  in the
# amplification efficiency plot.

effcalc(res[, 1], res[, 2])

dev.off()

#########################################################
pdf(file = "normalization.pdf", width = 11, height = 11)
par(mfrow = c(2,3), las = 0, bty = "n", cex.axis = 1.5, cex.lab = 1.5, 
    font = 2, cex.main = 1.1, oma = c(1,1,1,1))
tmp <- VIMCFX96_60

plot(NA, NA, xlim = c(1,40), ylim = c(0, 6000), xlab = "Cycle", 
     ylab = "RFU", main = "Raw data")
mtext("A", cex = 2, side = 3, adj = 0, font = 2) 
apply(tmp[, -1], 2, function(x) lines(tmp[, 1], x))
abline(lm(rowMeans(tmp[2:10, 2L:ncol(tmp)]) ~ tmp[2:10, 1]), col = 2)

plot(NA, NA, xlim = c(1,40), ylim = c(0, 3300), xlab = "Cycle", 
     ylab = "RFU", main = "Baselined data")
mtext("B", cex = 2, side = 3, adj = 0, font = 2) 
apply(tmp[, -1], 2, function(x) lines(tmp[, 1], CPP(tmp[, 1], x, 
	method.norm = "none")$y))
	

plot(NA, NA, xlim = c(1,40), ylim = c(0, 1.15), xlab = "Cycle", 
     ylab = "RFU", main = "MinMax-Normalization")
mtext("C", cex = 2, side = 3, adj = 0, font = 2) 
apply(tmp[, -1], 2, function(x) lines(tmp[, 1], CPP(tmp[, 1], x, 
	method.norm = "minmax")$y))
	
plot(NA, NA, xlim = c(1,40), ylim = c(0, 1.15), xlab = "Cycle", 
     ylab = "RFU", main = "Max-Normalization")
mtext("D", cex = 2, side = 3, adj = 0, font = 2) 
apply(tmp[, -1], 2, function(x) lines(tmp[, 1], CPP(tmp[, 1], x,, 
	method.norm = "max")$y))
     
plot(NA, NA, xlim = c(1,40), ylim = c(0, 1.15), xlab = "Cycle", 
     ylab = "RFU", main = "luqn-Normalization")
mtext("E", cex = 2, side = 3, adj = 0, font = 2) 
apply(tmp[, -1], 2, function(x) lines(tmp[, 1], CPP(tmp[, 1], x, 
	method.norm = "luqn", qnL = 0.03)$y))

plot(NA, NA, xlim = c(1,40), ylim = c(-1.5, 1.5), xlab = "Cycle", 
     ylab = "RFU", main = "zscore-Normalization")
mtext("F", cex = 2, side = 3, adj = 0, font = 2) 
apply(tmp[, -1], 2, function(x) lines(tmp[, 1], CPP(tmp[, 1], x, 
	method.norm = "zscore")$y))
    
dev.off()
#########################################################################
pdf(file = "HDA.pdf", width = 11, height = 8)
par(las = 0, bty = "n", cex.axis = 1.5, cex.lab = 1.5, 
     font = 2, cex.main = 1.1, oma = c(1,1,1,1))

data(C85)
# First example
plot(NA, NA, xlim = c(0,85), ylim = c(0,1), xlab = "Time [min]", 
     ylab = "refMFI", main = "HDA")
  points(C85[, 2]/60, C85[, 3], type = "b", col = 1, pch = 20)
  points(C85[, 4]/60, C85[, 5], type = "b", col = 2, pch = 20)
  points(C85[, 6]/60, C85[, 7], type = "b", col = 3, pch = 20)
legend(40, 0.5, c("D1, 1x", "D2, 1:10", "D3, 1:100"), col = c(1:3), 
	pch = rep(20,3))
dev.off()

##################################

colors <- rep(rainbow(7), each = 2)
par(mfrow = c(2,2))

plot(NA, NA, xlim = c(0,44), ylim = c(0, 6), xlab = "Cycles", ylab = "RFU")
legend(0, 6, colnames(C60.amp[, 4L:17]), ncol = 2, col = colors[1:14], pch = 19, bty = "n")
SDM.vim <- sapply(4L:17, function(i) {
	lines(C60.amp[, 1], C60.amp[, i], col = colors[i - 3])
	SDM <- summary(inder(C60.amp[, 1], C60.amp[, i]))[2]
	}
)

plot(NA, NA, xlim = c(0,44), ylim = c(0, 4), xlab = "Cycles", ylab = "RFU")
legend(0, 4, colnames(C60.amp[, 18L:31]), ncol = 2, col = colors[1:14], pch = 19, bty = "n")
SDM.mlc2v <- sapply(18L:31, function(i) {
	lines(C60.amp[, 1], C60.amp[, i], col = colors[i - 17])
	SDM <- summary(inder(C60.amp[, 1], C60.amp[, i]))[2]
	}
)

dil <- rep(sapply(-(0:6), function(i) {10^i}), each = 2)

res <- cbind(dil, SDM.vim, SDM.mlc2v)

effcalc(res[, 1], res[, 2])
effcalc(res[, 1], res[, 3])



##################################

library(reshape2)
library(ggplot2)
pdf(file = "fixNA_stats.pdf", width = 12, height = 10)
tester <- function(cycles, y.raw, y.na, bg.range = 1L:10) {  
  data <- cbind(cycles = cycles, raw = y.raw, na = fixNA(cycles, y.na))
  fit.raw <- pcrfit(data = data, 1, 2, l5)
  fit.na <- pcrfit(data = data, 1, 3, l5)
  
  bg.raw = mean(predict(fit.raw)[bg.range, ])
  cpD2.raw = efficiency(fit.raw, plot = FALSE)$cpD2
  Cy0.raw = Cy0(fit.raw)
  
  bg.fix = mean(predict(fit.na)[bg.range, ])
  cpD2.fix = efficiency(fit.na, plot = FALSE)$cpD2
  Cy0.fix = Cy0(fit.na)
  
  res.cor <- cor.test(data[, 3], data[, 2], method = "pearson")
  res.lm <- lm(data[, 3] ~ data[, 2])
  sum.lm <- summary(res.lm)
  mes <- res.lm[["model"]][, 1]
  RMSE <- sqrt(mean(residuals(res.lm)^2))
  NRMSE <- RMSE/(max(mes) - min(mes))
  NRMSE.warning <- ifelse(NRMSE > 0.08, TRUE, FALSE)
  
  c(R.squared = sum.lm$r.squared, 
    r = res.cor[["estimate"]], 
    r.p.value = res.cor[["p.value"]],
    bg.raw = bg.raw,
    bg.fix = bg.fix,
    cpD2.raw = cpD2.raw,
    cpD2.fix = cpD2.fix,
    Cy0.raw = Cy0.raw,
    Cy0.fix = Cy0.fix,
    NRMSE = NRMSE, 
    RMSE = RMSE, 
    NRMSE.warning = NRMSE.warning
    #     ,
    #     orig = data[is.na(y.na), "raw"],
    #     fixed = data[is.na(y.na), "na"]
  )
}


calc_fixna_error <- function(x, y, measure_range = 1L:10, three_points = FALSE) {
  t(apply(y, 2, function(i) {
    if (three_points) {
      removed_points <- sample(measure_range[-c(1, length(measure_range))], 1)
      removed_points <- c(removed_points - 1, removed_points, removed_points + 1)
    } else {
      removed_points <- sample(measure_range, 1)
    }
    with_na <- i
    with_na[removed_points] <- NA
    
    tester(cycles = x, y.raw = i, y.na = with_na, bg.range = measure_range)
  }))
}

set.seed(1988)

# Single point ------------------------------------------------------
#linear phase error
phasp_error <- calc_fixna_error(reps384[, 1], reps384[, -1], 1L:10)
#exponential phase error
expp_error <- calc_fixna_error(reps384[, 1], reps384[, -1], 11L:33)
#plateau phase error
platp_error <- calc_fixna_error(reps384[, 1], reps384[, -1], 34L:40)

# Three points ------------------------------------------------------
#linear phase error
phasp3_error <- calc_fixna_error(reps384[, 1], reps384[, -1], 1L:10, TRUE)
#exponential phase error
expp3_error <- calc_fixna_error(reps384[, 1], reps384[, -1], 11L:33, TRUE)
#plateau phase error
platp3_error <- calc_fixna_error(reps384[, 1], reps384[, -1], 34L:40, TRUE)


#plots
library(reshape2)
library(ggplot2)
library(gridExtra)

melt_test <- function(r1, r2, r3, number_of_points) {
  enhance_test <- function(i, region)
    cbind(i, region = rep(region, nrow(i)))
  
  melt_dat <- melt(rbind(enhance_test(data.frame(r1), "linear phase"),
                         enhance_test(data.frame(r2), "exponential phase"),
                         enhance_test(data.frame(r3), "plateau phase")))
  
  agg_dat <- aggregate(value ~ region + variable, melt_dat, mean)
  res <- enhance_test(agg_dat, number_of_points)
  colnames(res) <- c("region", "measure", "value", "NA")
  res
}

dat <- rbind(melt_test(phasp3_error, expp3_error, platp3_error, "3 NA"),
             melt_test(phasp_error, expp_error, platp_error, "1 NA"))

# ggplot(dat[!dat[["measure"]] %in% c("bg.raw", "bg.fix", "NRMSE.warning"), ], 
#        aes(x = region, y = value, fill = points)) +
#   geom_bar(stat = "identity", position = "dodge") + facet_wrap(~ measure, nrow = 3)


plot_measure <- function(what)
  ggplot(dat[dat[["measure"]] %in% what, ], 
         aes(x = region, y = value, fill = `NA`)) +
  geom_bar(stat = "identity", position = "dodge") + facet_grid(. ~ measure)

p1 <- plot_measure(c("R.squared", "R.cor"))
p2 <- plot_measure(c("cpD2.raw", "cpD2.fix"))
p3 <- plot_measure(c("Cy0.raw", "Cy0.fix"))
p4 <- plot_measure(c("r.p.value", "NRMSE"))

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp[["grobs"]], function(x) x[["name"]]) == "guide-box")
  tmp[["grobs"]][[leg]]
}

mylegend<-g_legend(p1)


grid.arrange(arrangeGrob(p1 + theme(legend.position="none"), 
                         p2 + theme(legend.position="none"), 
                         p3 + theme(legend.position="none"), 
                         p4 + theme(legend.position="none"), nrow = 2),
             mylegend, ncol=2,widths=c(10, 1))
dev.off()

head(phasp3_error)

data_sets_list <- list(phasp_error = phasp_error, 
                       expp_error = expp_error, 
                       platp_error = platp_error,
                       phasp3_error = phasp3_error, 
                       expp3_error = expp3_error, 
                       platp3_error = platp3_error)
fixNA_summ <- do.call(rbind, lapply(data_sets_list, function(data_set)
              apply(data_set[, c("bg.raw", "bg.fix", "cpD2.raw", 
                                 "cpD2.fix", "Cy0.raw", "Cy0.fix", "NRMSE")], 2,
                    function(i) paste0(format(mean(i), digits = 2), " \\pm ", 
                                       format(sd(i), digits = 3)))))
rownames(fixNA_summ) <- c("Linear phase - 1 NA", "Exponential phase - 1 NA",
                          "Plateau phase  - 1 NA", "Linear phase - 3 NA", 
                          "Exponential phase - 3 NA", "Plateau phase  - 3 NA")
colnames(fixNA_summ) <- c("Background - raw data", "Background - fixed NA",
                          "cpD2 - raw data", "cpD2 - fixed NA",
                          "Cy0 - raw data", "Cy0 - fixed NA",
                          "NRMSE")

imp_measurements <- c("cpD2.fix", "Cy0.fix")
kr_tests <- sapply(imp_measurements, function(data_column) {
           kr_dat <- melt(cbind(phasp_error = phasp_error[, data_column], 
                        expp_error = expp_error[, data_column], 
                        platp_error = platp_error[, data_column],
                        phasp3_error = phasp3_error[, data_column], 
                        expp3_error = expp3_error[, data_column], 
                        platp3_error = platp3_error[, data_column]))
           kruskal.test(value ~ Var2, kr_dat)[["p.value"]]
           })

kr_tests <- sapply(imp_measurements, function(data_column) {
  kr_dat <- melt(cbind(phasp_error = phasp_error[, paste0(strsplit(data_column, ".", 
                                                                   fixed = TRUE)[[1]][1], ".raw")],
                       phasp_error = phasp_error[, data_column], 
                       expp_error = expp_error[, data_column], 
                       platp_error = platp_error[, data_column],
                       phasp3_error = phasp3_error[, data_column], 
                       expp3_error = expp3_error[, data_column], 
                       platp3_error = platp3_error[, data_column]))
  kruskal.test(value ~ Var2, kr_dat)[["p.value"]]
})

