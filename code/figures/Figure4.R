###############
# DEPENDENCIES
###############
library(lattice)
library(latticeExtra)
library(wesanderson)
library(dplyr)


#######################
# AUXILLIARY FUNCTIONS
#######################

toPdf <- function(expr, filename, ...) {
  toDev(expr, pdf, filename, ...)
}

figPath  <-  function(name) {
  file.path('./', name)
}

toDev <- function(expr, dev, filename, ..., verbose=TRUE) {
  if ( verbose )
    cat(sprintf('Creating %s\n', filename))
  dev(filename, family='CM Roman', ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}



####################
# PLOTTING FUNCTIONS
####################

#' Plot text or points according to relative axis position.
#'
#' @title Plot text or points according to relative axis position
#' @param px Relative x-axis position (in proportion) where character is to be plotted.
#' @param py Relative y-axis position (in proportion) where character is to be plotted.
#' @param lab Plotted text. Works if argument \code{\link[graphics]{text}} is \code{TRUE}.
#' @param adj See argument of same name in R base function \code{\link[graphics]{par}}.
#' @param text Logical. Should text or points be plotted?
#' @param log Used if the original plot uses the argument log, e.g. \code{log='x'}, \code{log='y'} or \code{log='xy'}.
#' @param ... Additional arguments to R base function \code{\link[graphics]{text}}.
#' @export
proportionalLabel <- function(px, py, lab, adj=c(0, 1), text=TRUE, log=FALSE, ...) {
  usr  <-  par('usr')
  x.p  <-  usr[1] + px*(usr[2] - usr[1])
  y.p  <-  usr[3] + py*(usr[4] - usr[3])
  if(log=='x') {
    x.p<-10^(x.p)
  }
  if(log=='y') {
    y.p<-10^(y.p)
  }
  if(log=='xy') {
    x.p<-10^(x.p)
    y.p<-10^(y.p)
  }
  if(text){
    text(x.p, y.p, lab, adj=adj, ...)
  } else {
    points(x.p, y.p, ...)
  }
}



proportionalArrows <- function(px1, py1, px2, py2, adj=c(0, 1), log=FALSE, length=length, ...) {
  usr  <-  par('usr')
  x.p1  <-  usr[1] + px1*(usr[2] - usr[1])
  y.p1  <-  usr[3] + py1*(usr[4] - usr[3])
  x.p2  <-  usr[1] + px2*(usr[2] - usr[1])
  y.p2  <-  usr[3] + py2*(usr[4] - usr[3])
  if(log=='x') {
    x.p1  <-  10^(x.p1)
    x.p2  <-  10^(x.p2)
  }
  if(log=='y') {
    y.p1  <-  10^(y.p1)
    y.p2  <-  10^(y.p2)
  }
  if(log=='xy') {
    x.p1  <-  10^(x.p1)
    y.p1  <-  10^(y.p1)
    x.p2  <-  10^(x.p2)
    y.p2  <-  10^(y.p2)
  }
  arrows(x0=x.p1, y0=y.p1, x1=x.p2, y1=y.p2, length=length,...)
}

#' Draw equally-spaced white lines on plot window.
#'
#' @title Equally-spaced white lines on plot window
#' @param ... Additional arguments to internal function \code{\link{proportionalLabel}}.
#' @author Diego Barneche
#' @export
plotGrid  <-  function(lineCol='white',...) {
  proportionalLabel(rep(0.2, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
  proportionalLabel(rep(0.4, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
  proportionalLabel(rep(0.6, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
  proportionalLabel(rep(0.8, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
  proportionalLabel(c(0,1), rep(0.2, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
  proportionalLabel(c(0,1), rep(0.4, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
  proportionalLabel(c(0,1), rep(0.6, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
  proportionalLabel(c(0,1), rep(0.8, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
}


#' Internal. Create nice rounded numbers for plotting.
#'
#' @title Rounded numbers for plotting
#' @param value A numeric vector.
#' @param precision Number of rounding digits.
#' @return A character vector.
#' @author Diego Barneche.
rounded  <-  function(value, precision=1) {
  sprintf(paste0('%.', precision, 'f'), round(value, precision))
}


#' Creates transparent colours
#'
#' @title Creates transparent colours
#' @param col Colour.
#' @param opacity equivalent to alpha transparency parameter
#' @export
transparentColor <- function(col, opacity=0.5) {
  if (length(opacity) > 1 && any(is.na(opacity))) {
    n        <-  max(length(col), length(opacity))
    opacity  <-  rep(opacity, length.out=n)
    col      <-  rep(col, length.out=n)
    ok       <-  !is.na(opacity)
    ret      <-  rep(NA, length(col))
    ret[ok]  <-  Recall(col[ok], opacity[ok])
    ret
  } else {
    tmp  <-  col2rgb(col)/255
    rgb(tmp[1,], tmp[2,], tmp[3,], alpha=opacity)
  }
}

fibonacci.scale  <-  function(n) {
  fibs  <-  c(0,1)
  for(i in 2:n) {
    fibs  <-  c(fibs, (fibs[i] + fibs[i-1]))
  }
  (fibs/max(fibs))[-2]
}


############################
## 4 panel figure for paper
############################

# Import data for Full PAR
fullPARCoord  <- read.table(file="../../data/simulation_input/rho_r_Ne_table.txt", header=TRUE)
fullPAR2     <-  read.csv(file="../../simulation_output/fullPAR_Black_smooth_PlotData.csv", header=TRUE)
Z_pi <- read.table("../../data/diversity/genetic_variation/z/Z.200000.chr.coord.sfs.txt", header=F)
colnames(Z_pi) <- c("chr", "cstart", "cend", "scaffold",	"start",	"end", "swin_midpoint", "pi", "theta", "td",	"winsize",	"pi_resam_mean",	"pi_lowCI", "pi_upCI")
empPi_full <- Z_pi[Z_pi$cstart<52000000,]

fst <- read.table("../../data/fst/z/black_male_female_Z_200Kb.windowed.weir.Z.coord.fst")
colnames(fst) <- c("Chr", "Chr_start", "Chr_end", "Scaffold", "Window_start",
                   "Window_end", "N_VARIANTS", "WEIGHTED_FST", "MEAN_FST")

fst$WEIGHTED_FST[fst$WEIGHTED_FST<0] <- 0

empFst_full <- fst[fst$Window_start<52000000,]

fullPAR2$midpoint   <-  round(fullPAR2$start + (fullPAR2$end - fullPAR2$start)/2)
fullPAR2$midPosMB   <-  fullPAR2$midpoint/1000000
#    fullPAR2$midPosMB_fwd   <-  52-fullPAR2$midPosMB
empPi_full$midPosMB   <-  (empPi_full$cstart + empPi_full$cend)/2000000
empFst_full$midPosMB  <-  (empFst_full$Chr_start + empFst_full$Chr_end)/2000000
empPi_full  <-  empPi_full[empPi_full$midPosMB < 52,]
empFst_full  <-  empFst_full[empFst_full$midPosMB < 52,]


# Import data for PAR boundary
PARbound2    <-  read.table(file="../../simulation_output/PARboundary_smooth_PlotData_zcoord.csv", header=TRUE)

empPi_bound   <-  read.table(file="../../data/boundary/superscaffold36.1000.sfs.z.coord.txt", header=FALSE)
colnames(empPi_bound) <- c("chr", "cstart", "cend", "scaffold",	"start",	"end", "swin_midpoint", "pi", "theta", "td",	"winsize")

empFst_bound  <-  read.table(file="../../data/boundary/male_female_Z_1Kb.windowed.weir.zcoord.fst", header=FALSE)
colnames(empFst_bound) <- c("Chr", "Chr_start", "Chr_end", "Scaffold", "Window_start",
                   "Window_end", "N_VARIANTS", "WEIGHTED_FST", "MEAN_FST")

empFst_bound$WEIGHTED_FST[empFst_bound$WEIGHTED_FST<0] <- 0

PARbound2$pos.coord <- PARbound2$chr_end
PARbound2$pos.coord.Mb  <-   PARbound2$pos.coord/10^6

# Trim Empirical data where necessary, calculate positions in Mb

empPi_bound$midPos   <-   (empPi_bound$cstart+empPi_bound$cend)/2
empPi_bound$midPosMB   <-  empPi_bound$midPos/1000000
empFst_bound$midpoint  <-  (empFst_bound$Chr_start + empFst_bound$Chr_end)/2
empFst_bound$midPos   <-   empFst_bound$midpoint
empFst_bound$midPosMB  <-  empFst_bound$midPos/1000000
#empPi_bound            <-  empPi_bound[-c(1:4),]
#empFst_bound           <-  empFst_bound[-c(1:4),]

# Color scheme
COLS  <-  list(
  "Tt"        =  transparentColor('#252525', opacity=1),
  "FstXY"     =  transparentColor('dodgerblue', opacity=1),
  "FstFM"     =  transparentColor('tomato', opacity=1),
  "Tt_ci"     =  transparentColor('#252525', opacity=0.3),
  "FstXY_ci"  =  transparentColor('dodgerblue', opacity=0.3),
  "FstFM_ci"  =  transparentColor('tomato', opacity=0.3)                    )

#*******************************************************************************
# Save as PDF
#*******************************************************************************
pdf("../../figures/Figure4.pdf", height = 8.1, width = 9.8)

#*******************************************************************************
# Set plot layout
#*******************************************************************************
layout.mat  <- matrix(c(1,1,2,2,
                        3,3,4,4,
                        5,5,6,6), nrow=3, ncol=4, byrow=TRUE)
layout      <- layout(layout.mat,respect=TRUE)

par(mfrow=c(2,2))

#*******************************************************************************
# Autosomal average diversity
#*******************************************************************************
autoPi   <-  0.00165

#*******************************************************************************
#*******************************************************************************
## Panel A: Nucleotide Diversity - FULL PAR
#*******************************************************************************
#*******************************************************************************
# Make the plot
par(mar = c(3, 5, 3, 1))
plot(NA, axes=FALSE, type='n', main='', xlim = c((min(fullPAR2$midPosMB)), max(fullPAR2$midPosMB)), ylim = c(0,0.004), ylab='', xlab='', cex.lab=1.2)
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
#        plotGrid(lineCol='grey80')
box()
# Simulation Results
polygon(x = c(fullPAR2$midPosMB, rev(fullPAR2$midPosMB)), y=c((fullPAR2$Tt_ci_hi*autoPi), rev((fullPAR2$Tt_ci_lo*autoPi))), col=COLS$Tt_ci, border=COLS$Tt_ci)
lines(Tbar_t*autoPi ~ fullPAR2$midPosMB, lwd=2, col=COLS$Tt, data=fullPAR2)
abline(v=max(fullPAR2$midPosMB), lwd=2, lty=2)
abline(h=autoPi, lwd=1, lty=2, col=2)
# Empirical estimates
points((pi/rev(winsize))[empPi_full$midPosMB >= min(fullPAR2$midPosMB)] ~ empPi_full$midPosMB[empPi_full$midPosMB >= min(fullPAR2$midPosMB)], pch=21, col="black", bg=COLS$Tt_ci, data=empPi_full)
ss5  <-  smooth.spline(y=(empPi_full$pi/empPi_full$winsize)[empPi_full$midPosMB >= min(fullPAR2$midPosMB)], x=empPi_full$midPosMB[empPi_full$midPosMB >= min(fullPAR2$midPosMB)], df=5)
lines(ss5, lwd=1, col=4)
# axes
axis(1, las=1)
axis(2, las=1)
# Labels/annotations
proportionalLabel(0.03, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
proportionalLabel(-0.18, 0.5, expression(paste("Genetic diversity (", bar(italic(pi)), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
# legend
legend( x       =  24,
        y       =  usr[4],
        legend  =  c(
          expression(paste("Autosomal  ", bar(italic(pi)))),
          expression(paste("Expected ", bar(italic(pi))[PAR]))),
        cex     =  1.25,
        lty     =  c(2,1),
        lwd     =  c(1,2),
        col     =  c(2, COLS$Tt),
        xjust   =  1,
        yjust   =  1,
        bty     =  'n',
        border  =  NA)

#*******************************************************************************
#*******************************************************************************
## Panel B: Fst MF - FULL PAR
#*******************************************************************************
#*******************************************************************************
FstFM_ci_lo_zero  <-  fullPAR2$FstFM_ci_lo
FstFM_ci_lo_zero[FstFM_ci_lo_zero < 0]  <-  0

plotFST  <-  empFst_full$WEIGHTED_FST[empFst_full$midPosMB >= min(fullPAR2$midPosMB)]
plotFST[plotFST < 0]  <-  0

# Make the plot
par(mar = c(3, 5, 3, 1))
plot(NA, axes=FALSE, type='n', main='', xlim = c((min(fullPAR2$midPosMB)), max(fullPAR2$midPosMB)), ylim = c(0,max(fullPAR2$FstFM_ci_hi, 0.3)), ylab='', xlab='', cex.lab=1.2)
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
#        plotGrid(lineCol='grey80')
box()
# Benchmarks
abline(v=max(fullPAR2$midPosMB), lwd=2, lty=2)
abline(h=0, lwd=1, lty=1)
# Simulation Results
polygon(x = c(fullPAR2$midPosMB, rev(fullPAR2$midPosMB)), y=c(fullPAR2$FstFM_ci_hi, rev(FstFM_ci_lo_zero)), col=COLS$FstFM_ci, border=COLS$FstFM_ci)
lines(FstBarFM ~ fullPAR2$midPosMB, lwd=2, col=COLS$FstFM, data=fullPAR2)
# Empirical estimates
points(plotFST ~ empFst_full$midPosMB[empFst_full$midPosMB >= min(fullPAR2$midPosMB)], pch=21, col="black", bg=COLS$FstFM_ci, data=empFst_full)
ss5  <-  smooth.spline(y=plotFST, x = empFst_full$midPosMB[empFst_full$midPosMB >= min(fullPAR2$midPosMB)], df=5)
lines(ss5, lwd=1, col=4)
# axes
axis(1, las=1)
axis(2, las=1)
# Labels/annotations
proportionalLabel(0.03, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
proportionalLabel(-0.18, 0.5, expression(paste("Divergence (", bar(italic(F))[FM], ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
# legend
legend( x       =  24,
        y       =  usr[4],
        legend  =  c(
          expression(paste("Expected ", bar(italic(F))[FM]))),
        cex     =  1.25,
        lty     =  c(1),
        lwd     =  c(2),
        col     =  c(COLS$FstFM),
        xjust   =  1,
        yjust   =  1,
        bty     =  'n',
        border  =  NA)

#*******************************************************************************
#*******************************************************************************
## Panel C: Nucleotide Diversity - PAR BOUNDARY
#*******************************************************************************
#*******************************************************************************
# Make the plot
empPi_bound <- empPi_bound %>% filter(if_all(everything(), ~ !is.na(.x)))

par(mar = c(5, 5, 1, 1))
plot(NA, axes=FALSE, type='n', main='', xlim = c((min(PARbound2$pos.coord.Mb)), max(PARbound2$pos.coord.Mb)), ylim = c(0,0.006), ylab='', xlab='', cex.lab=1.2)
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
#        plotGrid(lineCol='grey80')
box()
# Simulation Results
polygon(x = c(PARbound2$pos.coord.Mb, rev(PARbound2$pos.coord.Mb)), y=c((PARbound2$Tt_ci_hi*autoPi), rev((PARbound2$Tt_ci_lo*autoPi))), col=COLS$Tt_ci, border=COLS$Tt_ci)
lines(Tbar_t*autoPi ~ PARbound2$pos.coord.Mb, lwd=2, col=COLS$Tt, data=PARbound2)
abline(v=max(PARbound2$pos.coord.Mb), lwd=2, lty=2)
abline(h=autoPi, lwd=1, lty=2, col=2)
# Empirical estimates
points((pi/winsize)[empPi_bound$midPosMB >= min(PARbound2$pos.coord.Mb)] ~ empPi_bound$midPosMB[empPi_bound$midPosMB >= min(PARbound2$pos.coord.Mb)], pch=21, col="black", bg=COLS$Tt_ci, data=empPi_bound)

ss5  <-  smooth.spline(y=(empPi_bound$pi/empPi_bound$winsize)[empPi_bound$midPosMB >= min(PARbound2$pos.coord.Mb)], x=empPi_bound$midPosMB[empPi_bound$midPosMB >= min(PARbound2$pos.coord.Mb)], df=3)
#lines(ss5, lwd=1, col=4)
# axes
axis(1, las=1)
axis(2, las=1)
# Labels/annotations
#proportionalLabel(1.125, 1.25, expression(paste("PAR boundary")), cex=2, adj=c(0.5, 0.5), xpd=NA, srt=0)
proportionalLabel(0.03, 1.05, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
proportionalLabel(-0.18, 0.5, expression(paste("Genetic diversity (", bar(italic(pi)), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
proportionalLabel(0.5, -0.17, expression(paste("Physical Position (Mb)")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)

#*******************************************************************************
#*******************************************************************************
## Panel D: Fst M-F - PAR BOUNDARY
#*******************************************************************************
#*******************************************************************************
FstFM_ci_lo_zero  <-  PARbound2$FstFM_ci_lo
FstFM_ci_lo_zero[FstFM_ci_lo_zero < 0]  <-  0

plotFST  <-  empFst_bound$MEAN_FST[empFst_bound$midPosMB >= min(PARbound2$pos.coord.Mb)]
plotFST[plotFST < 0]  <-  0
# Make the plot
par(mar = c(5, 5, 1, 1))
plot(NA, axes=FALSE, type='n', main='', xlim = c((min(PARbound2$pos.coord.Mb)), max(PARbound2$pos.coord.Mb)), ylim = c(0, 0.5), ylab='', xlab='', cex.lab=1.2)
usr  <-  par('usr')
rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
#        plotGrid(lineCol='grey80')
box()
# Benchmarks
abline(v=max(PARbound2$pos.coord.Mb), lwd=2, lty=2)
abline(h=0, lwd=1, lty=1)
# Simulation Results
polygon(x = c(PARbound2$pos.coord.Mb, rev(PARbound2$pos.coord.Mb)), y=c(PARbound2$FstFM_ci_hi, rev(FstFM_ci_lo_zero)), col=COLS$FstFM_ci, border=COLS$FstFM_ci)
lines(FstBarFM ~ PARbound2$pos.coord.Mb, lwd=2, col=COLS$FstFM, data=PARbound2)
# Empirical estimates
points(plotFST ~ empFst_bound$midPosMB[empFst_bound$midPosMB >= min(PARbound2$pos.coord.Mb)], pch=21, col="black", bg=COLS$FstFM_ci, data=empFst_bound)
ss5  <-  smooth.spline(y=plotFST, x=rev(empFst_bound$midPosMB[empFst_bound$midPosMB >= min(PARbound2$pos.coord.Mb)]), df=5)
#lines(ss5, lwd=1, col=4)
# axes
axis(1, las=1)
axis(2, las=1)
# Labels/annotations
proportionalLabel(0.03, 1.05, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
proportionalLabel(-0.15, 0.5, expression(paste("Divergence (", bar(italic(F))[FM], ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
proportionalLabel(0.5, -0.17, expression(paste("Physical Position (Mb)")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=0)

#*******************************************************************************
dev.off()
#*******************************************************************************







#### Smoothed r values
#par(omi=rep(0.5, 4), mar = c(5,5,2,3), bty='o', xaxt='s', yaxt='s')
#plot(NA, axes=FALSE, type='n', main='', xlim = c(min(rho.pos$positions), max(rho.pos$positions)), ylim = c(0,max(rho.pos$rho_rates_kb)), ylab='', xlab='', cex.lab=1.2)
#usr  <-  par('usr')
#rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
#        plotGrid(lineCol='grey80')
#box()
#abline(v=51.575, lwd=2, lty=2)
# Plot rho, cM
#lines(cumRho ~ positions, lwd=2, col = transparentColor('#252525', opacity=1), data=rho.pos)
# axes
#axis(1, las=1)
#axis(2, las=1)
#par(new=T)
#plot(NA, axes=FALSE, type='n', main='', xlim = c(min(cM.pos$positions), max(cM.pos$positions)), ylim = c(0,max(cM.pos$sex_averaged_cM)), ylab='', xlab='', cex.lab=1.2)
#lines(sex_averaged_cM_smooth_rev_predict ~ positions, lwd=2, col = transparentColor('darkorchid1', opacity=1), data=cM.pos)
#points(rev.cM ~ positions, pch=21, col = transparentColor('darkorchid1', opacity=1), bg = transparentColor('darkorchid1', opacity=0.5), data=cM.pos)
# 2nd y-axis
#axis(4, las=1)
# Labels/annotations
#proportionalLabel(1.12, 1.25, expression(paste("Full PAR")), cex=2, adj=c(0.5, 0.5), xpd=NA, srt=0)
#        proportionalLabel(0.5, 1.2, expression(paste("Cumulative recombination")), cex=2, adj=c(0.5, 0.5), xpd=NA, srt=0)
#proportionalLabel(0.03, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#proportionalLabel(-0.15, 0.5, expression(paste("Cumulative ", italic(rho))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
#proportionalLabel(1.15, 0.5, expression(cM), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=270)


#plot(NA, axes=FALSE, type='n', main='', xlim = c(min(rho.pos$rev.pos), max(rho.pos$rev.pos)), ylim = c(0,max(rho.pos$rho_rates_kb)), ylab='', xlab='', cex.lab=1.2)
