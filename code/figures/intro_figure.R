#****************************************************************************************************************************************
# Introductory figure
#****************************************************************************************************************************************
# Read gene coordinates
mRNA_coords <- read.table("../../../../sex_chr/ostrich_z/Z_mRNA_chro_coordinates.txt")
gameto_mRNA_coords <- read.table("/Users/homapapoli/Documents/projects/sex_chr/ostrich_z/gametolog_chro_coordinates.txt")

head(mRNA_coords)
plot(c(0, 83), c(0, 3), type= "n", xlab = "", ylab = "",  bty="l", axes=F)
mb = 10^6

for(i in seq(1,849)){
  rect(mRNA_coords[i,5]/mb, 0, mRNA_coords[i,6]/mb, 1, col = "grey", border = "transparent")
}

for(i in seq(1,42)){
  rect(gameto_mRNA_coords[i,5]/mb, 0, gameto_mRNA_coords[i,6]/mb, 1, col = "red", border = "transparent")
}

rect(76174922/mb, 0,  76230182/mb, 1, col = "blue", border = "black")
abline(v = 52193205/mb, col = "black", lty = 2)
rect(0, 0, 82, 1, border = "black")
axis(side = 1, at = seq(0, 82, 10))
text(x = 76174922/mb, y = 1.2, "DMRT1")
#****************************************************************************************************************************************
# Recombination frequency
female_PAR_length = seq(0, 80.628, 1)
male_PAR_length = seq(0, 42.641, 1)

female_PAR_r = round((0.5 * ((exp(4*(female_PAR_length/100))-1)/(exp(4*(female_PAR_length/100))+1))), 3)
male_PAR_r = round((0.5 * ((exp(4*(male_PAR_length/100))-1)/(exp(4*(male_PAR_length/100))+1))), 3)

pdf("../../figures/intro_rec_frequency.pdf", height = 3.71, width = 3.65)
plot(rev(seq(52.18, 0, -0.65225)), rev(female_PAR_r), col = "red", type = "l", lwd = 2, ylab = "", xlab = "", xaxt="n", yaxt="n", axes=F)
lines(rev(c(seq(52.18, 0, -1.213489)[-43], 0)), rev(c(male_PAR_r)), lwd = 2, col = "blue")
axis(1, seq(0, 60, 10), pos=0)
axis(2, c(0, 0.1, 0.2, 0.3, 0.4, 0.5), pos=0, las = 2)
title(ylab = "Recombination Frequency (r)", x = "Position (Mb)", line=2, cex.lab=1)
legend(32, 0.48, legend=c("Female", "Male"), col=c("red", "blue"), lty=1, lwd = 2, cex=0.8, bty = "n")
dev.off()

# I have used the total genetic map length of the PAR to obtain the recombination frequency for the PAR.
# To plot, I have divided the PAR length into equal distances and then plotted the according r for the given distance.
