# Scripts for the plots presented in the paper

Plots 2 to 4 and the supplementary plots S1 and S2 can be produced by running
the script:

`bash figures.sh` <br>

For figure 1, run the script `intro_figure.R` in Rstudio and produce the two plots separately. 
Using inkscape, you can assemble the two plots as figure 1.

The directory also contains **Z_mRNA_chro_coordinates.txt** used for gene coordinates
plotted in Figure 1.

For all figures, some edits to figure aesthetic might have been done in Inkscape with file saved as `.svg`.

For submission to **PLOS Genetics**, only "eps" or "tiff" formats are acceptable. To convert the figures from svg to tiff, I used the following code
in my macOS Catalina Version 10.15.7 after installing **imagemagick** using `brew install imagemagick`. 


```
#!/bin/sh

# Convert all arguments (assumed SVG) to a TIFF acceptable to PLOS
# Requires Inkscape and ImageMagick 6.8 (doesn't work with 6.6.9)

for i in $@; do
  BN=$(basename $i .svg)
  /Applications/Inkscape.app/Contents/MacOS/inkscape --without-gui --export-png="$BN.png" --export-dpi 300 $i
  convert -compress LZW -alpha remove $BN.png $BN.tiff
  mogrify -alpha off $BN.tiff
  rm $BN.png
done
```


