# Windows bases and genomic features

Many analyses in this paper are done in sliding window. Each sliding window summarises information in the region.
Choosing a size for the sliding window is not a trivial task. Smaller window sizes add noise and larger ones, migth
reduce the resolution in the data. Here we choose the following window sizes: 100 Kb,200 Kb and 1M. 

## Extract genomic coordinates for CDS, intron and intergenic regions
`bash extract_ostrich_coordinates.sh`

# Create window files
`bash generate_windows.sh`

