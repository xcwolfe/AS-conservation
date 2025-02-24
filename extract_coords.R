library("tidyverse")

# Extract start and end coordinates of each IR event:
## intron_sums_heatmap is the uniqueness index of intron retention events in C. elegans. This heatmap was published (alongisde similar heatmaps for other alternative splicing types) in the manuscript "Deep transcriptomics reveals cell-specific isoforms of pan-neuronal genes" (Wolfe et al. 2025) 
extracted_coords <- lapply(intron_sums_heatmap$AS_event_ID, function(id) {
  parts <- unlist(strsplit(id, "_"))
  chr <- as.character(parts[1])
  start <- as.numeric(parts[3])
  end <- as.numeric(parts[4])
  return(c(chr, start, end))
})

extracted_coords <- do.call(rbind, extracted_coords)
rownames(extracted_coords) <- intron_sums_heatmap$AS_event_ID
write.csv(extracted_coords, file = "extracted_coords_IR_events.csv")   ### move this file to HPCC and run average_phylop_score_for_IR_events.txt
