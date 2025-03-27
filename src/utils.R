
log <- function(text, log_file){
  text <- paste0("\n[", Sys.time(),"]  ", text)
  write(text, file = log_file, append = TRUE)
}

findOverapsAndSave <- function(grA, grB, save_file){
  overlaps <- GenomicRanges::findOverlaps(grA, grB)
  minCorrDF <- auxDF[queryHits(overlaps), ]
  write.table(minCorrDF, file = save_file, sep = "\t", row.names = F, quote = F, col.names = F)
}


# Returns mix-max range
minMaxRange <- function(gr){
  minPosition <- min(start(gr))
  maxPosition <- max(end(gr))
  limitRange <- GRanges(
    seqnames = Rle(seqnames(gr[1])),
    ranges = IRanges(start = minPosition, end = maxPosition)
  )
  return(limitRange)
}


# isoformStr: String defining isoform like "1_nE1q536, 2, 3, 4, 5, 6, 10, 11, 43082404-43082572, 13"
# exonsGR: all gene exons
# canonicalExonsGR: all canonical exons
# aberrantGaps: DF with known gaps coming from alternative exons like 4_dE4q22
# minMaxGR: Granges defining start and end of the region where abonormal events will be examined
# lengthTreshold: ignore events equal or lower than lengthTreshold, except for known deletions like "7_dE7p3"
# returns object containing:
#      events: abonormal events identified when compared the isoformStr against the canonical exons
#      GR: isoform in GR format 
getAbnormalities <- function(isoformStr, exonsGR, canonicalExonsGR, aberrantGaps,  minMaxGR = NULL, lengthTreshold = 5){
  
  # Get isoform array
  isoformStr <- strsplit(isoformStr, ", ")[[1]]
  
  # Define vars for later use
  exonsDF <- as.data.frame(exonsGR)
  
  # transform isoform into GRanges
  isoformRegionsGR <- isoformToGR(isoform = isoformStr, exonsGR)
  
  # Calculate min-max for later use
  if (is.null(minMaxGR)){
    minMax <- minMaxRange(isoformRegionsGR) 
  } else {
    minMax <- minMaxGR  
  }
  
  
  # Which events are present in isoform but not in the canonical form?
  
  # limit canonical regions to the range where the isoform exists
  auxCanonicalGR <- subsetByOverlaps(canonicalExonsGR,  minMax)
  
  # Get diffs
  insEventsDF <- as.data.frame(setdiff(isoformRegionsGR, auxCanonicalGR))
  delEventsGR <- setdiff(auxCanonicalGR, isoformRegionsGR)
  delEventsDF <- as.data.frame(delEventsGR)
  
  # Assign known labels if it is an known aberrant gap
  if (nrow(delEventsDF) > 0){
    delEventsDF$name <- NA  
  }
  for (i in seq_len(nrow(delEventsDF))){
    res <- aberrantGaps[aberrantGaps$start == delEventsDF$start[i] & aberrantGaps$end == delEventsDF$end[i], "name"]
    if (length(res) > 0){
      delEventsDF$name[i] <- aberrantGaps[aberrantGaps$start == delEventsDF$start[i] & aberrantGaps$end == delEventsDF$end[i], "name"]  
    }
  }
  
  
  # ignore events smaller or equal than threshold
  # Exception: those dels with known name (aberrant gaps)
  ignoredIns <- insEventsDF[insEventsDF$width <= lengthTreshold,] # to report
  ignoredDels <- delEventsDF[delEventsDF$width <= lengthTreshold & is.na(delEventsDF$name),] # to report
  insEventsDF <- insEventsDF[insEventsDF$width > lengthTreshold,]
  delEventsDF <- delEventsDF[delEventsDF$width > lengthTreshold | !is.na(delEventsDF$name),]
  
  
  # Build final composed event
  finalEvent <- ""
  
  # process ins
  for (i in seq_len(nrow(insEventsDF))){
    eventDef <- paste0("ins_", insEventsDF$seqnames[i], ":", insEventsDF$start[i], "-", insEventsDF$end[i])
    finalEvent <- paste(finalEvent, eventDef)
  }
  
  
  # process dels
  for (i in seq_len(nrow(delEventsDF))){
    
    # Find official exon number if the region completely spans 1 or more exons
    res <- exonsDF[exonsDF$start == delEventsDF$start[i] & exonsDF$end == delEventsDF$end[i], ]
    if(nrow(res) > 1){
      stop("unexpected")
    }
    
    if (nrow(res == 1)){  # it's a complete deletion
      
      eventDef <- paste0("del_", res$exon.[1])  
      
    } else { # it's a partial deletion
      
      # Does the deletion matches official exon gap?
      res <- aberrantGaps[aberrantGaps$start == delEventsDF$start[i] & aberrantGaps$end == delEventsDF$end[i], ]
      if (nrow(res == 1)){
        eventDef <- res$name
      } else {
        eventDef <- paste0("del_", delEventsDF$seqnames[i], ":", delEventsDF$start[i], "-", delEventsDF$end[i])
      }
    }
    
    
    finalEvent <- paste(finalEvent, eventDef)
  }
  
  
  finalEvent <- substring(finalEvent, 2) 
  
  return(list(
    events = finalEvent,
    GR = isoformRegionsGR,
    ignoredIns = ignoredIns,
    ignoredDels = ignoredDels
  ))  
}


# Give a GR of two non-adjacent regions, returns inner range
getInnerRange <- function(gr1, gr2){
  gr_ordered <- sort(c(gr1, gr2))
  
  # Calculate the "inner range"
  inner_range <- GRanges(
    seqnames = seqnames(gr_ordered)[1], # Use the same seqname
    ranges = IRanges(
      start = end(gr_ordered)[1] + 1,   # Start right after the first range
      end = start(gr_ordered)[2] - 1    # End right before the second range
    )
  )
  return(inner_range)
}


# transform isoform into GRanges
# isoform: isoform array like 2 3 4 5 43079334-43079399
# exonsGR: exons reference in GR format
isoformToGR <- function(isoform, exonsGR){
  chr <- seqnames(exonsGR)[1]
  isoformRegionsGR <- exonsGR[exonsGR$`exon#` %in% isoform, ]
  
  aberrantEvents <- setdiff(isoform, exonsGR$`exon#`)
  # Convert aberrantEvents to GRanges in a vectorized manner
  if (length(aberrantEvents) > 0) {
    splitRanges <- do.call(rbind, strsplit(aberrantEvents, ":|-"))
    aberrantGR <- GRanges(
      seqnames = splitRanges[, 1],
      ranges = IRanges(start = as.numeric(splitRanges[, 2]), end = as.numeric(splitRanges[, 3]))
    )
    
    # Combine GRanges
    isoformRegionsGR <- suppressWarnings(c(isoformRegionsGR, aberrantGR))
  }
  
  
  return(isoformRegionsGR)
}


# Converts isoform from GRanges to str like 1, 2, 3, 17:43079334-43079399
fromIsoformGRtoIsoformStr <- function(isoGR){
  
  # Compose isoform getting official exon names
  isoformStr <- paste(isoGR$`exon#`[!is.na(isoGR$`exon#`)], collapse=", ")
  
  # Add unknown regions
  abDF <- as.data.frame(isoGR[is.na(isoGR$`exon#`)], row.names = NULL)
  for (k in seq_len(nrow(abDF))){
    abEv <- paste0(abDF$seqnames[k], ":", abDF$start[k], "-", abDF$end[k])
    if (isoformStr == ""){
      isoformStr <- abEv
    } else {
      isoformStr <- paste(isoformStr, abEv, sep = ", ")
    }
  }
  return(isoformStr)
}


# Sorts event summary str, decreasing, using pct
# event: Example "Cont9 (473, 4.87%)  |  Cont10 (2810, 6.86%)" or "ins_17:43079334-43079399 10_d10q del_9 del_8 (91, 2.7%)"
sortEventStr <- function(eventSummaryStr){
  
  items <- unlist(strsplit(eventSummaryStr, "  \\|  "))
  
  # Extract the percentage values (second number inside parentheses)
  percentage_values <- as.numeric(gsub(".*\\(\\d+, ([0-9\\.]+)%\\).*", "\\1", items))
  
  # Sort items in decreasing order based on the percentage
  sorted_items <- items[order(percentage_values, decreasing = TRUE)]
  
  # Reconstruct the sorted string
  sorted_str <- paste(sorted_items, collapse = "  |  ")
  return(sorted_str)
}


auxFormatRows <- function(df) {
  apply(df, 1, function(row) paste0(row["seqnames"], ":", row["start"], "-", row["end"]))
}
