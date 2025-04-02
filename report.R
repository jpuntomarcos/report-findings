# Reports events within the defined ranges to explore

# CHECK INPUT ---

args <- commandArgs(trailingOnly = TRUE)

# Exactly one argument is provided
if (length(args) != 1) {
  stop("Error: Exactly one parameter (a .yaml file) must be provided")
}


# Check if the argument is a .yaml file
if (!grepl("\\.yaml$", args, ignore.case = TRUE)) {
  stop("Error: The provided parameter is not a .yaml file.")
}

paramsFileName <- args[1]


# INITIAL VARS ---

# load libs
suppressPackageStartupMessages(library(yaml))
library(GenomicRanges)
source("src/exons.R")
source("src/utils.R")


# load vars from yaml file
params <- yaml.load_file(paramsFileName)

# other vars
outputDir <- params$output_dir
dir.create(outputDir, recursive = T)


# SCRIPT ----

start.time <- Sys.time()

# Get ranges to evaluate for later use
rangesDF <- read.table(params$ranges_to_evaluate, header = T, sep = "\t", stringsAsFactors = F)

# Get counts (Ulf's) paths
countPathsDF <- read.table(params$ulf_count_paths, header = T, sep = "\t", stringsAsFactors = F)
genes <- unique(countPathsDF$gene) # for later use
samplesIds <- unique(countPathsDF$sample) # for later use

# Read events to evaluate. It will be the final df to report
evaluateDF <- read.table(params$events_to_evaluate, header = T, sep = "\t", stringsAsFactors = F)
evaluateDF$FoundInSample <- NA    
evaluateDF$FoundInControl <- NA   
evaluateDF$FoundInMut <- NA
evaluateDF$FIS_n <- 0
evaluateDF$FIS_pct <- 0
evaluateDF$FIS_isoforms <- ""
evaluateDF$FIS_d <- 0
evaluateDF$FIS_d_iso <- ""
evaluateDF$FIC_max <- 0
evaluateDF$FIM_max <- 0
evaluateDF$FIC <- ""
evaluateDF$FIM <- ""
evaluateDF$SampleFindings <- ""

# Create log data.frame to report isoforms supporting "n" for which some events where ignored due to size threshold
ignoredLogDF <- evaluateDF[,c(1,2)]
ignoredLogDF$FIS_n.isoformsWithIgnoredEvents <- ""
ignoredLogDF$count <- 0


for (gene in genes){
  
  # Build exons-data object and GR format
  exons <- Exons$new(gene,  params[[paste0(gene, "_exons")]])  
  exonsGR <- exons$getAllExonsGR()
  geneGR <- GRanges(seqnames = unique(seqnames(exonsGR)),
                    ranges = IRanges(start = min(start(exonsGR)),
                                     end = max(end(exonsGR))))
  canonicalExonsGR <- exons$getCanonicalExonsGR()
  aberrantGaps <- exons$getAberrantGaps()
  
  # Read all unique isoforms from samples
  allIsos <- NULL
  samples = list()  # list of samples. Keys: sample names. Values: data.frame containing counts (Ulf's 2nd script output)
  for (s in samplesIds){
    
    f <- countPathsDF[countPathsDF$sample == s & countPathsDF$gene == gene, "path"]
    isoDF <- read.table(f, sep = "\t", stringsAsFactors = F, header = F, col.names = c("id", "count", "regions"))
    
    # remove isoforms with NA
    idxToRemove <- which(isoDF$regions == "N/A")
    if (length(idxToRemove) > 0){
      cat("\nWarning: N/As were found in", basename(f), "and will not be included in the analysis")  
      isoDF <- isoDF[-idxToRemove,]
    }
    
    # Remove isoforms that do not overlap the gene
    isoGRList <- lapply(strsplit(isoDF$regions, ", "), isoformToGR, exonsGR) # Convert isoforms to GRanges once
    overlapGene <- vapply(isoGRList, function(isoGR) overlapsAny(geneGR, isoGR), logical(1))
    isoDF <- isoDF[which(overlapGene == T), ]  # Keep only those overlapping the gene
    
    # Assign
    samples[[s]] <- isoDF
    allIsos <- c(allIsos, samples[[s]]$regions)
  }
  
  
  # keep only unique isoforms
  allIsos <- unique(allIsos)
  
  # Build isoform dictionary. Items are objects with events and isoform in GR format
  allIsosDict <- sapply(allIsos,  function(isoform){
    getAbnormalities(isoformStr = isoform,
                     exonsGR =  exonsGR,
                     canonicalExonsGR = canonicalExonsGR,
                     aberrantGaps = aberrantGaps, 
                     lengthTreshold = params$events_to_ignore_threshold)
  }, simplify = F)  # Keep as a list
  
  
  
  # 2 - Report ----
  
  
  # Fill report data.frame
  for (i in seq_len(nrow(evaluateDF))){
    
    # Analyze only those corresponding with current gene
    if (evaluateDF$gene[i] != gene){
      next
    }
    
    sName <- evaluateDF$sample[i]
    
    # Split event if necessary
    expectedFullEvent <- evaluateDF$event[i]
    expectedEventItems <- strsplit(expectedFullEvent, "\\ ")[[1]]
    
    # Range for this event
    range <- rangesDF[rangesDF$sample == sName,] 
    startEndExonsGR <- exonsGR[exonsGR$`exon#` %in% c(range$start, range$end)] # start-end exons in the gene-exons file
    
    # GR to be examined: from the end of the 1st exon to the start of the last exon
    limGR <- getInnerRange(startEndExonsGR[1], startEndExonsGR[2])
    
    # Review all samples
    maxFICpct <- 0 # Maximum percentage across control samples
    maxFIMpct <- 0 # Maximum percentage across Mut samples (excluding the target sample)
    additionalFindings <- list() # Keys: unexpected events, values: n_additional_finding
    isosDict <- list()  # dict like allIsosDict, but containing only the observation for the desired range
    for (s in samplesIds){
      
      # Retrieve counts
      isoDF <- samples[[s]]
      
      # calculate n: find all sample isoforms containing the target event
      n <- 0  # isoform quantification
      suppIso <- "" # Concat iso ids supporting n
      d <- 0 # denominator: any isoform meeting cond 1 (start and end should be contained)
      suppIso_d <- ""  # Concat iso ids supporting d
      for (j in seq_len(nrow(isoDF))){
        
        # Convert whole isoform into GR and limit it to the evaluated range
        isoformGR <- allIsosDict[[isoDF$regions[j]]]$GR # try to find it in dict
        if (is.null(isoformGR)){
          stop("Unexpected")
        }
        isoformGrlim <- subsetByOverlaps(isoformGR, limGR)
        
        # Convert to str
        isoformLim <- fromIsoformGRtoIsoformStr(isoGR = isoformGrlim)
        
        # Cond 1: start and end should be at least partially contained
        if (!all(1:length(startEndExonsGR) %in% queryHits(findOverlaps(startEndExonsGR, isoformGR)))){ # both exons should be contained
          next
        }
        
        # If code reachs here: isoform is included in denominator
        d <- d + isoDF$count[j]
        suppIso_d <- paste(suppIso_d, isoDF$id[j])
        
        
        # Cond 2: contains target event and does not contain any new event in the regions
        if (isoformLim == ""){
          isoformKey <- "Empty" # We use "empty" because "" is not allowed
        } else {
          isoformKey <- isoformLim
        }
        isoObj <- isosDict[[isoformKey]]
        sEvents <- isoObj$events # try to find in the dict
        if (is.null(sEvents)){
          isoObj <- getAbnormalities(isoformStr = isoformLim,
                                     minMaxGR = limGR,
                                     exonsGR =  exonsGR,
                                     canonicalExonsGR = canonicalExonsGR,
                                     aberrantGaps = aberrantGaps,
                                     lengthTreshold = params$events_to_ignore_threshold)
          sEvents <- isoObj$events
          isosDict[[isoformKey]] <- isoObj  # save 
        }
        sEvents_parts <- strsplit(sEvents, "\\ ")[[1]]
        
        
        # Skip if cond 2 not satisfied
        if (!setequal(expectedEventItems, sEvents_parts)){
          
          # Before skipping save it if current sample is target sample. It's an additional finding
          if (s == sName & sEvents != "" & !sEvents %in% evaluateDF[evaluateDF$sample == sName, "event"] ){ # If target sample
            n_additionalFinding <- additionalFindings[[sEvents]]
            if (is.null(n_additionalFinding)){
              additionalFindings[[sEvents]] <- isoDF$count[j]
            } else {
              additionalFindings[[sEvents]] <- n_additionalFinding + isoDF$count[j]
            }
          }
          
          # skip
          next
        }
        
        # Success: add n
        n <- n + isoDF$count[j]
        suppIso <- paste(suppIso, isoDF$id[j])
        
        # Report ignored events for target sample
        if (s == sName){
          concat_str <- paste(c(auxFormatRows(isoObj$ignoredIns), auxFormatRows(isoObj$ignoredDels)), collapse = ", ")
          if (nchar(concat_str)> 0){
            ignoredLogDF$FIS_n.isoformsWithIgnoredEvents[i] <- paste(ignoredLogDF$FIS_n.isoformsWithIgnoredEvents[i], paste0("(", isoDF$id[j], ", ", concat_str, ") ")) 
            ignoredLogDF$count[i] <- ignoredLogDF$count[i] + isoDF$count[j]
          }  
        }
        
      }
      
      
      # Report additional findings for target sample if exists
      if (s == sName) {
        if (length(additionalFindings) > 0){
          evaluateDF$SampleFindings[i] <- paste(
            sapply(names(additionalFindings), function(key) {
              value <- additionalFindings[[key]]
              sprintf("%s (%s %.2f%%)", key, paste0(value, ","), (value/d*100.0))
            }),
            collapse = "  |  "
          )
        }
      }
      
      
      # Calculate pct
      pct <- round(n * 100.0 / d , 2) # round %
      
      
      # If iso has some quantification
      if (n > 0){
        
        # Final annotation 
        if (s == sName){ # If target sample
          
          evaluateDF$FoundInSample[i] <- "Yes"
          evaluateDF$FIS_n[i] <- n
          evaluateDF$FIS_pct[i] <- pct
          evaluateDF$FIS_isoforms[i] <- suppIso
          evaluateDF$FIS_d[i] <- d
          evaluateDF$FIS_d_iso[i] <- suppIso_d
          
        } else { # If Control or other Mut
          
          auxStr <- paste0(s, " (", n, ", ", pct, "%)")
          
          if (startsWith(s, "Con")){
            evaluateDF$FoundInControl[i] <- "Yes"
            evaluateDF$FIC[i] <- paste0(evaluateDF$FIC[i], "  |  " , auxStr)
            maxFICpct <- max(maxFICpct, pct)
          } else { # Mut
            evaluateDF$FoundInMut[i] <- "Yes"
            evaluateDF$FIM[i] <- paste0(evaluateDF$FIM[i], "  |  " , auxStr)
            maxFIMpct <- max(maxFIMpct, pct)
          }
        }
      }
      
    }
    
    
    # Set max Con/Mut pcts
    evaluateDF$FIC_max[i] <- maxFICpct
    evaluateDF$FIM_max[i] <- maxFIMpct
    
    # Fix initial chars if necessary
    evaluateDF$FIC[i] <- substring(evaluateDF$FIC[i], 5)
    evaluateDF$FIM[i] <- substring(evaluateDF$FIM[i], 5)
  }
}


# Clean incidental findings: show only 1 row per sample
for (sample in unique(evaluateDF$sample)){
  
  # Get the intersection of the incidental findings for each sample
  events <- evaluateDF[evaluateDF$sample == sample, "SampleFindings"]
  events <- lapply(events, function(x){
    strsplit(x, "  \\|  ")[[1]]
  })
  events <- Reduce(intersect, events)
  
  # Build final str
  eventStr <- paste0(sample, ": ", sortEventStr(paste0(events, collapse = "  |  ")))
  
  # Clean findings column and store new value
  evaluateDF[evaluateDF$sample == sample, "SampleFindings"] <- " "
  evaluateDF[which(evaluateDF$sample == sample)[1], "SampleFindings"] <- eventStr
}


# Sort FIC and FIM values
evaluateDF$FIC <- sapply(evaluateDF$FIC, sortEventStr)
evaluateDF$FIM <- sapply(evaluateDF$FIM, sortEventStr)


# Show total time
end.time <- Sys.time()
cat("\nTotal time: ", round(difftime(end.time, start.time, units = "h"), 2), "h")

# Write results
openxlsx::write.xlsx(evaluateDF, file.path(outputDir, paste0("EvaluatedMutations_", params$output_suffix, ".xlsx")))
openxlsx::write.xlsx(ignoredLogDF, file.path(outputDir, paste0("IgnoredEventsLog_", params$output_suffix, ".xlsx")))


