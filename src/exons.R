library(R6)

# Exons class: stores exon definition and provides them in multiple formats
Exons <- R6Class("Exons",
                  public = list(
                    name = NULL,
                    exonsDF = NULL,
                    
                    initialize = function(name, path) {
                      self$name <- name
                      self$exonsDF <- data.table::fread(path, select = 1:6, sep = "\t", header = T)
                    },
                    
                    # Returns canonical exons in GRanges format
                    getCanonicalExonsGR = function(){
                      canonicalExonsDF <- self$exonsDF[!grep(pattern = "_", x = self$exonsDF$`exon#`)]
                      return(regioneR::toGRanges(canonicalExonsDF[, c(1:5)]))  # we omit strand column to avoid compatibility problems
                    },
                    
                    
                    # Returns all exons definitions in GRanges format
                    getAllExonsGR = function(){
                      return(regioneR::toGRanges(self$exonsDF[, c(1:5)])) # we omit strand column to avoid compatibility problems
                    },
                     
                    
                    # Returns "aberrant gaps": those regions covered by a canonical exon (Ex: exon 1) that are not covered by aberrant exon (Ex: exon 1_dE1q6)
                    # DISCLAIMER: Does not reports "gaps" when the aberrant exon is greater than the canonical exon. It's not an issue with current BRCA1/BRCA2 exon definitions
                    getAberrantGaps = function(){
                      
                      # Get aberrant exons (those with _, i.e. 4_dE4q22)
                      aberrantExonsDF <- self$exonsDF[grep(pattern = "_", x = self$exonsDF$`exon#`)]
                    
                      # Get canonical exons
                      canonicalGR <- self$getCanonicalExonsGR()
                      
                      # Determine gaps
                      aberrantGaps <- data.frame()
                      for (i in seq_len(nrow(aberrantExonsDF))){
                        eName <- aberrantExonsDF$`exon#`[i]
                        parts <- strsplit(eName, "_")[[1]]
                        eNumber <- parts[1]
                        if (!startsWith(parts[2], "n")){ # disclaimer: ignore _n events
                          eGR <- regioneR::toGRanges(aberrantExonsDF[i])
                          gap <- setdiff(canonicalGR[canonicalGR$`exon#` == eNumber], eGR)
                          gap <- as.data.frame(gap)
                          gap$name <- eName
                          aberrantGaps <- rbind(aberrantGaps, as.data.frame(gap))
                        }
                      }
                      
                      return(aberrantGaps)
                    }
                    
                  ))
