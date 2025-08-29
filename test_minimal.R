# Test minimal CytoConverter functionality without optional packages
library(stringr)
library(stringi)
library(dplyr)

# Simple test
result <- data.frame(
  Sample_ID = c("ABC_1", "DEF_1"),
  Chr = c("chr1", "chr8"),
  Start = c(0, 0),
  End = c(125000000, 146364022),
  Type = c("Loss", "Gain"),
  Cells_Present = c("unknown", "unknown")
)

print("Test CNV data:")
print(result)