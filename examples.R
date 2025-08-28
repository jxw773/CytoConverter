#!/usr/bin/env Rscript

#' CytoConverter Usage Examples
#' 
#' This script demonstrates various ways to use CytoConverter for converting
#' cytogenetic nomenclature to genomic coordinates.

# Load required modules
if (!require("modules", quietly = TRUE)) {
    install.packages("modules")
    library(modules)
}

# Source the main CytoConverter function
source("modules/cytoconverter.R")

cat("=== CytoConverter Usage Examples ===\n\n")

# Example 1: Simple trisomy 21
cat("Example 1: Simple trisomy 21\n")
cat("Input: '47,XY,+21'\n")
result1 <- CytoConverter("47,XY,+21")
cat("Results:\n")
if (nrow(result1$Results) > 0) {
    print(result1$Results)
} else {
    cat("No results generated\n")
}
if (nrow(result1$Error_log) > 0) {
    cat("Errors/Warnings:\n")
    print(result1$Error_log)
}
cat("\n")

# Example 2: Deletion
cat("Example 2: Deletion of chromosome 5q\n")
cat("Input: '46,XY,del(5q13q33)'\n")
result2 <- CytoConverter("46,XY,del(5q13q33)")
cat("Results:\n")
if (nrow(result2$Results) > 0) {
    print(result2$Results)
} else {
    cat("No results generated\n")
}
if (nrow(result2$Error_log) > 0) {
    cat("Errors/Warnings:\n")
    print(result2$Error_log)
}
cat("\n")

# Example 3: Multiple samples using data frame
cat("Example 3: Multiple samples using data frame\n")
samples <- data.frame(
    Sample = c("Patient1", "Patient2", "Patient3"),
    Karyotype = c(
        "46,XX,del(5q13q33)",
        "47,XY,+21", 
        "46,XY,t(9;22)(q34;q11)"
    ),
    stringsAsFactors = FALSE
)
cat("Input data frame:\n")
print(samples)
cat("\n")

result3 <- CytoConverter(samples, build = "GRCh38")
cat("Results:\n")
if (nrow(result3$Results) > 0) {
    print(result3$Results)
} else {
    cat("No results generated\n")
}
if (nrow(result3$Error_log) > 0) {
    cat("Errors/Warnings:\n")
    print(result3$Error_log)
}
cat("\n")

# Example 4: Using different genome build
cat("Example 4: Using hg19 genome build\n")
cat("Input: '46,XY,+21' with build='hg19'\n")
result4 <- CytoConverter("46,XY,+21", build = "hg19")
cat("Results:\n")
if (nrow(result4$Results) > 0) {
    print(result4$Results)
} else {
    cat("No results generated\n")
}
if (nrow(result4$Error_log) > 0) {
    cat("Errors/Warnings:\n")
    print(result4$Error_log)
}
cat("\n")

# Example 5: Error handling
cat("Example 5: Error handling with invalid karyotype\n")
cat("Input: 'invalid_karyotype'\n")
result5 <- CytoConverter("invalid_karyotype")
cat("Results:\n")
if (nrow(result5$Results) > 0) {
    print(result5$Results)
} else {
    cat("No results generated\n")
}
if (nrow(result5$Error_log) > 0) {
    cat("Errors/Warnings:\n")
    print(result5$Error_log)
}
cat("\n")

cat("=== End of Examples ===\n")