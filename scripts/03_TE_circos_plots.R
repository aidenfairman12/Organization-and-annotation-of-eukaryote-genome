#!/usr/bin/env Rscript

# TE Distribution Circos Plots
# This script generates circular plots showing the distribution of transposable elements
# across the genome using the circlize package.
# 
# Inputs:
# - assembly.fasta.mod.EDTA.TEanno.gff3: TE annotations from EDTA
# - assembly.fasta.fai: Fasta index file with scaffold lengths
# - assembly.fasta.mod.LTR.intact.raw.fa.rexdb-plant.cls.tsv: TEsorter clade classifications
#
# Outputs:
# - 03_TE_density_all_superfamilies.pdf: Circos plot with all major TE superfamilies
# - 03_TE_density_centromeric.pdf: Circos plot highlighting Athila and CRM clades

library(circlize)
library(grid)  # For legend drawing

# Helper functions to replace tidyverse functionality
filter <- function(df, condition) {
    df[condition, ]
}

select <- function(df, ...) {
    cols <- as.character(substitute(list(...)))[-1]
    df[, cols, drop = FALSE]
}

mutate <- function(df, ...) {
    args <- list(...)
    for (name in names(args)) {
        df[[name]] <- eval(args[[name]], envir = df, enclos = parent.frame())
    }
    df
}

# Set working directory to EDTA output folder
setwd("/data/users/afairman/euk_org_ann/results/edta_output")

#-------------------------------------------------
# 1. Load the TE annotation GFF3 file
#-------------------------------------------------
message("Loading TE annotation GFF3 file...")
gff_file <- "assembly.fasta.mod.EDTA.TEanno.gff3"
gff_data <- read.table(gff_file, header = FALSE, sep = "\t", 
                       stringsAsFactors = FALSE, comment.char = "#")

# Assign column names
colnames(gff_data) <- c("seqid", "source", "type", "start", "end", 
                        "score", "strand", "phase", "attributes")

# Check the superfamilies present in the GFF3 file
message("\nTE superfamily counts in GFF3:")
superfamily_counts <- table(gff_data$type)
print(sort(superfamily_counts, decreasing = TRUE))

#-------------------------------------------------
# 2. Create custom ideogram from FAI file
#-------------------------------------------------
message("\nCreating custom ideogram from scaffold lengths...")
custom_ideogram <- read.table("assembly.fasta.fai", header = FALSE, 
                              stringsAsFactors = FALSE)
custom_ideogram$chr <- custom_ideogram$V1
custom_ideogram$start <- 1
custom_ideogram$end <- custom_ideogram$V2
custom_ideogram <- custom_ideogram[, c("chr", "start", "end")]

# Sort by length (descending)
custom_ideogram <- custom_ideogram[order(custom_ideogram$end, decreasing = TRUE), ]

# Report total length of top scaffolds
message("Total length of top 30 scaffolds: ", 
        format(sum(custom_ideogram$end[1:30]), big.mark = ","), " bp")

# Select only the first 30 longest scaffolds for visualization
# This includes scaffolds with centromeric TEs (Athila/CRM)
custom_ideogram <- custom_ideogram[1:30, ]

message("Selected scaffolds for circos plot:")
print(custom_ideogram)

#-------------------------------------------------
# 2b. Find genes GFF (optional) and load gene coordinates
#-------------------------------------------------
message("\nSearching for gene annotation GFF files...")
gene_candidates <- c(
    file.path("/data/users/afairman/euk_org_ann/results/Maker/final/filtered.genes.renamed2.gff3"),
    file.path("/data/users/afairman/filtered.genes.renamed2.gff3"),
    file.path("..","..","filtered.genes.renamed2.gff3"),
    file.path("..","results","Maker","final","filtered.genes.renamed2.gff3"),
    file.path("..","..","..","filtered.genes.renamed2.gff3")
)
gene_file <- NULL
for (g in gene_candidates) {
    if (file.exists(g)) {
        # skip empty files
        finfo <- file.info(g)
        if (!is.na(finfo$size) && finfo$size > 0) {
            gene_file <- g
            break
        } else {
            message("Found gene GFF but it's empty, skipping: ", g)
        }
    }
}
if (!is.null(gene_file)) {
    message("Found gene GFF: ", gene_file)
    gene_gff <- tryCatch(
        read.table(gene_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#"),
        error = function(e) {
            message("Error reading gene GFF: ", conditionMessage(e))
            return(NULL)
        }
    )
    if (!is.null(gene_gff) && nrow(gene_gff) > 0) {
        if (ncol(gene_gff) >= 9) {
            colnames(gene_gff)[1:9] <- c("seqid","source","type","start","end","score","strand","phase","attributes")
        }
        # keep features that look like genes
        gene_rows <- gene_gff[grepl("gene|mRNA|transcript", gene_gff[,3], ignore.case = TRUE), ]
        if (nrow(gene_rows) > 0) {
            gene_data <- data.frame(chrom = gene_rows[,1], start = as.numeric(gene_rows[,4]), end = as.numeric(gene_rows[,5]), strand = gene_rows[,7], stringsAsFactors = FALSE)
            gene_data <- gene_data[gene_data$chrom %in% custom_ideogram$chr, ]
            message("  - Genes on selected scaffolds: ", nrow(gene_data))
        } else {
            gene_data <- data.frame(chrom = character(), start = numeric(), end = numeric(), strand = character(), stringsAsFactors = FALSE)
            message("  - No gene features found in GFF (gene/mRNA/transcript)")
        }
    } else {
        message("Failed to parse gene GFF; skipping gene track")
        gene_data <- data.frame(chrom = character(), start = numeric(), end = numeric(), strand = character(), stringsAsFactors = FALSE)
    }
} else {
    message("No non-empty gene GFF found in candidate locations; skipping gene track")
    gene_data <- data.frame(chrom = character(), start = numeric(), end = numeric(), strand = character(), stringsAsFactors = FALSE)
}

#-------------------------------------------------
# 3. Function to filter GFF3 data by superfamily
#-------------------------------------------------
filter_superfamily <- function(gff_data, superfamily, custom_ideogram) {
    # Filter by type
    filtered_data <- gff_data[gff_data$type == superfamily, ]
    # Filter by chromosome
    filtered_data <- filtered_data[filtered_data$seqid %in% custom_ideogram$chr, ]
    # Select and rename columns
    result <- data.frame(
        chrom = filtered_data$seqid,
        start = filtered_data$start,
        end = filtered_data$end,
        strand = filtered_data$strand,
        stringsAsFactors = FALSE
    )
    
    message("  - ", superfamily, ": ", nrow(result), " elements")
    return(result)
}

#-------------------------------------------------
# 4. Identify most abundant TE superfamilies
#-------------------------------------------------
message("\nIdentifying most abundant TE superfamilies...")

# Get top superfamilies (excluding generic categories and Helitrons)
exclude_types <- c("repeat_region", "target_site_duplication", "long_terminal_repeat", "Helitron")
te_data_filtered <- gff_data[!gff_data$type %in% exclude_types & 
                              gff_data$seqid %in% custom_ideogram$chr, ]
te_counts <- as.data.frame(table(te_data_filtered$type), stringsAsFactors = FALSE)
colnames(te_counts) <- c("type", "n")
te_counts <- te_counts[order(te_counts$n, decreasing = TRUE), ]

message("\nTop TE superfamilies (on selected scaffolds):")
print(te_counts)

# Select top superfamilies for plotting (adjust as needed)
# Limit to top 8 to ensure they fit in the plot
top_superfamilies <- head(te_counts$type[te_counts$n >= 400], 8)

message("\nSuperfamilies selected for plotting: ", 
        paste(top_superfamilies, collapse = ", "))

#-------------------------------------------------
# 5. Define color palette for superfamilies
#-------------------------------------------------
# Create a color palette that distinguishes between TE types
te_colors <- c(
    "Gypsy_LTR_retrotransposon" = "#d95f02",      # orange
    "Copia_LTR_retrotransposon" = "#1b9e77",      # teal
    "Ty3_LTR_retrotransposon" = "#7570b3",        # purple
    "LINE_element" = "#e7298a",                    # magenta
    "SINE_element" = "#66a61e",                    # green
    "TIR_transposon" = "#e6ab02",                  # gold
    "DNA_transposon" = "#a6761d",                  # brown
    "Helitron" = "#666666",                        # gray
    "Mutator_TIR_transposon" = "#377eb8",         # blue
    "hAT_TIR_transposon" = "#4daf4a",             # light green
    "CACTA_TIR_transposon" = "#984ea3",           # dark purple
    "PIF_Harbinger_TIR_transposon" = "#ff7f00",   # orange-red
    "Tc1_Mariner_TIR_transposon" = "#ffff33"      # yellow
)

# Assign colors to selected superfamilies
plot_colors <- te_colors[top_superfamilies]
# For any superfamily without a predefined color, assign from rainbow palette
missing_colors <- is.na(plot_colors)
if (any(missing_colors)) {
    n_missing <- sum(missing_colors)
    plot_colors[missing_colors] <- rainbow(n_missing)
}

#-------------------------------------------------
# 6. PLOT 1: All abundant TE superfamilies
#-------------------------------------------------
message("\n========================================")
message("Creating circos plot with all TE superfamilies...")
message("========================================\n")

# Create plots directory if it doesn't exist
plots_dir <- "/data/users/afairman/euk_org_ann/plots"
if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
}

pdf(file.path(plots_dir, "03_TE_density_all_superfamilies.pdf"), width = 12, height = 12)

# Set circos parameters
gaps <- c(rep(1, length(custom_ideogram$chr) - 1), 5)
circos.par(start.degree = 90, gap.after = 1, track.margin = c(0.01, 0.01), 
           gap.degree = gaps, cell.padding = c(0, 0, 0, 0))

# Initialize the circos plot
circos.genomicInitialize(custom_ideogram)

# Plot TE density for each superfamily
message("Plotting TE density for each superfamily:")
for (i in seq_along(top_superfamilies)) {
    sf <- top_superfamilies[i]
    color <- plot_colors[i]
    
    te_data <- filter_superfamily(gff_data, sf, custom_ideogram)
    
    if (nrow(te_data) > 0) {
        circos.genomicDensity(te_data, 
                             count_by = "number", 
                             col = color, 
                             track.height = 0.08,  # Slightly larger tracks
                             window.size = 1e5)
    }
}

# Add gene density track (if gene_data present)
if (exists("gene_data") && nrow(gene_data) > 0) {
    message("Adding gene density track to circos plot")
    circos.genomicDensity(gene_data, count_by = "number", col = "#000000", track.height = 0.04, window.size = 1e5)
}

# Add scaffold labels on a new outer track
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[2], 
               CELL_META$sector.index, 
               facing = "clockwise", niceFacing = TRUE,
               adj = c(0, 0.5), cex = 0.5)
}, track.height = 0.03, bg.border = NA)

circos.clear()

# Add legend using base R (include Genes if present)
legend_entries <- names(plot_colors)
legend_colors <- as.character(plot_colors)
if (exists("gene_data") && nrow(gene_data) > 0) {
    legend_entries <- c(legend_entries, "Genes")
    legend_colors <- c(legend_colors, "#000000")
}
legend("center", legend = legend_entries, 
    fill = legend_colors, 
    title = "TE Superfamily",
    ncol = 2, cex = 0.8, bty = "n")

dev.off()

message("\n✓ Plot saved: ", file.path(plots_dir, "03_TE_density_all_superfamilies.pdf"))

#-------------------------------------------------
# 7. Load TEsorter classification for clade-level analysis
#-------------------------------------------------
message("\n========================================")
message("Loading TEsorter clade classifications...")
message("========================================\n")

cls_file <- "/data/users/afairman/euk_org_ann/results/edta_output/assembly.fasta.mod.EDTA.raw/assembly.fasta.mod.LTR.intact.raw.fa.rexdb-plant.cls.tsv"

if (file.exists(cls_file)) {
    # Read classification file
    # Note: First line starts with '#TE' which is the header, not a comment
    cls_data <- read.table(cls_file, sep = "\t", header = TRUE, 
                          stringsAsFactors = FALSE, comment.char = "")
    
    # Rename first column from #TE to TE_id
    colnames(cls_data)[1] <- "TE_id"
    
    message("TEsorter classification loaded successfully")
    message("Number of classified TEs: ", nrow(cls_data))
    
    # Check clade distribution
    if ("Clade" %in% colnames(cls_data)) {
        message("\nClade distribution:")
        print(table(cls_data$Clade))
        
        # (previously extracted Athila/Ale clades here; removed per user request)
        
        #-------------------------------------------------
        # 8. Parse coordinates from TE IDs and create GFF-like data
        #-------------------------------------------------
        # TEsorter TE IDs are formatted like: contig_1:123456..789012#LTR/Gypsy
        parse_te_coords <- function(te_id) {
            # Extract scaffold:start..end using regex
            if (grepl("^(.+):(\\d+)\\.\\.(\\d+)", te_id)) {
                # Remove everything after # if present
                te_id_clean <- sub("#.*$", "", te_id)
                # Split by colon
                parts <- strsplit(te_id_clean, ":")[[1]]
                if (length(parts) == 2) {
                    seqid <- parts[1]
                    # Split coordinates by ..
                    coords <- strsplit(parts[2], "\\.\\.")[[1]]
                    if (length(coords) == 2) {
                        coord1 <- as.numeric(coords[1])
                        coord2 <- as.numeric(coords[2])
                        # Ensure start < end (coordinates might be reversed)
                        return(data.frame(
                            seqid = seqid,
                            start = min(coord1, coord2),
                            end = max(coord1, coord2),
                            stringsAsFactors = FALSE
                        ))
                    }
                }
            }
            return(NULL)
        }
        
        # Athila/Ale extraction removed; centromeric plot will be generated based on presence of Gypsy/Copia data
        
        #-------------------------------------------------
        # 9. PLOT 2: Centromeric TEs (context plot with Gypsy/Copia; genes on outer track)
        # Ensure gypsy/copia data frames are available for the check and plotting
        gypsy_data <- filter_superfamily(gff_data, "Gypsy_LTR_retrotransposon", custom_ideogram)
        copia_data <- filter_superfamily(gff_data, "Copia_LTR_retrotransposon", custom_ideogram)

        if (nrow(gypsy_data) > 0 || nrow(copia_data) > 0) {
            message("\n========================================")
            message("Creating circos plot for centromeric TEs...")
            message("========================================\n")
            
            pdf(file.path(plots_dir, "03_TE_density_centromeric3.pdf"), width = 12, height = 12)
            
            # Set circos parameters
            circos.par(start.degree = 90, gap.after = 1, track.margin = c(0.01, 0.01), 
                      gap.degree = gaps, cell.padding = c(0, 0, 0, 0))
            
            # Initialize the circos plot
            circos.genomicInitialize(custom_ideogram)
            
            # Plot major LTR retrotransposon superfamilies for context
            gypsy_data <- filter_superfamily(gff_data, "Gypsy_LTR_retrotransposon", custom_ideogram)
            copia_data <- filter_superfamily(gff_data, "Copia_LTR_retrotransposon", custom_ideogram)
            
            if (nrow(gypsy_data) > 0) {
                circos.genomicDensity(gypsy_data, count_by = "number", 
                                     col = "#d95f02", track.height = 0.08, 
                                     window.size = 1e5)
            }
            
            if (nrow(copia_data) > 0) {
                circos.genomicDensity(copia_data, count_by = "number", 
                                     col = "#1b9e77", track.height = 0.08, 
                                     window.size = 1e5)
            }
            
            # Add scaffold labels on a new inner track (labels placed inside gene track)
            circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
                circos.text(CELL_META$xcenter, CELL_META$ylim[2], 
                           CELL_META$sector.index, 
                           facing = "clockwise", niceFacing = TRUE,
                           adj = c(0, 0.5), cex = 0.5)
            }, track.height = 0.03, bg.border = NA)

            # Plot gene positions as an outermost point track (if gene_data present)
            if (exists("gene_data") && nrow(gene_data) > 0) {
                circos.genomicTrackPlotRegion(gene_data, track.height = 0.03, ylim = c(0,1), panel.fun = function(region, value, ...) {
                    if (nrow(region) > 0) {
                        centers <- (region$start + region$end) / 2
                        circos.points(centers, rep(0.5, length(centers)), pch = 16, cex = 0.3, col = "#000000")
                    }
                })
            }
            
            circos.clear()
            
            # Add legend using base R
                # Build legend entries for centromeric plot (Athila/Ale removed; Genes shown if present)
                cen_legend <- c("Gypsy (all)", "Copia (all)")
                cen_colors <- c("#d95f02", "#1b9e77")
                if (exists("gene_data") && nrow(gene_data) > 0) {
                    cen_legend <- c(cen_legend, "Genes")
                    cen_colors <- c(cen_colors, "#000000")
                }
                legend("center", legend = cen_legend, fill = cen_colors, title = "TE Type", cex = 0.9, bty = "n")
            
            dev.off()
            
            message("\n✓ Plot saved: ", file.path(plots_dir, "03_TE_density_centromeric3.pdf"))
        } else {
            message("\nWARNING: No Athila or CRM elements found on selected scaffolds")
            message("Centromeric TE plot will not be generated")
        }
        
    } else {
        message("WARNING: 'Clade' column not found in TEsorter output")
    }
    
} else {
    message("WARNING: TEsorter classification file not found: ", cls_file)
    message("Skipping centromeric TE analysis")
    message("Run TEsorter first with: sbatch 02_TEsorter_clade_classification.sh")
}

#-------------------------------------------------
# 10. Summary statistics
#-------------------------------------------------
message("\n========================================")
message("SUMMARY")
message("========================================")
message("Total scaffolds analyzed: ", nrow(custom_ideogram))
message("Total genome length analyzed: ", 
        format(sum(custom_ideogram$end), big.mark = ","), " bp")
message("Number of TE superfamilies plotted: ", length(top_superfamilies))
message("\nPlots generated:")
message("  1. ", file.path(plots_dir, "03_TE_density_all_superfamilies.pdf"))
if (exists("gypsy_data") && exists("copia_data")) {
    if (nrow(gypsy_data) > 0 || nrow(copia_data) > 0) {
        message("  2. ", file.path(plots_dir, "03_TE_density_centromeric.pdf"))
    }
}
message("\n✓ Analysis complete!")
