# In this script, we analyze the identity of full-length LTR retrotransposons (LTR-RTs)
# annotated in a genome. We read a GFF file with LTR-RT annotations,
# extract the LTR identity values, merge with classification data from TEsorter,
# and generate plots showing the distribution of LTR identities per clade within
# the Copia and Gypsy superfamilies.

library(tidyverse)
library(data.table)
library(cowplot)

#-------------------------------------------------
# Input files (edit paths if needed)
#-------------------------------------------------
# default input files (EDTA + TEsorter outputs)
gff_file <- "/data/users/afairman/euk_org_ann/results/edta_output/assembly.fasta.mod.EDTA.raw/assembly.fasta.mod.LTR.intact.raw.gff3"
cls_file <- "/data/users/afairman/euk_org_ann/results/edta_output/assembly.fasta.mod.EDTA.raw/assembly.fasta.mod.LTR.intact.raw.fa.rexdb-plant.cls.tsv"
# cls_file is the output from TEsorter on the raw LTR-RT fasta file
#-------------------------------------------------
# Read and preprocess input data
#-------------------------------------------------
message("Reading GFF: ", gff_file)
# fread on this cluster version doesn't accept comment.char in all builds, so strip '#' lines first
gff_lines <- readLines(gff_file)
gff_lines <- gff_lines[!grepl("^\\s*#", gff_lines)]
if(length(gff_lines) == 0) stop("No GFF data found in: ", gff_file)
anno_dt <- fread(text = paste(gff_lines, collapse = "\n"), sep = "\t", header = FALSE, data.table = FALSE)
colnames(anno_dt)[1:9] <- c("seqid","source","type","start","end","score","strand","phase","attributes")

# convert to tibble for tidyverse processing
anno <- as_tibble(anno_dt)

# Remove subfeatures (Terminal repeats, TSDs) so we keep top-level TE annotations
exclude_feats <- c("long_terminal_repeat", "repeat_region", "target_site_duplication")
anno <- anno %>% filter(!type %in% exclude_feats)

# Extract Name and ltr_identity from the ninth column (attributes). This uses regex.

anno <- anno %>%
  rowwise() %>%
  mutate(
    # extract Name=... from attributes
    Name = str_extract(attributes, "(?<=Name=)[^;]+"),
    # extract ltr_identity=... (EDTA stores it 0-1)
    Identity = as.numeric(str_extract(attributes, "(?<=ltr_identity=)[^;]+")),
    # compute length as end - start
    length = as.numeric(end) - as.numeric(start),
    # extract classification from EDTA gff if present (superfamily)
    EDTA_class = str_extract(attributes, "(?<=classification=)[^;]+")
  ) %>%
  # keep only the columns used downstream and rename for clarity
  select(seqid, start, end, type, Name, Identity, length, EDTA_class)

message("Reading classification (TEsorter): ", cls_file)
# Read cls.tsv but strip leading comment lines (fread on this system/version doesn't accept comment.char)
cls_lines <- readLines(cls_file)
cls_lines <- cls_lines[!grepl("^\\s*#", cls_lines)]
if(length(cls_lines) == 0) stop("No data found in cls file: ", cls_file)
cls_df <- fread(text = paste(cls_lines, collapse = "\n"), sep = "\t", header = TRUE)

# Expected columns from TEsorter .cls.tsv: first column is TE id (matching fasta header), other columns contain order/superfamily/clade
# Ensure first column is named TE (if not, set it)
if(ncol(cls_df) >= 1) setnames(cls_df, 1, "TE")

# TEsorter TE IDs in this dataset are formatted like: contig:start..end#LTR/Copia
# Keep TE as-is and parse the Clade column (column named 'Clade' or 4th column). We'll try to be flexible.
if("Clade" %in% colnames(cls_df)) {
  cls_proc <- cls_df %>% select(TE, Clade)
} else if (ncol(cls_df) >= 4) {
  cls_proc <- cls_df %>% select(TE, Clade = 4)
} else {
  stop("Unexpected cls.tsv format; expected a 'Clade' column or >=4 columns")
}

# simplify TE ids: remove any trailing descriptors - keep the Name that matches EDTA 'Name' field
# In our files, TE names look like: contig_1:326079..321038#LTR/Copia
cls_proc <- cls_proc %>% mutate(Name = TE) %>% select(Name, Clade)

## Merge annotation with classification table
# Use a left join so all annotated TEs are kept even if they have no classification match
# join EDTA GFF annotations (by Name) with TEsorter clade assignments
anno_cls <- anno %>% left_join(cls_proc, by = "Name") %>%
  # If TEsorter has no Clade but EDTA_class is available, use EDTA_class as a fallback
  mutate(Clade = ifelse(is.na(Clade) | Clade == "", EDTA_class, Clade),
         Superfamily = case_when(
           str_detect(Clade, "Copia") ~ "Copia",
           str_detect(Clade, "Gypsy") ~ "Gypsy",
           TRUE ~ ifelse(str_detect(EDTA_class, "Copia"), "Copia",
                         ifelse(str_detect(EDTA_class, "Gypsy"), "Gypsy", NA_character_))
         ))

# Quick checks: how many per Superfamily/Clade (may be NA if classification missing)
message("Counts per Superfamily")
print(table(anno_cls$Superfamily, useNA = "ifany"))
message("Counts per Clade")
print(table(anno_cls$Clade, useNA = "ifany"))

#-------------------------------------------------
# Plot setup
#-------------------------------------------------
# binwidth controls histogram resolution around identity values 
binwidth <- 0.005
# x axis limits (EDTA stores identity 0-1). We'll plot as percent (0-100) so adjust accordingly.
xlims <- c(80, 100)

# Compute a single y-max across ALL Copia and Gypsy clades. This ensures consistent y-axis scaling.
## Convert Identity to percent for plotting
anno_cls <- anno_cls %>% mutate(Identity_pct = Identity * 100)

global_ymax <- anno_cls %>%
  filter(Superfamily %in% c("Copia", "Gypsy"), !is.na(Identity_pct)) %>%
  # bin Identity_pct into consistent breaks and count occurrences
  count(Superfamily, Clade, Identity_bin = cut(Identity_pct, seq(xlims[1], xlims[2], by = binwidth*100))) %>%
  pull(n) %>%
  max(na.rm = TRUE)

message("Global y-limit (shared for overview plots): ", global_ymax)

#-------------------------------------------------
# Plot function for one superfamily
#-------------------------------------------------
plot_by_clade <- function(df, sf, ymax) {
  df %>%
    filter(Superfamily == sf, !is.na(Identity_pct), !is.na(Clade)) %>%
    ggplot(aes(x = Identity_pct)) +
    # histogram with color coding per superfamily
    geom_histogram(binwidth = binwidth*100,
                   color = "black",
                   fill = ifelse(sf == "Copia", "#1b9e77", "#d95f02")) +
    # one facet per Clade with fixed y scale so bars are comparable
    facet_wrap(~Clade, ncol = 1, scales = "fixed") +
    # x axis focused around xlims values
    scale_x_continuous(limits = xlims, breaks = seq(xlims[1], xlims[2], 2)) +
    # set y limit to the provided ymax (useful for consistent overview plots)
    scale_y_continuous(limits = c(0, ymax), expand = c(0, 0)) +  
    theme_cowplot() +
    theme(strip.background = element_rect(fill = "#f0f0f0"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold", hjust = 0.5)) +
    labs(title = sf, x = "LTR identity (%)", y = "Count")
}

#-------------------------------------------------
# Generate Copia and Gypsy plots
#-------------------------------------------------
p_copia <- plot_by_clade(anno_cls, "Copia", global_ymax)
p_gypsy <- plot_by_clade(anno_cls, "Gypsy", global_ymax)

# Combine with cowplot side-by-side 
combined <- plot_grid(p_copia, p_gypsy, ncol = 2, rel_widths = c(1, 1))

dir.create("plots", showWarnings = FALSE)
ggsave("plots/02_LTR_Copia_Gypsy_cladelevel.png", combined, width = 12, height = 10, dpi = 300)
ggsave("plots/02_LTR_Copia_Gypsy_cladelevel.pdf", combined, width = 12, height = 10)

message("Plots saved to ./plots/")