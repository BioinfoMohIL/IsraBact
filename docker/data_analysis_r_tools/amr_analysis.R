library(tidyverse)
library(janitor)
library(pheatmap)
library(patchwork)
library(argparse)

# Argument parser
parser <- ArgumentParser(description = 'AMR Data Analysis Script')
parser$add_argument('--amr_input', required=TRUE, help='Path of arm phenotype file (csv, xlsx)')
parser$add_argument('--remove_null', action = 'store_true', help = 'If set, rows with NA resistance will be removed')
parser$add_argument('--output', required=TRUE, help='Dir path to save the plots')

args <- parser$parse_args()
output_dir <- args$output


# Ex: Load amr_phenotype_eu.csv
phenotype_raw <- read.csv(args$amr_input)

#Clean and standardize the phenotype dataset
phenotype <- phenotype_raw %>%  # Start from the raw dataset
  clean_names()  # Convert column names to lowercase and snake_case for consistency (from janitor package)

phenotype <- phenotype %>%
  mutate(
    # Capitalize the first letter of each word in species and country (e.g., "france" → "France")
    species = str_to_title(species),
    country = str_to_title(country),
    
    # Convert antibiotic test result values to uppercase ("r" → "R", "s" → "S")
    # The tilde (~) introduces an anonymous function, and the dot (.) refers to each value in the columns
    across(c(ciprofloxacin:colistin), ~str_to_upper(.))
  )


## Convert to Long Format and handle missing values

#Reshape the phenotype data to long format and handle missing values
phenotype_long <- phenotype %>%
  # Convert from wide to long format: one row per antibiotic result
  pivot_longer(
    cols = ciprofloxacin:colistin,          # All antibiotic test columns
    names_to = "antibiotic",                # New column for antibiotic names
    values_to = "resistance"                # New column for resistance results (S, R, or NA)
  ) %>%
  mutate(
    # Capitalize antibiotic names for clean plotting
    antibiotic = str_to_title(antibiotic),
    
    # Convert resistance column to a factor with levels in logical order
    resistance = factor(resistance, levels = c("S", "R"))
  )


# Option: Remove rows with missing resistance data (uncomment to use)
if (args$remove_null) {
  phenotype_long <- phenotype_long %>%
    filter(!is.na(resistance))
}


## Resistance Proportion by Country and Year
# Summarize resistance proportions
res_prop <- phenotype_long %>%
  # Count the number of isolates by country, year, antibiotic, and resistance status ("S" or "R")
  group_by(country, year, antibiotic, resistance) %>%
  summarise(count = n(), .groups = 'drop') %>%  # Drop groups to prevent grouped behavior downstream

  # Re-group by country, year, and antibiotic (excluding resistance category)
  group_by(country, year, antibiotic) %>%
  # Calculate the proportion of each resistance type (e.g., R or S) within each group
  mutate(proportion = count / sum(count))

# Preview resistant proportions only
res_prop %>%
  filter(resistance == "R")


## Ciprofloxacin Resistance by Country (Styled Plot)

cipro_plot <- phenotype_long %>%
  filter(antibiotic == "Ciprofloxacin") %>% #add !is.na(resistance) to remove NAs 
  ggplot(aes(x = country, fill = resistance)) +
  geom_bar(position = "fill") +
  labs(title = "Ciprofloxacin Resistance by Country",
       y = "Proportion of Isolates", x = "Country") +
  theme_classic() +
  scale_fill_manual(values = c("S" = "steelblue", "R" = "firebrick"))


# Example: Visualizing resistance to Ampicillin
amp_plot <- phenotype_long %>%
  filter(antibiotic == "Ampicillin", !is.na(resistance)) %>%
  ggplot(aes(x = country, fill = resistance)) +
  geom_bar(position = "fill") +
  labs(title = "Ampicillin Resistance by Country",
       y = "Proportion of Isolates", x = "Country") +
  theme_minimal() +
  scale_fill_manual(values = c("S" = "steelblue", "R" = "firebrick"))



## Faceted Plot: Resistance by Antibiotic and Year
# Create a faceted bar plot to visualize resistance trends by year (excluding missing resistance data)

facet_plot <- phenotype_long %>%
  filter(!is.na(resistance)) %>%                    # Remove rows with missing resistance values
  ggplot(aes(x = antibiotic, fill = resistance)) +  # X-axis: antibiotics; fill by resistance status (S or R)
  geom_bar(position = "fill") +                     # Use "fill" to show proportions within each antibiotic group
  facet_wrap(~year) +                               # Create a separate panel (facet) for each year
  labs(
    title = "Resistance by Antibiotic and Year",    # Main plot title
    y = "Proportion",                               # Y-axis label showing proportions (stacked bars)
    x = "Antibiotic"                                # X-axis label
  ) +
  theme_bw() +                                      # Use a clean black-and-white theme
  scale_fill_manual(                                # Set custom fill colors for "S" and "R"
    values = c("S" = "seagreen", "R" = "tomato")
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for readability
  )
  

# We start by ensuring that all isolate IDs are unique.
# This is important because we'll use these IDs as row names in our matrix.
# If any duplicates exist, this step will make them unique by adding suffixes like ".1", ".2", etc.
phenotype_unique <- phenotype %>%
  mutate(isolate_id = make.unique(as.character(isolate_id)))


# Step 1: Ensure isolate IDs are unique
phenotype_unique <- phenotype %>%
  mutate(isolate_id = make.unique(as.character(isolate_id)))


##Step 2: Create Binary Resistance Matrix

# Next, we create a resistance matrix where:
# - Rows are isolates
# - Columns are antibiotics
# - "R" is converted to 1 (resistant), and everything else to 0 (not resistant)
# This binary matrix is what weparser$add_argument('--output', required=TRUE, help='Path to save the heatmap image (PNG, PDF, etc.)')
# will use in our heatmap.
phenotype_matrix <- phenotype_unique %>%
  select(isolate_id, ciprofloxacin:colistin) %>%
  mutate(across(c(ciprofloxacin:colistin), ~ ifelse(. == "R", 1, 0))) %>%
  column_to_rownames("isolate_id") %>%
  as.matrix()


##Step 3a: Create Annotations for Heatmap Rows

# To make our heatmap easier to interpret, we add annotations:
# - Species and country will be shown alongside the rows of the heatmap.
# - We use str_to_title() to clean and standardize the formatting (e.g., "e.coli" → "E. Coli")
annotation_row <- phenotype_unique %>%
  select(isolate_id, species, country) %>%
  mutate(
    Species = str_to_title(species),
    Country = str_to_title(country)
  ) %>%
  select(isolate_id, Species, Country) %>%
  column_to_rownames("isolate_id")



##Step 3b: Order the Matrix and Annotations by Country 

# Optionally, we can sort the heatmap by country (and implicitly by species within each country).
# This makes the heatmap visually cleaner and easier to interpret geographically.
annotation_row <- annotation_row %>%
  arrange(Country)

# We also need to reorder the resistance matrix to match the sorted annotation rows.
phenotype_matrix <- phenotype_matrix[rownames(annotation_row), ]

##Step 4: Define Color Palettes for Annotations 
# Now we assign specific colors to species and countries for the heatmap annotations.
# This ensures consistent and interpretable color coding across the plot.
species_colors <- setNames(
  RColorBrewer::brewer.pal(length(unique(annotation_row$Species)), "Set2"),
  unique(annotation_row$Species)
)

country_colors <- setNames(
  RColorBrewer::brewer.pal(length(unique(annotation_row$Country)), "Set3"),
  unique(annotation_row$Country)
)

# Combine the palettes into one list to pass to the heatmap
annotation_colors <- list(
  Species = species_colors,
  Country = country_colors
)


##Step 4: Preview Heatmap Without Annotations
# If you want to preview the heatmap without any annotations, you can use this simpler version.
# This shows resistance patterns only, without the context of species or country.
pdf(file.path(output_dir, "heatmap_non_annotated.pdf"), width = 8, height = 6)
pheatmap(
  phenotype_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  # annotation_row = annotation_row,         # <- comment out to remove annotations
  # annotation_colors = annotation_colors,   # <- comment out to remove annotation coloring
  display_numbers = TRUE,
  color = colorRampPalette(c("green", "red"))(100),
  main = "AMR Heatmap with Species and Antibiotics"
)
dev.off()


##Step 5: Final Heatmap with Annotations
# Finally, we generate the full heatmap:
# - Each cell shows
# - Row annotations for species and country add biological and geographic context
# - No clustering is applied — rows and columns are displayed in the order we specified

pdf(file.path(output_dir, "heatmap_annotated.pdf"), width = 8, height = 6)
pheatmap(
  phenotype_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
  display_numbers = TRUE,
  color = colorRampPalette(c("green", "red"))(100),
  main = "AMR Heatmap with Species and Country Annotations"
)
dev.off()


## Saving
# Plots
plots <- list(
  facet_plot = facet_plot,
  cipro_plot = cipro_plot,
  amp_plot = amp_plot
)

for (plot_name in names(plots)) {
  ggsave(
    filename = file.path(output_dir, paste0(plot_name, ".pdf")),
    plot = plots[[plot_name]],
    width = 10,
    height = 7
  )
}
