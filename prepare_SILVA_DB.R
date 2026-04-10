# R script used to download and modify the SILVA to species database for the class
# Kept only one sequence per species
# Remove sequences for which the classification did not include all 7 ranks
# Modify E.Coli name

install.packages("BiocManager")
BiocManager::install("Biostrings")

library(Biostrings)
library(dplyr)
library(tidyr)
library(seqinr)

dna_sequences = Biostrings::readDNAStringSet("silva_nr99_v138.2_toSpecies_trainset.fa")
df = data.frame(name = names(dna_sequences),sequence = as.character(dna_sequences), Sequence_lenght = width(dna_sequences))

df_unique = df %>%
  group_by(name) %>%
  arrange(Sequence_lenght) %>%
  slice_head(n = 1) %>%
  mutate(split_count = str_count(name, ";"))

df_tospecies = df_unique %>%
  filter(split_count == 7) %>%
  separate(name, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
  data.frame()

df_tospecies$Genus = gsub("Escherichia-Shigella", "Escherichia", df_tospecies$Genus)

df_tospecies2 = df_tospecies %>%
  mutate(sequence_name = paste(Domain, Phylum, Class, Order, Family, Genus, Species, sep = ";"))

write.fasta(as.list(df_tospecies2$sequence), names = df_tospecies2$sequence_name, nbchar = 80, as.string = TRUE, file.out = "silva_nr99_v138_2_BIA1401.fasta")
