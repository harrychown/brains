library(ggplot2)
library(dplyr)
library(readr)
library(ggsci)
library(forcats)
library(rstatix)

metadata <-  read.csv("your/filepath/sample_metadata.csv")

raw <- list.files(path="your/filepath/taxa_calc", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows 

raw$BarcodeID <- sapply(raw$BarcodeID,function(x) gsub("barcode_","",as.character(x)))
#raw$BarcodeID <- factor(raw$BarcodeID, levels=c("0", "1", "2", "3", "4", "5","6", 
#                                                "7", "8", "9", "10", "11", "12", 
#                                                "13","14", "15", "16", "17", 
#                                                "18", "19"))


raw$GenusLevel[raw$GenusLevel == "Filobasidiaceae"] <- "Unclassified Filobasidiaceae"
GenusLevel <- raw %>%
  select(BarcodeID, GenusLevel, ReadCount) %>%
  group_by(BarcodeID, GenusLevel) %>%
  summarize(TotalCount = sum(ReadCount)) %>%
  group_by(BarcodeID) %>%
  mutate(SumCount = sum(TotalCount)) %>%
  mutate(RelativePercentage = (TotalCount / SumCount) * 100)

Top20Genera <- GenusLevel %>%
  group_by(GenusLevel) %>%
  summarise(TotalRelativePercentage = sum(RelativePercentage)) %>%
  arrange(desc(TotalRelativePercentage)) %>%
  head(19)

Top20GeneraDF <- GenusLevel %>%
  filter(GenusLevel %in% Top20Genera$GenusLevel)

# Step 3: Sum the RelativePercentage of all other genera into a new GenusLevel called "other"
OtherDF <- GenusLevel %>%
  filter(!GenusLevel %in% Top20Genera$GenusLevel) %>%
  group_by(BarcodeID) %>%
  summarise(GenusLevel = "Other",
            TotalCount = sum(TotalCount),
            RelativePercentage = sum(RelativePercentage))

# Step 4: Combine the top 20 genera with the "other" category
FinalGenusDF <- bind_rows(Top20GeneraDF, OtherDF)
FinalGenusDF$GenusLevel <- fct_relevel(FinalGenusDF$GenusLevel, c("Other","Unknown fungal classification", "Homo"), after = Inf)


FinalGenusDF <- merge(FinalGenusDF, metadata, by.x = "BarcodeID", by.y = "Code", all.x = TRUE)


# Calculate reads per sample 
SampleCounts <- GenusLevel %>%
  group_by(BarcodeID) %>%
  summarise(TotalCount = sum(TotalCount))

FinalSampleCounts <- merge(SampleCounts, metadata, by.x = "BarcodeID", by.y = "Code", all.x = TRUE)


SpeciesLevel <- raw %>%
  select(BarcodeID, SpeciesLevel, ReadCount) %>%
  group_by(BarcodeID, SpeciesLevel) %>%
  summarize(TotalCount = sum(ReadCount)) %>%
  group_by(BarcodeID) %>%
  mutate(SumCount = sum(TotalCount)) %>%
  mutate(RelativePercentage = (TotalCount / SumCount) * 100)



# Arrange by condition type
FinalGenusDF <- FinalGenusDF %>%
  arrange(Condition) %>%  
  mutate(BarcodeID = factor(BarcodeID, levels = unique(BarcodeID[order(Condition)])))

relAbu <- ggplot(FinalGenusDF, aes(x = BarcodeID, y = RelativePercentage, fill = GenusLevel)) + 
  geom_bar(stat = "identity", width=0.8, colour="black") +
  #facet_wrap(~ Condition) +  # Faceting based on BarcodeID
  scale_fill_d3(palette="category20") +
  ylab("Relative Abundance (%)") +
  #facet_wrap(phase~part*treatment, scales = "free_x", nrow = 6, ncol=4) +
  labs(title = "Relative Abundance of Genus") +
  scale_y_continuous(sec.axis = sec_axis(~. * 600, name = "Read Count")) +
  geom_point(inherit.aes = FALSE, data = FinalSampleCounts, 
            aes(x = BarcodeID, y = TotalCount/600), size = 2) +
  xlab("Sample") +
  theme_minimal()

# Extract colourscheme
gg_build <- ggplot_build(relAbu)
colour_scheme <- gg_build[["plot"]][["scales"]][["scales"]][[1]][["palette.cache"]]
colour_names <- gg_build[["plot"]][["scales"]][["scales"]][[1]][["range"]][["range"]]
names(colour_scheme) <- colour_names

pdf(file = "your/filepath/relative_abundance.pdf", width=12, height=6)
print(relAbu)
dev.off()

readCount <- ggplot(FinalSampleCounts, aes(x=Condition, y=TotalCount)) +
  geom_boxplot(aes(fill=Condition)) +
  ylab("Total Read Number") +
  theme_minimal()
pdf(file = "your/filepath/read_count.pdf", width=4, height=6)
print(readCount)
dev.off()

write.csv(FinalSampleCounts, "your/filepath/read_count.csv")
write.csv(GenusLevel, "your/filepath/genus_data.csv")
# Sample size too small for significance
#t.test(TotalCount ~ Condition, data=FinalSampleCounts)
#stat.test <- FinalSampleCounts %>%
#  t_test(TotalCount ~ Condition) %>%
#  adjust_pvalue(method = "BH") %>%
#  add_significance()
#stat.test

# Focus on samples 5,6 and 8
focus_group <- FinalGenusDF[FinalGenusDF$BarcodeID %in% c("5", "6", "8") & FinalGenusDF$GenusLevel != "Homo", ]
focusGroup <- ggplot(focus_group, aes(x = BarcodeID, y = RelativePercentage, fill = GenusLevel)) + 
  geom_bar(stat = "identity", width=0.8, colour="black") +
  #facet_wrap(~ Condition) +  # Faceting based on BarcodeID
  scale_fill_manual(values = colour_scheme)  +
  ylab("Relative Abundance (%)") +
  #facet_wrap(phase~part*treatment, scales = "free_x", nrow = 6, ncol=4) +
  labs(title = "Relative Abundance of Genus") +
  xlab("Sample") +
  theme_minimal()

pdf(file = "your/filepath/focus_samples.pdf", width=6.7, height=6)
print(focusGroup)
dev.off()
