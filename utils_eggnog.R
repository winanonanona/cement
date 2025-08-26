# Script for helper functions
# date: 20241225

library(ggrepel)
library(tidyverse)
## Auxilary function to get COG counts per dataframe
get_cog_count <- function(eggnog_annot) {
  cog_count <- eggnog_annot %>% 
    mutate(COG = str_extract(eggNOG_OGs, "(?<=,)[^@,]+(?=@\\d+\\|(?:Bacteria))")) %>% 
    dplyr::select(eggNOG_OGs, COG) %>% 
    dplyr::count(COG)
  return(cog_count)
  
}

## Auxilary function to merge two profiles
merge2_cog <- function(x, y){
  aux <- function(x) read_tsv(file.path(cog_list_path, x),
                              col_names=c(names(test)),
                              show_col_types = FALSE, skip = 5)
  # extract COG column
  if(is.character(x)){
    dat.x <- get_cog_count(aux(x)) 
    names(dat.x)[2] <- strsplit(x, '[.]')[[1]][1]
  }else{
    dat.x <- x
  }
  dat.y <- get_cog_count(aux(y)) 
  names(dat.y)[2] <- strsplit(y, '[.]')[[1]][1]
  full_join(dat.x, dat.y)
}

## Auxilary function to get preferred_name counts per dataframe
get_gene_count <- function(eggnog_annot) {
  cog_count <- eggnog_annot %>% 
    # if preferred_name is empty, use COG
    mutate(COG = case_when(Preferred_name == '-' ~ 
                             str_extract(eggNOG_OGs, "(?<=,)[^@,]+(?=@\\d+\\|(?:Bacteria|Archaea))"),
                           TRUE ~ Preferred_name)) %>% 
    dplyr::select(eggNOG_OGs, COG) %>% 
    dplyr::count(COG)
  return(cog_count)
  
}

## Auxilary function to merge two profiles
merge2_gene <- function(x, y){
  aux <- function(x) read_tsv(file.path(cog_list_path, x),
                              col_names=c(names(test)),
                              show_col_types = FALSE, skip = 5)
  # extract COG column
  if(is.character(x)){
    dat.x <- get_gene_count(aux(x)) 
    names(dat.x)[2] <- strsplit(x, '[.]')[[1]][1]
  }else{
    dat.x <- x
  }
  dat.y <- get_gene_count(aux(y)) 
  names(dat.y)[2] <- strsplit(y, '[.]')[[1]][1]
  full_join(dat.x, dat.y)
}

# Consider non-top features as 'Other', just for genus & species
top_features_other <- function(cog_profiles) {
  means <- rowMeans(cog_profiles)
  cog_profiles.fil <- cog_profiles[means>=sort(means, decreasing = TRUE)[n_features], ]
  
  # Add 'Other' to the dataframe
  others <- cog_profiles %>% 
    filter(!(row.names(.) %in% rownames(cog_profiles.fil)))
  sums <- t(data.frame(colSums(others)))
  rownames(sums) <- 'Other'
  cog_profiles.fil <- rbind(cog_profiles.fil, sums)
  
  ## convert to long
  cog_profiles.dat <- tibble::rownames_to_column(cog_profiles.fil) %>% 
    reshape2::melt()
  
  # Arrange rowname based on mean
  cog_profiles.dat <- cog_profiles.dat %>%
    mutate(rowname = forcats::fct_reorder(rowname, value, .fun = mean))
  return(cog_profiles.dat)
}  

top_features <- function(cog_profiles) {
  if(dim(cog_profiles)[[1]] < n_features) {
    cog_profiles.fil <- cog_profiles
    cog_profiles.dat <- tibble::rownames_to_column(cog_profiles.fil) %>% 
      reshape2::melt()
  }
  else {
    means <- rowMeans(cog_profiles)
    cog_profiles.fil <- cog_profiles[means>=sort(means, decreasing = TRUE)[n_features], ]
    
    ## renormalize to 100%
    cog_profiles.fil <- data.frame(apply(cog_profiles.fil, 2, function(x) x/sum(x)) * 100)
    cog_profiles.dat <- tibble::rownames_to_column(cog_profiles.fil) %>% 
      reshape2::melt()  
  }
  
  # Arrange rowname based on mean
  cog_profiles.dat <- cog_profiles.dat %>%
    mutate(rowname = forcats::fct_reorder(rowname, value, .fun = mean))
  return(cog_profiles.dat)
}

c24 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2",
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

plot_overview <- function(top_cog_features) {
  ggplot(top_cog_features, aes(x=variable, y=value, fill=rowname)) + 
    geom_bar(stat = 'identity') + 
    scale_fill_manual(values=c24) + 
    labs(x="Samples", y="Relative abudances (%)") + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1))
}

## similar to above but faceted according to category
plot_overview_stratified <- function(cog_profiles.dat) {
  tmp <- merge(cog_profiles.dat, metadata, by.x=2, by.y=5, all.x=TRUE)
  p <- ggplot(tmp, aes(x=variable, y=value, fill=rowname)) + 
    geom_bar(stat = 'identity') + 
    #scale_fill_manual(values=color.list, na.value='grey') + 
    theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
    labs(x="Samples", y="Relative abudances (%)", fill = 'COG') + 
    #facet_grid(vars(Time), vars(Group), scales = "free_x", space="free")
    facet_grid(col = vars(Side), scales = 'free_x',
               space = 'free') + 
    scale_fill_manual(values = c24) +
    theme_bw() +
    theme(text = element_text(size = 28),
          strip.background = element_rect(fill = 'white'),
          strip.text = element_text(size = 32))
  return(p)
  
}

# run wilcoxon
run_wilcox <- function(cog_profiles, group1, group2) {
  cog_profiles.filt <- cog_profiles %>% 
    filter(row.names(.) %in% rownames(cog_profiles))
  # renormalize
  cog_profiles.filt <- data.frame(apply(cog_profiles.filt, 2, 
                                         function(x) x/sum(x)) * 100)
  #print(colSums(cog_profiles.filt))
  group1_libraries <- metadata %>% 
    filter(Side==group1) %>% pull(SampleID)
  print(group1_libraries)
  group2_libraries <- metadata %>% 
    filter(Side==group2)%>% pull(SampleID)
  print(group2_libraries)
  
  p.values <- apply(cog_profiles, 1, function(x)
    wilcox.test(x[group1_libraries],
                x[group2_libraries])$p.value)
  q.values <- p.adjust(p.values, method="fdr")
  
  return(data.frame(p.values) %>% rownames_to_column('COG') %>% 
           inner_join(data.frame(q.values) %>% rownames_to_column('COG'), 
                      by = 'COG'))
}

## Plot DA results from Wilcoxon
plot_da_results <- function(combined_res, profiles, Group1, Group2, comparison) {
  # get specified column of comparison
  combined_res <- combined_res %>% dplyr::select(COG, comparison) %>% 
    na.omit()
  da_relabund <- (profiles * 100) %>% 
    filter(row.names(.) %in% combined_res$Taxa) %>% 
    tibble::rownames_to_column('Taxa') %>% 
    pivot_longer(-Taxa, names_to ='LibraryID') %>%
    left_join(metadata, by = 'LibraryID')
  p <- da_relabund %>% 
    filter(Group %in% c(Group1, Group2)) %>% 
    ggplot(aes(x = Time, y = value, col = Group)) + 
    geom_boxplot(outlier.shape = NA) + geom_jitter(position = position_jitterdodge()) +
    #stat_compare_means(method='wilcox', label='p.signif') +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    facet_wrap(~Taxa, ncol = 5, scales='free') + 
    labs(y = '% Relative abundance') +
    theme(strip.text = element_text(size = 12),
          legend.position = 'top')
  return(p)
}

# get deseq2 results
get_de_res <- function(deseq2_res, padj_thresh) {
  signif_res <- deseq2_res %>% 
    as.data.frame() %>% 
    filter(padj < padj_thresh) %>% 
    tibble::rownames_to_column('COG') %>% 
    left_join(bacteria_og_annot, by = 'COG')
  return(signif_res)
}

# plot volcano
plot_volcano <- function(deseq2_res, padj_thresh, plot_title) {
  add_signif_res <- deseq2_res %>% 
    as.data.frame() %>% 
    na.omit() %>% 
    #filter(padj < 0.1) %>% 
    tibble::rownames_to_column('COG') %>% 
    mutate(signif = case_when(padj < 0.1 ~ paste('padj < ', 0.1, sep = ''),
                              TRUE ~ 'n.s.')) %>% 
    mutate(signif = case_when(padj < 0.05 ~ paste('padj < ', 0.05, sep = ''),
                              TRUE ~ signif))
  
  add_signif_res %>% 
    ggplot(aes(x = log2FoldChange, y = -log10(padj), col = signif)) +
    geom_point(size = 4, alpha = 0.75) +
    ggrepel::geom_text_repel(data = add_signif_res %>% 
                      filter(padj < 0.1),
                             aes(label = COG),
                    show.legend = F,
                    size = 7, max.overlaps = 20) +
    geom_hline(yintercept = -log10(0.1), linetype = 'dashed', size = 1.5) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed', size = 1.5) +
    ggtitle(plot_title) +
    scale_color_manual(values= c('grey', 'dodgerblue2', '#E31A1C')) +
    theme_bw() + 
    theme(text = element_text(size = 32))
}

# Run GSEA
GSEA_from_DESeq <- function(genelist, gene_names, term2gene, term2name = KEGG_term2name){
  
  names(genelist) <- as.character(gene_names)
  # omit any NA values 
  genelist<-na.omit(genelist)
  #sort the genelist in decreasing order.
  genelist <- sort(genelist, decreasing = TRUE)
  
  set.seed(123)
  GSEA_out <- clusterProfiler::GSEA(geneList = genelist, TERM2GENE = term2gene, 
                   TERM2NAME = term2name, eps= 0, nPermSimple=10000, seed=TRUE)
  
  return(GSEA_out)
  
}

#https://www.biostars.org/p/467197/ for plotting pathway GSEA normalized enrichment scores in a bar plot

make_NES_barplot <- function(GSEA_df, title, x_axis_lab="Function/Process"){
  
  GSEA_df$direction <- ifelse(GSEA_df$NES > 0, "upregulated", "downregulated")
  cols <- c("downregulated" = "darkblue", "upregulated" = "red")
  
  GSEA_df$ID_and_Desc <- paste0(GSEA_df$ID," ",GSEA_df$Description)
  
  
  ggplot(GSEA_df, aes(reorder(ID_and_Desc, NES), NES, fill = direction)) +
    geom_col() +
    scale_fill_manual(values = cols) +
    theme(axis.text.x = element_text(vjust = 0.5, hjust=1)) +
    coord_flip() +
    labs(x=x_axis_lab, y="Normalized Enrichment Score") + ggtitle(title) + theme_classic()
  
}

make_NES_barplot_combined <- function(GSEA_df, x_axis_lab="Function/Process"){
  
  GSEA_df$direction <- ifelse(GSEA_df$NES > 0, "upregulated", "downregulated")
  cols <- c("downregulated" = "darkblue", "upregulated" = "darkred")
  
  GSEA_df$ID_and_Desc <- paste0(GSEA_df$ID," ",GSEA_df$Description)
  
  
  ggplot(GSEA_df, aes(reorder(ID_and_Desc, NES), NES, fill = direction)) +
    geom_col() +
    scale_fill_manual(values = cols) +
    theme(axis.text.x = element_text(vjust = 0.5, hjust=1), legend.position = "none") +
    guides(fill = "none") +
    coord_flip() +
    labs(x=x_axis_lab, y="Normalized Enrichment Score") +
    theme_classic() +
    facet_wrap(~comparison, ncol = 1, scales = "free_y")
  
}

extract_archaea_cogs <- function(filepath) {
  
  # Check if the file exists to provide a clear error message.
  if (!file.exists(filepath)) {
    stop("Error: File not found at the specified path: ", filepath)
  }
  
  # Read the eggNOG annotation file.
  # read_tsv is used because the file is tab-separated.
  # The 'comment' argument automatically skips the header lines starting with '##'.
  # The first non-commented line is correctly treated as the column headers.
  annotations <- read_tsv(filepath, comment = "##", col_types = cols(.default = "c"))
  
  # Check if the data frame is empty after reading
  if (nrow(annotations) == 0) {
    message("File is empty or contains only header comments: ", filepath)
    # Return an empty tibble with the expected columns
    return(tibble(arCOG = character(), COG_category = character(), Description = character()))
  }
  
  # Process the data to find and extract Archaea information.
  archaea_data <- annotations %>%
    # Step 1: Filter rows where 'eggNOG_OGs' contains the text "|Archaea".
    # We need to escape the pipe character with '\\' because it has a special
    # meaning in regular expressions.
    filter(str_detect(eggNOG_OGs, "\\|Archaea")) %>%
    # Step 2: Create a new column 'arCOG' by extracting it from 'eggNOG_OGs'.
    # The pattern "arCOG\\d+" matches the literal string "arCOG" followed by
    # one or more digits (\\d+).
    mutate(arCOG = str_extract(eggNOG_OGs, "arCOG\\d+")) %>%
    # Step 3: Keep only the columns of interest.
    select(arCOG, COG_category, Description)
  
  return(archaea_data)
}
