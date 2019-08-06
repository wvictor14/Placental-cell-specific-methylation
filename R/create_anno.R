library(annotatr)
library(tidyverse)

b4 <- read_csv('../../data/external/MethylationEPIC_v-1-0_B4.csv', skip =7)
b4_missing <- read_csv('../../data/external/MethylationEPIC Missing Legacy CpG (v1.0_B3 vs. v1.0_B2) Annotations.csv')
anno <- bind_rows(b4, b4_missing %>% mutate_at(c('AddressA_ID', 'AddressB_ID'), as.character))

anno <- anno %>% mutate(Strand = ifelse(Strand == 'F', '+', '-'),
                Start = MAPINFO,
                End = MAPINFO,
                CHR = paste0('chr', CHR))%>%
  dplyr::select(Strand, Start, End, CHR, Name) %>%
  
  filter(CHR != 'chrNA')



# select annotations
annots <- c('hg19_enhancers_fantom',
            'hg19_cpgs', 
            'hg19_genes_1to5kb',
            'hg19_genes_promoters', 
            'hg19_genes_5UTRs',
            'hg19_genes_exons',
            'hg19_genes_intronexonboundaries',
            'hg19_genes_introns',
            'hg19_genes_3UTRs',
            'hg19_genes_intergenic')

annots <- build_annotations(genome = 'hg19', annotations = annots)

# coerce to granges
anno_GR <- makeGRangesFromDataFrame(anno,
                                    start.field = 'Start', end.field = 'End', 
                                    seqnames.field = 'CHR', keep.extra.columns = T)

annotated <- annotate_regions(regions = anno_GR, annotations = annots, ignore.strand = T) %>% as_tibble

# clean it up
x <- annotated
x1 <- length(unique(x$Name))

system.time(
x <- x %>% 
  
  # 1. remove identifiers from annot.id, i.e "intron.9999" should be just called "intron"
  mutate(annot.id = str_extract(annot.id, '^[A-z0-9]+'),
         annot.type = str_extract(annot.type, '(?<=hg19_)[^_]+')
         ) %>%
  dplyr::select(-annot.seqnames, -annot.start, -annot.end, -annot.strand) %>%
  
  # 2. paste together all elements, grouped by cpg-related and gene-related elements
  ### I'm not sure but this step might be able to skip
  group_by(Name,annot.type) %>%
  summarize(
    # paste(unique(.)) ensures all nonunique mappings are retained
    cpg = paste(unique(Name)),
    chr = paste(unique(seqnames)),
    start = paste(unique(start)),
    end = paste(unique(end)),
    
    id = paste(annot.id, collapse = ", "),
    width = paste(annot.width, collapse = ", "),
    tx_id = paste(annot.tx_id, collapse = ", "),
    gene_id = paste(annot.gene_id, collapse = ", "),
    symbol = paste(annot.symbol, collapse = ", "))  %>%
  
  # 3. now spread multiple values across columns
  # first create a temporary variable to spread by, using gather() and unite()
  ungroup() %>%
  dplyr::select(-Name)  %>%
  gather(variable, value, -(annot.type:end)) %>%
  unite(temp, annot.type, variable, sep = '_') %>%
  spread(temp, value) %>%
  
  # 4. clean up
  # remove uninformative columns
  dplyr::select(-enhancers_gene_id , -enhancers_symbol, -enhancers_tx_id, 
                -cpg_gene_id, -cpg_symbol, -cpg_tx_id) %>%
  
  # cpg_id == inter replace with "sea"
  mutate(cpg_id = ifelse(cpg_id == 'inter', 'sea', cpg_id)) %>%
  mutate_all(~replace(., . == 'NA', NA_character_))
) # takes 10 hours
# ensure the number of cpgs ended up with equals the number of unique cpgs starting out with
x1 == nrow(x) # must be T 
  

saveRDS(x, '../../data/main/interim/annotation.rds')
write_csv(x, '../../data/main/processed/annotation.csv')

## column descriptors
colnames(x)
# cpg: cpg identifier
# chr: chromosome position
# start, end: cpg position in bp 
# cpg_id: relation to cpg island features. uses UCSC cpg islands as reference, then definitions for shelves, shores, and seas are used
# cpg_width: width in bp of the cpg island -related feature
# enhancers_id: fantom5 enhancer
# enhancers_width: width in bp of the fantom5 enhancer
# genes_: columns are retrieved from txdb.hsapiens.ucsc.hg19.knowngene
# genes_gene_id: I'm not sure what this is TBH
# genes_id: gene-related elemtn, e.g. promoters, introns/exons, based off of UCSC
# genes_symbol: Gene symbol
# genes_tx_id: transcript identifier
# genes_width: width in bp of gene element
