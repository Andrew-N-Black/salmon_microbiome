library(phyloseqCompanion)

sample_data_df = phyloseqCompanion::sample.data.frame(ps_filtered)
sample_data_df %>% kruskal_test(TotalReads ~ hatchery)
# A tibble: 1 × 6
#  .y.            n statistic    df       p method        
#* <chr>      <int>     <dbl> <int>   <dbl> <chr>         
#1 TotalReads    61      16.5     5 0.00555 Kruskal-Wallis
