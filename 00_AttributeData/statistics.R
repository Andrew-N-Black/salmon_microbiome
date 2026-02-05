library(rstatix)
library(dplyr)


#Kruskal test for epithelium remaining by hatchery
my_data <- read.table("~/SMB_n61/input/metadata61_ext.txt",header=TRUE)
my_data %>% kruskal_test(epithelium_remaining ~ hatchery)

# A tibble: 1 × 6
#  .y.                      n statistic    df             p method        
#* <chr>                <int>     <dbl> <int>         <dbl> <chr>         
#1 epithelium_remaining    61      49.5     5      0.00000000177 Kruskal-Wallis



#Kruskal test for Enteritis Score by hatchery
my_data %>% kruskal_test(enteritis ~ hatchery)
# A tibble: 1 × 6
  .y.           n statistic    df      p method        
* <chr>     <int>     <dbl> <int>  <dbl> <chr>         
1 enteritis    61      15.2     2 0.0005 Kruskal-Wallis

# Epithelium integrity vs Enteritis score
cor.test(my_data$epithelium_remaining, my_data$enteritis, method = "spearman",exact=FALSE)

	#Spearman's rank correlation rho

#data:  my_data$epithelium_remaining and my_data$enteritis
#S = 23904, p-value = 0.00353
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#      rho 
#0.3679424 


 my_data %>% kruskal_test(es ~ hatchery)
# A tibble: 1 × 6
  .y.       n statistic    df       p method        
* <chr> <int>     <dbl> <int>   <dbl> <chr>         
1 es       61      10.1     2 0.00644 Kruskal-Wallis

# Cshasta
 my_data %>% kruskal_test(cshasta ~ hatchery)
# A tibble: 1 × 6
  .y.         n statistic    df        p method        
* <chr>   <int>     <dbl> <int>    <dbl> <chr>         
1 cshasta    61      14.5     2 0.000714 Kruskal-Wallis
