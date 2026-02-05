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

library(dplyr)
library(betareg)
library(ordinal)

#Formulate a dataframe with squeezed epithelium integrity values.
df2 <- metadata %>%
    mutate(
        # outcome as proportion
        epithelium_prop = percent_epithelium / 100,
        
        # scores as ordered predictors (0 < 1 < 2 < 3)
        Cshasta_score   = ordered(cshasta, levels = c(0, 1, 2, 3)),
        Esherekii_score = ordered(es, levels = c(0, 1, 2)),Enteritis_ordered   = ordered(enteritis, levels = c("E0", "E1", "E2", "E3")))

## "Squeeze" 0 and 1 into (0,1) so beta regression is valid. Beta regression only works on 0 - 1 non-inclusive. 
So we add a tiny think to those values which are exactly 1 or 0. 
This is an adjustment that is commonly done for beta regression which takes into account the sample size and the y value.

n <- nrow(df2)
df2 <- df2 %>%
    mutate(epithelium_prop_squeezed = (epithelium_prop * (n - 1) + 0.5) / n)

#Run a betaregression
m_beta <- betareg(
  epithelium_prop_squeezed ~ Cshasta_score + Esherekii_score + hatchery,
  data = df2,
  link = "logit"
)

#Summarize model
summary(m_beta)

Quantile residuals:
    Min      1Q  Median      3Q     Max 
-2.1493 -0.4646  0.2667  0.2667  2.5471 

Coefficients (mean model with logit link):
                        Estimate Std. Error z value Pr(>|z|)    
(Intercept)            5.268e+00  5.957e-01   8.844  < 2e-16 ***
Cshasta_score.L        4.308e-01  3.497e-01   1.232   0.2181    
Cshasta_score.Q       -1.540e+00  3.229e-01  -4.770 1.84e-06 ***
Cshasta_score.C        1.072e+00  2.519e-01   4.256 2.08e-05 ***
Esherekii_score.L      1.343e+00  5.482e-01   2.449   0.0143 *  
Esherekii_score.Q     -1.218e-01  3.453e-01  -0.353   0.7243    
hatcheryround_butte   -7.664e+00  7.358e-01 -10.415  < 2e-16 ***
hatcherysandy         -4.617e+00  6.278e-01  -7.355 1.91e-13 ***
hatcherysouth_santiam -3.900e+00  6.080e-01  -6.415 1.41e-10 ***
hatcherywhite_river   -5.430e-16  4.399e-01   0.000   1.0000    
hatcherywillamette    -5.651e+00  6.450e-01  -8.762  < 2e-16 ***

Phi coefficients (precision model with identity link):
      Estimate Std. Error z value Pr(>|z|)    
(phi)    7.248      1.413   5.131 2.89e-07 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 

Type of estimator: ML (maximum likelihood)
Log-likelihood: 82.56 on 12 Df
Pseudo R-squared: 0.9038
Number of iterations: 42 (BFGS) + 3 (Fisher scoring) 


#Enteritis###
mm<-clmm(
    Enteritis_ordered ~ Cshasta_score + Esherekii_score + (1 | hatchery),
    data = df2)

summary(mm)
Cumulative Link Mixed Model fitted with the Laplace approximation

formula: Enteritis_ordered ~ Cshasta_score + Esherekii_score + (1 | hatchery)
data:    df2

 link  threshold nobs logLik AIC    niter     max.grad cond.H 
 logit flexible  61   -45.34 108.68 432(2903) 6.78e-06 1.5e+02

Random effects:
 Groups   Name        Variance Std.Dev.
 hatchery (Intercept) 22.92    4.788   
Number of groups:  hatchery 6 

Coefficients:
                  Estimate Std. Error z value Pr(>|z|)
Cshasta_score.L     0.6889     1.0903   0.632    0.527
Cshasta_score.Q     0.2927     1.0278   0.285    0.776
Cshasta_score.C    -0.1049     0.8134  -0.129    0.897
Esherekii_score.L  -2.5201     1.9524  -1.291    0.197
Esherekii_score.Q  -1.3148     1.2054  -1.091    0.275

Threshold coefficients:
      Estimate Std. Error z value
E0|E1   -6.212      2.594  -2.395
E1|E2   -1.933      2.383  -0.811
E2|E3    2.697      2.422   1.113
