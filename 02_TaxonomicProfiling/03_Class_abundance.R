library(phyloseq)
library(qiime2R)
library(microViz)
library(ggplot2)
library(tibble)
library(dplyr)
library(reshape2)



###Class Abundance plots##

glom <- tax_glom(ps_rarefied, taxrank = 'Class', NArm = FALSE)
#Melt and merge dataframe to work with ggplot2
ps.melt <- psmelt(glom)
#Set as character
ps.melt$Class <- as.character(ps.melt$Class)

ps.melt <- ps.melt %>%group_by(ASE, Class) %>% mutate(median=median(Abundance))
keep <- unique(ps.melt$Class[ps.melt$median > 2.5])
ps.melt$Class[!(ps.melt$Class %in% keep)] <- "< 2.5%"
ps.melt_sum <- ps.melt %>% group_by(Sample,ASE,Class) %>% summarise(Abundance=sum(Abundance))

ggplot(ps.melt_sum, aes(x = Sample, y = Abundance, fill = Class)) + 
    geom_bar(stat = "identity", aes(fill=Class)) + 
    labs(x="", y="%") +
    facet_wrap(~ASE, scales= "free_x", nrow=1) +
    theme_q2r()+scale_fill_brewer(palette = "Dark2",name = "Class")+theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())
ggsave("~/Figure_3.svg")
