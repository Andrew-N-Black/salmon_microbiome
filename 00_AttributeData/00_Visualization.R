# =============================================================================
# Purpose:  Visualize host pathology attributes across hatcheries and in
#           relation to parasite load. Produces Figures 1b–1e for the manuscript.
# Inputs:   ps — filtered phyloseq object from 02_BioinformaticProcessingFiltering/01_phyloseq.R
# Outputs:  ~/Figure_1b.svg  — epithelial integrity by hatchery and ASE
#           ~/Figure_1c.svg  — enteritis score distribution by hatchery
#           ~/Figure_1d.svg  — epithelium vs enteritis score scatter
#           ~/Figure_1e.svg  — enteritis score by hatchery, faceted by parasite (es, cshasta)
# Key parameters:
#   ASE    — acute systemic enteritis status (binary: NEG vs ++)
#   enteritis — histological enteritis score (E0–E3)
#   epithelium_remaining — proportion of gut epithelium intact
#   es     — Enterocytozoon schreckii load (ordinal 0–2)
#   cshasta — Ceratonova shasta load (ordinal 0–3)
# =============================================================================

library(ggplot2)
library(dplyr)
library(phyloseqCompanion)

# --- Load metadata ---
#Load in metadata from filtered phyloseq object (See 02_TaxonomicProfiling/01_phyloseq.R)
meta = phyloseqCompanion::sample.data.frame(ps)

# Order hatcheries geographically (north to south, then inland)
desired_order <- c("minter_creek","white_river", "south_santiam", "sandy", "willamette","round_butte")
meta$hatchery <- factor(meta$hatchery, levels = desired_order)

# --- Figure 1b: Epithelial integrity by hatchery and ASE status ---
# Violin + boxplot + jitter shows distribution of epithelium_remaining
# across hatcheries, split by ASE (NEG vs positive)
ggplot(meta, aes(y =epithelium_remaining, x=hatchery, fill = ASE,color=ASE)) +
    geom_violin(alpha=0.2) +
    labs(x = "",
         y = "Epithelial Integrity",
         fill = "ASE") +
    theme_classic(base_size = 14)+scale_fill_brewer(palette = "Dark2",name = "ASE")+geom_boxplot(aes(x=hatchery, y=epithelium_remaining, fill=ASE), width=0.9,alpha=0.2)+geom_jitter(aes(x=hatchery, y=epithelium_remaining), width=0.1)+ theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+scale_color_brewer(palette = "Dark2",name = "ASE")
ggsave("~/Figure_1b.svg")

# --- Figure 1c: Enteritis score distribution by hatchery ---
# Stacked bar shows proportion of fish at each histological enteritis grade (E0–E3)
ggplot(meta, aes(x=hatchery,fill=enteritis)) +geom_bar(position = "stack")+labs(x = "", y = "Count",fill = "hatchery") +theme_q2r()+scale_fill_brewer(palette = "Dark2",name = "Enteritis Score")+theme_classic(base_size = 14)+
 theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave("~/Figure_1c.svg")

# --- Figure 1d: Epithelial integrity vs enteritis score scatter ---
# Confirms expected negative correlation: higher enteritis grade = lower epithelium remaining
# Shape encodes hatchery to reveal whether hatchery confounds the relationship
ggplot(meta, aes(x=enteritis,y=epithelium_remaining,color=ASE,shape=hatchery)) +geom_point(size=5,position = "jitter",aes(color=ASE))+labs(x = "Enteritis Score", y = "Epithelial Integrity",color = "ASE",shape="Hatchery") +theme_classic(base_size = 14)+scale_color_brewer(palette = "Dark2",name = "ASE")
ggsave("~/Figure_1d.svg")


# --- Figure 1e: Enteritis score by hatchery, faceted by parasite species ---
# Melt to long format so es and cshasta can be displayed as separate facets
meta$ID <- rownames(meta)
melted_scores <- melt(meta,id.vars = c("ID","hatchery","enteritis"), measure.vars = c("es", "cshasta"))

# Color = ordinal parasite load score; facet rows = parasite species
ggplot(melted_scores, aes(x=hatchery,y=enteritis,color=as.factor(value))) +geom_point(size=5,position="jitter")+labs(x = "", y = "Enteritis Score",color = "variable") +facet_wrap(~variable,ncol=1)+scale_color_brewer(palette = "Dark2",name = "Parasite Load")+theme_classic(base_size = 14)+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave("~/Figure_1e.svg")



###OLD####

#Figure 1e. Es parasite by hatchery
#Set rownames as ID before melting dataframe
meta$ID <- rownames(meta)
melted_scores <- melt(meta,id.vars = c("ID","hatchery"), measure.vars = c("es", "cshasta"))

ggplot(melted_scores, aes(x=hatchery,fill=as.factor(value))) +geom_bar(position = "stack")+labs(x = "", y = "Count",fill = "variable") +facet_wrap(~variable,ncol=1)+scale_fill_brewer(palette = "Dark2",name = "Parasite Load")+theme_classic(base_size = 14)+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave("~/Figure_1e.svg")
