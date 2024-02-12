# Plot PCs for all populations, then for each subpopulation queried (AFR,AMR,EUR)

pcs = vroom("regenie_covar_20pcs.txt",show_col_types = F)
anc = vroom("ancestry_preds.tsv",show_col_types = F) %>% rename(IID = research_id)
pcs_anc = pcs %>% inner_join(anc,by = "IID")

# Plot results
plt_all = ggplot(pcs_anc,aes(x=PC1_C,y=PC2_C,color=ancestry_pred)) + geom_point() + 
    theme(text = element_text(size=20)) + xlab("PC1") + ylab("PC2") +
    scale_color_viridis(option = "inferno",discrete=T) + xlim(c(-0.002,0.005)) + ylim(c(-0.0025,0.01))
ggsave(plot=plt_all,filename="Output/FigureS1A.png",width=6,height=4)

plt_eur = ggplot(pcs_anc %>% filter(ancestry_pred == "eur"),aes(x=PC1_C,y=PC2_C,color=ancestry_pred)) + geom_point() + 
    theme(text = element_text(size=20)) + xlab("PC1") + ylab("PC2") +
    scale_color_viridis(option = "inferno",discrete=T) + xlim(c(-0.002,0.005)) + ylim(c(-0.0025,0.01))
ggsave(plot=plt_eur,filename="Output/FigureS1D.png",width=6,height=4)

plt_afr = ggplot(pcs_anc %>% filter(ancestry_pred == "afr"),aes(x=PC1_C,y=PC2_C,color=ancestry_pred)) + geom_point() + 
    theme(text = element_text(size=20)) + xlab("PC1") + ylab("PC2") +
    scale_color_viridis(option = "inferno",discrete=T) + xlim(c(-0.002,0.005)) + ylim(c(-0.0025,0.01))
ggsave(plot=plt_afr,filename="Output/FigureS1B.png",width=6,height=4)

plt_amr = ggplot(pcs_anc %>% filter(ancestry_pred == "amr"),aes(x=PC1_C,y=PC2_C,color=ancestry_pred)) + geom_point() + 
    theme(text = element_text(size=20)) + xlab("PC1") + ylab("PC2") +
    scale_color_viridis(option = "inferno",discrete=T) + xlim(c(-0.002,0.005)) + ylim(c(-0.0025,0.01))
ggsave(plot=plt_amr,filename="Output/FigureS1C.png",width=6,height=4)
