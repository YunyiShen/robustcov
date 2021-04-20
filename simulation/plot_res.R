library(ggplot2)
res_files_all <- list.files("./simulation/Res", full.names = TRUE) 
res_all <- lapply(res_files_all, read.csv)
res_all <- Reduce(rbind, res_all)

res_all_p120 <- res_all[res_all$p==120 & !is.na(res_all$distribution),]

ggplot(res_all_p120, aes(x=cov_measure, y = log(norm_Omega)))+
  geom_boxplot() + 
  geom_point(alpha=0.1, size=1)+
  facet_grid(graph~distribution, scales = "free_y")+
  xlab("cov estimation") + 
  ylab("log L_inf loss of precision") + 
  theme(text = element_text(size=14), 
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.margin = margin(.15, .15, .15, .15, "cm"))
ggsave("./simulation/Figs/loss_Omega_120.pdf", width = 14, height = 9, scale = 0.8)

ggplot(res_all_p120, aes(x=cov_measure, y = log(norm_Sigma)))+
  geom_boxplot() + 
  geom_point(alpha=0.1, size=1)+
  facet_grid(graph~distribution, scales = "free_y")+
  xlab("cov estimation") + 
  ylab("log L_inf loss of covariance") + 
  theme(text = element_text(size=14), 
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.margin = margin(.15, .15, .15, .15, "cm"))
ggsave("./simulation/Figs/loss_Sigma_120.pdf", width = 14, height = 9, scale = 0.8)


ggplot(res_all_p120[!is.na(res_all_p120$FP), ], aes(x=cov_measure, y = FP))+
  geom_boxplot() + 
  geom_point(alpha=0.1, size=1)+
  facet_grid(graph~distribution, scales = "free_y")+
  xlab("cov estimation") + 
  ylab("False positive rate") + 
  theme(text = element_text(size=14), 
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.margin = margin(.15, .15, .15, .15, "cm"))

ggsave("./simulation/Figs/FP_120.pdf", width = 14, height = 6, scale = 0.8)


ggplot(res_all_p120[!is.na(res_all_p120$FN), ], aes(x=cov_measure, y = FN))+
  geom_boxplot() + 
  geom_point(alpha=0.1, size=1)+
  facet_grid(graph~distribution, scales = "free_y")+
  xlab("cov estimation") + 
  ylab("False negative rate") + 
  theme(text = element_text(size=14), 
        axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
        plot.margin = margin(.15, .15, .15, .15, "cm"))
ggsave("./simulation/Figs/FN_120.pdf", width = 14, height = 9, scale = 0.8)
