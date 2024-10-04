Z.slp.cats$z.cat<-NA

Z.slp.cats[which(Z.slp.cats$z>0 & Z.slp.cats$z<=3),8] <-"poor"
Z.slp.cats[which(Z.slp.cats$z>3 & Z.slp.cats$z<=5),8] <-"low"
Z.slp.cats[which(Z.slp.cats$z>5 & Z.slp.cats$z<=10),8] <-"intermediate"
Z.slp.cats[which(Z.slp.cats$z>10),8] <-"rich"

Z.slp.cats |> 
  filter(cat=="C") |> 
  ggplot(aes(x=slope, y=factor(z.cat, levels=c("poor", "low", "intermediate", "rich"))))+
  geom_boxplot()+
  labs(y="Species richness", x="Terrain slope") +
  scale_x_continuous(trans="pseudo_log", limits=c(0, 35))

names(Z.slp.cats)
