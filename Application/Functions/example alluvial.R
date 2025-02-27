library(ggplot2)
#install.packages("remotes")
#remotes::install_github("davidsjoberg/ggsankey")
library(ggsankey)
library(ggpubr)

nsim <- 10000
N<- 100
scenario <- "B"
p1 <- load(paste0("Simulation/Simulation Outputs/pop_study_",scenario,"_",N,"_",nsim,".RData"))
pop1<- get(p1)


pop_ms <- pop1 %>% filter(dataset_id==1)
source("Functions//Apply_LCA.R")

pop_ms<- apply_LCA(pop_ms,sim_obj, scenario)

pop_ms <- pop_ms %>%
  mutate(death = ifelse(Age_exit < Age_death, 0, 1))
pop_ms_death_1 <- pop_ms %>%
  dplyr::filter(death == 1) %>%
  group_by(dataset_id, patient_id) %>%
  group_split() %>%
  map_df(~ bind_rows(
    .x,
    slice_tail(.x, n = 1) %>% mutate(MP = dim(sim_obj$tmat)[1],MP_sim=dim(sim_obj$tmat)[1], age = Age_death)
  )) %>%
  ungroup()

pop_ms_death_0 <- pop_ms %>%
  dplyr::filter(death == 0)

pop_ms <- bind_rows(pop_ms_death_0, pop_ms_death_1)

pop_ms <- pop_ms %>%
  group_by(dataset_id,patient_id) %>%
  filter(n() > 1) %>%
  ungroup()

if (scenario=="B"){
  pop_ms <- pop_ms %>%
    group_by(patient_id) %>%
    mutate(flag = cumsum(MP == 2),
           # Adjust MP based on the flag
           MP = case_when(
             flag > 0 & MP == 1 ~ 2, # If flag is active and MP is 1, change it to 2
             TRUE ~ MP           # Otherwise, keep MP as it is
           )
    ) %>%
    ungroup() %>%
    dplyr::select(-flag) 
}

prova <- expand_grid(patient_id=unique(pop_ms$patient_id),
                     age=c(60:100))
n<- dim(sim_obj$tmat)[1]
pop_ms_all <- pop_ms %>% 
  dplyr::select(patient_id, dataset_id, MP, MP_sim, age, Age_exit, Age_entry) %>%
  mutate(age=round(age), Age_entry=round(Age_entry))%>%
  full_join(prova)%>%
  arrange(patient_id, age)%>%
  fill(Age_entry, Age_exit) %>%
  filter(age>=Age_entry) %>%
  group_by( patient_id) %>%
  fill(MP, MP_sim, .direction = "down") %>%
  mutate(MP=ifelse(age>Age_exit & MP!=n, NA_integer_, MP),
         MP_sim=ifelse(age>Age_exit & MP!=n, NA_integer_, MP_sim),
         next_age = lead(age), next_MP = lead(MP), next_MP_sim = lead(MP_sim)) %>%
  ungroup() %>%
  drop_na(MP)
mycolor2 <- c("#440154FF","#3B528BFF" ,'grey90' 
              #,"#FDE725FF", "#21908CFF" ,'grey90'
              #,'grey30'
              )

# collonna age, collonna next_age, MP  e next MP. provare a mettere NA per i next non disponibili
ggalluvial <- ggplot(pop_ms_all, aes(x = age,
                                       next_x = next_age, 
                                       node = MP, 
                                       next_node = next_MP,
                                       fill = factor(MP),
                                       label=MP,
                                       node.fill=MP)) +
  geom_sankey(flow.alpha=0.8,node.col=1) +
  scale_fill_manual("",values=mycolor2,labels=rownames(sim_obj$tmat)) +
  theme_sankey(base_size = 30)+
  scale_x_continuous("Age"
                     #breaks = unique(dat_long_ext$Age_x),
                     #labels = unique(dat_long_ext$Age_cut)
                     )+
  ggtitle("Assigned by LCA") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 18), 
        axis.text = element_text(size = 18),  
        legend.text = element_text(size = 15))
ggalluvial
gg_sim  <- ggplot(pop_ms_all, aes(x = age,
                                  next_x = next_age, 
                                  node = MP_sim, 
                                  next_node = next_MP_sim,
                                  fill = factor(MP_sim),
                                  label=MP_sim,
                                  node.fill=MP_sim)) +
  geom_sankey(flow.alpha=0.8,node.col=1) +
  scale_fill_manual("",values=mycolor2,labels=c('Mild', 'Complex', "Death")) +
  theme_sankey(base_size = 30)+
  scale_x_continuous("Age"
                     #breaks = unique(dat_long_ext$Age_x),
                     #labels = unique(dat_long_ext$Age_cut)
  )+
  ggtitle("Ground truth") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 20),
        axis.title = element_text(size = 18), 
        axis.text = element_text(size = 18),  
        legend.text = element_text(size = 15))
gg_sim
ggarrange(ggalluvial, gg_sim, nrow=2, common.legend=TRUE,heights = c(1.2, 1.2))

