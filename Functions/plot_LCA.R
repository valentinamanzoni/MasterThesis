library(ggplot2)
library(dplyr)
library(survminer)
library(tidyverse) 
library(poLCA)
library(magrittr)
library(ggprism) 
ggOE_ex <- function(obj,nclass,cutoff_OE=2,cutoff_Ex=0.25){
  E <-apply(obj$y-1,2,mean)
  n <- list()
  for (j in 1:nclass){
    n[[j]] <- apply(obj$y[obj$predclass==j,]-1,2,mean)
  }
  
  
  O <- do.call("cbind",n)
  R <- O/E
  R %<>% as.data.frame() %>% rownames_to_column("Disease") %>% pivot_longer(2:(nclass+1),
                                                                            names_to = "Multimorbidity profile",
                                                                            values_to = "O/E") %>% 
    mutate(`Multimorbidity profile`=as.numeric (gsub("\\D", "", `Multimorbidity profile`)))%>% 
    mutate(label=ifelse(`O/E`<cutoff_OE,NA_integer_,Disease))
  
  
  ######## Exclusivity ########
  
  N <-apply(obj$y-1,2,sum)
  
  n <- list()
  for (j in 1:nclass){
    n[[j]] <- apply(obj$y[obj$predclass==j,]-1,2,sum)
  }
  
  Ex <- do.call("cbind",n)
  
  Ex <- Ex/N
  Ex %<>% as.data.frame() %>% rownames_to_column("Disease") %>% pivot_longer(2:(nclass+1),
                                                                             names_to = "Multimorbidity profile",
                                                                             values_to = "Exclusivity") %>% 
    mutate(`Multimorbidity profile`=as.numeric (gsub("\\D", "", `Multimorbidity profile`)))%>% 
    mutate(label2=ifelse(`Exclusivity`<cutoff_Ex,NA_integer_,Disease))
  
  
  Char_MP <- R %>% left_join(Ex) 
  Char_MP %<>% mutate(char=ifelse(!is.na(label) & !is.na(label2),1,NA_integer_)) 
  
  ggOE <- ggplot(Char_MP)+
    geom_bar(aes(`O/E`,Disease,fill=Disease),stat = "identity")+
    geom_vline(aes(xintercept=2),linetype="dashed")+
    # geom_text(aes(5,Disease[char==1],label=label[char==1],vjust="topright"),check_overlap = T)+
    facet_grid(.~`Multimorbidity profile`)+
    # scale_color_manual(values = c("black","indianred"))+
    scale_y_discrete("Chronic conditions")+
    theme_prism()+
    theme(legend.position = "null",
          axis.text.y =element_text( hjust = 1),
          strip.text.x.top =element_text(size=16),
          axis.ticks.y = element_blank())+
    ggtitle("Multimorbidity Profiles")
  
  ggex <- ggplot(Char_MP)+
    geom_bar(aes(Exclusivity,Disease,fill=Disease),stat = "identity")+
    geom_vline(aes(xintercept=0.25),linetype="dashed")+
    # geom_text(aes(0.6,Disease[char==1],label=label2[char==1],vjust="topright"),check_overlap = T)+
    facet_grid(.~`Multimorbidity profile`)+
    # scale_color_manual(values = c("black","indianred"))+
    scale_y_discrete("Chronic conditions")+
    scale_x_continuous(limits = c(0,1))+
    theme_prism()+
    theme(legend.position = "null",
          axis.text.y = element_text( hjust = 1),
          strip.text.x.top = element_blank(),
          axis.ticks.y = element_blank())
  
  Char_MP %<>% mutate(`Multimorbidity profile`=as.factor(`Multimorbidity profile`))%>% filter(char==1) %>% 
    group_by(`Multimorbidity profile`) %>% 
    mutate(index=row_number()) 
  
  # %>% 
  #   select(MP,label2) %>% flextable::flextable() %>% flextable::merge_at(j=1)
  
  # Char_MP$x <- rep(c(0.1,0.2),floor(nrow(Char_MP)/2))[1:nrow(Char_MP)]
  
  require(ggrepel)
  ggnames <- ggplot(Char_MP)+
    geom_text(aes(0.1,index,label=label2,hjust="left"),size=8)+
    facet_grid(.~`Multimorbidity profile`,drop = F)+
    # scale_color_manual(values = c("black","indianred"))+
    scale_y_reverse("Chronic conditions")+
    scale_x_continuous(limits = c(0,1))+
    theme_void()+
    theme(legend.position = "null",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text.x.top = element_blank(),
          text = element_text(face = "bold",size = 16))
  
  datn <- data.frame(MP=1:nclass,
                     N=as.numeric(table(obj$predclass)),
                     P=round(as.numeric(table(obj$predclass))/length(obj$predclass)*100,0))
  print(datn)
  
  ggn <- ggplot(datn) +
    geom_point(aes(1,1,size=N*100),alpha=0.2)+
    geom_text(aes(1,1,size=N,label=paste(P,"%"),vjust="middle",hjust="center"),size=12)+
    facet_grid(.~MP)+
    theme_void()+
    scale_y_continuous(limits = c(0.9,1.1))+
    scale_x_continuous(limits = c(0.9,1.1))+
    theme(legend.position = "null",
          strip.text.x.top = element_blank(),
          text = element_text(face = "bold",size = 24))
  
  library(patchwork)
  # gg <- ggpubr::ggarrange(ggOE,ggex,ggnames,ggn,nrow=4,heights = c(1,1,0.5,0.5),align = "v")
  gg <- ggOE/ggex/ggnames/ggn +plot_layout(heights = c(1,1,1,0.5))
  print(gg)
  return(gg)
}

ggOE_ex(sim_obj$pattern_obj$obj,2)


ggOE <-function(obj,nclass,cutoff_OE=2,cutoff_Ex=0.5){
  E <-apply(obj$y-1,2,mean)
  n <- list()
  for (j in 1:nclass){
    n[[j]] <- apply(obj$y[obj$predclass==j,]-1,2,mean)
  }
  
  
  O <- do.call("cbind",n)
  R <- O/E
  R %<>% as.data.frame() %>% rownames_to_column("Disease") %>% pivot_longer(2:(nclass+1),
                                                                            names_to = "Multimorbidity profile",
                                                                            values_to = "O/E") %>% 
    mutate(`Multimorbidity profile`=as.numeric (gsub("\\D", "", `Multimorbidity profile`)))%>% 
    mutate(label=ifelse(`O/E`<cutoff_OE,NA_integer_,Disease))
  
  
  ######## Exclusivity ########
  
  N <-apply(obj$y-1,2,sum)
  
  n <- list()
  for (j in 1:nclass){
    n[[j]] <- apply(obj$y[obj$predclass==j,]-1,2,sum)
  }
  
  Ex <- do.call("cbind",n)
  
  Ex <- Ex/N
  Ex %<>% as.data.frame() %>% rownames_to_column("Disease") %>% pivot_longer(2:(nclass+1),
                                                                             names_to = "Multimorbidity profile",
                                                                             values_to = "Exclusivity") %>% 
    mutate(`Multimorbidity profile`=as.numeric (gsub("\\D", "", `Multimorbidity profile`)))%>% 
    mutate(label2=ifelse(`Exclusivity`<cutoff_Ex,NA_integer_,Disease))
  
  
  Char_MP <- R %>% left_join(Ex) 
  Char_MP %<>% mutate(char=ifelse(!is.na(label) & !is.na(label2),1,NA_integer_)) 
  
  ggOE <- ggplot(Char_MP)+
    geom_bar(aes(`O/E`,Disease,fill=Disease),stat = "identity")+
    geom_vline(aes(xintercept=cutoff_OE),linetype="dashed")+
    # geom_text(aes(5,Disease[char==1],label=label[char==1],vjust="topright"),check_overlap = T)+
    facet_grid(.~`Multimorbidity profile`,
               labeller = as_labeller(c("1" = "Mild", "2" = "Complex")))+
    # scale_color_manual(values = c("black","indianred"))+
    scale_y_discrete("Chronic conditions")+
    theme_prism()+
    theme(legend.position = "null",
          axis.text.y = element_text(hjust = 1, size = 8),  # Correct usage
          strip.text.x.top = element_text(size = 16),
          axis.ticks.y = element_line(size = 1))+
    ggtitle("Multimorbidity Profiles")
  
  Char_MP %<>% mutate(`Multimorbidity profile`=as.factor(`Multimorbidity profile`))%>% filter(char==1) %>% 
    group_by(`Multimorbidity profile`) %>% 
    mutate(index=row_number()) 
  
  # %>% 
  #   select(MP,label2) %>% flextable::flextable() %>% flextable::merge_at(j=1)
  
  # Char_MP$x <- rep(c(0.1,0.2),floor(nrow(Char_MP)/2))[1:nrow(Char_MP)]
  
  require(ggrepel)
  ggnames <- ggplot(Char_MP)+
    geom_text(aes(0.1,index,label=label2,hjust="left"),size=8)+
    facet_grid(.~`Multimorbidity profile`,labeller = as_labeller(c("1" = "Mild", "2" = "Complex")),
                 drop = F)+
    # scale_color_manual(values = c("black","indianred"))+
    scale_y_reverse("Chronic conditions")+
    scale_x_continuous(limits = c(0,1))+
    theme_void()+
    theme(legend.position = "null",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text.x.top = element_blank(),
          text = element_text(face = "bold",size = 16))
  
  datn <- data.frame(MP=1:nclass,
                     N=as.numeric(table(obj$predclass)),
                     P=round(as.numeric(table(obj$predclass))/length(obj$predclass)*100,0))
  print(datn)
  
  ggn <- ggplot(datn) +
    geom_point(aes(1,1,size=N*100),alpha=0.2)+
    geom_text(aes(1,1,size=N,label=paste(P,"%"),vjust="middle",hjust="center"),size=12)+
    facet_grid(.~MP)+
    theme_void()+
    scale_y_continuous(limits = c(0.9,1.1))+
    scale_x_continuous(limits = c(0.9,1.1))+
    theme(legend.position = "null",
          strip.text.x.top = element_blank(),
          text = element_text(face = "bold",size = 24))
  
  
  return(ggOE)
}
ggEx <-function(obj,nclass,cutoff_OE=2,cutoff_Ex=0.5){
  E <-apply(obj$y-1,2,mean)
  n <- list()
  for (j in 1:nclass){
    n[[j]] <- apply(obj$y[obj$predclass==j,]-1,2,mean)
  }
  
  
  O <- do.call("cbind",n)
  R <- O/E
  R %<>% as.data.frame() %>% rownames_to_column("Disease") %>% pivot_longer(2:(nclass+1),
                                                                            names_to = "Multimorbidity profile",
                                                                            values_to = "O/E") %>% 
    mutate(`Multimorbidity profile`=as.numeric (gsub("\\D", "", `Multimorbidity profile`)))%>% 
    mutate(label=ifelse(`O/E`<cutoff_OE,NA_integer_,Disease))
  
  
  ######## Exclusivity ########
  
  N <-apply(obj$y-1,2,sum)
  
  n <- list()
  for (j in 1:nclass){
    n[[j]] <- apply(obj$y[obj$predclass==j,]-1,2,sum)
  }
  
  Ex <- do.call("cbind",n)
  
  Ex <- Ex/N
  Ex %<>% as.data.frame() %>% rownames_to_column("Disease") %>% pivot_longer(2:(nclass+1),
                                                                             names_to = "Multimorbidity profile",
                                                                             values_to = "Exclusivity") %>% 
    mutate(`Multimorbidity profile`=as.numeric (gsub("\\D", "", `Multimorbidity profile`)))%>% 
    mutate(label2=ifelse(`Exclusivity`<cutoff_Ex,NA_integer_,Disease))
  
  
  Char_MP <- R %>% left_join(Ex) 
  Char_MP %<>% mutate(char=ifelse(!is.na(label) & !is.na(label2),1,NA_integer_)) 
  
  ggex <- ggplot(Char_MP)+
    geom_bar(aes(Exclusivity,Disease,fill=Disease),stat = "identity")+
    geom_vline(aes(xintercept=cutoff_Ex),linetype="dashed")+
    # geom_text(aes(0.6,Disease[char==1],label=label2[char==1],vjust="topright"),check_overlap = T)+
    facet_grid(.~`Multimorbidity profile`,labeller = as_labeller(c("1" = "Mild", "2" = "Complex")))+
    # scale_color_manual(values = c("black","indianred"))+
    scale_y_discrete("Chronic conditions")+
    scale_x_continuous(limits = c(0,1))+
    theme_prism()+
    theme(legend.position = "null",
          axis.text.y = element_text(hjust = 1, size = 8),
          strip.text.x.top = element_text(size = 16),
          axis.ticks.y = element_blank()) +
    ggtitle("Multimorbidity Profiles")
  
  Char_MP %<>% mutate(`Multimorbidity profile`=as.factor(`Multimorbidity profile`))%>% filter(char==1) %>% 
    group_by(`Multimorbidity profile`) %>% 
    mutate(index=row_number()) 
  
  # %>% 
  #   select(MP,label2) %>% flextable::flextable() %>% flextable::merge_at(j=1)
  
  # Char_MP$x <- rep(c(0.1,0.2),floor(nrow(Char_MP)/2))[1:nrow(Char_MP)]
  
  require(ggrepel)
  ggnames <- ggplot(Char_MP)+
    geom_text(aes(0.1,index,label=label2,hjust="left"),size=8)+
    facet_grid(.~`Multimorbidity profile`,drop = F)+
    # scale_color_manual(values = c("black","indianred"))+
    scale_y_reverse("Chronic conditions")+
    scale_x_continuous(limits = c(0,1))+
    theme_void()+
    theme(legend.position = "null",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text.x.top = element_blank(),
          text = element_text(face = "bold",size = 16))
  
  datn <- data.frame(MP=1:nclass,
                     N=as.numeric(table(obj$predclass)),
                     P=round(as.numeric(table(obj$predclass))/length(obj$predclass)*100,0))
  print(datn)
  return(ggex)
}

ggOE(sim_obj$pattern_obj$obj,2)
ggEx(sim_obj$pattern_obj$obj,2)
