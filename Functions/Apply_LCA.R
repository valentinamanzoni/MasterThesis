
apply_LCA <- function(pop,scenario_obj, scenario){
  
  colnames(scenario_obj$pattern_obj$obj$y)<- tolower(colnames(scenario_obj$pattern_obj$obj$y))
  X <- pop %>%  dplyr::select(any_of(colnames(scenario_obj$pattern_obj$obj$y))) %>% dplyr::mutate_all(function(x)x+1)
  if (scenario %in% c("A", "Av2")){
    post <- poLCA.posterior(scenario_obj$pattern_obj$obj,y=X)
  } else if (scenario %in% c("B")){
    post <- poLCA.posterior(scenario_obj$pattern_obj$obj,y=X,x=as.matrix(cbind(pop$age), ncol=1))
  }
  pop$MP <-as.numeric(apply(post,1,which.max))# MODE
  
  # creare tante versioni dello stesso dataset (30-50) tutto uguale tranne MP 
  # mp_list <- list()
  # # array 5x5x50 per salvare table
  # table_array <- array(data=NA, dim= c(5,5,200))
  # table_ind <- array(data=NA, dim =c(5,5,200))
  # pop2 <- pop
  # for (i in 1:200){
  #     mp_list[[i]]<- unlist(lapply(1:nrow(post), function(x) which(rmultinom(n=1, size=1, prob=post[x,])==TRUE)))
  #     pop2$MP <- mp_list[[i]]
  #     table_array[,,i] <- statetable.msm(MP, patient_id, data=pop2)
  #     table_ind[,,i] <- table_array[,,i]>0
  #     
  # } 
  # table_array
  # table_mean <- apply(table_array, 1:2, mean)
  # table_mean
  
  # renaming MP with right order:
  #scenario_obj$right_order # for renaming MP
  n <- dim(scenario_obj$tmat)[1]
  recode_map <- setNames(1:length(scenario_obj$right_order), scenario_obj$right_order)
  pop <- pop %>%
    mutate(MP = recode(MP, !!!recode_map))
  # move MP col next to MP_sim
  pop <- pop %>%
    dplyr::select(
      1:which(names(.) == "MP_sim"), "MP", setdiff(names(.), c("MP", "MP_sim"))   
    )
  
  return(pop)
  
}



