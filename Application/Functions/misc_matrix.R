
# get_internal_validation_matrix <- function(fit){
#   pClY = poLCA.posterior(fit,y=X)
#   ng <- ncol(pClY)
#   pred <-as.numeric(apply(pClY,1,function(x)which.max(x)))# MODE
#   
#   
#   
#   Ptable <- cbind(pClY,pred)
#   Pmatrix<-matrix(0,ng,ng)
#   Npmatrix <- matrix(0,ng,ng)
#   modclass <- pred
#   
#   
#   for (i in 1:ng){
#     for (j in 1:ng){
#       Pmatrix[i,j]<-sum(subset(Ptable,modclass==i)[,j])/table(modclass)[i]
#       Npmatrix[i,j]<-Pmatrix[i,j]*table(modclass)[i]
#     }}
#   
#   
#   denom<-colSums(Npmatrix)
#   Qmatrix<-matrix(0,ng,ng)
#   
#   
#   for (i in 1:ng){
#     for (j in 1:ng){
#       Qmatrix[j,i]<-Npmatrix[i,j]/denom[j]
#     }}
#   
#   Qmatrix
#   
# }
get_internal_validation_matrix <- function(fit,X,covs=NULL){
  pClY = poLCA.posterior(fit,y=X,x=covs)
  
  
  ng <- ncol(pClY)
  pred <-as.numeric(apply(pClY,1,function(x)which.max(x)))# MODE
  
  
  
  Ptable <- cbind(pClY,pred)
  Pmatrix<-matrix(0,ng,ng)
  Npmatrix <- matrix(0,ng,ng)
  modclass <- pred
  
  
  for (i in 1:ng){
    for (j in 1:ng){
      Pmatrix[i,j]<-sum(subset(Ptable,modclass==i)[,j])/table(modclass)[i]
      Npmatrix[i,j]<-Pmatrix[i,j]*table(modclass)[i]
    }}
  
  
  denom<-colSums(Npmatrix)
  Qmatrix<-matrix(0,ng,ng)
  
  
  for (i in 1:ng){
    for (j in 1:ng){
      Qmatrix[j,i]<-Npmatrix[i,j]/denom[j]
    }}
  
  Qmatrix
  
}
add_death<- function(matrix_misc, n){
  matrix_nxn <- as.matrix(matrix_misc)
  matrix_extended <- matrix(0, nrow = n + 1, ncol = n + 1)
  matrix_extended[1:n, 1:n] <- matrix_nxn
  matrix_extended[n + 1, n + 1] <- 1
  return(matrix_extended)
}
calculate_m_matrix <- function(pop, index){
  pop_mp_base <- pop %>% 
    filter(visit_number==index) %>%
    dplyr::select(MP, MP_sim)
  
  misc<- table(class.true=pop_mp_base$MP_sim, class.assigned=pop_mp_base$MP)
  misc <- misc/rowSums(misc)
  n <- dim(sim_obj$tmat)[1]
  misc <- add_death(misc, n-1)
  rownames(misc) <- rownames(sim_obj$tmat)
  colnames(misc) <- rownames(sim_obj$tmat)
  return(misc)
}

calc_approx_m_matrix <- function(pop, index){
  pop_mp_base <- pop %>%
    group_by(dataset_id, patient_id) %>%
    slice(index) %>%
    #dplyr::select(MP, MP_base)%>%
    ungroup()
  X <- pop_mp_base %>% dplyr::select(any_of(colnames(sim_obj$pattern_obj$obj$y))) %>% mutate_all(function(x)x+1)
  misc_approx <- NULL
  if (scenario %in% c("A", "Av2")){
    misc_approx <- get_internal_validation_matrix(sim_obj$pattern_obj$obj, X)
  } else if (scenario %in% c("B")){
    misc_approx <- get_internal_validation_matrix(sim_obj$pattern_obj$obj, X, covs =as.matrix(cbind(pop_mp_base$age), ncol=1) )
  } else {
    print("Insert valid scenario")
  }
  n <- dim(sim_obj$tmat)[1]
  misc_approx <- add_death(misc_approx, n-1)
  rownames(misc_approx) <- rownames(sim_obj$tmat)
  colnames(misc_approx) <- rownames(sim_obj$tmat)
  return(misc_approx)
}

perform_LCA = FALSE
if (perform_LCA){
  q_1 <- calculate_m_matrix(pop,1)
  q_2 <- calculate_m_matrix(pop,2)
  q_3 <- calculate_m_matrix(pop,3)

  print(q_1)
  print(q_2)
  print(q_3)


  pop_mp_base <- pop %>%
    group_by(dataset_id, patient_id) %>%
    slice(1) %>%
    #dplyr::select(MP, MP_base)%>%
    ungroup()
  X <- pop_mp_base %>% dplyr::select(any_of(colnames(sim_obj$pattern_obj$obj$y))) %>% mutate_all(function(x)x+1)
  if (scenario %in% c("A", "Av2")){
    misc_approx <- get_internal_validation_matrix(sim_obj$pattern_obj$obj, X)
  } else if (scenario %in% c("B")){
    misc_approx <- get_internal_validation_matrix(sim_obj$pattern_obj$obj, X, covs =as.matrix(cbind(pop_mp_base$age), ncol=1) )
  } else {
      print("Insert valid scenario")
    }
  
  print(misc_approx)
  #matrix_reordered
  #mean(diag(matrix_reordered))

}
