##Modelos de nicho
##instalación de KUENM

# Installing and loading packages
# if(!require(devtools)){
#   install.packages("devtools")
# }
# 
# if(!require(kuenm)){
#   devtools::install_github("marlonecobos/kuenm")
# }

library(kuenm)
#install.packages("betareg")
#sass, hier.part


# Bucle KUENM -------------------------------------------------------------

library(kuenm)
spp.carp.list<-list.files("D:Mauro/modelos/M_spec.vbles_DD/", full.names = T)

DATOS <- spp.carp.list
#DATOS <- paste("C:/DOCUMENTOS/Project folder/", list.files("C:/DOCUMENTOS/Project folder/"), sep="")
#DATOS <- paste("D:/Maestria_DD/spp. records_DD/specialists_DD/spec. shapes_DD/M_spec.vbles_DD/", list.files("D:/Maestria_DD/spp. records_DD/specialists_DD/spec. shapes_DD/M_spec.vbles_DD/"), sep="")

for(i in 1:length(DATOS)){
  
  setwd(DATOS[i])
  
  
  ##The next chunk of code is for preparing the arguments for using the function following the modularity principle. These variables can be changed according to each case.
  occ_joint <- "Sp_joint.csv"
  occ_tra <- "Sp_train.csv"
  M_var_dir <- "M_variables"
  batch_cal <- "Candidate_models"
  out_dir <- "Candidate_Models"
  reg_mult <- c(seq(0.1, 1, 0.1), seq(2, 6, 1), 8, 10)
  f_clas <- c("l","lq","lqp", "lp")
  background <- 10000
  maxent_path <- "C:/maxent"
  wait <- FALSE
  run <- TRUE
  
  ##The following is the code for using the function.
  
  kuenm_cal(occ.joint = occ_joint, occ.tra = occ_tra, M.var.dir = M_var_dir, batch = batch_cal,
            out.dir = out_dir, reg.mult = reg_mult, f.clas = f_clas, 
            maxent.path = maxent_path, wait = wait, run = run)
  
  
  ## Evaluation and selection of best models
  
  occ_test <- "Sp_test.csv"
  out_eval <- "Calibration_results"
  threshold <- 5
  rand_percent <- 50
  iterations <- 500
  kept <- TRUE #eliminar modelos candidatos
  selection <- "OR_AICc"
  paral_proc <- FALSE 
  
  #This code also allows evaluating candidate models that were created previously, selecting those with best performance based on the three criteria.
  
  cal_eval <- kuenm_ceval(path = out_dir, occ.joint = occ_joint, occ.tra = occ_tra, occ.test = occ_test,
                          batch = batch_cal, out.eval = out_eval, threshold = threshold,
                          rand.percent = rand_percent, iterations = iterations, kept = kept,
                          selection = selection, parallel.proc = paral_proc)
  
  ### Final model creation
  
  #For preparing the arguments for this function use the following chunk of code.
  
  batch_fin <- "Final_models"
  mod_dir <- "Final_Models"
  rep_n <- 10
  rep_type <- "Bootstrap"
  jackknife <- FALSE
  out_format <- "logistic"
  project <- FALSE
  G_var_dir <- "G_variables"
  ext_type <- "all"
  write_mess <- FALSE
  write_clamp <- FALSE
  wait1 <- FALSE
  run1 <- TRUE
  args <- NULL
  
  #The kuenm_mod function has the following syntax:
  kuenm_mod(occ.joint = occ_joint, M.var.dir = M_var_dir, out.eval = out_eval, batch = batch_fin,
            rep.n = rep_n, rep.type = rep_type, jackknife = jackknife, out.dir = mod_dir,
            out.format = out_format, project = project, G.var.dir = G_var_dir, ext.type = ext_type,
            write.mess = write_mess, write.clamp = write_clamp, maxent.path = maxent_path,
            args = args, wait = wait1, run = run1)
}




