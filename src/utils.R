input_ini <- function(X, inl, i){
    input_param <- c(X$yTi[1] ,  inl[["Ti_lb"]][i],  inl[["Ti_ub"]][i],  # Ti
                   inl[["Tm_ini"]][i] ,  inl[["Tm_lb"]][i],  inl[["Tm_ub"]][i],  # Tm
                   inl[["Te_ini"]][i] ,  inl[["Te_lb"]][i],  inl[["Te_ub"]][i],  # Te
                   inl[["Th_ini"]][i],    -1E-5,  inl[["Th_ub"]][i],  # Th
                   X$yTi[1] ,  inl[["Ti_lb"]][i],  inl[["Ti_ub"]][i],  # Ts
                   inl[["Ci_ini"]][i] , 1E-1,  inl[["Ci_ub"]][i],  # Ci
                   2  , 1E-1,  20,  # Cm
                   2  , 1E-1,  20,  # Ce
                   1  , 1E-1,  10,  # Ch
                   1  , 1E-1,  20,  # Cs
                   inl[["Rie_ini"]][i] , 1E-2,  50,  # Rie
                   inl[["Ria_ini"]][i] , 1E-2,  50,  # Ria
                   inl[["Rea_ini"]][i] , 1E-2,  50,  # Rea
                   inl[["Rim_ini"]][i] , 1E-2,  50,  # Rim
                   inl[["Rih_ini"]][i] , 1E-2,  50,  # Rih
                   inl[["Ris_ini"]][i] , 1E-2,  50,  # Ris
                   inl[["Rh_ini"]][i] , 1E-2,  10,  # Rh
                   5  ,  0.1, 200,  # Aw
                   5  ,  0.1, 200,  # Ae
                   inl[["p_ini"]][i] ,  inl[["p_lb"]][i],  inl[["p_ub"]][i],  # P11
                   inl[["p_ini"]][i] ,  inl[["p_lb"]][i],  inl[["p_ub"]][i],  # P22
                   inl[["p_ini"]][i] ,  inl[["p_lb"]][i],  inl[["p_ub"]][i],  # P33
                   inl[["p_ini"]][i] ,  inl[["p_lb"]][i],  inl[["p_ub"]][i],  # P44
                   inl[["p_ini"]][i] ,  inl[["p_lb"]][i],  inl[["p_ub"]][i],  # P55
                   inl[["e_ini"]][i] ,  inl[["e_lb"]][i],  inl[["e_ub"]][i])  # e11
  names(input_param) <- c("Ti_ini" , "Ti_lb" , "Ti_ub" ,
                          "Tm_ini" , "Tm_lb" , "Tm_ub" , 
                          "Te_ini" , "Te_lb" , "Te_ub" ,
                          "Th_ini" , "Th_lb" , "Th_ub" ,
                          "Ts_ini" , "Ts_lb" , "Ts_ub" ,
                          "Ci_ini" , "Ci_lb" , "Ci_ub" ,
                          "Cm_ini" , "Cm_lb" , "Cm_ub" ,
                          "Ce_ini" , "Ce_lb" , "Ce_ub" ,
                          "Ch_ini" , "Ch_lb" , "Ch_ub" ,
                          "Cs_ini" , "Cs_lb" , "Cs_ub" ,
                          "Rie_ini", "Rie_lb", "Rie_ub",
                          "Ria_ini", "Ria_lb", "Ria_ub",
                          "Rea_ini", "Rea_lb", "Rea_ub",
                          "Rim_ini", "Rim_lb", "Rim_ub",
                          "Rih_ini", "Rih_lb", "Rih_ub",
                          "Ris_ini", "Ris_lb", "Ris_ub",
                          "Rh_ini" , "Rh_lb" , "Rh_ub" ,
                          "Aw_ini" , "Aw_lb" , "Aw_ub" ,
                          "Ae_ini" , "Ae_lb" , "Ae_ub" ,
                          "p11_ini", "p11_lb", "p11_ub",
                          "p22_ini", "p22_lb", "p22_ub",
                          "p33_ini", "p33_lb", "p33_ub",
                          "p44_ini", "p44_lb", "p44_ub",
                          "p55_ini", "p55_lb", "p55_ub",
                          "e11_ini", "e11_lb", "e11_ub")
  return(input_param)
}