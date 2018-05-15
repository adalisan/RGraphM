# //*********************GRAPHS**********************************
#   //graph_1,graph_2 are graph adjacency matrices, 
# //C_matrix is the matrix of local similarities  between vertices of graph_1 and graph_2. 
# //If graph_1 is NxN and graph_2 is MxM then C_matrix should be NxM
# graph_1=../test_simple/g s
# graph_2=../test_simple/h s
# C_matrix=../test_simple/null s
# //*******************ALGORITHMS********************************
#   //used algorithms and what should be used as initial solution in corresponding algorithms
# algo=I U RANK QCV rand PATH s
# algo_init_sol=unif unif unif unif unif unif s
# solution_file=solution_im.txt s
# //coeficient of linear combination between (1-alpha_ldh)*||graph_1-P*graph_2*P^T||^2_F +alpha_ldh*C_matrix 
# alpha_ldh=0 d
# cdesc_matrix=A c
# cscore_matrix=A c
# //**************PARAMETERS SECTION*****************************
#   hungarian_max=10000 d
# algo_fw_xeps=0.01 d
# algo_fw_feps=0.01 d
# //0 - just add a set of isolated nodes to the smallest graph, 1 - double size 
# dummy_nodes=0 i
# // fill for dummy nodes (0.5 - these nodes will be connected with all other by edges of weight 0.5(min_weight+max_weight))
# dummy_nodes_fill=0 d
# // fill for linear matrix C, usually that's the minimum (dummy_nodes_c_coef=0),
# // but may be the maximum (dummy_nodes_c_coef=1)
# dummy_nodes_c_coef=0.01 d
# 
# qcvqcc_lambda_M=10 d
# qcvqcc_lambda_min=1e-5 d
# 
# 
# //0 - all matching are possible, 1-only matching with positive local similarity are possible
# blast_match=1 i
# blast_match_proj=0 i
# 
# 
# //****************OUTPUT***************************************
# //output file and its format 
# exp_out_file=../test_simple/exp_out_file s
# exp_out_format=Parameters Compact Permutation s
# //other
# debugprint=0 i
# debugprint_file=debug.txt s
# verbose_mode=1 i
# //verbose file may be a file or just a screen:cout
# verbose_file=cout s

test_algo_pars<-list(
  # Already provided as A and B matrices
  # *******************ALGORITHMS********************************
  # used algorithms and what should be used as
  #initial solution in corresponding algorithms
  algo="I QCV PATH",
  algo_init_sol="unif unif unif",
  solution_file="solution_im.txt",
  # coeficient of linear combination between
  # (1-alpha_ldh)*||graph_1-P*graph_2*P^T||^2_F +alpha_ldh*C_matrix
  alpha_ldh=0 ,
  cdesc_matrix="A" ,
  cscore_matrix="A" ,
  C_matrix = "none",
  # **************PARAMETERS SECTION*****************************
  hungarian_max=10000 ,
  algo_fw_xeps=0.01 ,
  algo_fw_feps=0.01 ,
  # 0 - just add a set of isolated nodes to the smallest graph
  # 1 - double size
  dummy_nodes=as.integer(0),
  # fill for dummy nodes (0.5 - these nodes will be connected with all other
  # by edges of weight 0.5(min_weight+max_weight))
  dummy_nodes_fill=as.integer(0),
  # fill for linear matrix C, usually that's the minimum (dummy_nodes_c_coef=0),
  # but may be the maximum (dummy_nodes_c_coef=1)
  dummy_nodes_c_coef=0.01,
  
  qcvqcc_lambda_M=10,
  qcvqcc_lambda_min=1E-5,
  
  
  # 0 - all matching are possible, 1-only matching with positive local similarity are possible
  blast_match=as.integer(1) ,
  blast_match_proj=as.integer(0) ,
  
  
  
  #****************OUTPUT***************************************
  #output file and its format
  exp_out_file="exp_out_file" ,
  exp_out_format="Parameters Compact Permutation",
  #other
  debugprint=as.integer(1) ,
  debugprint_file="debug.txt",
  verbose_mode=as.integer(1) ,
  # verbose file may be a file or just a screen:cout
  verbose_file="verbose_debug.txt",
  graph_dot_print = as.integer(1)
)


test_graph_match <- function(A,B,algorithm_params) {
  writeMat(A,"./src/graphm/Rpkgtest/testA")
  writeMat(B,"./src/graphm/Rpkgtest/testB")
  write_config("./src/graphm/Rpkgtest/config.txt",A,B,algorithm_params)
  system("./src/graphm/bin/graphm.exe ./src/graphm/Rpkgtest/config.txt",show.output.on.console = T)
}

writeMat <- function (A,filepath) {
  
  write.table(A,file = filepath,row.names= F, col.names=F)
}
write_config <- function  (filepath,A,B,algo_params){
  new_f = file(filepath,open = 'wt')
  #open(new_f)
  writeLines(text= "graph_1=./src/graphm/Rpkgtest/testA s",new_f)
  writeLines(text= "graph_2=./src/graphm/Rpkgtest/testB s",new_f)
  close(new_f)
  for ( par_name in names(algo_params)) {
    par  = algo_params[[par_name]]
    print (par)
    if (is.numeric(par)) {
      if (is.double(par)) {
      type = 'd'
      }
      else if  (is.integer(par)) {
        type = 'i'
      }
    } else if (is.character(par)){
      if (nchar(par)==0) {
        print ("invalid parameter value for " + names(par))
      }
      if (nchar(par)==1) {
      type = 'c'
      } else {
        type = 's'
      }
      
    }
    write(paste0(par_name,'=', par ," " ,type ,collapse = ""),file = filepath, append = TRUE)
    
  }
  
}
