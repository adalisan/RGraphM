\name{run_graph_match}
\alias{run_graph_match}
\title{Matches two graphs with adjacency matrices}
\usage{
run_graph_match(A,B,algorithm_params)
}
\description{
Returns an \R \code{list} containing the character vector
\code{c("foo", "bar")} and the numeric vector \code{c(0, 1)}.
}
\examples{
A<-matrix( sample(c(0,1),100, replace = TRUE),10,10)
perm <- sample.int(10)
B<-A [perm,perm]

# parameters copied from test_simple example or graphm
algorithm_params<-list(
# Already provided as A and B matrices
# *******************ALGORITHMS********************************
# used algorithms and what should be used as
#initial solution in corresponding algorithms
algo="PATH",
algo_init_sol="unif",
solution_file="solution_im.txt",
# coeficient of linear combination between
# (1-alpha_ldh)*||graph_1-P*graph_2*P^T||^2_F +alpha_ldh*C_matrix
alpha_ldh=0 ,
cdesc_matrix="A" ,
cscore_matrix="A" ,
# **************PARAMETERS SECTION*****************************
hungarian_max=10000 ,
algo_fw_xeps=0.01 ,
algo_fw_feps=0.01 ,
# 0 - just add a set of isolated nodes to the smallest graph
# 1 - double size
dummy_nodes=0,
# fill for dummy nodes (0.5 - these nodes will be connected with all other
# by edges of weight 0.5(min_weight+max_weight))
dummy_nodes_fill=0,
# fill for linear matrix C, usually that's the minimum (dummy_nodes_c_coef=0),
# but may be the maximum (dummy_nodes_c_coef=1)
dummy_nodes_c_coef=0.01,

qcvqcc_lambda_M=10,
qcvqcc_lambda_min=1E-5,


# 0 - all matching are possible, 1-only matching with positive local similarity are possible
blast_match=1 ,
blast_match_proj=0 ,


#****************OUTPUT***************************************
#output file and its format
exp_out_file="exp_out_file" ,
exp_out_format="Parameters Compact Permutation",
#other
debugprint=0 ,
debugprint_file="debug.txt",
verbose_mode=1 ,
# verbose file may be a file or just a screen:cout
verbose_file="verbose_debug.txt"
)
print( algorithm_params)
 result <- run_graph_match(A, B, algorithm_params)
 print(result$Pmat)
}
