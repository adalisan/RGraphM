\name{run_graph_match}
\alias{run_graph_match}
\title{Matches two graphs with adjacency matrices}
\usage{
run_graph_match(A,B,algorithm_params)
}
\arguments{
\item{A}{ the adjacency matrix for first graph }
\item{B}{ the adjacency matrix for second graph}
\item{algorithm_params}{ A named list containing the parameters of the graph matching algorithm(s). The parameters can be integer,numeric or strings. See example and the source code for graphm package.}
}
\value{
A named list containing
\code{debugprint_file} Name of debug output  file
\code{algo_names} The names of GM algorithms tested in experiment
\code{Pmat} The list of  permutation matrices as found by the graph matching algorithms for matching graph \code{B} to \code{A}.
\code{exp_count} The number of GM algorithms tested in experiment
\code{Pvec} The  list of  index permutation vectors as found by the graph matching algorithms for matching graph \code{B} to \code{A}.
}
\description{
Runs the graph matching algorithm for graphs represented by \code{A} and \code{B}. Returns an \R \code{list} containing the permutation matrix \code{Pmat} and other computed variables.
}
\examples{
B<-matrix( sample(c(0,1),225, replace = TRUE,prob=c(0.75,0.25)),15,15)
diag(B)<- 1
perm <- sample.int(15)
B<-(B+t(B))/2
B[B>0.1] <- 1
A<-B [perm,][,perm]

# parameters copied from test_simple example or graphm
algorithm_params<-list(
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
debugprint=1 ,
debugprint_file="debug.txt",
verbose_mode=1 ,
# verbose file may be a file or just a screen:cout
verbose_file="verbose_debug.txt",
graph_dot_print = 1
)
print( algorithm_params)
print ("Dimensions of A")
print (dim(A))
print ("Dimensions of B")
print (dim(B))
if ( dim(A)[1]>0 && dim(B)[1]>0) {
 result <- run_graph_match(A, B, algorithm_params)

 print ("Original permutation")
 print(perm)
 if (!is.null(result) && !is.null(result$algo_names)) {
 		for (algo_i in 1:length(result$algo_names) ) {
 		algo = result$algo_names[algo_i]
 		perm_vector <- result$Pvec[[algo_i]]
 		print (paste(algo ,perm_vector , sep = " ",collapse= " "))
 		}
 }
}
#Example for matching weighted graphs. Disabled to reduce elapsed time for R check
if (FALSE){



 algorithm_params["algo"]<- "IDEN UMEY RANK QCV RAND PATH  CA"
 algorithm_params["algo_init_sol"] <- "unif unif unif unif unif unif"
 algorithm_params["alpha_ldh"] <- 0.44
 anotherresult <- run_graph_match(a,b, algorithm_params)

C<-matrix( runif(	min = 0, max = 1,n= 225),15,15)
diag(C)<- 1
B <- C*B
perm <- sample.int(15)
A<-B [perm,perm]
weightedmatrixresult <- run_graph_match(A,B, algorithm_params)
}
}
