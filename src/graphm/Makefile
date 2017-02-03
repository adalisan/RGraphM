CXX ?= g++
CFLAGS = -Wno-deprecated -O3 -fPIC
LDFLAGS= -lgsl -lgslcblas
SHVER = 1
all: ./bin/graphm
./bin/graphm: ./algorithm_ca.o ./algorithm_path.o ./algorithm_ext.o ./algorithm_fsol.o  ./algorithm_iden.o ./algorithm_lp.o ./algorithm.o ./algorithm_path.o ./algorithm_qcv.o ./algorithm_rand.o ./algorithm_rank.o ./algorithm_umeyama.o ./algorithm_unif.o ./experiment.o ./graph.o ./hungarian.o  ./rpc.o main.cpp
	$(CXX) $(CFLAGS)  main.cpp ./*o $(LDFLAGS) -o ./bin/graphm
#interface classes
./experiment.o:   experiment.cpp experiment.h
		$(CXX) $(CFLAGS) -c experiment.cpp  -o ./experiment.o
./rpc.o:   rpc.cpp rpc.h
		$(CXX) $(CFLAGS) -c rpc.cpp  -o ./rpc.o
./graph.o:   graph.cpp graph.h
		$(CXX) $(CFLAGS) -c graph.cpp  -o ./graph.o
./hungarian.o:   hungarian.cpp hungarian.h
		$(CXX) $(CFLAGS) -c hungarian.cpp  -o ./hungarian.o
./algorithm.o:   algorithm.cpp algorithm.h
		$(CXX) $(CFLAGS) -c algorithm.cpp  -o ./algorithm.o
#graph matching algorithms
./algorithm_ca.o: algorithm_ca.cpp algorithm_ca.h
		$(CXX) $(CFLAGS) -c algorithm_ca.cpp -o ./algorithm_ca.o
./algorithm_path.o: algorithm_path.cpp algorithm_path.h
		$(CXX) $(CFLAGS) -c algorithm_path.cpp -o ./algorithm_path.o
./algorithm_ext.o:   algorithm_ext.cpp algorithm_ext.h
		$(CXX) $(CFLAGS) -c algorithm_ext.cpp  -o ./algorithm_ext.o
./algorithm_fsol.o:   algorithm_fsol.cpp algorithm_fsol.h
		$(CXX) $(CFLAGS) -c algorithm_fsol.cpp  -o ./algorithm_fsol.o
./algorithm_iden.o:   algorithm_iden.cpp algorithm_iden.h
		$(CXX) $(CFLAGS) -c algorithm_iden.cpp  -o ./algorithm_iden.o
./algorithm_lp.o:   algorithm_lp.cpp algorithm_lp.h
		$(CXX) $(CFLAGS) -c algorithm_lp.cpp  -o ./algorithm_lp.o
./algorithm_qcv.o:   algorithm_qcv.cpp algorithm_qcv.h
		$(CXX) $(CFLAGS) -c algorithm_qcv.cpp  -o ./algorithm_qcv.o
./algorithm_rand.o:   algorithm_rand.cpp algorithm_rand.h
		$(CXX) $(CFLAGS) -c algorithm_rand.cpp  -o ./algorithm_rand.o
./algorithm_rank.o:   algorithm_rank.cpp algorithm_rank.h
		$(CXX) $(CFLAGS) -c algorithm_rank.cpp  -o ./algorithm_rank.o
./algorithm_umeyama.o:   algorithm_umeyama.cpp algorithm_umeyama.h
		$(CXX) $(CFLAGS) -c algorithm_umeyama.cpp  -o ./algorithm_umeyama.o
./algorithm_unif.o:   algorithm_unif.cpp algorithm_unif.h
		$(CXX) $(CFLAGS) -c algorithm_unif.cpp  -o ./algorithm_unif.o

clean:
	rm -f ./*o
