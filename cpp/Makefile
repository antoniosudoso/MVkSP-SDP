MATLAB_PATH = /usr/local/MATLAB/R2022b
GUROBI_PATH = /home/antonio/gurobi1300

sources = main.cpp matlab_util.cpp JobQueue.cpp sdp_branch_and_bound.cpp ThreadPool.cpp
objects = $(sources:.cpp=.o)

all: bb

bb: $(objects)
	g++ -L $(MATLAB_PATH)/extern/bin/glnxa64/ -L $(GUROBI_PATH)/lib/ -o $@ $^ -larmadillo -pthread -lMatlabDataArray -lMatlabEngine -Wl,-rpath,$(MATLAB_PATH)/extern/bin/glnxa64/

%.o: %.cpp
		g++ -g -Wall -O2 -o $@ -c $< -I $(MATLAB_PATH)/extern/include/ -I $(GUROBI_PATH)/include/


.PHONY: clean
		
clean:
		rm -f bb *.o
