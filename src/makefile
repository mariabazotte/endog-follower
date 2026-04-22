CXX	:= c++
CXXFLAGS := -Wall -Wextra -MMD -O2 -m64 -g -std=c++17
# 
# -O0 -g -fsanitize=address -fno-omit-frame-pointer
INCLUDES := /home/gurobi/12.0.0/linux64/include/ ../auxiliar/eigen-5.0.0 
LIBSINCLUDES := /home/gurobi/12.0.0/linux64/lib/ 
LIBS = gurobi_c++ gurobi120 m

CPPFLAGS += $(foreach includedir,$(INCLUDES),-I$(includedir))
LDFLAGS += $(foreach librarydir,$(LIBSINCLUDES),-L$(librarydir))
LDFLAGS += $(foreach library,$(LIBS),-l$(library))

PYBIND11_INCLUDE := /home/bazmar/.local/lib/python3.12/site-packages/pybind11/include
PYTHON_INCLUDE := /cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/python/3.12.4/include/python3.12
PYTHON_LIB_DIR := /cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v3/Compiler/gcccore/python/3.12.4/lib
PYTHON_LIBS := -lpython3.12 -lpthread -ldl -lutil -lm

CPPFLAGS += -I$(PYBIND11_INCLUDE) -I$(PYTHON_INCLUDE)
LDFLAGS += -L$(PYTHON_LIB_DIR) $(PYTHON_LIBS)

JL_SHARE = $(shell julia -e 'print(joinpath(Sys.BINDIR, Base.DATAROOTDIR, "julia"))')
CPPFLAGS += $(shell $(JL_SHARE)/julia-config.jl --cflags)
LDFLAGS  += $(shell $(JL_SHARE)/julia-config.jl --ldflags)
LDFLAGS  += $(shell $(JL_SHARE)/julia-config.jl --ldlibs)

OBJECTS := main.o \
		input/input.o \
		instance/instance.o \
		instance/bilevelmodel/bilevelmodel.o \
		instance/performance/performance.o \
		leadersolver/status/status.o \
		leadersolver/leadersolver.o \
		followersolver/followersolver.o \
		followersolver/fixstrongweak/fixstrongweakfollower.o \
		followersolver/depstrongweak/saadepstrongweakfollower.o \
		followersolver/depstrongweak/dedepstrongweakfollower.o \
		followersolver/depgeneral/depgeneralfollower.o

DIRS = $(OBJECTS:.o=.d)

EXE:=exe

$(EXE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $(EXE) $(OBJECTS) $(LDFLAGS) 

%.o: %.cxx
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@ $(LDFLAGS)

clean:
	$(RM) *.o *~ $(EXE)
	rm -rf $(OBJECTS)
	rm -rf $(DIRS)