CXX	:= c++
CXXFLAGS := -Wall -Wextra -MMD -O2 -m64 -g 
INCLUDES := /home/gurobi/12.0.0/linux64/include/ ../auxiliar/eigen-5.0.0
LIBSINCLUDES := /home/gurobi/12.0.0/linux64/lib/ 
LIBS = gurobi_c++ gurobi120 m

CPPFLAGS += $(foreach includedir,$(INCLUDES),-I$(includedir))
LDFLAGS += $(foreach librarydir,$(LIBSINCLUDES),-L$(librarydir))
LDFLAGS += $(foreach library,$(LIBS),-l$(library))

OBJECTS := main.o input/input.o instance/instance.o instance/bilevelmodel/bilevelmodel.o leadersolver/status/status.o leadersolver/leadersolver.o followersolver/followersolver.o followersolver/fixstrongweak/fixstrongweakfollower.o followersolver/depstrongweak/saadepstrongweakfollower.o followersolver/depstrongweak/dedepstrongweakfollower.o followersolver/depgeneral/depgeneralfollower.o

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