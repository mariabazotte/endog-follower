#include <string>
#include <iostream>
#include <fstream>
#include <exception>
#include <iterator>  
#include <pybind11/embed.h> // for py::scoped_interpreter

#include "input/input.hpp"
#include "instance/instance.hpp"
#include "leadersolver/leadersolverfactory.hpp"
#include "leadersolver/leadersolver.hpp"

int main(int argc, char *argv[]) {
	// Create
	py::scoped_interpreter guard{};  // start Python

	/* Data information. */
	try{
        // Read input.
		Input input(argc,argv);
		input.display();

        // Read instance.
        Instance instance(input);
		instance.display();
		input.setInstAlign(instance.getModel()->alignment);

        // Define and solve leader problem.
		LeaderSolverFactory factory;
		AbstractLeaderSolver * leader_solver = factory.createLeaderSolver(input,instance,"Leader");
        leader_solver->solve();

        // Write solution files.
		input.writeHead(leader_solver->writeHead());
		input.write(leader_solver->write());

        input.writeHeadComp(leader_solver->writeHeadComp());
		input.writeComp(leader_solver->writeComp());

		delete leader_solver;
	}catch (GRBException e){
    	std::cout << "Error code = " << e.getErrorCode() << std::endl;
    	std::cout << e.getMessage() << std::endl;
  	}
	catch (const std::string & e) { std::cout << "EXCEPTION | " << e << std::endl; }
	catch (const std::exception & e) { std::cout << "EXCEPTION | " << e.what() << std::endl; }

    return 0; 
}