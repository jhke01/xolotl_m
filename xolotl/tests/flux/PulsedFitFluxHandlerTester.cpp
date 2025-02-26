#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include "PulsedFitFluxHandler.h"
#include <mpi.h>
#include <HDF5NetworkLoader.h>
#include <DummyHandlerRegistry.h>
#include <XolotlConfig.h>
#include <Options.h>
#include <iostream>
#include <fstream>

using namespace std;
using namespace xolotlCore;

/**
 * The test suite is responsible for testing the PulsedFitFluxHandler.
 */
BOOST_AUTO_TEST_SUITE (PulsedFitFluxHandlerTester_testSuite)

BOOST_AUTO_TEST_CASE(checkComputeIncidentFlux) {
	// Create the option to create a network
	xolotlCore::Options opts;
	// Create a good parameter file
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=0 0 0 2 2" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	int argc = 2;
	char **argv = new char*[3];
	std::string appName = "fakeXolotlAppNameForTests";
	argv[0] = new char[appName.length() + 1];
	strcpy(argv[0], appName.c_str());
	std::string parameterFile = "param.txt";
	argv[1] = new char[parameterFile.length() + 1];
	strcpy(argv[1], parameterFile.c_str());
	argv[2] = 0; // null-terminate the array
	// Initialize MPI for HDF5
	MPI_Init(&argc, &argv);
	opts.readParams(argc, argv);

	// Create the network loader
	HDF5NetworkLoader loader = HDF5NetworkLoader(
			make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Create the network
	auto network = loader.generate(opts);
	// Get its size
	const int dof = network->getDOF();

	// Create a grid
	std::vector<double> grid;
	for (int l = 0; l < 100; l++) {
		grid.push_back((double) l * 50.0);
	}
	// Specify the surface position
	int surfacePos = 0;

	// Create the ulsed flux handler
	auto testFitFlux = make_shared<PulsedFitFluxHandler>();
	// Set the flux amplitude and pulse parameters
	testFitFlux->setFluxAmplitude(1.0);
	testFitFlux->setPulseTime(1.0e-3);
	testFitFlux->setProportion(0.2);
	// Initialize the flux handler
	testFitFlux->initializeFluxHandler(*network, surfacePos, grid);

	// Create a time
	double currTime = 1.0e-4;

	// The array of concentration
	double newConcentration[100 * dof];

	// Initialize their values
	for (int i = 0; i < 100 * dof; i++) {
		newConcentration[i] = 0.0;
	}

	// The pointer to the grid point we want
	double *updatedConc = &newConcentration[0];
	double *updatedConcOffset = updatedConc + 50 * dof;

	// Update the concentrations at some grid point
	testFitFlux->computeIncidentFlux(currTime, updatedConcOffset, 50,
			surfacePos);

	// Check the value at some grid point
	BOOST_REQUIRE_CLOSE(newConcentration[250], 5.0295072885924443e-08, 0.01);
	BOOST_REQUIRE_CLOSE(newConcentration[252], 5.0295072885924443e-08, 0.01);

	// Check that the flux is 0.0 at later time
	currTime = 8.0e-4;
	updatedConcOffset = updatedConc + 22 * dof;

	// Update the concentrations at some grid point
	testFitFlux->computeIncidentFlux(currTime, updatedConcOffset, 22,
			surfacePos);

	// Check the value at some grid point
	BOOST_REQUIRE_CLOSE(newConcentration[110], 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(newConcentration[112], 0.0, 0.01);

	// Finalize MPI
	MPI_Finalize();

	return;
}

BOOST_AUTO_TEST_SUITE_END()
