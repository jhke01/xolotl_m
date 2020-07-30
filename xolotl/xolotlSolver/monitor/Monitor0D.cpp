// Includes
#include "PetscSolver.h"
#include <xolotlPerf.h>
#include <VizHandlerRegistryFactory.h>
#include <PlotType.h>
#include <CvsXDataProvider.h>
#include <CvsXYDataProvider.h>
#include <LabelProvider.h>
#include <Constants.h>
#include <petscts.h>
#include <petscsys.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <memory>
#include <NESuperCluster.h>
#include <PSISuperCluster.h>
#include <FeSuperCluster.h>
#include <UZrCluster.h>
#include <NEClusterReactionNetwork.h>
#include <PSIClusterReactionNetwork.h>
#include <FeClusterReactionNetwork.h>
#include <AlloyClusterReactionNetwork.h>
#include <UZrClusterReactionNetwork.h>
#include <AlloySuperCluster.h>
#include "xolotlCore/io/XFile.h"
#include "xolotlSolver/monitor/Monitor.h"

namespace xolotlSolver {

// Declaration of the functions defined in Monitor.cpp
extern PetscErrorCode checkTimeStep(TS ts);
extern PetscErrorCode monitorTime(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx);
extern PetscErrorCode computeFluence(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx);
extern PetscErrorCode monitorPerf(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx);

// Declaration of the variables defined in Monitor.cpp
extern std::shared_ptr<xolotlViz::IPlot> perfPlot;
extern double previousTime;
extern double timeStepThreshold;

//! The pointer to the plot used in monitorScatter0D.
std::shared_ptr<xolotlViz::IPlot> scatterPlot0D;
//! How often HDF5 file is written
PetscReal hdf5Stride0D = 0.0;
//! Previous time for HDF5
PetscInt hdf5Previous0D = 0;
//! HDF5 output file name
std::string hdf5OutputName0D = "xolotlStop.h5";
// Declare the vector that will store the Id of the helium clusters
std::vector<int> indices0D;
// Declare the vector that will store the weight of the helium clusters
// (their He composition)
std::vector<int> weights0D;
// Declare the vector that will store the radii of bubbles
std::vector<double> radii0D;

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "startStop0D")
/**
 * This is a monitoring method that update an hdf5 file at each time step.
 */
PetscErrorCode startStop0D(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void*) {
	// Initial declaration
	PetscErrorCode ierr;
	const double **solutionArray, *gridPointSolution;

	PetscFunctionBeginUser;

	// Compute the dt
	double dt = time - previousTime;

	// Don't do anything if it is not on the stride
	if (((int) ((time + dt / 10.0) / hdf5Stride0D) <= hdf5Previous0D)
			&& timestep > 0)
		PetscFunctionReturn(0);

	// Update the previous time
	if ((int) ((time + dt / 10.0) / hdf5Stride0D) > hdf5Previous0D)
		hdf5Previous0D++;

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Get the solver handler
	auto &solverHandler = PetscSolver::getSolverHandler();

	// Get the network and dof
	auto &network = solverHandler.getNetwork();
	const int dof = network.getDOF();

	// Create an array for the concentration
	double concArray[dof][2];

	// Open the existing HDF5 file
	xolotlCore::XFile checkpointFile(hdf5OutputName0D, PETSC_COMM_WORLD,
			xolotlCore::XFile::AccessMode::OpenReadWrite);

	// Get the current time step
	double currentTimeStep;
	ierr = TSGetTimeStep(ts, &currentTimeStep);
	CHKERRQ(ierr);

	// Add a concentration time step group for the current time step.
	auto concGroup = checkpointFile.getGroup<
			xolotlCore::XFile::ConcentrationGroup>();
	assert(concGroup);
	auto tsGroup = concGroup->addTimestepGroup(timestep, time, previousTime,
			currentTimeStep);

	// Determine the concentration values we will write.
	// We only examine and collect the grid points we own.
	// TODO measure impact of us building the flattened representation
	// rather than a ragged 2D representation.
	XFile::TimestepGroup::Concs1DType concs(1);

	// Access the solution data for the current grid point.
	gridPointSolution = solutionArray[0];

	for (auto l = 0; l < dof; ++l) {
		if (std::fabs(gridPointSolution[l]) > 1.0e-16) {
			concs[0].emplace_back(l, gridPointSolution[l]);
		}
	}

	// Write our concentration data to the current timestep group
	// in the HDF5 file.
	// We only write the data for the grid points we own.
	tsGroup->writeConcentrations(checkpointFile, 0, concs);

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "computeXenonContent0D")
/**
 * This is a monitoring method that will compute the xenon content in UZr
 * TODO: Let's discuss what type of output you will need
 */
PetscErrorCode computeXenonContent0D(TS ts, PetscInt, PetscReal time,
		Vec solution, void*) {
	// Initial declarations
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	// Get the solver handler
	auto &solverHandler = PetscSolver::getSolverHandler();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the network
	auto &network = solverHandler.getNetwork();

	// Get the array of concentration
	PetscReal **solutionArray;
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);

	// Declare the pointer for the concentrations at a specific grid point
	PetscReal *gridPointSolution;

	// Get the pointer to the beginning of the solution data for this grid point
	gridPointSolution = solutionArray[0];

	// Update the concentration in the network
	network.updateConcentrationsFromArray(gridPointSolution);
	CHKERRQ(ierr);

	// Get the minimum size for the radius
	auto minSizes = solverHandler.getMinSizes();

	// Store the concentration and other values over the grid
	double bubbleConcentration = 0.0, radii = 0.0,
			partialBubbleConcentration = 0.0, partialRadii = 0.0, partialSize =
					0.0, Xe_Density = 0.0, Xe_Size_t = 0.0, Xe_Par_Density = 0.0, Xe_Par_Size_t = 0.0;

	// Loop on all the indices
	for (unsigned int i = 0; i < indices0D.size(); i++) {
		// Add the current concentration times the number of xenon in the cluster
		// (from the weight vector)
		double conc = gridPointSolution[indices0D[i]];
		bubbleConcentration += conc;
		radii += conc * radii0D[i];
	}

	Xe_Par_Density = 0.0;
	Xe_Par_Size_t = 0.0;
	Xe_Density = 0.0;
	Xe_Size_t = 0.0;

	for (auto const &currMapItem : network.getAll(ReactantType::Xe)) {
		// Get the cluster
    auto const &cluster = *(currMapItem.second);
		double Xe_conc = cluster.getConcentration();
		Xe_Density += Xe_conc;
		Xe_Size_t += Xe_conc * cluster.getSize();
		if (cluster.getSize() >= minSizes[0] && Xe_conc > 1.0e-16) {
		  Xe_Par_Density += Xe_conc;
			Xe_Par_Size_t += Xe_conc * cluster.getSize();
			//Xe_Par_Radius += gridPointSolution[cluster.getId() - 1] * cluster.getReactionRadius();
		}
	}

	double Xe_Size = Xe_Size_t / Xe_Density;
	double Xe_Par_Size = (Xe_Par_Size_t) / (Xe_Par_Density);

	double xeConcentration = network.getTotalAtomConcentration();
	double vConcentration = network.getTotalVConcentration();

	// Print the result
	std::cout << "\nTime: " << time << std::endl;
	std::cout << "Xenon concentration = " << xeConcentration << std::endl;
	std::cout << "Vacancy concentration = " << vConcentration << std::endl
			<< std::endl;


	// Write more data in a file
	std::ofstream outputFile;
	outputFile.open("contentOut.txt", ios::app);
	outputFile << time << " " << xeConcentration << " " << vConcentration << " "
			<< xeConcentration / bubbleConcentration << " "
			<< vConcentration / bubbleConcentration << " "
			//<< bubbleConcentration << " "
			<< Xe_Density << " "
			<< Xe_Size << " "
			<< Xe_Par_Density << " "
			<< Xe_Par_Size << " "
			//<< radii / bubbleConcentration
			<< std::endl;
			outputFile.close();

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "computeXenonRetention0D")
/**
 * This is a monitoring method that will compute the xenon retention
 */
PetscErrorCode computeXenonRetention0D(TS ts, PetscInt, PetscReal time,
		Vec solution, void*) {
	// Initial declarations
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	// Get the solver handler
	auto &solverHandler = PetscSolver::getSolverHandler();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the network
	auto &network = solverHandler.getNetwork();

	// Get the array of concentration
	PetscReal **solutionArray;
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Store the concentration and other values over the grid
	double xeConcentration = 0.0, bubbleConcentration = 0.0, radii = 0.0,
			partialBubbleConcentration = 0.0, partialRadii = 0.0, partialSize =
					0.0;

	// Declare the pointer for the concentrations at a specific grid point
	PetscReal *gridPointSolution;

	// Get the pointer to the beginning of the solution data for this grid point
	gridPointSolution = solutionArray[0];

	// Update the concentration in the network
	network.updateConcentrationsFromArray(gridPointSolution);

	// Get the minimum size for the radius
	auto minSizes = solverHandler.getMinSizes();

	// Loop on all the indices
	for (unsigned int i = 0; i < indices0D.size(); i++) {
		// Add the current concentration times the number of xenon in the cluster
		// (from the weight vector)
		double conc = gridPointSolution[indices0D[i]];
		xeConcentration += conc * weights0D[i];
		bubbleConcentration += conc;
		radii += conc * radii0D[i];
		if (weights0D[i] >= minSizes[0] && conc > 1.0e-16) {
			partialBubbleConcentration += conc;
			partialRadii += conc * radii0D[i];
			partialSize += conc * weights0D[i];
		}
	}

	// Loop on all the super clusters
	for (auto const &superMapItem : network.getAll(ReactantType::NESuper)) {
		auto const &cluster =
				static_cast<NESuperCluster&>(*(superMapItem.second));
		double conc = cluster.getTotalConcentration();
		xeConcentration += cluster.getTotalXenonConcentration();
		bubbleConcentration += conc;
		radii += conc * cluster.getReactionRadius();
		if (cluster.getSize() >= minSizes[0] && conc > 1.0e-16) {
			partialBubbleConcentration += conc;
			partialRadii += conc * cluster.getReactionRadius();
			partialSize += cluster.getTotalXenonConcentration();
		}
	}

	// Print the result
	std::cout << "\nTime: " << time << std::endl;
	std::cout << "Xenon concentration = " << xeConcentration << std::endl
			<< std::endl;

	// Make sure the average partial radius makes sense
	double averagePartialRadius = partialRadii / partialBubbleConcentration;
	double minRadius = pow(
			(3.0 * (double) minSizes[0])
					/ (4.0 * xolotlCore::pi * network.getDensity()),
			(1.0 / 3.0));

	if (partialBubbleConcentration < 1.e-16 || averagePartialRadius < minRadius)
		averagePartialRadius = minRadius;

	// Uncomment to write the retention and the fluence in a file
	std::ofstream outputFile;
	outputFile.open("retentionOut.txt", ios::app);
	outputFile << time << " " << xeConcentration << " "
			<< radii / bubbleConcentration << " " << averagePartialRadius << " "
			<< partialBubbleConcentration << " "
			<< partialSize / partialBubbleConcentration << std::endl;
	outputFile.close();

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "computeAlloy0D")
/**
 * This is a monitoring method that will compute average density and diameter
 */
PetscErrorCode computeAlloy0D(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {

	// Initial declarations
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	// Get the solver handler
	auto &solverHandler = PetscSolver::getSolverHandler();

	// Get the physical grid and its length
	auto grid = solverHandler.getXGrid();
	int xSize = grid.size();

	// Get the position of the surface
	int surfacePos = solverHandler.getSurfacePosition();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the array of concentration
	PetscReal **solutionArray;
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Get the network
	auto &network = solverHandler.getNetwork();

	// Get degrees of freedom
	auto dof = network.getDOF();

	// Initial declarations for the density and diameter
	double iDensity = 0.0, vDensity = 0.0, voidDensity = 0.0,
			frankDensity = 0.0, faultedDensity = 0.0, perfectDensity = 0.0,
			voidPartialDensity = 0.0, frankPartialDensity = 0.0,
			faultedPartialDensity = 0.0, perfectPartialDensity = 0.0,
			iDiameter = 0.0, vDiameter = 0.0, voidDiameter = 0.0,
			frankDiameter = 0.0, faultedDiameter = 0.0, perfectDiameter = 0.0,
			voidPartialDiameter = 0.0, frankPartialDiameter = 0.0,
			faultedPartialDiameter = 0.0, perfectPartialDiameter = 0.0;

	// Get the minimum size for the loop densities and diameters
	auto minSizes = solverHandler.getMinSizes();

	// Declare the pointer for the concentrations at a specific grid point
	PetscReal *gridPointSolution;

	// Get the pointer to the beginning of the solution data for this grid point
	gridPointSolution = solutionArray[0];

	// Update the concentration in the network
	network.updateConcentrationsFromArray(gridPointSolution);

	// Loop on I
	for (auto const &iMapItem : network.getAll(ReactantType::I)) {
		// Get the cluster
		auto const &cluster = *(iMapItem.second);
		iDensity += gridPointSolution[cluster.getId() - 1];
		iDiameter += gridPointSolution[cluster.getId() - 1]
				* cluster.getReactionRadius() * 2.0;
	}

	// Loop on V
	for (auto const &vMapItem : network.getAll(ReactantType::V)) {
		// Get the cluster
		auto const &cluster = *(vMapItem.second);
		vDensity += gridPointSolution[cluster.getId() - 1];
		vDiameter += gridPointSolution[cluster.getId() - 1]
				* cluster.getReactionRadius() * 2.0;
	}

	// Loop on Void
	for (auto const &voidMapItem : network.getAll(ReactantType::Void)) {
		// Get the cluster
		auto const &cluster = *(voidMapItem.second);
		voidDensity += gridPointSolution[cluster.getId() - 1];
		voidDiameter += gridPointSolution[cluster.getId() - 1]
				* cluster.getReactionRadius() * 2.0;
		if (cluster.getSize() >= minSizes[0]) {
			voidPartialDensity += gridPointSolution[cluster.getId() - 1];
			voidPartialDiameter += gridPointSolution[cluster.getId() - 1]
					* cluster.getReactionRadius() * 2.0;
		}
	}
	for (auto const &voidMapItem : network.getAll(ReactantType::VoidSuper)) {
		// Get the cluster
		auto const &cluster =
				static_cast<AlloySuperCluster&>(*(voidMapItem.second));
		voidDensity += cluster.getTotalConcentration();
		voidDiameter += cluster.getTotalConcentration()
				* cluster.getReactionRadius() * 2.0;
		if (cluster.getSize() >= minSizes[0]) {
			voidPartialDensity += cluster.getTotalConcentration();
			voidPartialDiameter += cluster.getTotalConcentration()
					* cluster.getReactionRadius() * 2.0;
		}
	}

	// Loop on Faulted
	for (auto const &faultedMapItem : network.getAll(ReactantType::Faulted)) {
		// Get the cluster
		auto const &cluster = *(faultedMapItem.second);
		faultedDensity += gridPointSolution[cluster.getId() - 1];
		faultedDiameter += gridPointSolution[cluster.getId() - 1]
				* cluster.getReactionRadius() * 2.0;
		if (cluster.getSize() >= minSizes[1]) {
			faultedPartialDensity += gridPointSolution[cluster.getId() - 1];
			faultedPartialDiameter += gridPointSolution[cluster.getId() - 1]
					* cluster.getReactionRadius() * 2.0;
		}
	}
	for (auto const &faultedMapItem : network.getAll(ReactantType::FaultedSuper)) {
		// Get the cluster
		auto const &cluster =
				static_cast<AlloySuperCluster&>(*(faultedMapItem.second));
		faultedDensity += cluster.getTotalConcentration();
		faultedDiameter += cluster.getTotalConcentration()
				* cluster.getReactionRadius() * 2.0;
		if (cluster.getSize() >= minSizes[1]) {
			faultedPartialDensity += cluster.getTotalConcentration();
			faultedPartialDiameter += cluster.getTotalConcentration()
					* cluster.getReactionRadius() * 2.0;
		}
	}

	// Loop on Perfect
	for (auto const &perfectMapItem : network.getAll(ReactantType::Perfect)) {
		// Get the cluster
		auto const &cluster = *(perfectMapItem.second);
		perfectDensity += gridPointSolution[cluster.getId() - 1];
		perfectDiameter += gridPointSolution[cluster.getId() - 1]
				* cluster.getReactionRadius() * 2.0;
		if (cluster.getSize() >= minSizes[2]) {
			perfectPartialDensity += gridPointSolution[cluster.getId() - 1];
			perfectPartialDiameter += gridPointSolution[cluster.getId() - 1]
					* cluster.getReactionRadius() * 2.0;
		}
	}
	for (auto const &perfectMapItem : network.getAll(ReactantType::PerfectSuper)) {
		// Get the cluster
		auto const &cluster =
				static_cast<AlloySuperCluster&>(*(perfectMapItem.second));
		perfectDensity += cluster.getTotalConcentration();
		perfectDiameter += cluster.getTotalConcentration()
				* cluster.getReactionRadius() * 2.0;
		if (cluster.getSize() >= minSizes[2]) {
			perfectPartialDensity += cluster.getTotalConcentration();
			perfectPartialDiameter += cluster.getTotalConcentration()
					* cluster.getReactionRadius() * 2.0;
		}
	}

	// Loop on Frank
	for (auto const &frankMapItem : network.getAll(ReactantType::Frank)) {
		// Get the cluster
		auto const &cluster = *(frankMapItem.second);
		frankDensity += gridPointSolution[cluster.getId() - 1];
		frankDiameter += gridPointSolution[cluster.getId() - 1]
				* cluster.getReactionRadius() * 2.0;
		if (cluster.getSize() >= minSizes[3]) {
			frankPartialDensity += gridPointSolution[cluster.getId() - 1];
			frankPartialDiameter += gridPointSolution[cluster.getId() - 1]
					* cluster.getReactionRadius() * 2.0;
		}
	}
	for (auto const &frankMapItem : network.getAll(ReactantType::FrankSuper)) {
		// Get the cluster
		auto const &cluster =
				static_cast<AlloySuperCluster&>(*(frankMapItem.second));
		frankDensity += cluster.getTotalConcentration();
		frankDiameter += cluster.getTotalConcentration()
				* cluster.getReactionRadius() * 2.0;
		if (cluster.getSize() >= minSizes[3]) {
			frankPartialDensity += cluster.getTotalConcentration();
			frankPartialDiameter += cluster.getTotalConcentration()
					* cluster.getReactionRadius() * 2.0;
		}
	}

	// Set the output precision
	const int outputPrecision = 5;

	// Average the diameters
	iDiameter = iDiameter / iDensity;
	vDiameter = vDiameter / vDensity;
	voidDiameter = voidDiameter / voidDensity;
	perfectDiameter = perfectDiameter / perfectDensity;
	faultedDiameter = faultedDiameter / faultedDensity;
	frankDiameter = frankDiameter / frankDensity;
	voidPartialDiameter = voidPartialDiameter / voidPartialDensity;
	perfectPartialDiameter = perfectPartialDiameter / perfectPartialDensity;
	faultedPartialDiameter = faultedPartialDiameter / faultedPartialDensity;
	frankPartialDiameter = frankPartialDiameter / frankPartialDensity;

	// Open the output file
	std::fstream outputFile;
	outputFile.open("Alloy.dat", std::fstream::out | std::fstream::app);
	outputFile << std::setprecision(outputPrecision);

	// Output the data
	outputFile << timestep << " " << time << " " << iDensity << " " << iDiameter
			<< " " << vDensity << " " << vDiameter << " " << voidDensity << " "
			<< voidDiameter << " " << faultedDensity << " " << faultedDiameter
			<< " " << perfectDensity << " " << perfectDiameter << " "
			<< frankDensity << " " << frankDiameter << " " << voidPartialDensity
			<< " " << voidPartialDiameter << " " << faultedPartialDensity << " "
			<< faultedPartialDiameter << " " << perfectPartialDensity << " "
			<< perfectPartialDiameter << " " << frankPartialDensity << " "
			<< frankPartialDiameter << std::endl;

	// Close the output file
	outputFile.close();

	// Restore the PETSC solution array
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);

}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "monitorScatter0D")
/**
 * This is a monitoring method that will save 1D plots of the xenon concentration
 * distribution.
 */
PetscErrorCode monitorScatter0D(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void*) {
	// Initial declarations
	PetscErrorCode ierr;
	double **solutionArray, *gridPointSolution;

	PetscFunctionBeginUser;

	// Don't do anything if it is not on the stride
	if (timestep % 10 != 0)
		PetscFunctionReturn(0);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Get the solver handler
	auto &solverHandler = PetscSolver::getSolverHandler();

	// Get the network and its size
	auto &network = solverHandler.getNetwork();
	int networkSize = network.size();
	auto &superClusters = network.getAll(ReactantType::NESuper);

	// Create a Point vector to store the data to give to the data provider
	// for the visualization
	auto myPoints = std::make_shared<std::vector<xolotlViz::Point> >();

	// Get the pointer to the beginning of the solution data for this grid point
	gridPointSolution = solutionArray[0];

	// Update the concentration in the network
	network.updateConcentrationsFromArray(gridPointSolution);

	for (int i = 0; i < networkSize - superClusters.size(); i++) {
		// Create a Point with the concentration[i] as the value
		// and add it to myPoints
		xolotlViz::Point aPoint;
		aPoint.value = gridPointSolution[i];
		aPoint.t = time;
		aPoint.x = (double) i + 1.0;
		myPoints->push_back(aPoint);
	}
	int nXe = networkSize - superClusters.size() + 1;
	// Loop on the super clusters
	for (auto const &superMapItem : superClusters) {
		// Get the cluster
		auto const &cluster =
				static_cast<NESuperCluster&>(*(superMapItem.second));
		// Get the width
		int width = cluster.getSectionWidth();
		// Loop on the width
		for (int k = 0; k < width; k++) {
			// Compute the distance
			double dist = cluster.getDistance(nXe + k);
			// Create a Point with the concentration[i] as the value
			// and add it to myPoints
			xolotlViz::Point aPoint;
			aPoint.value = cluster.getConcentration(dist);
			aPoint.t = time;
			aPoint.x = (double) nXe + k;
			myPoints->push_back(aPoint);
		}

		// update nXe
		nXe += width;
	}

	// Get the data provider and give it the points
	scatterPlot0D->getDataProvider()->setPoints(myPoints);

	// Change the title of the plot and the name of the data
	std::stringstream title;
	title << "Size Distribution";
	scatterPlot0D->getDataProvider()->setDataName(title.str());
	scatterPlot0D->plotLabelProvider->titleLabel = title.str();
	// Give the time to the label provider
	std::stringstream timeLabel;
	timeLabel << "time: " << std::setprecision(4) << time << "s";
	scatterPlot0D->plotLabelProvider->timeLabel = timeLabel.str();
	// Get the current time step
	PetscReal currentTimeStep;
	ierr = TSGetTimeStep(ts, &currentTimeStep);
	CHKERRQ(ierr);
	// Give the timestep to the label provider
	std::stringstream timeStepLabel;
	timeStepLabel << "dt: " << std::setprecision(4) << currentTimeStep << "s";
	scatterPlot0D->plotLabelProvider->timeStepLabel = timeStepLabel.str();

	// Render and save in file
	std::stringstream fileName;
	fileName << "Scatter_TS" << timestep << ".png";
	scatterPlot0D->write(fileName.str());

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "monitorBubble0D")
/**
 * This is a monitoring method that will create files with the mean
 * concentration of each bubble at each time step.
 */
PetscErrorCode monitorBubble0D(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	// Initial declaration
	PetscErrorCode ierr;
	double **solutionArray, *gridPointSolution;

	PetscFunctionBeginUser;

//	// Don't do anything if it is not on the stride
//	if (timestep % 10 != 0)
//		PetscFunctionReturn(0);

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the solutionArray
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Get the solver handler
	auto &solverHandler = PetscSolver::getSolverHandler();

	// Get the network
	auto &network = solverHandler.getNetwork();
	int dof = network.getDOF();

	// Create the output file
	std::ofstream outputFile;
	std::stringstream name;
	name << "bubble_" << timestep << ".dat";
	outputFile.open(name.str());

	// Get the pointer to the beginning of the solution data for this grid point
	gridPointSolution = solutionArray[0];

	// Update the concentration in the network
	network.updateConcentrationsFromArray(gridPointSolution);

	// Initialize the total helium and concentration before looping
	double concTot = 0.0, heliumTot = 0.0;

	// Consider each super cluster.
	for (auto const &superMapItem : network.getAll(ReactantType::FeSuper)) {
		// Get the super cluster
		auto const &superCluster =
				static_cast<FeSuperCluster&>(*(superMapItem.second));
		// Get its boundaries
		auto const &heBounds = superCluster.getHeBounds();
		auto const &vBounds = superCluster.getVBounds();
		// Get its diameter
		double diam = 2.0 * superCluster.getReactionRadius();
		// Get its concentration
		double conc = superCluster.getConcentration(0.0, 0.0);

		// For compatibility with previous versions, we output
		// the value of a closed upper bound of the He and V intervals.
		outputFile << *(heBounds.begin()) << " " << *(heBounds.end()) - 1 << " "
				<< *(vBounds.begin()) << " " << *(vBounds.end()) - 1 << " "
				<< conc << std::endl;
	}

	// Close the file
	outputFile.close();

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

/**
 * This operation sets up different monitors
 *  depending on the options.
 * @param ts The time stepper
 * @return A standard PETSc error code
 */
PetscErrorCode setupPetsc0DMonitor(TS ts) {
	PetscErrorCode ierr;

	// Get xolotlViz handler registry
	auto vizHandlerRegistry = xolotlFactory::getVizHandlerRegistry();

	// Flags to launch the monitors or not
	PetscBool flagCheck, flag1DPlot, flagBubble, flagPerf, flagStatus,
			flagAlloy, flagXeRetention, flagXeContent;

	// Check the option -check_collapse
	ierr = PetscOptionsHasName(NULL, NULL, "-check_collapse", &flagCheck);
	checkPetscError(ierr,
			"setupPetsc0DMonitor: PetscOptionsHasName (-check_collapse) failed.");

	// Check the option -plot_perf
	ierr = PetscOptionsHasName(NULL, NULL, "-plot_perf", &flagPerf);
	checkPetscError(ierr,
			"setupPetsc0DMonitor: PetscOptionsHasName (-plot_perf) failed.");

	// Check the option -plot_1d
	ierr = PetscOptionsHasName(NULL, NULL, "-plot_1d", &flag1DPlot);
	checkPetscError(ierr,
			"setupPetsc0DMonitor: PetscOptionsHasName (-plot_1d) failed.");

	// Check the option -start_stop
	ierr = PetscOptionsHasName(NULL, NULL, "-start_stop", &flagStatus);
	checkPetscError(ierr,
			"setupPetsc0DMonitor: PetscOptionsHasName (-start_stop) failed.");

	// Check the option -bubble
	ierr = PetscOptionsHasName(NULL, NULL, "-bubble", &flagBubble);
	checkPetscError(ierr,
			"setupPetsc0DMonitor: PetscOptionsHasName (-bubble) failed.");

	// Check the option -alloy
	ierr = PetscOptionsHasName(NULL, NULL, "-alloy", &flagAlloy);
	checkPetscError(ierr,
			"setupPetsc0DMonitor: PetscOptionsHasName (-alloy) failed.");

	// Check the option -xenon_retention
	ierr = PetscOptionsHasName(NULL, NULL, "-xenon_retention",
			&flagXeRetention);
	checkPetscError(ierr,
			"setupPetsc0DMonitor: PetscOptionsHasName (-xenon_retention) failed.");

	// Check the option -xenon_content
	ierr = PetscOptionsHasName(NULL, NULL, "-xenon_content", &flagXeContent);
	checkPetscError(ierr,
			"setupPetsc0DMonitor: PetscOptionsHasName (-xenon_content) failed.");

	// Get the solver handler
	auto &solverHandler = PetscSolver::getSolverHandler();

	// Get the network and its size
	auto &network = solverHandler.getNetwork();
	const int networkSize = network.size();

	// Determine if we have an existing restart file,
	// and if so, it it has had timesteps written to it.
	std::unique_ptr<xolotlCore::XFile> networkFile;
	std::unique_ptr<xolotlCore::XFile::TimestepGroup> lastTsGroup;
	std::string networkName = solverHandler.getNetworkName();
	bool hasConcentrations = false;
	if (not networkName.empty()) {
		networkFile.reset(new xolotlCore::XFile(networkName));
		auto concGroup = networkFile->getGroup<
				xolotlCore::XFile::ConcentrationGroup>();
		hasConcentrations = (concGroup and concGroup->hasTimesteps());
		if (hasConcentrations) {
			lastTsGroup = concGroup->getLastTimestepGroup();
		}
	}

	// Set the post step processing to stop the solver if the time step collapses
	if (flagCheck) {
		// Find the threshold
		PetscBool flag;
		ierr = PetscOptionsGetReal(NULL, NULL, "-check_collapse",
				&timeStepThreshold, &flag);
		checkPetscError(ierr,
				"setupPetsc0DMonitor: PetscOptionsGetInt (-check_collapse) failed.");
		if (!flag)
			timeStepThreshold = 1.0e-16;

		// Set the post step process that tells the solver when to stop if the time step collapse
		ierr = TSSetPostStep(ts, checkTimeStep);
		checkPetscError(ierr,
				"setupPetsc0DMonitor: TSSetPostStep (checkTimeStep) failed.");
	}

	// Set the monitor to save the status of the simulation in hdf5 file
	if (flagStatus) {
		// Find the stride to know how often the HDF5 file has to be written
		PetscBool flag;
		ierr = PetscOptionsGetReal(NULL, NULL, "-start_stop", &hdf5Stride0D,
				&flag);
		checkPetscError(ierr,
				"setupPetsc0DMonitor: PetscOptionsGetInt (-start_stop) failed.");
		if (!flag)
			hdf5Stride0D = 1.0;

		// Compute the correct hdf5Previous0D for a restart
		// Get the last time step written in the HDF5 file
		if (hasConcentrations) {

			assert(lastTsGroup);

			// Get the previous time from the HDF5 file
			previousTime = lastTsGroup->readPreviousTime();
			hdf5Previous0D = (int) (previousTime / hdf5Stride0D);
		}

		// Don't do anything if both files have the same name
		if (hdf5OutputName0D != solverHandler.getNetworkName()) {

			PetscInt Mx;
			PetscErrorCode ierr;

			// Get the da from ts
			DM da;
			ierr = TSGetDM(ts, &da);
			checkPetscError(ierr, "setupPetsc0DMonitor: TSGetDM failed.");

			// Get the size of the total grid
			ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE);
			checkPetscError(ierr, "setupPetsc0DMonitor: DMDAGetInfo failed.");

			// Get the solver handler
			auto &solverHandler = PetscSolver::getSolverHandler();

			// Get the physical grid (which is empty)
			auto grid = solverHandler.getXGrid();

			// Get the compostion list and save it
			auto compList = network.getCompositionList();

			// Create and initialize a checkpoint file.
			// We do this in its own scope so that the file
			// is closed when the file object goes out of scope.
			// We want it to close before we (potentially) copy
			// the network from another file using a single-process
			// MPI communicator.
			{
				xolotlCore::XFile checkpointFile(hdf5OutputName0D, grid,
						compList, PETSC_COMM_WORLD);
			}

			// Copy the network group from the given file (if it has one).
			// We open the files using a single-process MPI communicator
			// because it is faster for a single process to do the
			// copy with HDF5's H5Ocopy implementation than it is
			// when all processes call the copy function.
			// The checkpoint file must be closed before doing this.
			writeNetwork(PETSC_COMM_WORLD, solverHandler.getNetworkName(),
					hdf5OutputName0D, network);
		}

		// startStop0D will be called at each timestep
		ierr = TSMonitorSet(ts, startStop0D, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc0DMonitor: TSMonitorSet (startStop0D) failed.");
	}

	// Set the monitor to save 1D plot of xenon distribution
	if (flag1DPlot) {
		// Create a ScatterPlot
		scatterPlot0D = vizHandlerRegistry->getPlot("scatterPlot0D",
				xolotlViz::PlotType::SCATTER);

//		scatterPlot0D->setLogScale();

		// Create and set the label provider
		auto labelProvider = std::make_shared<xolotlViz::LabelProvider>(
				"labelProvider");
		labelProvider->axis1Label = "Xenon Size";
		labelProvider->axis2Label = "Concentration";

		// Give it to the plot
		scatterPlot0D->setLabelProvider(labelProvider);

		// Create the data provider
		auto dataProvider = std::make_shared<xolotlViz::CvsXDataProvider>(
				"dataProvider");

		// Give it to the plot
		scatterPlot0D->setDataProvider(dataProvider);

		// monitorScatter0D will be called at each timestep
		ierr = TSMonitorSet(ts, monitorScatter0D, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc0DMonitor: TSMonitorSet (monitorScatter0D) failed.");
	}

	// Set the monitor to save performance plots (has to be in parallel)
	if (flagPerf) {
		// Create a ScatterPlot
		perfPlot = vizHandlerRegistry->getPlot("perfPlot",
				xolotlViz::PlotType::SCATTER);

		// Create and set the label provider
		auto labelProvider = std::make_shared<xolotlViz::LabelProvider>(
				"labelProvider");
		labelProvider->axis1Label = "Process ID";
		labelProvider->axis2Label = "Solver Time";

		// Give it to the plot
		perfPlot->setLabelProvider(labelProvider);

		// Create the data provider
		auto dataProvider = std::make_shared<xolotlViz::CvsXDataProvider>(
				"dataProvider");

		// Give it to the plot
		perfPlot->setDataProvider(dataProvider);

		// monitorPerf will be called at each timestep
		ierr = TSMonitorSet(ts, monitorPerf, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc0DMonitor: TSMonitorSet (monitorPerf) failed.");
	}

	// Set the monitor to save text file of the mean concentration of bubbles
	if (flagBubble) {
		// monitorBubble0D will be called at each timestep
		ierr = TSMonitorSet(ts, monitorBubble0D, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc0DMonitor: TSMonitorSet (monitorBubble0D) failed.");
	}

	// Set the monitor to output data for Alloy
	if (flagAlloy) {
		// Create/open the output files
		std::fstream outputFile;
		outputFile.open("Alloy.dat", std::fstream::out);
		outputFile.close();

		// computeAlloy0D will be called at each timestep
		ierr = TSMonitorSet(ts, computeAlloy0D, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc0DMonitor: TSMonitorSet (computeAlloy0D) failed.");
	}

	// Set the monitor to compute the xenon fluence and the retention
	// for the retention calculation
	if (flagXeRetention) {
		// Loop on the xenon clusters
		for (auto const &xeMapItem : network.getAll(ReactantType::Xe)) {
			auto const &cluster = *(xeMapItem.second);

			int id = cluster.getId() - 1;
			// Add the Id to the vector
			indices0D.push_back(id);
			// Add the number of xenon of this cluster to the weight
			weights0D.push_back(cluster.getSize());
			radii0D.push_back(cluster.getReactionRadius());
		}

		// Get the previous time if concentrations were stored and initialize the fluence
		if (hasConcentrations) {

			assert(lastTsGroup);

			// Get the previous time from the HDF5 file
			double time = lastTsGroup->readPreviousTime();
			// Initialize the fluence
			auto fluxHandler = solverHandler.getFluxHandler();
			// The length of the time step
			double dt = time;
			// Increment the fluence with the value at this current timestep
			fluxHandler->incrementFluence(dt);
			// Get the previous time from the HDF5 file
			previousTime = time;
		}

		// computeFluence will be called at each timestep
		ierr = TSMonitorSet(ts, computeFluence, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc0DMonitor: TSMonitorSet (computeFluence) failed.");

		// computeXenonRetention0D will be called at each timestep
		ierr = TSMonitorSet(ts, computeXenonRetention0D, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc0DMonitor: TSMonitorSet (computeXenonRetention0D) failed.");

		// Uncomment to clear the file where the retention will be written
		std::ofstream outputFile;
		outputFile.open("retentionOut.txt");
		outputFile.close();
	}

	// Set the monitor to compute the content in UZr material
	if (flagXeContent) {
		// Loop on all clusters
		for (IReactant &currReactant : network.getAll()) {
			auto& cluster = static_cast<UZrCluster&>(currReactant);
			int id = cluster.getId() - 1;
			// Add the Id to the vector
			indices0D.push_back(id);
			// Add the radius
			radii0D.push_back(cluster.getReactionRadius());
		}

		// Get the previous time if concentrations were stored and initialize the fluence
		if (hasConcentrations) {

			assert(lastTsGroup);

			// Get the previous time from the HDF5 file
			double time = lastTsGroup->readPreviousTime();
			// Initialize the fluence
			auto fluxHandler = solverHandler.getFluxHandler();
			// The length of the time step
			double dt = time;
			// Increment the fluence with the value at this current timestep
			fluxHandler->incrementFluence(dt);
			// Get the previous time from the HDF5 file
			previousTime = time;
		}

		// computeXenonContent0D will be called at each timestep
		ierr = TSMonitorSet(ts, computeXenonContent0D, NULL, NULL);
		checkPetscError(ierr,
				"setupPetsc0DMonitor: TSMonitorSet (computeXenonContent0D) failed.");

		// Reset the output file
		std::ofstream outputFile;
		outputFile.open("contentOut.txt");
		outputFile.close();
	}

	// Set the monitor to simply change the previous time to the new time
	// monitorTime will be called at each timestep
	ierr = TSMonitorSet(ts, monitorTime, NULL, NULL);
	checkPetscError(ierr,
			"setupPetsc0DMonitor: TSMonitorSet (monitorTime) failed.");

	PetscFunctionReturn(0);
}

}

/* end namespace xolotlSolver */
