#ifndef UZRXECLUSTER_H
#define UZRXECLUSTER_H

// Includes
#include "UZrCluster.h"
#include <xolotlPerf.h>
#include <Constants.h>

namespace xolotlCore {

/**
 * This class represents a cluster composed entirely of xenon.
 */
class UZrXeCluster: public UZrCluster {

private:
	static std::string buildName(IReactant::SizeType nXe) {
		std::stringstream nameStream;
		nameStream << "Xe_" << nXe;
		return nameStream.str();
	}

public:

	/**
	 * Default constructor, deleted because we require info to construct.
	 */
	UZrXeCluster() = delete;

	/**
	 * The constructor. All UZrXeClusters must be initialized with a size.
	 *
	 * @param nXe the number of xenon atoms in the cluster
	 * @param registry The performance handler registry
	 */
	UZrXeCluster(int nXe, IReactionNetwork &_network,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
			UZrCluster(_network, registry, buildName(nXe)) {

		// Set the size
		size = nXe;
		// Update the composition map
		composition[toCompIdx(Species::Xe)] = size;

		// Set the typename appropriately
		type = ReactantType::Xe;

		// TODO: Compute the reaction radius
		// You can look at NEXeCluster.h

		double FourPi = 4.0 * xolotlCore::pi;
		//reactionRadius = network.getLatticeParameter()
			//	* pow((3.0 * nXe) / xolotlCore::pi, (1.0 / 3.0)) * 0.5;
		reactionRadius = pow(
				(3.0 * (double) size) / (FourPi * 10.162795276841),
				(1.0 / 3.0));

		// cout << "jvar = " << network.getDensity() << endl;
		// seems that network.getDensity() will always return zero under UZr


		if (size == 1)
			reactionRadius = network.getImpurityRadius();

		return;
	}

	/**
	 * Copy constructor, deleted to prevent use.
	 */
	UZrXeCluster(const UZrXeCluster &other) = delete;

	/**
	 * Destructor
	 */
	~UZrXeCluster() {
	}

	/**
	 * Add grid points to the vector of diffusion coefficients or remove
	 * them if the value is negative.
	 *
	 * @param i The number of grid point to add or remove
	 */
	void addGridPoints(int i) override {
		if (diffusionFactor > 0.0) {
			Reactant::addGridPoints(i);
		}

		// Don't do anything
		return;
	}

	/**
	 * This operation sets the temperature at which the reactant currently
	 * exists. Temperature-dependent quantities are recomputed when this
	 * operation is called, so the temperature should always be set first.
	 *
	 * @param temp The new cluster temperature
	 * @param i The location on the grid
	 */
	void setTemperature(double temp, int i) override {
		if (diffusionFactor > 0.0) {
			Reactant::setTemperature(temp, i);
		}

		// Don't do anything
		return;
	}

	/**
	 * This operation returns the diffusion coefficient for this reactant and is
	 * calculated from the diffusion factor.
	 *
	 * @param i The position on the grid
	 * @return The diffusion coefficient
	 */
	double getDiffusionCoefficient(int i) const override {
		if (diffusionFactor > 0.0) {
			return Reactant::getDiffusionCoefficient(i);
		}

		return 0.0;
	}

};
//end class UZrXeCluster

} /* namespace xolotlCore */
#endif
