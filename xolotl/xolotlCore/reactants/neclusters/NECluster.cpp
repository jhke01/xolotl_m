#include <algorithm>
#include "NECluster.h"
#include <xolotlPerf.h>
#include <Constants.h>
#include <MathUtils.h>

using namespace xolotlCore;

void NECluster::recomputeDiffusionCoefficient(double temp, int i) {
	// Return zero if the diffusion factor is zero.
	if (xolotlCore::equal(diffusionFactor, 0.0))
		return;

	// If the diffusivity is given
	if (migrationEnergy < 0.0) {
		diffusionCoefficient[i] = diffusionFactor;
		return;
	}
  /**
	// Intrinsic diffusion
	double kernel = -3.04 / (xolotlCore::kBoltzmann * temp);
	double D3 = 7.6e8 * exp(kernel); // nm2/s

	// We need the fission rate now
	double fissionRate = network.getFissionRate() * 1.0e27; // #/m3/s

	// Athermal diffusion
	double D1 = (8e-40 * fissionRate) * 1.0e18; // nm2/s

	// Radiation-enhanced diffusion
	kernel = -1.2 / (xolotlCore::kBoltzmann * temp);
	double D2 = (5.6e-25 * sqrt(fissionRate) * exp(kernel)) * 1.0e18; // nm2/s
	
	diffusionCoefficient[i] = 1;
	**/

	return;
}

void NECluster::resultFrom(ProductionReaction& reaction, IReactant& product) {
	// Cast the reacting clusters
	auto& cluster1 = static_cast<NECluster&>(reaction.first);
	auto& cluster2 = static_cast<NECluster&>(reaction.second);
	auto& prodCluster = static_cast<NECluster&>(product);

	// Compute the overlap
	int width1 = cluster1.getSectionWidth();
	int size1 = cluster1.getSize();
	int width2 = cluster2.getSectionWidth();
	int size2 = cluster2.getSize();
	int prodWidth = prodCluster.getSectionWidth(), prodSize =
			prodCluster.getSize();
	int lo1 = ((int) ((double) size1 - (double) width1 / 2.0) + 1), lo2 =
			((int) ((double) size2 - (double) width2 / 2.0) + 1), hi1 =
			((int) ((double) size1 + (double) width1 / 2.0)), hi2 =
			((int) ((double) size2 + (double) width2 / 2.0));
	int prodLo = ((int) ((double) prodSize - (double) prodWidth / 2.0) + 1),
			prodHi = ((int) ((double) prodSize + (double) prodWidth / 2.0));

	int overlap = std::min(prodHi, hi1 + hi2) - std::max(prodLo, lo1 + lo2) + 1;

	// Skip if the reaction doesn't overlap
	if (overlap < 1)
		return;

	// Create the pair
	ClusterPair pair(reaction, &cluster1, &cluster2);
	// Compute the coefficients
	pair.a0 = overlap;
	if (width1 > 1) {
		pair.a1 = 2.0
				* firstOrderSum(std::max(prodLo - lo2, lo1),
						std::min(prodHi - hi2, hi1), (double) (lo1 + hi1) / 2.0)
				/ (double) (hi1 - lo1);
	}
	if (width1 > 1 && width2 > 1) {
		// Should never happen for now
		std::cout << "Both reactants are super: " << cluster1.getName() << " + "
				<< cluster2.getName() << " -> " << name << std::endl;
	}

	// Add the pair
	reactingPairs.emplace_back(pair);

	return;
}

void NECluster::resultFrom(ProductionReaction& reaction, int[4], int[4]) {

	// Add a cluster pair for given reaction
	reactingPairs.emplace_back(
			reaction,  // TODO verify this is correct
			&static_cast<NECluster&>(reaction.first),
			&static_cast<NECluster&>(reaction.second));

	return;
}

void NECluster::resultFrom(ProductionReaction& reaction, double *coef) {

	// Add a cluster pair for given reaction
	reactingPairs.emplace_back(
			reaction,  // TODO verify this is correct
			&static_cast<NECluster&>(reaction.first),
			&static_cast<NECluster&>(reaction.second));
	auto& newPair = reactingPairs.back();

	// Update the coefs
	newPair.a0 = coef[0];
	newPair.a1 = coef[1];

	return;
}

void NECluster::participateIn(ProductionReaction& reaction,
		IReactant& product) {
	// Look for the other cluster
	auto& otherCluster = static_cast<NECluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);
	auto& prodCluster = static_cast<NECluster&>(product);

	// Compute the overlap
	int width1 = otherCluster.getSectionWidth();
	int size1 = otherCluster.getSize();
	int width2 = getSectionWidth();
	int size2 = getSize();
	int prodWidth = prodCluster.getSectionWidth(), prodSize =
			prodCluster.getSize();
	int lo1 = ((int) ((double) size1 - (double) width1 / 2.0) + 1), lo2 =
			((int) ((double) size2 - (double) width2 / 2.0) + 1), hi1 =
			((int) ((double) size1 + (double) width1 / 2.0)), hi2 =
			((int) ((double) size2 + (double) width2 / 2.0));
	int prodLo = ((int) ((double) prodSize - (double) prodWidth / 2.0) + 1),
			prodHi = ((int) ((double) prodSize + (double) prodWidth / 2.0));

	int overlap = std::min(prodHi, hi1 + hi2) - std::max(prodLo, lo1 + lo2) + 1;

	// Skip if the reaction doesn't overlap
	if (overlap < 1)
		return;

	// Create the pair
	CombiningCluster pair(reaction, &otherCluster);

	// Compute the coefficients
	pair.a0 = overlap;
	if (width1 > 1) {
		pair.a1 = 2.0
				* firstOrderSum(std::max(prodLo - lo2, lo1),
						std::min(prodHi - hi2, hi1), (double) (lo1 + hi1) / 2.0)
				/ (double) (hi1 - lo1);
	}

	// Add the combining cluster to list of clusters that combine with us
	combiningReactants.emplace_back(pair);

	return;
}

void NECluster::participateIn(ProductionReaction& reaction, int[4]) {
	// Look for the other cluster
	auto& otherCluster = static_cast<NECluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);
//	// Build a production reaction for it
//	std::unique_ptr<ProductionReaction> newReaction(
//			new ProductionReaction(otherCluster, *this));
//	// Add it to the network
//	auto& prref = network.add(std::move(newReaction));

	// Add the combining cluster to list of clusters that combine with us
	combiningReactants.emplace_back(reaction, &otherCluster);

	return;
}

void NECluster::participateIn(ProductionReaction& reaction, double *coef) {
	// Look for the other cluster
	auto& otherCluster = static_cast<NECluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);
//	// Build a production reaction for it
//	std::unique_ptr<ProductionReaction> newReaction(
//			new ProductionReaction(otherCluster, *this));
//	// Add it to the network
//	auto& prref = network.add(std::move(newReaction));

	// Add the combining cluster to list of clusters that combine with us
	combiningReactants.emplace_back(reaction, &otherCluster);
	auto& newComb = combiningReactants.back();

	// Update the coefs
	newComb.a0 = coef[0];
	newComb.a1 = coef[1];

	return;
}

void NECluster::participateIn(DissociationReaction& reaction,
		IReactant& disso) {
	// Look for the other cluster
	auto& emittedCluster = static_cast<NECluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);
	auto& dissoCluster = static_cast<NECluster&>(disso);

	// Compute the overlap
	int width1 = emittedCluster.getSectionWidth();
	int size1 = emittedCluster.getSize();
	int width2 = getSectionWidth();
	int size2 = getSize();
	int dissoWidth = dissoCluster.getSectionWidth(), dissoSize =
			dissoCluster.getSize();
	int lo1 = ((int) ((double) size1 - (double) width1 / 2.0) + 1), lo2 =
			((int) ((double) size2 - (double) width2 / 2.0) + 1), hi1 =
			((int) ((double) size1 + (double) width1 / 2.0)), hi2 =
			((int) ((double) size2 + (double) width2 / 2.0));
	int dissoLo = ((int) ((double) dissoSize - (double) dissoWidth / 2.0) + 1),
			dissoHi = ((int) ((double) dissoSize + (double) dissoWidth / 2.0));

	int overlap = std::min(dissoHi, hi1 + hi2) - std::max(dissoLo, lo1 + lo2)
			+ 1;

	// Skip if the reaction doesn't overlap
	if (overlap < 1)
		return;

	// Create the pair
	ClusterPair pair(reaction, &dissoCluster, &emittedCluster);
	// Compute the coefficients
	pair.a0 = overlap;
	if (dissoWidth > 1) {
		pair.a1 = 2.0
				* firstOrderSum(std::max(dissoLo, lo1 + lo2),
						std::min(dissoHi, hi1 + hi2),
						(double) (dissoLo + dissoHi) / 2.0)
				/ (double) (dissoHi - dissoLo);
	}

	// Add the pair to the vector
	dissociatingPairs.emplace_back(pair);

	return;
}

void NECluster::participateIn(DissociationReaction& reaction, int[4], int[4]) {
	// Look for the other cluster
	auto& emittedCluster = static_cast<NECluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);

	// Add a pair where it is important that the
	// dissociating cluster is the first one
	dissociatingPairs.emplace_back(
			reaction,  // TODO is this correct?
			&static_cast<NECluster&>(reaction.dissociating),
			&static_cast<NECluster&>(emittedCluster));

	return;
}

void NECluster::participateIn(DissociationReaction& reaction, double *coef) {
	// Look for the other cluster
	auto& emittedCluster = static_cast<NECluster&>(
			(reaction.first.getId() == id) ? reaction.second : reaction.first);

	// Add a pair where it is important that the
	// dissociating cluster is the first one
	dissociatingPairs.emplace_back(
			reaction,  // TODO is this correct?
			&static_cast<NECluster&>(reaction.dissociating),
			&static_cast<NECluster&>(emittedCluster));
	auto& newPair = dissociatingPairs.back();

	// Set the coefs
	newPair.a0 = coef[0];
	newPair.a1 = coef[1];

	return;
}

void NECluster::emitFrom(DissociationReaction& reaction, IReactant& disso) {
	// Cast the reacting clusters
	auto& cluster1 = static_cast<NECluster&>(reaction.first);
	auto& cluster2 = static_cast<NECluster&>(reaction.second);
	auto& dissoCluster = static_cast<NECluster&>(disso);

	// Compute the overlap
	int width1 = cluster1.getSectionWidth();
	int size1 = cluster1.getSize();
	int width2 = cluster2.getSectionWidth();
	int size2 = cluster2.getSize();
	int dissoWidth = dissoCluster.getSectionWidth(), dissoSize =
			dissoCluster.getSize();
	int lo1 = ((int) ((double) size1 - (double) width1 / 2.0) + 1), lo2 =
			((int) ((double) size2 - (double) width2 / 2.0) + 1), hi1 =
			((int) ((double) size1 + (double) width1 / 2.0)), hi2 =
			((int) ((double) size2 + (double) width2 / 2.0));
	int dissoLo = ((int) ((double) dissoSize - (double) dissoWidth / 2.0) + 1),
			dissoHi = ((int) ((double) dissoSize + (double) dissoWidth / 2.0));

	int overlap = std::min(dissoHi, hi1 + hi2) - std::max(dissoLo, lo1 + lo2)
			+ 1;

	// Skip if the reaction doesn't overlap
	if (overlap < 1)
		return;

	// Create the pair
	ClusterPair pair(reaction, &cluster1, &cluster2);

	// Compute the coefficients
	pair.a0 = overlap;

	// Add the pair
	emissionPairs.emplace_back(pair);

	return;
}

void NECluster::emitFrom(DissociationReaction& reaction, int[4]) {

	// Add the pair of emitted clusters.
	emissionPairs.emplace_back(
			reaction, // TODO is this correct?
			&static_cast<NECluster&>(reaction.first),
			&static_cast<NECluster&>(reaction.second));

	return;
}

void NECluster::emitFrom(DissociationReaction& reaction, double *coef) {

	// Add the pair of emitted clusters.
	emissionPairs.emplace_back(
			reaction, // TODO is this correct?
			&static_cast<NECluster&>(reaction.first),
			&static_cast<NECluster&>(reaction.second));
	auto& newPair = emissionPairs.back();

	// Set the coefs
	newPair.a0 = coef[0];

	return;
}

static std::vector<int> getFullConnectivityVector(std::set<int> connectivitySet,
		int size) {
	// Create a vector of zeroes with size equal to the network size
	std::vector<int> connectivity(size);

	// Set the value of the connectivity array to one for each element that is
	// in the set.
	for (auto it = connectivitySet.begin(); it != connectivitySet.end(); ++it) {
		connectivity[(*it) - 1] = 1;
	}

	return connectivity;
}

std::vector<int> NECluster::getReactionConnectivity() const {
	// Create the full vector from the set and return it
	return getFullConnectivityVector(reactionConnectivitySet, network.getDOF());
}

std::vector<int> NECluster::getDissociationConnectivity() const {
	// Create the full vector from the set and return it
	return getFullConnectivityVector(dissociationConnectivitySet,
			network.getDOF());
}

void NECluster::resetConnectivities() {
	// Shrink the arrays to save some space. (About 10% or so.)
	reactingPairs.shrink_to_fit();
	combiningReactants.shrink_to_fit();
	dissociatingPairs.shrink_to_fit();
	emissionPairs.shrink_to_fit();

	// Clear both sets
	reactionConnectivitySet.clear();
	dissociationConnectivitySet.clear();

	// Connect this cluster to itself since any reaction will affect it
	setReactionConnectivity(id);
	setDissociationConnectivity(id);
	setReactionConnectivity(momId[0]);
	setDissociationConnectivity(momId[0]);

	// Apply to each reacting pair.
	std::for_each(reactingPairs.begin(), reactingPairs.end(),
			[this](const ClusterPair& currPair) {
				// The cluster is connecting to both clusters in the pair
				setReactionConnectivity(currPair.first->id);
				setReactionConnectivity(currPair.first->momId[0]);
				setReactionConnectivity(currPair.second->id);
				setReactionConnectivity(currPair.second->momId[0]);
			});

	// Apply to each combining cluster.
	std::for_each(combiningReactants.begin(), combiningReactants.end(),
			[this](const CombiningCluster& cc) {
				// The cluster is connecting to the combining cluster
				NECluster const& combCluster = *cc.combining;
				setReactionConnectivity(combCluster.id);
				setReactionConnectivity(combCluster.momId[0]);
			});

	// Apply to each effective dissociating pair
	std::for_each(dissociatingPairs.begin(), dissociatingPairs.end(),
			[this](const ClusterPair& currPair) {
				// The cluster is connecting to the dissociating cluster which
				// is the first one by definition
				setDissociationConnectivity(currPair.first->id);
				setDissociationConnectivity(currPair.first->momId[0]);
			});

	// Don't apply to the emission pairs because
	// this cluster is not connected to them

	return;
}

double NECluster::getMoment() const {
	return 0.0;
}

double NECluster::getTotalFlux(int i) {
	// Get the fluxes
	double prodFlux = getProductionFlux(i);
	double dissFlux = getDissociationFlux(i);
	double combFlux = getCombinationFlux(i);
	double emissFlux = getEmissionFlux(i);

	return prodFlux - combFlux + dissFlux - emissFlux;
}

double NECluster::getDissociationFlux(int xi) const {

	// Sum dissociation flux over all pairs that dissociate to form this one.
	double flux =
			std::accumulate(dissociatingPairs.begin(), dissociatingPairs.end(),
					0.0, [&xi](double running, const ClusterPair& currPair) {
						// Get the dissociating clusters
					NECluster* dissociatingCluster = currPair.first;
					double l0A = dissociatingCluster->getConcentration();
					double l1A = dissociatingCluster->getMoment();
					// Update the flux
					return running + (currPair.reaction.kConstant[xi] * (currPair.a0 * l0A + currPair.a1 * l1A));
				});

	// Return the flux
	return flux;
}

double NECluster::getEmissionFlux(int xi) const {

	// Sum reaction rate constants over all emission pair reactions.
	double flux = std::accumulate(emissionPairs.begin(), emissionPairs.end(),
			0.0, [&xi](double running, const ClusterPair& currPair) {
				return running + currPair.reaction.kConstant[xi] * currPair.a0;
			});

	return flux * concentration;
}

double NECluster::getProductionFlux(int xi) const {
	// Local declarations
	double flux = 0.0;

	// Sum over all the reacting pairs
	std::for_each(reactingPairs.begin(), reactingPairs.end(),
			[&flux,&xi](ClusterPair const& currPair) {
				// Get the two reacting clusters
				NECluster* firstReactant = currPair.first;
				NECluster* secondReactant = currPair.second;
				// We know the first one is always the single one
				double l0A = firstReactant->getConcentration();
				double l0B = secondReactant->getConcentration();
				double l1B = secondReactant->getMoment();
				// Update the flux
				flux += currPair.reaction.kConstant[xi] * l0A * (currPair.a0 * l0B + currPair.a1 * l1B);
			});

	// Return the production flux
	return flux;
}

double NECluster::getCombinationFlux(int xi) const {

	double flux = std::accumulate(combiningReactants.begin(),
			combiningReactants.end(), 0.0,
			[xi](double running, const CombiningCluster& currPair) {
				// Get the cluster that combines with this one
				NECluster const& combiningCluster = *currPair.combining;
				double l0A = combiningCluster.getConcentration();
				double l1A = combiningCluster.getMoment();
				// Update the flux
				return running +
				(currPair.reaction.kConstant[xi] *
						(currPair.a0 * l0A + currPair.a1 * l1A));
			});

	return flux * concentration;
}

std::vector<double> NECluster::getPartialDerivatives(int i) const {
	// Local Declarations
	std::vector<double> partials(network.getDOF(), 0.0);

	// Get the partial derivatives for each reaction type
	getProductionPartialDerivatives(partials, i);
	getCombinationPartialDerivatives(partials, i);
	getDissociationPartialDerivatives(partials, i);
	getEmissionPartialDerivatives(partials, i);

	return partials;
}

void NECluster::getPartialDerivatives(std::vector<double> & partials,
		int i) const {
	// Get the partial derivatives for each reaction type
	getProductionPartialDerivatives(partials, i);
	getCombinationPartialDerivatives(partials, i);
	getDissociationPartialDerivatives(partials, i);
	getEmissionPartialDerivatives(partials, i);

	return;
}

void NECluster::getProductionPartialDerivatives(std::vector<double> & partials,
		int xi) const {

	// Production
	// A + B --> D, D being this cluster
	// The flux for D is
	// F(C_D) = k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_D)/dC_A = k+_(A,B)*C_B
	// dF(C_D)/dC_B = k+_(A,B)*C_A
	std::for_each(reactingPairs.begin(), reactingPairs.end(),
			[&partials,&xi](ClusterPair const& currPair) {
				// Get the two reacting clusters
				NECluster* firstReactant = currPair.first;
				NECluster* secondReactant = currPair.second;
				double l0A = firstReactant->getConcentration();
				double l0B = secondReactant->getConcentration();
				double l1B = secondReactant->getMoment();

				// Compute the contribution from the first part of the reacting pair
				auto value = currPair.reaction.kConstant[xi];
				auto index = currPair.first->id - 1;
				partials[index] += value * (currPair.a0 * l0B + currPair.a1 * l1B);
				// Compute the contribution from the second part of the reacting pair
				index = currPair.second->id - 1;
				partials[index] += value * currPair.a0 * l0A;
				index = currPair.second->momId[0] - 1;
				partials[index] += value * currPair.a1 * l0A;
			});

	return;
}

void NECluster::getCombinationPartialDerivatives(std::vector<double> & partials,
		int xi) const {

	// Combination
	// A + B --> D, A being this cluster
	// The flux for A is outgoing
	// F(C_A) = - k+_(A,B)*C_A*C_B
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k+_(A,B)*C_B
	// dF(C_A)/dC_B = - k+_(A,B)*C_A
	std::for_each(combiningReactants.begin(), combiningReactants.end(),
			[this,&partials,&xi](const CombiningCluster& cc) {
				NECluster const& cluster = *cc.combining;
				Reaction const& currReaction = cc.reaction;
				double l0A = cluster.getConcentration();
				double l1A = cluster.getMoment();

				// Remember that the flux due to combinations is OUTGOING (-=)!
				// Compute the contribution from this cluster
				partials[id - 1] -= currReaction.kConstant[xi] * (cc.a0 * l0A + cc.a1 * l1A);
				// Compute the contribution from the combining cluster
				double value = currReaction.kConstant[xi] * concentration;

				partials[cluster.id - 1] -= value * cc.a0;
				partials[cluster.momId[0] - 1] -= value * cc.a1;
			});

	return;
}

void NECluster::getDissociationPartialDerivatives(
		std::vector<double> & partials, int xi) const {

	// Dissociation
	// A --> B + D, B being this cluster
	// The flux for B is
	// F(C_B) = k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_B)/dC_A = k-_(B,D)
	std::for_each(dissociatingPairs.begin(), dissociatingPairs.end(),
			[&partials,&xi](const ClusterPair& currPair) {
				// Get the dissociating cluster
				NECluster* cluster = currPair.first;
				Reaction const& currReaction = currPair.reaction;
				partials[cluster->id - 1] += currReaction.kConstant[xi] * currPair.a0;
				partials[cluster->momId[0] - 1] += currReaction.kConstant[xi] *
				currPair.a1;
			});

	return;
}

void NECluster::getEmissionPartialDerivatives(std::vector<double> & partials,
		int xi) const {

	// Emission
	// A --> B + D, A being this cluster
	// The flux for A is
	// F(C_A) = - k-_(B,D)*C_A
	// Thus, the partial derivatives
	// dF(C_A)/dC_A = - k-_(B,D)
	double emissionFlux = std::accumulate(emissionPairs.begin(),
			emissionPairs.end(), 0.0,
			[&xi](double running, const ClusterPair& currPair) {
				return running + currPair.reaction.kConstant[xi] * currPair.a0;
			});

	// Recall emission flux is OUTGOING
	partials[id - 1] -= emissionFlux;

	return;
}

double NECluster::getLeftSideRate(int i) const {

	// Sum reaction rate contributions over all combining clusters.
	double combiningRateTotal = std::accumulate(combiningReactants.begin(),
			combiningReactants.end(), 0.0,
			[&i](double running, const CombiningCluster& currPair) {
				NECluster const& cluster = *currPair.combining;

				return running + (currPair.reaction.kConstant[i] *
						cluster.concentration);
			});

	// Sum reaction rate constants over all emission pairs.
	double emissionRateTotal = std::accumulate(emissionPairs.begin(),
			emissionPairs.end(), 0.0,
			[&i](double running, const ClusterPair& currPair) {
				return running + currPair.reaction.kConstant[i];
			});

	return combiningRateTotal + emissionRateTotal;
}

std::vector<std::vector<double> > NECluster::getProdVector() const {
	// Initial declarations
	std::vector<std::vector<double> > toReturn;

	// Loop on the reacting pairs
	std::for_each(reactingPairs.begin(), reactingPairs.end(),
			[&toReturn](const ClusterPair& currPair) {
				// Build the vector containing ids and rates
				std::vector<double> tempVec;
				tempVec.push_back(currPair.first->getId() - 1);
				tempVec.push_back(currPair.second->getId() - 1);
				tempVec.push_back(currPair.a0);
				tempVec.push_back(currPair.a1);

				// Add it to the main vector
				toReturn.push_back(tempVec);
			});

	return toReturn;
}

std::vector<std::vector<double> > NECluster::getCombVector() const {
	// Initial declarations
	std::vector<std::vector<double> > toReturn;

	// Loop on the combining reactants
	std::for_each(combiningReactants.begin(), combiningReactants.end(),
			[&toReturn](const CombiningCluster& cc) {
				// Build the vector containing ids and rates
				std::vector<double> tempVec;
				tempVec.push_back(cc.combining->getId() - 1);
				tempVec.push_back(cc.a0);
				tempVec.push_back(cc.a1);

				// Add it to the main vector
				toReturn.push_back(tempVec);
			});

	return toReturn;
}

std::vector<std::vector<double> > NECluster::getDissoVector() const {
	// Initial declarations
	std::vector<std::vector<double> > toReturn;

	// Loop on the dissociating pairs
	std::for_each(dissociatingPairs.begin(), dissociatingPairs.end(),
			[&toReturn](const ClusterPair& currPair) {
				// Build the vector containing ids and rates
				std::vector<double> tempVec;
				tempVec.push_back(currPair.first->getId() - 1);
				tempVec.push_back(currPair.second->getId() - 1);
				tempVec.push_back(currPair.a0);
				tempVec.push_back(currPair.a1);

				// Add it to the main vector
				toReturn.push_back(tempVec);
			});

	return toReturn;
}

std::vector<std::vector<double> > NECluster::getEmitVector() const {
	// Initial declarations
	std::vector<std::vector<double> > toReturn;

	// Loop on the emitting pairs
	std::for_each(emissionPairs.begin(), emissionPairs.end(),
			[&toReturn](const ClusterPair& currPair) {
				// Build the vector containing ids and rates
				std::vector<double> tempVec;
				tempVec.push_back(currPair.first->getId() - 1);
				tempVec.push_back(currPair.second->getId() - 1);
				tempVec.push_back(currPair.a0);

				// Add it to the main vector
				toReturn.push_back(tempVec);
			});

	return toReturn;
}

std::vector<int> NECluster::getConnectivity() const {
	int connectivityLength = network.getDOF();
	std::vector<int> connectivity = std::vector<int>(connectivityLength, 0);
	auto reactionConnVector = getReactionConnectivity();
	auto dissociationConnVector = getDissociationConnectivity();

	// The reaction and dissociation vectors must have a length equal to the
	// number of clusters
	if (reactionConnVector.size() != (unsigned int) connectivityLength) {
		throw std::string("The reaction vector has an incorrect length");
	}
	if (dissociationConnVector.size() != (unsigned int) connectivityLength) {
		throw std::string("The dissociation vector has an incorrect length");
	}

	// Merge the two vectors such that the final vector contains
	// a 1 at a position if either of the connectivity arrays
	// have a 1
	for (int i = 0; i < connectivityLength; i++) {
		// Consider each connectivity array only if its type is enabled
		connectivity[i] = reactionConnVector[i] || dissociationConnVector[i];
	}

	return connectivity;
}
