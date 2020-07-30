#ifndef UZRMATERIALHANDLERFACTORY_H
#define UZRMATERIALHANDLERFACTORY_H

#include <memory>
#include <MaterialFactory.h>
#include <FuelFitFluxHandler.h>
#include <DummyAdvectionHandler.h>
#include <DummyTrapMutationHandler.h>

namespace xolotlFactory {

/**
 * Subclass of MaterialFactory for an UZr material.
 */
class UZrMaterialFactory: public MaterialFactory {
public:

	/**
	 * The constructor creates the handlers.
	 *
	 * @param dim The number of dimensions for the problem
	 */
	UZrMaterialFactory(const xolotlCore::Options &options) :
			MaterialFactory(options) {
		theFluxHandler = std::make_shared<xolotlCore::FuelFitFluxHandler>();
		theAdvectionHandler.push_back(
				std::make_shared<xolotlCore::DummyAdvectionHandler>());
		theTrapMutationHandler = std::make_shared<
				xolotlCore::DummyTrapMutationHandler>();
		theReSolutionHandler = std::make_shared<
				xolotlCore::DummyReSolutionHandler>();
		theNucleationHandler = std::make_shared<
				xolotlCore::DummyNucleationHandler>();

		return;
	}

	/**
	 * The destructor
	 */
	 void initializeMaterial(const xolotlCore::Options &options) {
 		// Call the mother method
 		MaterialFactory::initializeMaterial(options);

 		// Change the flux amplitude because we have to take into account
 		// that there are one xenon created every 4 fissions.
 		theFluxHandler->setFluxAmplitude(
 				options.getFluxAmplitude() * options.getFissionYield());

 		// Pass the fission yield to the re-solution and heterogenenous nucletation handlers
 		theReSolutionHandler->setFissionYield(options.getFissionYield());
 		theNucleationHandler->setFissionYield(options.getFissionYield());

 		return;
 	}



	~UZrMaterialFactory() {
	}
};

} // end namespace xolotlFactory

#endif // UZRMATERIALHANDLERFACTORY_H
