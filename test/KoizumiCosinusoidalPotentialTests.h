#include<memory>
#include<fstream>
#include"gtest/gtest.h"
#include"HarmonicInteractions.h"
#include"KoizumiInteractions.h"
#include"SquareUnitCell.h"
#include"PerfectCrystal.h"

TEST (KoizumiCosinusoidalPotentialTests, BindingEnergyKoizumiSquare){
	double harmLatConst=1.;
	double forceConst=1.;
	double mass=1.;
	std::shared_ptr<Interactions> interactions(KoizumiInteractions::Create(forceConst, harmLatConst, mass).release());
	std::shared_ptr<RectangularUnitCell> lattice(SquareUnitCell::Create(harmLatConst).release());
	std::unique_ptr<PerfectCrystal>  crystal=std::unique_ptr<PerfectCrystal>(new PerfectCrystal(lattice, interactions));
	EXPECT_FLOAT_EQ (0, crystal -> EvaluatePotentialEnergyPerAtom());
}

TEST (KoizumiCosinusoidalPotentialTests, AnalyticalNumericalForcesDifferenceKoizumiSquare){
	double eqmSep=1.9123;
	double forceConst=.492;
	double mass=4.30;
	std::shared_ptr<Interactions> interactions(KoizumiInteractions::Create(forceConst, eqmSep, mass, 1e-5).release());
	std::shared_ptr<RectangularUnitCell> lattice(SquareUnitCell::Create(eqmSep).release());
	std::unique_ptr<PerfectCrystal>  crystal=std::unique_ptr<PerfectCrystal>(new PerfectCrystal(lattice, interactions));
	RectangularSimulationBox simBox;
	Eigen::Vector3d minimumLengths(2*interactions->GetInteractionCutoff(), 2*interactions->GetInteractionCutoff(), 0);
	simBox.GenerateSimulationBox(lattice, minimumLengths);
	simBox.MoveParticleByAmount(2, 0, .18492);
	simBox.MoveParticleByAmount(2, 1, .164020);
	simBox.MoveParticleByAmount(2, 2, .104947);
	Eigen::Vector3d analyticalForce = interactions->EvaluateForceOnParticle(simBox, 2);
	Eigen::Vector3d numForce = interactions->EvaluateForceOnParticleNumerically(simBox, 2);
	for(int i=0; i<=2; i++){
		EXPECT_NEAR(analyticalForce[i], numForce[i], 1e-7);
	}
}

TEST (KoizumiCosinusoidalPotentialTests, DispersionKoizumiSquare){
	double eqmSep=1.9123;
	double forceConst=.492;
	double mass=4.30;
	std::shared_ptr<Interactions> interactions(KoizumiInteractions::Create(forceConst, eqmSep, mass, 1e-5).release());
	std::shared_ptr<RectangularUnitCell> lattice(SquareUnitCell::Create(eqmSep).release());
	std::unique_ptr<PerfectCrystal>  crystal=std::unique_ptr<PerfectCrystal>(new PerfectCrystal(lattice, interactions));
	EXPECT_NEAR(crystal -> CheckIfAsrIsObeyed(), 0, 1e-9);
	std::vector<Eigen::Vector3d> kvects;
	std::vector<double>freqs;
	Eigen::Vector3d k(1, 1, 1);
	for(double kmag=0; kmag<5; kmag+=.01){
		Eigen::Vector3d phononFrequencies = crystal->EvaluatePhononFrequencySquared(k*kmag).real();
		kvects.push_back(k*kmag);
		freqs.push_back(phononFrequencies[2]);
	}
	/* fit expected function to dispersion */
	for(int i=0; i<kvects.size(); i++){
		Eigen::Vector3d k = kvects.at(i);
		EXPECT_NEAR(freqs.at(i), 16*M_PI*M_PI*forceConst/(mass*eqmSep*eqmSep)*(pow(sin(M_PI*k[0]),2) + pow(sin(M_PI*k[1]),2)), 1e-8);
	}
}
