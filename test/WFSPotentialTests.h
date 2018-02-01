// Unit Test to make sure potential energy of Finnis-Sinclair tungsten is what it should be

#include<memory>
#include"gtest/gtest.h"
#include"FsInteractions.h"
#include"BccUnitCell.h"
#include"RectangularSimulationBox.h"
#include"PerfectCrystal.h"
#include "Eigen/Dense"

class WPotentialTests : public ::testing::Test{
	protected:
		static std::shared_ptr<Interactions> interactions;
		static std::shared_ptr<RectangularUnitCell> unitCell;
		static std::shared_ptr<RectangularUnitCell> rotatedUnitCell;
		static std::unique_ptr<PerfectCrystal> tungstenCrystal;
		static std::unique_ptr<PerfectCrystal> rotatedTungstenCrystal;
		static void SetUpTestCase(){
			double latConst=3.1652;
			double tungstenMass=183.84;
			std::string fsParamsFile ="W87.eam.fs";
			interactions = std::shared_ptr<Interactions>(FsInteractions::Create(fsParamsFile).release());
			unitCell = std::shared_ptr<RectangularUnitCell>(BccUnitCell::Create(latConst).release());
			rotatedUnitCell = std::shared_ptr<RectangularUnitCell>(BccUnitCell::Create(latConst, true).release());
			tungstenCrystal=std::unique_ptr<PerfectCrystal>(new PerfectCrystal(unitCell, interactions));
			rotatedTungstenCrystal=std::unique_ptr<PerfectCrystal>(new PerfectCrystal(rotatedUnitCell, interactions));
		}
};

std::shared_ptr<Interactions> WPotentialTests::interactions = NULL;
std::shared_ptr<RectangularUnitCell> WPotentialTests::unitCell = NULL;
std::shared_ptr<RectangularUnitCell> WPotentialTests::rotatedUnitCell = NULL;
std::unique_ptr<PerfectCrystal> WPotentialTests::tungstenCrystal = NULL;
std::unique_ptr<PerfectCrystal> WPotentialTests::rotatedTungstenCrystal = NULL;

TEST_F (WPotentialTests, BindingEnergy){
	EXPECT_FLOAT_EQ (-8.90, tungstenCrystal -> EvaluatePotentialEnergyPerAtom());
	EXPECT_FLOAT_EQ (-8.90, rotatedTungstenCrystal -> EvaluatePotentialEnergyPerAtom());
}

TEST_F (WPotentialTests, AnalyticalNumericalForcesDifference){
	RectangularSimulationBox simBox;
	simBox.GenerateSimulationBox(unitCell, 3*interactions->GetInteractionCutoff());
	simBox.MoveParticleByAmount(2, 0, .18492);
	simBox.MoveParticleByAmount(2, 1, .164020);
	simBox.MoveParticleByAmount(2, 2, .104947);
	Eigen::Vector3d analyticalForce = interactions->EvaluateForceOnParticle(simBox, 2);
	Eigen::Vector3d numForce = interactions->EvaluateForceOnParticleNumerically(simBox, 2);
	for(int i=0; i<=2; i++){
		EXPECT_NEAR(analyticalForce[i], numForce[i], 1e-7);
	}
}

TEST_F (WPotentialTests, ForceConstantsSumRule){
	EXPECT_NEAR(tungstenCrystal -> CheckIfAsrIsObeyed(), 0, 1e-9);
}

TEST_F (WPotentialTests, DispersionZeroAtGamma){
	Eigen::Vector3d k(0, 0, 0);
	Eigen::Vector3d freqs = tungstenCrystal -> EvaluatePhononFrequencySquared(k);
	for(int i=0; i<freqs.size(); i++){
		EXPECT_DOUBLE_EQ(freqs[i], 0);
	}
}

TEST_F (WPotentialTests, DispersionZeroAtBZBoundaries){
	Eigen::Vector3d k(2, 5, 9);
	Eigen::Vector3d freqs = tungstenCrystal -> EvaluatePhononFrequencySquared(k);
	for(int i=0; i<freqs.size(); i++){
		EXPECT_NEAR(freqs[i], 0, 1e-12);
	}
}
