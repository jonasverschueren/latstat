#include<memory>
#include"gtest/gtest.h"
#include"EamSetflInteractions.h"
#include"FccUnitCell.h"
#include"RectangularSimulationBox.h"
#include"PerfectCrystal.h"

class AlPotentialTests : public ::testing::Test{
	protected:
		static std::shared_ptr<Interactions> interactions;
		static std::shared_ptr<RectangularUnitCell> unitCell;
		static std::shared_ptr<RectangularUnitCell> rotatedUnitCell;
		static std::unique_ptr<PerfectCrystal> aluCrystal;
		static std::unique_ptr<PerfectCrystal> rotatedAluCrystal;
		static void SetUpTestCase(){
			double latConst=4.05;
			std::string setflFile="Al99.eam.alloy";
			interactions = std::shared_ptr<Interactions>(EamSetflInteractions::Create(setflFile, 1e-5).release());
			unitCell = std::shared_ptr<RectangularUnitCell>(FccUnitCell::Create(latConst, false).release());
			rotatedUnitCell = std::shared_ptr<RectangularUnitCell>(FccUnitCell::Create(latConst, true).release());
			aluCrystal=std::unique_ptr<PerfectCrystal>(new PerfectCrystal(unitCell, interactions));
			rotatedAluCrystal=std::unique_ptr<PerfectCrystal>(new PerfectCrystal(rotatedUnitCell, interactions));
		}
};

std::shared_ptr<Interactions> AlPotentialTests::interactions = NULL;
std::shared_ptr<RectangularUnitCell> AlPotentialTests::unitCell = NULL;
std::shared_ptr<RectangularUnitCell> AlPotentialTests::rotatedUnitCell = NULL;
std::unique_ptr<PerfectCrystal> AlPotentialTests::aluCrystal = NULL;
std::unique_ptr<PerfectCrystal> AlPotentialTests::rotatedAluCrystal = NULL;

TEST_F (AlPotentialTests, BindingEnergy){
	EXPECT_FLOAT_EQ (-3.36, aluCrystal -> EvaluatePotentialEnergyPerAtom());
	EXPECT_FLOAT_EQ (-3.36, rotatedAluCrystal -> EvaluatePotentialEnergyPerAtom());
}

TEST_F (AlPotentialTests, ForceConstantsSumRule){
	EXPECT_NEAR(aluCrystal -> CheckIfAsrIsObeyed(), 0, 1e-9);
}
