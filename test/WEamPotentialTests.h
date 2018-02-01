#include<memory>
#include"gtest/gtest.h"
#include"EamSetflInteractions.h"
#include"BccUnitCell.h"
#include"RectangularSimulationBox.h"
#include"PerfectCrystal.h"

class WEamPotentialTests : public ::testing::Test{
	protected:
		static std::shared_ptr<Interactions> interactions;
		static std::shared_ptr<RectangularUnitCell> unitCell;
		static std::shared_ptr<RectangularUnitCell> rotatedUnitCell;
		static std::unique_ptr<PerfectCrystal> tungstenCrystal;
		static std::unique_ptr<PerfectCrystal> rotatedTungstenCrystal;
		static void SetUpTestCase(){
			double latConst=3.14;
			std::string setflFile="w_eam2.fs";
			interactions = std::shared_ptr<Interactions>(EamSetflInteractions::Create(setflFile, 1e-5).release());
			unitCell = std::shared_ptr<RectangularUnitCell>(BccUnitCell::Create(latConst, false).release());
			rotatedUnitCell = std::shared_ptr<RectangularUnitCell>(BccUnitCell::Create(latConst, true).release());
			tungstenCrystal=std::unique_ptr<PerfectCrystal>(new PerfectCrystal(unitCell, interactions));
			rotatedTungstenCrystal=std::unique_ptr<PerfectCrystal>(new PerfectCrystal(rotatedUnitCell, interactions));
		}
};

std::shared_ptr<Interactions> WEamPotentialTests::interactions = NULL;
std::shared_ptr<RectangularUnitCell> WEamPotentialTests::unitCell = NULL;
std::shared_ptr<RectangularUnitCell> WEamPotentialTests::rotatedUnitCell = NULL;
std::unique_ptr<PerfectCrystal> WEamPotentialTests::tungstenCrystal = NULL;
std::unique_ptr<PerfectCrystal> WEamPotentialTests::rotatedTungstenCrystal = NULL;

TEST_F (WEamPotentialTests, BindingEnergy){
	EXPECT_FLOAT_EQ (-8.89998, tungstenCrystal -> EvaluatePotentialEnergyPerAtom());
	EXPECT_FLOAT_EQ (-8.89998, rotatedTungstenCrystal -> EvaluatePotentialEnergyPerAtom());
}

TEST_F (WEamPotentialTests, ForceConstantsSumRule){
	EXPECT_NEAR(tungstenCrystal -> CheckIfAsrIsObeyed(), 0, 1e-9);
}
