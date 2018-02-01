#include<memory>
#include"gtest/gtest.h"
#include"HarmonicInteractions.h"
#include"ScUnitCell.h"
#include"PerfectCrystal.h"

TEST (NnScPotentialTests, ForceConstantsSumRuleNNSc){
	double harmLatConst=.491238;
	double forceConst=1.0912384;
	double mass=1.;
	std::shared_ptr<Interactions> interactions(HarmonicInteractions::Create(forceConst, harmLatConst, mass).release());
	std::shared_ptr<RectangularUnitCell> lattice(ScUnitCell::Create(harmLatConst).release());
	std::unique_ptr<PerfectCrystal>  crystal=std::unique_ptr<PerfectCrystal>(new PerfectCrystal(lattice, interactions));
	EXPECT_NEAR(crystal -> CheckIfAsrIsObeyed(), 0, 1e-9);
}


