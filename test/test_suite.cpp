#include"gtest/gtest.h"
#include"NnScPotentialTests.h"
#include"KoizumiCosinusoidalPotentialTests.h"
#include"WFSPotentialTests.h"
#include"WEamPotentialTests.h"
#include"AlPotentialTests.h"

int main(int argc, char* argv[]){
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
