#include"ImperfectCrystal.h"

ImperfectCrystal::ImperfectCrystal(const char* xyz_filename, std::shared_ptr<Interactions> read_interactions) : Crystal(read_interactions){
	RectangularSimulationBox::ReadFromXYZDump(xyz_filename, this->simBox);
}

ImperfectCrystal::ImperfectCrystal(RectangularSimulationBox& read_simBox, std::shared_ptr<Interactions> read_interactions) : Crystal(read_interactions), simBox(read_simBox){
	;
}

double ImperfectCrystal::EvaluateTotalPotentialEnergy(){
	return ptr_interactions->EvaluatePotentialEnergy(this->simBox);
}
