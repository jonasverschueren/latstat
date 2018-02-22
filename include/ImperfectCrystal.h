// Loosely speaking: ImperfectCrystal = RectangularSimulationBox + Interactions
#pragma once

#include<memory>
#include"Crystal.h"
#include"RectangularSimulationBox.h"
#include"Interactions.h"

class ImperfectCrystal : public Crystal{
	protected:
		ImperfectCrystal(){;}
	public:
		RectangularSimulationBox simBox;
		ImperfectCrystal(const char* xyz_filename, std::shared_ptr<Interactions>);
		ImperfectCrystal(RectangularSimulationBox&, std::shared_ptr<Interactions>);
		double EvaluateTotalPotentialEnergy();
};
