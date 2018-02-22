// Mostly virtual class for set of particles with an interaction potential in periodic simulation box
#pragma once

#include<memory>
#include"RectangularSimulationBox.h"
#include"Interactions.h"

class Crystal{
	protected:
		Crystal(){;}
		Crystal(std::shared_ptr<Interactions> interactions);
		std::shared_ptr<Interactions> ptr_interactions;
	public:
		double GetMass();
};
