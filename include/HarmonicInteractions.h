// Class Harmonic nearest-neighbour type interactions
#pragma once

#include<string>
#include<vector>
#include<map>
#include<cmath>
#include<memory>
#include"Interactions.h"
#include"RectangularSimulationBox.h"

class HarmonicInteractions: public Interactions{
	private:
		double forceConst;						// force constant
		double eqmSep;					// eqm seperation
		HarmonicInteractions(double forceConst, double eqmDistance, double mass, double numericalStep=1e-4, int dim=-1);
	public:
		HarmonicInteractions();
		static std::unique_ptr<Interactions> Create(double forceConst, double eqmDistance, double mass, double numericalStep=1e-4, int dim=-1);
		double EvaluatePotentialEnergy(RectangularSimulationBox&);
		Eigen::Vector3d EvaluateForceOnParticle(RectangularSimulationBox& rectangularSimulationBox, int pid);
		Eigen::MatrixXd EvaluateForceConstants(RectangularSimulationBox&, int pid);
		double GetForceConst();
		double GetEqmSep();
};

