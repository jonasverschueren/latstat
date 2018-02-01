// Class transverse harmonic nearest-neighbour type interactions
// U = forceConst * (1 - cos(2*pi*u_ij/eqmDistance))
// taken from Koizumi et al., PRB 65(21)
#pragma once

#include<string>
#include<vector>
#include<map>
#include<cmath>
#include<memory>
#include"Interactions.h"
#include"RectangularSimulationBox.h"

class KoizumiInteractions: public Interactions{
	private:
		double forceConst;						// force constant
		double eqmSep;					// eqm seperation
		KoizumiInteractions(double forceConst, double eqmDistance, double mass, double numericalStep=1e-4);
	public:
		KoizumiInteractions();
		static std::unique_ptr<Interactions> Create(double forceConst, double eqmDistance, double mass, double numericalStep=1e-4);
		double EvaluatePotentialEnergy(RectangularSimulationBox&);
		Eigen::Vector3d EvaluateForceOnParticle(RectangularSimulationBox& rectangularSimulationBox, int pid);
		Eigen::MatrixXd EvaluateForceConstants(RectangularSimulationBox&, int pid);
		double GetForceConst();
		double GetEqmSep();
};
