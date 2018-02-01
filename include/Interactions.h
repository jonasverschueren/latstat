// Calculate Interactions between particles
#pragma once

#include<vector>
#include"Eigen/Dense"
#include"RectangularSimulationBox.h"
#include<memory>

class Interactions{
	private:
		double interactionCutoff;
		double numericalStep;
		int dim=-1;								//sentinel value equivalent to 3d forces. Set to 0,1 for 1D and 2D force constant evaluation
	protected:
		double mass = 0;
	public:
		Interactions();
		Interactions(double numericalStep, int dim);
		Interactions(double interactionCutoff, double numericalStep=1e-4, int dim=-1);
		virtual double EvaluatePotentialEnergy(RectangularSimulationBox&)=0;
		Eigen::Vector3d EvaluateForceOnParticleNumerically(RectangularSimulationBox& rectangularSimulationBox, int pid);
		virtual Eigen::Vector3d EvaluateForceOnParticle(RectangularSimulationBox& rectangularSimulationBox, int pid)=0;
		virtual Eigen::MatrixXd EvaluateForceConstants(RectangularSimulationBox&, int pid);
		double GetInteractionCutoff();
		void SetNumericalStep(double step);
		void SetCutoff(double cutoff);
		void SetDim(int dim);
		int GetDim();
		virtual double GetMass();
};
