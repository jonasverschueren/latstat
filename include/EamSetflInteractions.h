// Class EAM potentials fom setfl tabulated values
#pragma once

#include<string>
#include<sstream>
#include<stdlib.h>
#include<fstream>
#include<vector>
#include<map>
#include<cmath>
#include<memory>
#include"Eigen/Dense"
#include"Interactions.h"
#include"RectangularSimulationBox.h"

class EamSetflInteractions: public Interactions{
	private:
		int Nrho;					// number of density points in embedding function tabulation
		int Nr;						// number of potential and electron density points
		double dr;					//distance between different pair and density terms in tabulated values
		double drho;					//distance between different embedding function in tabulated values
		double cutoff;					// potential cutoff
		Eigen::VectorXd densityValues;			//electron density terms from file
		Eigen::VectorXd embeddingValues;		// embedding function values from file
		Eigen::VectorXd pairValues;			// pair interaction values from file
		Eigen::VectorXd y2densityValues;		// 2nd derivatives necessary for cubic spline interpolation
		Eigen::VectorXd y2embeddingValues;
		Eigen::VectorXd y2pairValues;
		void ReadSetflFile(std::string setflFileName);
		void ReadSetflFileHeader(std::string setflFileName);
		void Evaluate2ndDerivatives();
		Eigen::VectorXd EvaluateSpline2ndDerivative(Eigen::VectorXd&, double dx);
		double EvaluateCubicSpline(double x, Eigen::VectorXd& yvals, Eigen::VectorXd& y2ndDerivs, double dx);
		double EvaluateCubicSpline1stDerivative(double x, Eigen::VectorXd& yvals, Eigen::VectorXd& y2ndDerivs, double dx);
		EamSetflInteractions(std::string setflFileName, double numericalStep=1e-4, int dim=-1);						// constructor private because of Factory method Create
	public:
		EamSetflInteractions();
		static std::unique_ptr<Interactions> Create(std::string setflFileName, double numericalStep=1e-4, int dim=-1);
		double EvaluatePotentialEnergy(RectangularSimulationBox&);
		Eigen::Vector3d EvaluateForceOnParticle(RectangularSimulationBox& rectangularSimulationBox, int pid);
		double GetMass();
};

