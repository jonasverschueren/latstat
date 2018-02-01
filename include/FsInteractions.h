// Class Finnis-Sinclair type interactions
#pragma once

#include<string>
#include<sstream>
#include<vector>
#include<map>
#include<cmath>
#include<fstream>
#include"Interactions.h"
#include"RectangularSimulationBox.h"

class FsInteractions: public Interactions{
	private:
		std::map<std::string, double> FsParameters;
		FsInteractions(std::map<std::string, double>FsParameters, double numericalStep=1e-4, int dim=-1);			// because of factory method create
		FsInteractions(const std::string& paramFileName, double numericalStep=1e-4, int dim=-1);					// because of factory method create
		void ReadFsParamsFile(const std::string& paramFileName);
	public:
		FsInteractions();
		static std::unique_ptr<Interactions> Create(std::map<std::string, double>&, double numericalStep=1e-4, int dim=-1);
		static std::unique_ptr<Interactions> Create(const std::string& fname, double numericalStep=1e-4, int dim=-1);
		double EvaluatePotentialEnergy(RectangularSimulationBox&);
		Eigen::Vector3d EvaluateForceOnParticle(RectangularSimulationBox& rectangularSimulationBox, int pid);
		double GetMass();
};
