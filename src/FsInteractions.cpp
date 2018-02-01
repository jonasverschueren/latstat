#include"FsInteractions.h"

FsInteractions::FsInteractions():Interactions(){
	;
}

std::unique_ptr<Interactions> FsInteractions::Create(std::map<std::string, double>& interactionParams, double numericalStep, int dim){
	return std::unique_ptr<Interactions>(new FsInteractions(interactionParams, numericalStep, dim));
}

FsInteractions::FsInteractions(std::map<std::string, double> interactionParameters, double numericalStep, int dim): Interactions(interactionParameters["d"], numericalStep, dim){
	this->FsParameters=interactionParameters;
}

std::unique_ptr<Interactions> FsInteractions::Create(const std::string& paramFileName, double numericalStep, int dim){
	return std::unique_ptr<Interactions>(new FsInteractions(paramFileName, numericalStep, dim));
}

FsInteractions::FsInteractions(const std::string& paramFileName, double numericalStep, int dim){
	this->ReadFsParamsFile(paramFileName);
	this->SetCutoff(this->FsParameters["d"]);
	this->SetNumericalStep(numericalStep);
	this->SetDim(dim);
}

void FsInteractions::ReadFsParamsFile(const std::string& paramsFileName){
	this->mass = 0;
	std::map<std::string, double> fsParams={	//FS parameters (M.W. Finnis and J.E. Sinclair, Phil. Mag. A 50, 45 (1984))
		{"d",0},			
		{"A",0},
		{"beta",0},
		{"c",0},
		{"c0",0},
		{"c1",0},
		{"c2",0},
		{"B",0},											//Ackland and Thetford corrections	
		{"alpha",0},
		{"b0",0}
	};
	int numLinesIgnore=2;
	std::ifstream inputFileStream(paramsFileName);
	for(int i=0; i<numLinesIgnore; i++)
		inputFileStream.ignore(1000, '\n');
	// read in FS parameters
	std::string line, readString, paramName;
	while(std::getline(inputFileStream, line)){
		std::istringstream lineStream(line);
		std::getline(lineStream, paramName, '\t' );
		std::getline(lineStream, readString, '\t' );
		if (paramName == "m")
			this->mass = std::stod(readString);
		else
			fsParams[paramName] = std::stod(readString);
	}
	// check that all values have been set
	// not for beta since it is only a correction used for Iron and Chromium
	if (this->mass == 0)
		std::cout<<"WARNING: FS parameter mass has not been set."<<std::endl;
	for (auto const& el : fsParams)
		if (el.first != "beta" && el.second == 0)	
			std::cout<<"WARNING: FS parameter "<<el.first<<" has not been set."<<std::endl;
	inputFileStream.close();
	this->FsParameters = fsParams;
}

double FsInteractions::EvaluatePotentialEnergy(RectangularSimulationBox& rectangularSimulationBox){
	double potentialEnergy=0;
	rectangularSimulationBox.EvaluateNeighbourListsIfNecessary(this->GetInteractionCutoff()+1);
	std::vector<Eigen::Vector3d> simulationParticles=rectangularSimulationBox.GetSimulationParticlesReference();
	std::map<std::string, double>& params=this->FsParameters;
	for(int pid1=0; pid1<simulationParticles.size(); pid1++){
		Eigen::Vector3d p1=simulationParticles.at(pid1);
		double dens=0, quart=0, ATcorr=0;
		for(int pid2: rectangularSimulationBox.GetNeighbourListsReference().at(pid1)){
			Eigen::Vector3d p2=simulationParticles.at(pid2);
			double p12=rectangularSimulationBox.CalculatePeriodicDistanceMagnitude(p1, p2);
			if(p12 < params["d"]){
				dens+=(p12-params["d"])*(p12-params["d"]) + params["beta"]/params["d"]*(p12-params["d"])*(p12-params["d"])*(p12-params["d"]);
				if (p12 < params["c"]){
					quart+=(p12-params["c"])*(p12-params["c"])*(params["c0"]+params["c1"]*p12+params["c2"]*p12*p12);
				}
				if (p12 > 0.4*params["b0"] and p12 < params["b0"]){
					ATcorr+=(params["b0"]-p12)*(params["b0"]-p12)*(params["b0"]-p12)*exp(-params["alpha"]*p12);
				}
			}
		}
		potentialEnergy+=-params["A"]*sqrt(dens)+.5*(quart+params["B"]*ATcorr);
	}
	return potentialEnergy;

}

Eigen::Vector3d FsInteractions::EvaluateForceOnParticle(RectangularSimulationBox& rectangularSimulationBox, int pid){
	Eigen::Vector3d force=Eigen::Vector3d::Zero();
	rectangularSimulationBox.EvaluateNeighbourListsIfNecessary(this->GetInteractionCutoff()+.1);
	std::vector<Eigen::Vector3d> simulationParticles=rectangularSimulationBox.GetSimulationParticlesReference();
	std::map<std::string, double>& params=this->FsParameters;
	Eigen::Vector3d p1=simulationParticles.at(pid);
	double dens=0;
	for(int pid2 : rectangularSimulationBox.GetNeighbourListsReference().at(pid)){
		Eigen::Vector3d p2=simulationParticles.at(pid2);
		double p12=rectangularSimulationBox.CalculatePeriodicDistanceMagnitude(p1, p2);
		if(p12 < params["d"]){
			dens+=(p12-params["d"])*(p12-params["d"]) + params["beta"]/params["d"]*(p12-params["d"])*(p12-params["d"])*(p12-params["d"]);
		}
	}
	for (int pid2 : rectangularSimulationBox.GetNeighbourListsReference().at(pid)){
		Eigen::Vector3d p2=simulationParticles.at(pid2);
		double p12=rectangularSimulationBox.CalculatePeriodicDistanceMagnitude(p1, p2);
		Eigen::Vector3d p12_diff=rectangularSimulationBox.CalculatePeriodicDistanceVector(p1, p2);
		double dens_num=0;
		double invdens2=0;
		double pair=0;
		double ATcorr=0;
		if (p12 < params["d"]){
			dens_num=p12-params["d"] + 1.5*params["beta"]/params["d"]*(p12-params["d"])*(p12-params["d"]);
			double dens2=0;
			for (int pid3: rectangularSimulationBox.GetNeighbourListsReference().at(pid2)){
				Eigen::Vector3d p3=simulationParticles.at(pid3);
				double p23=rectangularSimulationBox.CalculatePeriodicDistanceMagnitude(p2, p3);
				if(p23 < params["d"]){
					dens2+=(p23-params["d"])*(p23-params["d"]) + params["beta"]/params["d"]*(p23-params["d"])*(p23-params["d"])*(p23-params["d"]);
				}
			}
			invdens2=1./sqrt(dens2);
		}
		if (p12 < params["c"]){
			pair=2*(p12-params["c"])*(params["c0"]+params["c1"]*p12+params["c2"]*p12*p12)+(p12-params["c"])*(p12-params["c"])*(params["c1"]+2*params["c2"]*p12);
		}
		if (p12>0.4*params["b0"] and p12<params["b0"]){
			ATcorr=(-3*(params["b0"]-p12)*(params["b0"]-p12)-params["alpha"]*(params["b0"]-p12)*(params["b0"]-p12)*(params["b0"]-p12))*exp(-params["alpha"]*p12);
		}
		double fmag = (-params["A"]*dens_num*(1./sqrt(dens)+invdens2)+pair+params["B"]*ATcorr) * 1./p12;
		force-=fmag*p12_diff;
	}
	return force;
}

double FsInteractions::GetMass(){
	return mass*1.6605/1.602*1e-4;					// to get time units of ps
}
