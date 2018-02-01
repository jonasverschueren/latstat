#include"HarmonicInteractions.h"

HarmonicInteractions::HarmonicInteractions():Interactions(){
	;
}

HarmonicInteractions::HarmonicInteractions(double forceConst, double eqmDistance, double mass, double numericalStep, int dim): Interactions(eqmDistance+.1, numericalStep, dim), forceConst(forceConst), eqmSep(eqmSep){
	this->mass = mass;;
}

std::unique_ptr<Interactions> HarmonicInteractions::Create(double forceConst, double eqmDistance, double mass, double numericalStep, int dim){
	return std::unique_ptr<Interactions>(new HarmonicInteractions(forceConst, eqmDistance, mass, numericalStep, dim));
}

Eigen::MatrixXd HarmonicInteractions::EvaluateForceConstants(RectangularSimulationBox& simulationBox, int pid1){
	int numParts=simulationBox.GetNumberOfSimulationParticles();
	Eigen::MatrixXd forceConsts=Eigen::MatrixXd::Zero(3, 3*numParts);
	for (int dim1=0; dim1<=2; dim1++){
		for(int pid2=0; pid2<numParts; pid2++){
			if (pid1 != pid2){
				Eigen::Vector3d positionDifference = simulationBox.CalculatePeriodicDistanceVector(simulationBox.GetSimulationParticlesReference().at(pid1), simulationBox.GetSimulationParticlesReference().at(pid2));
				if(positionDifference.norm() < GetInteractionCutoff()){
					if(GetDim() == -1 or dim1 == GetDim()){
						int dim2 = dim1;
						forceConsts(dim1, 3*pid1+dim2)+=forceConst;						// pid2=0
						forceConsts(dim1, 3*pid2+dim2)-=forceConst;
					}
				}
			}
		}
	}
	return forceConsts;
}

double HarmonicInteractions::EvaluatePotentialEnergy(RectangularSimulationBox& rectangularSimulationBox){
	double potentialEnergy=0;
//	rectangularSimulationBox.EvaluateNeighbourListsIfNecessary(this->GetInteractionCutoff()+.25*latConst);
//	std::vector<Eigen::Vector3d> simulationParticles=rectangularSimulationBox.GetSimulationParticlesReference();
//	for(int pid1=0; pid1<simulationParticles.size(); pid1++){
//		Eigen::Vector3d p1=simulationParticles.at(pid1);
//		for(int pid2: rectangularSimulationBox.GetNeighbourListsReference().at(pid1)){
//			Eigen::Vector3d p2=simulationParticles.at(pid2);
//			double p12=rectangularSimulationBox.CalculatePeriodicDistanceMagnitude(p1, p2);
//			if(p12 < params["d"]){
//				dens+=(p12-params["d"])*(p12-params["d"]);
//				if (p12 < params["c"]){
//					quart+=(p12-params["c"])*(p12-params["c"])*(params["c0"]+params["c1"]*p12+params["c2"]*p12*p12);
//				}
//				if (p12 > 0.4*params["b0"] and p12 < params["b0"]){
//					ATcorr+=(params["b0"]-p12)*(params["b0"]-p12)*(params["b0"]-p12)*exp(-params["alpha"]*p12);
//				}
//			}
//		}
//		potentialEnergy+=-params["A"]*sqrt(dens)+.5*(quart+params["B"]*ATcorr);
//	}
	return potentialEnergy;
}

Eigen::Vector3d HarmonicInteractions::EvaluateForceOnParticle(RectangularSimulationBox& rectangularSimulationBox, int pid){
	Eigen::Vector3d force=Eigen::Vector3d::Zero();
//	rectangularSimulationBox.EvaluateNeighbourListsIfNecessary(this->GetInteractionCutoff()+.1);
//	std::vector<Eigen::Vector3d> simulationParticles=rectangularSimulationBox.GetSimulationParticlesReference();
//	std::map<std::string, double>& params=this->FsParameters;
//	Eigen::Vector3d p1=simulationParticles.at(pid);
//	double dens=0;
//	for(int pid2 : rectangularSimulationBox.GetNeighbourListsReference().at(pid)){
//		Eigen::Vector3d p2=simulationParticles.at(pid2);
//		double p12=rectangularSimulationBox.CalculatePeriodicDistanceMagnitude(p1, p2);
//		if(p12 < params["d"]){
//			dens+=(p12-params["d"])*(p12-params["d"]);
//		}
//	}
//	for (int pid2 : rectangularSimulationBox.GetNeighbourListsReference().at(pid)){
//		Eigen::Vector3d p2=simulationParticles.at(pid2);
//		double p12=rectangularSimulationBox.CalculatePeriodicDistanceMagnitude(p1, p2);
//		Eigen::Vector3d p12_diff=rectangularSimulationBox.CalculatePeriodicDistanceVector(p1, p2);
//		double dens_num=0;
//		double invdens2=0;
//		double pair=0;
//		double ATcorr=0;
//		if (p12 < params["d"]){
//			dens_num=p12-params["d"];
//			double dens2=0;
//			for (int pid3: rectangularSimulationBox.GetNeighbourListsReference().at(pid2)){
//				Eigen::Vector3d p3=simulationParticles.at(pid3);
//				double p23=rectangularSimulationBox.CalculatePeriodicDistanceMagnitude(p2, p3);
//				if(p23 < params["d"]){
//					dens2+=(p23-params["d"])*(p23-params["d"]);
//				}
//			}
//			invdens2=1./sqrt(dens2);
//		}
//		if (p12 < params["c"]){
//			pair=2*(p12-params["c"])*(params["c0"]+params["c1"]*p12+params["c2"]*p12*p12)+(p12-params["c"])*(p12-params["c"])*(params["c1"]+2*params["c2"]*p12);
//		}
//		if (p12>0.4*params["b0"] and p12<params["b0"]){
//			ATcorr=(-3*(params["b0"]-p12)*(params["b0"]-p12)-params["alpha"]*(params["b0"]-p12)*(params["b0"]-p12)*(params["b0"]-p12))*exp(-params["alpha"]*p12);
//		}
//		double fmag = (-params["A"]*dens_num*(1./sqrt(dens)+invdens2)+pair+params["B"]*ATcorr) * 1./p12;
//		force-=fmag*p12_diff;
//	}
	return force;
}

double HarmonicInteractions::GetForceConst(){
	return forceConst;
}

double HarmonicInteractions::GetEqmSep(){
	return eqmSep;
}
