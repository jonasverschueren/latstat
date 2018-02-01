#include"KoizumiInteractions.h"

KoizumiInteractions::KoizumiInteractions():Interactions(){
	;
}

KoizumiInteractions::KoizumiInteractions(double forceConst, double eqmDistance, double mass, double numericalStep): Interactions(eqmDistance+.1, numericalStep, 2), forceConst(forceConst), eqmSep(eqmDistance){
	this->mass = mass;
}

std::unique_ptr<Interactions> KoizumiInteractions::Create(double forceConst, double eqmDistance, double mass, double numericalStep){
	return std::unique_ptr<Interactions>(new KoizumiInteractions(forceConst, eqmDistance, mass, numericalStep));
}

Eigen::MatrixXd KoizumiInteractions::EvaluateForceConstants(RectangularSimulationBox& simulationBox, int pid1){
	int numParts=simulationBox.GetNumberOfSimulationParticles();
	Eigen::MatrixXd forceConsts=Eigen::MatrixXd::Zero(3, 3*numParts);
	for (int dim1=0; dim1<=2; dim1++){
		for(int pid2=1; pid2<numParts; pid2++){
			Eigen::Vector3d positionDifference = simulationBox.CalculatePeriodicDistanceVector(simulationBox.GetSimulationParticlesReference().at(pid1), simulationBox.GetSimulationParticlesReference().at(pid2));
			if(positionDifference.norm() < GetInteractionCutoff()){
				if(dim1 == GetDim()){
					int dim2 = dim1;
					forceConsts(dim1, dim2)+=4*M_PI*M_PI*forceConst/(eqmSep*eqmSep);						// pid2=0
					forceConsts(dim1, 3*pid2+dim2)-=4*M_PI*M_PI*forceConst/(eqmSep*eqmSep);
				}
			}
		}
	}
	return forceConsts;
}

double KoizumiInteractions::EvaluatePotentialEnergy(RectangularSimulationBox& rectangularSimulationBox){
	double potentialEnergy=0;
	rectangularSimulationBox.EvaluateNeighbourListsIfNecessary(this->GetInteractionCutoff());
	std::vector<Eigen::Vector3d> simulationParticles=rectangularSimulationBox.GetSimulationParticlesReference();
	for(int pid1=0; pid1<simulationParticles.size(); pid1++){
		Eigen::Vector3d p1=simulationParticles.at(pid1);
		for(int pid2: rectangularSimulationBox.GetNeighbourListsReference().at(pid1)){
			double u_ij=simulationParticles.at(pid2)[this->GetDim()] - p1[this->GetDim()];
			potentialEnergy+= forceConst*(1-cos(2*M_PI*u_ij/eqmSep));
		}
	}
	return potentialEnergy;
}

Eigen::Vector3d KoizumiInteractions::EvaluateForceOnParticle(RectangularSimulationBox& rectangularSimulationBox, int pid){
	Eigen::Vector3d force=Eigen::Vector3d::Zero();
	rectangularSimulationBox.EvaluateNeighbourListsIfNecessary(this->GetInteractionCutoff()+.1);
	std::vector<Eigen::Vector3d> simulationParticles=rectangularSimulationBox.GetSimulationParticlesReference();
	Eigen::Vector3d p1=simulationParticles.at(pid);
	for(int pid2: rectangularSimulationBox.GetNeighbourListsReference().at(pid)){
			double u_ij=simulationParticles.at(pid2)[this->GetDim()] - p1[this->GetDim()];
			force[this->GetDim()] += 4*M_PI*forceConst/eqmSep*sin(2*M_PI*u_ij/eqmSep);
	}
	return force;
}

double KoizumiInteractions::GetForceConst(){
	return forceConst;
}

double KoizumiInteractions::GetEqmSep(){
	return eqmSep;
}
