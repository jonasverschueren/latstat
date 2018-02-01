#include "Interactions.h"

Interactions::Interactions(){
	;
}

Interactions::Interactions(double interactionCutoff, double numericalStep, int dim): interactionCutoff(interactionCutoff), dim(dim), numericalStep(numericalStep){
	;
}

Interactions::Interactions(double numericalStep, int dim): dim(dim), numericalStep(numericalStep){
	;
}

Eigen::Vector3d Interactions::EvaluateForceOnParticleNumerically(RectangularSimulationBox& rectangularSimulationBox, int pid){
	Eigen::Vector3d force=Eigen::Vector3d::Zero();
	for(int dim=0; dim<=2; dim++){
		rectangularSimulationBox.MoveParticleByAmount(pid, dim, numericalStep);
		double vPosStep=EvaluatePotentialEnergy(rectangularSimulationBox);
		rectangularSimulationBox.MoveParticleByAmount(pid, dim, -2*numericalStep);
		double vNegStep=EvaluatePotentialEnergy(rectangularSimulationBox);
		rectangularSimulationBox.MoveParticleByAmount(pid, dim, numericalStep);
		force(dim)=(vNegStep-vPosStep)/(2.*numericalStep);
	}
	return force;

}

Eigen::MatrixXd Interactions::EvaluateForceConstants(RectangularSimulationBox& rectangularSimulationBox, int pid1){
	int numParts=rectangularSimulationBox.GetNumberOfSimulationParticles();
	Eigen::MatrixXd forceConsts=Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>::Zero(3, 3*numParts);
	for (int dim1=0; dim1<=2; dim1++){
		for(int pid2=0; pid2<numParts; pid2++){
			for(int dim2=0; dim2<=2; dim2++){
				if(dim==-1 or dim1==dim and dim2==dim){
					rectangularSimulationBox.MoveParticleByAmount(pid1, dim1, this->numericalStep);
					double F2_p1PosStep=this->EvaluateForceOnParticle(rectangularSimulationBox, pid2)(dim2);
					rectangularSimulationBox.MoveParticleByAmount(pid1, dim1, -2*this->numericalStep);
					double F2_p1NegStep=this->EvaluateForceOnParticle(rectangularSimulationBox, pid2)(dim2);
					rectangularSimulationBox.MoveParticleByAmount(pid1, dim1, this->numericalStep);
					forceConsts(dim1, 3*pid2+dim2)=(F2_p1NegStep-F2_p1PosStep)/(2.*this->numericalStep);
				}
			}
		}
	}
	return forceConsts;
}

double Interactions::GetInteractionCutoff(){
	return interactionCutoff;
}

void Interactions::SetNumericalStep(double step){
	numericalStep=step;
}

int Interactions::GetDim(){
	return dim;
}

void Interactions::SetCutoff(double cutoff){
	interactionCutoff=cutoff;
}

void Interactions::SetDim(int new_dim){
	dim = new_dim;
}

double Interactions::GetMass(){
	return mass;
}
