#include"PerfectCrystal.h"
#include<fstream>

PerfectCrystal::PerfectCrystal(){
	;
}

PerfectCrystal::PerfectCrystal(std::shared_ptr<RectangularUnitCell> read_rectangularUnitCell, std::shared_ptr<Interactions> read_interactions): Crystal(read_interactions), ptr_rectangularUnitCell(read_rectangularUnitCell){
	this->EvaluateForceConstantsIfNecessary();
}

double PerfectCrystal::EvaluatePotentialEnergyPerAtom(){
	RectangularSimulationBox rectangularSimulationBox;
	rectangularSimulationBox.GenerateSimulationBox(this->ptr_rectangularUnitCell, 2*ptr_interactions->GetInteractionCutoff(), ptr_interactions->GetDim());
	double potE=ptr_interactions->EvaluatePotentialEnergy(rectangularSimulationBox);
	return potE/double(rectangularSimulationBox.GetNumberOfSimulationParticles());
}

Eigen::Vector3d PerfectCrystal::EvaluatePhononFrequencySquared(Eigen::Vector3d kfrac_vect){
	Eigen::Matrix3cd DynMat_k=Eigen::Matrix3cd::Zero();
	this->EvaluateDynamicalMatrix(DynMat_k, kfrac_vect);
	Eigen::ComplexEigenSolver<Eigen::Matrix3cd> eigenValueSolver(DynMat_k);
	Eigen::Vector3cd eigenvalues = eigenValueSolver.eigenvalues();
	for (int i = 0; i <=2; i++)
		if(abs(eigenvalues.imag()[0]) > 1e-5)
			std::cout<<"imag component: "<<eigenvalues.imag()<<std::endl;
	return 1./(this->GetMass())*eigenvalues.real();
}

Eigen::Matrix3d PerfectCrystal::EvaluatePhononPolarisations(Eigen::Vector3d kfrac_vect){
	Eigen::Matrix3cd DynMat_k=Eigen::Matrix3cd::Zero();
	this->EvaluateDynamicalMatrix(DynMat_k, kfrac_vect);
	Eigen::ComplexEigenSolver<Eigen::Matrix3cd> eigenValueSolver(DynMat_k);
	Eigen::Matrix3cd polarisations = eigenValueSolver.eigenvectors();
	for (int i = 0; i <=2; i++)
		if(abs(polarisations.imag()(i,i)) > 1e-5)
			std::cout<<"imag component: "<<polarisations.imag()<<std::endl;
	return 1./(this->GetMass())*polarisations.real();
}

void PerfectCrystal::EvaluateDynamicalMatrix(Eigen::Matrix3cd& DynMat_k, Eigen::Vector3d kfrac_vect){
	this->EvaluateForceConstantsIfNecessary();
	Eigen::Vector3d k = this->EvaluateKpointFromFractional(kfrac_vect);
	Eigen::Vector3d p0 = this->forceConstantsSimulationBox.GetSimulationParticlesReference().at(0);
	for(int pid1=0; pid1 < this->forceConstantsSimulationBox.GetNumberOfSimulationParticles(); pid1++){
		Eigen::Vector3d p1 = this->forceConstantsSimulationBox.GetSimulationParticlesReference().at(pid1);
		Eigen::Vector3d rdiff = this->forceConstantsSimulationBox.CalculatePeriodicDistanceVector(p0, p1);
		double rdotk = rdiff.dot(k);
		for(int dim0=0; dim0<=2; dim0++){
			for(int dim1=0; dim1<=2; dim1++){
				DynMat_k(dim0, dim1)-= 2*this->forceConstants(dim0, 3*pid1+dim1)*sin(.5*rdotk)*sin(.5*rdotk);
			}
		}
	}
}

Eigen::Vector3d PerfectCrystal::EvaluateKpointFromFractional(Eigen::Vector3d& kfrac_vect){
	Eigen::Vector3d k(0,0,0);
	for(int i=0; i<=2; i++){
		k += kfrac_vect(i)*(this->ptr_rectangularUnitCell->GetReciprocalLatticeVectorsReference()).at(i);
	}
	return k;
}

void PerfectCrystal::EvaluateForceConstantsIfNecessary(){
	if(this->forceConstantsSimulationBox.GetNumberOfSimulationParticles()==0 or this->forceConstants.cols()!=3*this->forceConstantsSimulationBox.GetNumberOfSimulationParticles()){
		this->EvaluateForceConstants();
	}
}

void PerfectCrystal::EvaluateForceConstants(){
	this->forceConstantsSimulationBox.GenerateSimulationBox(this->ptr_rectangularUnitCell, 6*ptr_interactions->GetInteractionCutoff(), ptr_interactions->GetDim());
	this-> forceConstants=ptr_interactions->EvaluateForceConstants(this->forceConstantsSimulationBox, 0);
}

double PerfectCrystal::CheckIfAsrIsObeyed(){
	double sum = 0;
	bool allZero = true;
	this->EvaluateForceConstantsIfNecessary();
	for(int i = 0; i < forceConstants.rows(); i++){
		for(int j = 0; j < forceConstants.cols(); j++){
			sum += forceConstants(i, j);
			if(allZero == true && forceConstants(i, j) !=0)
				allZero = false;
		}
	}
	if (allZero)
		std::cout<<"WARNING: All force constants are 0."<<std::endl;
	return sum;
}
