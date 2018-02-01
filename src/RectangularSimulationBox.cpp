#include"RectangularSimulationBox.h"

RectangularSimulationBox::RectangularSimulationBox(){
	this->boxLengths=std::vector<double>(3);
	this->simulationParticles=std::vector<Eigen::Vector3d>(0);
}

Eigen::Vector3d RectangularSimulationBox::CalculatePeriodicDistanceVector(Eigen::Vector3d v1, Eigen::Vector3d v2){
	Eigen::Vector3d distanceVector=v1-v2;
	for(int i=0; i<3; i++){
		if (distanceVector[i] > this->boxLengths[i]/2.){
			distanceVector[i]-=this->boxLengths[i];
		}
		else if(distanceVector[i] < -this->boxLengths[i]/2.){
			distanceVector[i]+=this->boxLengths[i];
		}
	}
	return distanceVector;
}

double RectangularSimulationBox::CalculatePeriodicDistanceMagnitude(Eigen::Vector3d v1, Eigen::Vector3d v2){
	return this->CalculatePeriodicDistanceVector(v1, v2).norm();
}

std::vector<Eigen::Vector3d>& RectangularSimulationBox::GetSimulationParticlesReference(){
	return this->simulationParticles;
}

void RectangularSimulationBox::GenerateSimulationBox(std::shared_ptr<RectangularUnitCell> rectangularUnitCell, double minimumLength, int dim){
	double latConst=rectangularUnitCell->GetLatticeConstant();
	Eigen::Vector3d minimumLengths(minimumLength, minimumLength, minimumLength);
	if(dim != -1){
		minimumLengths[dim]=0;
	}
	std::vector<int> numberOfUnitCellsPerDirection=this->CalculateNumberUnitCellsPerDirection(rectangularUnitCell, minimumLengths);
	this->SetBoxLengths(rectangularUnitCell, numberOfUnitCellsPerDirection);
	return this -> GenerateSimulationBox(rectangularUnitCell, numberOfUnitCellsPerDirection);
}

void RectangularSimulationBox::GenerateSimulationBox(std::shared_ptr<RectangularUnitCell> rectangularUnitCell, Eigen::Vector3d minimumLengths){	
	double latConst=rectangularUnitCell->GetLatticeConstant();
	std::vector<int> numberOfUnitCellsPerDirection=this->CalculateNumberUnitCellsPerDirection(rectangularUnitCell, minimumLengths);
	this->SetBoxLengths(rectangularUnitCell, numberOfUnitCellsPerDirection);
	return this -> GenerateSimulationBox(rectangularUnitCell, numberOfUnitCellsPerDirection);
}

std::vector<int> RectangularSimulationBox::CalculateNumberUnitCellsPerDirection(std::shared_ptr<RectangularUnitCell> rectangularUnitCell, Eigen::Vector3d minimumLengths){
	std::vector<int> lengths(3);
	for(int dim=0; dim<=2; dim++){
		lengths.at(dim)=int(minimumLengths[dim]/rectangularUnitCell->GetCellVectors().at(dim).norm())+1;
		if(minimumLengths[dim] == 0){
			lengths.at(dim)=0;
		}
	}
	return lengths;
}

void RectangularSimulationBox::SetBoxLengths(std::vector<double> readBoxLengths){
	this->boxLengths=readBoxLengths;
}

void RectangularSimulationBox::SetBoxLengths(std::shared_ptr<RectangularUnitCell> rectangularUnitCell, std::vector<int> numberOfUnitCellsPerDirection){
	for(int i=0; i<=2; i++){
		boxLengths[i]=numberOfUnitCellsPerDirection[i]*(rectangularUnitCell->GetCellVectors()).at(i).norm();
	}
}

void RectangularSimulationBox::GenerateSimulationBox(std::shared_ptr<RectangularUnitCell> rectangularUnitCell, std::vector<int> numberOfUnitCellsPerDirection){
	this->simulationParticles=rectangularUnitCell->GetGenerationParticles();
	for (int dim=0; dim<=2; dim++){
		std::vector<Eigen::Vector3d> existingSimulationParticles=this->simulationParticles;
		for(int i=1; i<numberOfUnitCellsPerDirection[dim]; i++){
			for(auto generationParticle:existingSimulationParticles){
				this->simulationParticles.push_back(generationParticle+rectangularUnitCell->GetCellVectors()[dim]*i);
			}
		}
	}
}

std::vector<double> RectangularSimulationBox::GetBoxLengths(){
	return this->boxLengths;
}

int RectangularSimulationBox::GetNumberOfSimulationParticles(){
	return this->simulationParticles.size();
}

std::vector<std::vector<int>>& RectangularSimulationBox::GetNeighbourListsReference(){
	return this->neighbourLists;
}

void RectangularSimulationBox::EvaluateNeighbourListsIfNecessary(double cutoff){
	if (this->GetNeighbourListsReference().size()==0){
		this->EvaluateNeighbourLists(cutoff);
	}
}

void RectangularSimulationBox::EvaluateNeighbourLists(double cutoff){
	this->neighbourLists.resize(simulationParticles.size());
	int i=0;
	for (Eigen::Vector3d p1: this->simulationParticles){
		int j=0;
		for (Eigen::Vector3d p2: this->simulationParticles){
			double distance_p12=this->CalculatePeriodicDistanceMagnitude(p1, p2);
			if(distance_p12< cutoff and i>j){
				this->neighbourLists.at(i).push_back(j);
				this->neighbourLists.at(j).push_back(i);
			}
			j++;
		}
		i++;
	}
}

void RectangularSimulationBox::MoveParticleByAmount(int pid, int dim, double amount){
	this->simulationParticles.at(pid)(dim)+=amount;
}
