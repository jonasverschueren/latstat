#include"ScUnitCell.h"

ScUnitCell::ScUnitCell(){
	;
}

ScUnitCell::ScUnitCell(double latConst): RectangularUnitCell(latConst){
	this->primitiveVectors={Eigen::Vector3d(latConst,0,0),Eigen::Vector3d(0,latConst,0), Eigen::Vector3d(0,0,latConst)};
	this ->generationParticles={Eigen::Vector3d(0,0,0)};
	this->cellVectors=primitiveVectors;
	this->GenerateReciprocalLatticeVectors();
}

std::unique_ptr<RectangularUnitCell> ScUnitCell::Create(double latConst){
	return std::unique_ptr<RectangularUnitCell>(new ScUnitCell(latConst));
}
