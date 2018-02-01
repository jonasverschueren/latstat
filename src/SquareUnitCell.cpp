#include"SquareUnitCell.h"

SquareUnitCell::SquareUnitCell(){
	;
}

SquareUnitCell::SquareUnitCell(double latConst): RectangularUnitCell(latConst){
	this->primitiveVectors={Eigen::Vector3d(latConst,0,0), Eigen::Vector3d(0,latConst,0), Eigen::Vector3d(0,0,latConst)};
	this ->generationParticles={Eigen::Vector3d(0,0,0)};
	this->cellVectors=primitiveVectors;
	this->GenerateReciprocalLatticeVectors();
}

std::unique_ptr<RectangularUnitCell> SquareUnitCell::Create(double latConst){
	return std::unique_ptr<RectangularUnitCell>(new SquareUnitCell(latConst));
}
