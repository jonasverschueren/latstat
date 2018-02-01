#include"RectangularUnitCell.h"

RectangularUnitCell::RectangularUnitCell(){
	;
}

RectangularUnitCell::RectangularUnitCell(double readLatConst): latConst(readLatConst){
	this->primitiveVectors=std::vector<Eigen::Vector3d>(3);
}

void RectangularUnitCell::SetLatticeConstant(double readLatConst){
	this->latConst=readLatConst;
}

double RectangularUnitCell::GetLatticeConstant(){
	return this->latConst;
}

std::vector<Eigen::Vector3d> RectangularUnitCell::GetGenerationParticles(){
	return this->generationParticles;
}

std::vector<Eigen::Vector3d> RectangularUnitCell::GetPrimitiveVectors(){
	return (this->primitiveVectors);
}

std::vector<Eigen::Vector3d> RectangularUnitCell::GetCellVectors(){
	return cellVectors;
}

void RectangularUnitCell::GenerateReciprocalLatticeVectors(){
	double primitiveCellVolume=this->CalculatePrimitiveCellVolume();
	this->reciprocalLatticeVectors.resize(3);
	this->reciprocalLatticeVectors.at(0)=2*M_PI/primitiveCellVolume*(this->primitiveVectors.at(1).cross(this->primitiveVectors.at(2)));
	this->reciprocalLatticeVectors.at(1)=2*M_PI/primitiveCellVolume*this->primitiveVectors.at(2).cross(this->primitiveVectors.at(0));
	this->reciprocalLatticeVectors.at(2)=2*M_PI/primitiveCellVolume*this->primitiveVectors.at(0).cross(this->primitiveVectors.at(1));
}

double RectangularUnitCell::CalculatePrimitiveCellVolume(){
	return primitiveVectors.at(0).dot(primitiveVectors.at(1).cross(primitiveVectors.at(2)));
}

double RectangularUnitCell::Calculate1stBZVolume(){
	return reciprocalLatticeVectors.at(0).dot(reciprocalLatticeVectors.at(1).cross(reciprocalLatticeVectors.at(2)));
}

std::vector<Eigen::Vector3d>& RectangularUnitCell::GetReciprocalLatticeVectorsReference(){
	return (this->reciprocalLatticeVectors);
}
