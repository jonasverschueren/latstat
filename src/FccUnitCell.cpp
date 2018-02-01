#include"FccUnitCell.h"

FccUnitCell::FccUnitCell(){
	;
}

FccUnitCell::FccUnitCell(double latConst, bool rot): RectangularUnitCell(latConst){
	this->primitiveVectors={.5*latConst*Eigen::Vector3d(1,1,0),.5*latConst*Eigen::Vector3d(0,1,1),.5*latConst*Eigen::Vector3d(1,0,1)};
	this ->generationParticles={Eigen::Vector3d(0,0,0), .5*latConst*Eigen::Vector3d(1,1,0),.5*latConst*Eigen::Vector3d(0,1,1),.5*latConst*Eigen::Vector3d(1,0,1)};
	this->cellVectors={Eigen::Vector3d(latConst, 0, 0), Eigen::Vector3d(0, latConst, 0), Eigen::Vector3d(0, 0, latConst)};
	if(rot){
		RotateFccUnitCell2();
	}
	else{
		std::cout<<"Generating FCC unit cell with cartesian axes along [1 0 0], [0 1 0], [0 0 1] respectively"<<std::endl;
	}
	this->GenerateReciprocalLatticeVectors();
}

std::unique_ptr<RectangularUnitCell> FccUnitCell::Create(double latConst, bool rot){
	return std::unique_ptr<RectangularUnitCell>(new FccUnitCell(latConst, rot));
}

void FccUnitCell::RotateFccUnitCell(){
	std::cout<<"Generating FCC unit cell with cartesian axes along [-1 1 0], [1 1 0], [0 0 1] respectively"<<std::endl;
	Eigen::Matrix3d rotationMatrix;							// rotationMatrix*v rotates to rotated coordinates for v in cartesian coordinates
	rotationMatrix << -1./sqrt(2), 1./sqrt(2), 0,
			1./sqrt(2), 1./sqrt(2), 0,
			0, 0, 1; 
	for(auto& primVec : primitiveVectors){
		primVec=rotationMatrix*primVec;
	}
	cellVectors={Eigen::Vector3d(latConst/sqrt(2), 0, 0), Eigen::Vector3d(0, latConst/sqrt(2), 0), Eigen::Vector3d(0, 0, latConst)};
	generationParticles={0*cellVectors[0], .5*cellVectors[0]+.5*cellVectors[1]+.5*cellVectors[2]};
}

void FccUnitCell::RotateFccUnitCell2(){
	std::cout<<"Generating FCC unit cell with cartesian axes along [1 -2 1], [1 1 1], [-1 0 1] respectively"<<std::endl;
	Eigen::Matrix3d rotationMatrix;							// rotationMatrix*v rotates to rotated coordinates for v in cartesian coordinates
	rotationMatrix << 1./sqrt(6),  -2./sqrt(6), 1./sqrt(6),
					1./sqrt(3), 1./sqrt(3), 1./sqrt(3),
					-1./sqrt(2), 0, 1./sqrt(2);
	for(auto& primVec : primitiveVectors){
		primVec=rotationMatrix*primVec;
	}
	cellVectors={Eigen::Vector3d(.5*sqrt(6)*latConst, 0, 0), Eigen::Vector3d(0, sqrt(3)*latConst, 0), Eigen::Vector3d(0, 0, .5*sqrt(2)*latConst)};
	generationParticles={0*cellVectors[0], 1./3*cellVectors[0]+2./3*cellVectors[1], 2./3*cellVectors[0]+1./3*cellVectors[1], .5*cellVectors[2]+1./6*cellVectors[0]+1./3*cellVectors[1], .5*cellVectors[2]+.5*cellVectors[0], .5*cellVectors[2]+5./6*cellVectors[0]+2./3*cellVectors[1]};
}
