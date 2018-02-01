#include "EamSetflInteractions.h"

EamSetflInteractions::EamSetflInteractions(){
	;
}

EamSetflInteractions::EamSetflInteractions(std::string setflFileName, double numericalStep, int dim):Interactions(numericalStep, dim){
	ReadSetflFile(setflFileName);
	Evaluate2ndDerivatives();
	SetCutoff(cutoff);
}

std::unique_ptr<Interactions> EamSetflInteractions::Create(std::string setflFileName, double numericalStep, int dim){
	return std::unique_ptr<Interactions>(new EamSetflInteractions(setflFileName, numericalStep, dim));
}

void EamSetflInteractions::ReadSetflFile(std::string fileName){
	ReadSetflFileHeader(fileName);
	std::ifstream inputFileStream(fileName);
	// Ignore the first 6 lines
	int numLinesIgnore=6;
	for(int i=0; i<numLinesIgnore; i++){
		inputFileStream.ignore(1000, '\n');
	}
	std::string line;
	// Read in embedding values for a given density rho
	embeddingValues=Eigen::VectorXd(Nrho);
	for(int i=0; i<Nrho; i++){
		std::getline(inputFileStream, line);
		embeddingValues(i)=std::stod(line);
	}
	// Read in density values for a given displacement magnitude r
	densityValues=Eigen::VectorXd(Nr);
	for(int i=0; i<Nr; i++){
		std::getline(inputFileStream, line);
		densityValues(i)=std::stod(line);
	}
	// Read in pair interaction values for a given displacement magnitude r
	pairValues=Eigen::VectorXd(Nr);
	for(int i=0; i<Nr; i++){
		std::getline(inputFileStream, line);
		pairValues(i)=std::stod(line);										
	}
	inputFileStream.close();
}

void EamSetflInteractions::ReadSetflFileHeader(std::string fileName){
	std::ifstream inputFileStream(fileName);
	if(inputFileStream.is_open()==false){
		std::cout<<"No Potential file read in, check the file name."<<std::endl;
		exit(EXIT_FAILURE);
	}
	// Ignore the first 4 lines
	int numLinesIgnore=4;
	for(int i=0; i<numLinesIgnore; i++){
		inputFileStream.ignore(1000, '\n');
	}
	// read in line 5 and get Nrho, drho, Nr, dr, cutoff
	std::string line, readString;
	std::getline(inputFileStream, line);
	std::istringstream lineStream(line);
	std::getline(lineStream, readString, '\t' );
	Nrho=std::stoi(readString);
	std::getline(lineStream, readString, '\t');
	drho=std::stod(readString);
	std::getline(lineStream, readString, '\t');
	Nr=std::stoi(readString);
	std::getline(lineStream, readString, '\t');
	dr=std::stod(readString);
	std::getline(lineStream, readString, '\t');
	cutoff=std::stod(readString);
	// read in line 6 and get mass 
	std::getline(inputFileStream, line);
	lineStream=std::istringstream(line);
	std::getline(lineStream, readString, '\t' );
	std::getline(lineStream, readString, '\t' );
	mass=std::stod(readString);
	if (this->mass == 0)
		std::cout<<"WARNING: EAM parameter mass has not been set."<<std::endl;
	inputFileStream.close();
}

// Evaluate 2nd derivatives necessary for cubic spline interpolation
void EamSetflInteractions::Evaluate2ndDerivatives(){
	y2embeddingValues=EvaluateSpline2ndDerivative(embeddingValues, drho);			// 2nd derivative vectors
	y2densityValues=EvaluateSpline2ndDerivative(densityValues, dr);
	y2pairValues=EvaluateSpline2ndDerivative(pairValues, dr);
}

// Solve A*x=y for A tridiagonal Matrix consisting of diagonals a, b and c
// Algorithm taken from wikipedia
void ThomasAlgorithm(Eigen::VectorXd& a, Eigen::VectorXd& b, Eigen::VectorXd& c, Eigen::VectorXd& y, Eigen::VectorXd& x){
	int numvals=b.size();
	// forward step
	Eigen::VectorXd c2(numvals-1), y2(numvals);
	c2(0)=c(0)/b(0);
	y2(0)=y(0)/b(0);
	for(int i=1; i<numvals-1; i++){
		c2(i)=c(i)/(b(i)-a(i-1)*c2(i-1));
		y2(i)=(y(i)-a(i-1)*y2(i-1))/(b(i)-a(i-1)*c2(i-1));
	}
	y2(numvals-1)=(y(numvals-1)-a(numvals-2)*y2(numvals-2))/(b(numvals-1)-a(numvals-2)*c2(numvals-2));
	// back substitution
	x(numvals-1)=y2(numvals-1);
	for(int i=numvals-2; i>0; i--){
		x(i)=y2(i)-c2(i)*x(i+1);
	}
}


Eigen::VectorXd EamSetflInteractions::EvaluateSpline2ndDerivative(Eigen::VectorXd& yvals, double dx){
	int numvals=yvals.size();								
	Eigen::VectorXd a(numvals-1), b(numvals), c(numvals-1);			// coefficients of tridiagonal Matrix
	Eigen::VectorXd F=Eigen::VectorXd::Zero(numvals);				// right hand-sides
	b(0)=2.;														// values set to get natural splines
	b(numvals-1)=b(0);
	c(0)=0.;
	a(numvals-2)=0.;
	F(0)=0;
	F(numvals-1)=0.;
	for(int i=1; i<numvals-1; i++){
		a(i-1)=dx/6.;
		c(i)=dx/6.;
		b(i)=2.*dx/3.;
		F(i)=(yvals(i+1)-yvals(i))/dx-(yvals(i)-yvals(i-1))/dx ;
	}
	Eigen::VectorXd solution(numvals);
	ThomasAlgorithm(a, b, c, F, solution);
	return solution;
}

double EamSetflInteractions::EvaluateCubicSpline(double x, Eigen::VectorXd& yvals, Eigen::VectorXd& y2ndDerivs, double dx){
	int lowerIndex=int(x/dx);					// lower index of the interval in which x lies
	double A=((lowerIndex+1)*dx-x)/dx;
	double B=1.-A;
	double C=(A*A*A-A)*dx*dx/6.;
	double D=(B*B*B-B)*dx*dx/6.;
	return A*yvals(lowerIndex)+B*yvals(lowerIndex+1)+C*y2ndDerivs(lowerIndex)+D*y2ndDerivs(lowerIndex+1);
}

double EamSetflInteractions::EvaluateCubicSpline1stDerivative(double x, Eigen::VectorXd& yvals, Eigen::VectorXd& y2ndDerivs, double dx){
	int lowerIndex=int(x/dx);					// lower index of the interval in which x lies
	double A=((lowerIndex+1)*dx-x)/dx;
	double dA=-1./dx;
	double B=1.-A;
	double dB=-dA;
	double dC=(3*A*A-1)*dx*dx*dA/6.;
	double dD=(3*B*B-1)*dx*dx*dB/6.;
	return dA*yvals(lowerIndex)+dB*yvals(lowerIndex+1)+dC*y2ndDerivs(lowerIndex)+dD*y2ndDerivs(lowerIndex+1);
}

double EamSetflInteractions::EvaluatePotentialEnergy(RectangularSimulationBox& rectangularSimulationBox){
	double potentialEnergy=0;
	rectangularSimulationBox.EvaluateNeighbourListsIfNecessary(this->GetInteractionCutoff()+1);
	std::vector<Eigen::Vector3d> simulationParticles=rectangularSimulationBox.GetSimulationParticlesReference();
	for(int pid1=0; pid1<simulationParticles.size(); pid1++){
		Eigen::Vector3d p1=simulationParticles.at(pid1);
		double dens=0, pairPot=0;
		for(int pid2: rectangularSimulationBox.GetNeighbourListsReference().at(pid1)){
			Eigen::Vector3d p2=simulationParticles.at(pid2);
			double p12=rectangularSimulationBox.CalculatePeriodicDistanceMagnitude(p1, p2);
			if(p12 < cutoff){
				dens+=EvaluateCubicSpline(p12, densityValues, y2densityValues, dr);
				pairPot+=EvaluateCubicSpline(p12, pairValues, y2pairValues, dr)/p12;						// divide by p12 because pair values are tabulated as r*pairInteractions(r)
			}
		}
		potentialEnergy+=EvaluateCubicSpline(dens, embeddingValues, y2embeddingValues, drho)+.5*pairPot;
	}
	return potentialEnergy;
}

Eigen::Vector3d EamSetflInteractions::EvaluateForceOnParticle(RectangularSimulationBox& rectangularSimulationBox, int pid){
	Eigen::Vector3d force=Eigen::Vector3d::Zero();
	rectangularSimulationBox.EvaluateNeighbourListsIfNecessary(this->GetInteractionCutoff()+.1);
	std::vector<Eigen::Vector3d> simulationParticles=rectangularSimulationBox.GetSimulationParticlesReference();
	Eigen::Vector3d p1=simulationParticles.at(pid);
	double dens=0;
	for(int pid2 : rectangularSimulationBox.GetNeighbourListsReference().at(pid)){
		Eigen::Vector3d p2=simulationParticles.at(pid2);
		double p12=rectangularSimulationBox.CalculatePeriodicDistanceMagnitude(p1, p2);
		if(p12 < cutoff){
			dens+=EvaluateCubicSpline(p12, densityValues, y2densityValues, dr);
		}
	}
	for (int pid2 : rectangularSimulationBox.GetNeighbourListsReference().at(pid)){
		Eigen::Vector3d p2=simulationParticles.at(pid2);
		double p12=rectangularSimulationBox.CalculatePeriodicDistanceMagnitude(p1, p2);
		Eigen::Vector3d p12_diff=rectangularSimulationBox.CalculatePeriodicDistanceVector(p1, p2);
		double pair=0;
		double dens2=0;
		double densDeriv=0;
		if (p12 < cutoff){
			densDeriv=EvaluateCubicSpline1stDerivative(p12, densityValues, y2densityValues, dr);
			pair=EvaluateCubicSpline1stDerivative(p12, pairValues, y2pairValues, dr)-EvaluateCubicSpline(p12, pairValues, y2pairValues, dr)/p12;
			for (int pid3: rectangularSimulationBox.GetNeighbourListsReference().at(pid2)){
				Eigen::Vector3d p3=simulationParticles.at(pid3);
				double p23=rectangularSimulationBox.CalculatePeriodicDistanceMagnitude(p2, p3);
				if(p23 < cutoff){
					dens2+=EvaluateCubicSpline(p23, densityValues, y2densityValues, dr);
				}
			}
		}
		double embeddingContribution=EvaluateCubicSpline1stDerivative(dens, embeddingValues, y2embeddingValues, drho)+EvaluateCubicSpline1stDerivative(dens2, embeddingValues, y2embeddingValues, drho);
		force-=densDeriv*embeddingContribution*p12_diff/p12+pair*p12_diff/(p12*p12);
	}
	return force;
}

double EamSetflInteractions::GetMass(){
	return mass*1.6605/1.602*1e-4;					// to get time units of ps
}
