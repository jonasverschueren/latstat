#include<memory>
#include<fstream>
#include<cmath>
#include"EamSetflInteractions.h"
#include"FccUnitCell.h"
#include"PerfectCrystal.h"
#include "Eigen/Dense"

int main(){
	double latConst=4.05;
	std::string setflFile="Al99.eam.alloy";
	std::shared_ptr<Interactions> interactions = std::shared_ptr<Interactions>(EamSetflInteractions::Create(setflFile, 1e-5).release());
	std::shared_ptr<RectangularUnitCell> rotatedUnitCell = std::shared_ptr<RectangularUnitCell>(FccUnitCell::Create(latConst, false).release());
	std::unique_ptr<PerfectCrystal> aluCrystal=std::unique_ptr<PerfectCrystal>(new PerfectCrystal(rotatedUnitCell, interactions));
	
	Eigen::Vector3d kfrac(0, 1, 0);
	std::ofstream ex_phonon_file;
	ex_phonon_file.open("ex_phonon_file.txt");
	for (double kmag = 0; kmag <= 2; kmag+=.01){
		Eigen::Vector3d freqs2 = aluCrystal -> EvaluatePhononFrequencySquared(kmag*kfrac);
		for (int i = 0; i <=2; i++)
			ex_phonon_file<<kmag<<"\t"<<sqrt(freqs2[i])/(2*M_PI)<<std::endl;
	}
	ex_phonon_file.close();
	return 0;
}
