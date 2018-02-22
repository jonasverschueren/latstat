// Loosely speaking: PerfectCrystal = RectangularUnitCell + Interactions
#pragma once

#include<vector>
#include<complex>
#include<memory>
#include"Eigen/Dense"
#include"Eigen/Eigenvalues"
#include"Crystal.h"
#include"RectangularSimulationBox.h"
#include"RectangularUnitCell.h"
#include"Interactions.h"

class PerfectCrystal : public Crystal{
	private:
		Eigen::MatrixXd forceConstants;
		void EvaluateDynamicalMatrix(Eigen::Matrix3cd& DynMat_k, Eigen::Vector3d k);
		void EvaluateForceConstantsIfNecessary();
		void EvaluateForceConstants();
	public:
		PerfectCrystal();
		PerfectCrystal(std::shared_ptr<RectangularUnitCell>, std::shared_ptr<Interactions>);
		std::shared_ptr<RectangularUnitCell> ptr_rectangularUnitCell;
		RectangularSimulationBox forceConstantsSimulationBox;
		double EvaluatePotentialEnergyPerAtom();
		Eigen::Vector3d EvaluatePhononFrequencySquared(Eigen::Vector3d k);
		Eigen::Matrix3d EvaluatePhononPolarisations(Eigen::Vector3d k);
		Eigen::Vector3d EvaluateKpointFromFractional(Eigen::Vector3d& kfrac_vect);
		double CheckIfAsrIsObeyed();
};
