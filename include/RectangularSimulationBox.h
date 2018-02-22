#pragma once

#include<vector>
#include<memory>
#include"Eigen/Dense"
#include"RectangularUnitCell.h"

class RectangularSimulationBox{
	private:
		std::vector<double> boxLengths;
		std::vector<Eigen::Vector3d> simulationParticles;
		std::vector<std::vector<int>> neighbourLists;
		std::vector<int> CalculateNumberUnitCellsPerDirection(std::shared_ptr<RectangularUnitCell> rectangularUnitCell, Eigen::Vector3d minimumLengths);
		void GenerateSimulationBox(std::shared_ptr<RectangularUnitCell>, std::vector<int> numberOfUnitCellsPerDirection);
		void SetBoxLengths(std::vector<double> boxLengths);
		void SetBoxLengths(std::shared_ptr<RectangularUnitCell>, std::vector<int> numberOfUnitCellsPerDirection);
	public:
		RectangularSimulationBox();
		Eigen::Vector3d CalculatePeriodicDistanceVector(Eigen::Vector3d v1, Eigen::Vector3d v2);		// v1-v2
		double CalculatePeriodicDistanceMagnitude(Eigen::Vector3d v1, Eigen::Vector3d v2 );
		std::vector<Eigen::Vector3d>& GetSimulationParticlesReference();
		void GenerateSimulationBox(std::shared_ptr<RectangularUnitCell>, double minimumLength, int dim=-1);
		void GenerateSimulationBox(std::shared_ptr<RectangularUnitCell>, Eigen::Vector3d minimumLengths);
		std::vector<double> GetBoxLengths();
		int GetNumberOfSimulationParticles();
		void EvaluateNeighbourListsIfNecessary(double cutoff);
		void EvaluateNeighbourLists(double cutoff);
		std::vector<std::vector<int>>& GetNeighbourListsReference();
		void MoveParticleByAmount(int pid, int dim, double amount);
};
