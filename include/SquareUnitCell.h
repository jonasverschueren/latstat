// square lattice unit cell
#pragma once

#include<vector>
#include<memory>
#include"RectangularUnitCell.h"

class SquareUnitCell: public RectangularUnitCell{
	private:
		SquareUnitCell(double latConst);
	public:
		SquareUnitCell();
		static std::unique_ptr<RectangularUnitCell> Create(double latConst);
};
