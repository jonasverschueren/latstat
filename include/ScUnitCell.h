// simple cubic unit cell class
#pragma once

#include<vector>
#include<memory>
#include"RectangularUnitCell.h"

class ScUnitCell: public RectangularUnitCell{
	private:
		ScUnitCell(double latConst);
	public:
		ScUnitCell();
		static std::unique_ptr<RectangularUnitCell> Create(double latConst);
};
