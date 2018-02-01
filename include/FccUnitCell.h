// BCC unit cell
#pragma once

#include<vector>
#include<memory>
#include"RectangularUnitCell.h"

class FccUnitCell: public RectangularUnitCell{
	private:
		FccUnitCell(double latConst, bool rot=false);
	public:
		FccUnitCell();
		static std::unique_ptr<RectangularUnitCell> Create(double latConst, bool rot=false);
		void RotateFccUnitCell();
		void RotateFccUnitCell2();
};
