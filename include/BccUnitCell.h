// BCC unit cell
#pragma once

#include<vector>
#include<memory>
#include"RectangularUnitCell.h"

class BccUnitCell: public RectangularUnitCell{
	private:
		BccUnitCell(double latConst, bool rot=false);
	public:
		BccUnitCell();
		static std::unique_ptr<RectangularUnitCell> Create(double latConst, bool rot=false);
		void RotateBccUnitCell();
};
