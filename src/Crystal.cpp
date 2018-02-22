#include"Crystal.h"

Crystal::Crystal(std::shared_ptr<Interactions> read_interactions): ptr_interactions(read_interactions){
	;
}

double Crystal::GetMass(){
	return ptr_interactions->GetMass();
}
