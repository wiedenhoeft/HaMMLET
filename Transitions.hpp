//     Copyright 2014-2015 John Wiedenhoeft, Eric Brugel
//
//     This file is part of HaMMLET.
// 
//     HaMMLET is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.
// 
//     HaMMLET is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
// 
//     You should have received a copy of the GNU General Public License
//     along with HaMMLET.  If not, see <http://www.gnu.org/licenses/>.


#ifndef TRANSITIONS_H
#define TRANSITIONS_H

#include "Categorical.hpp"
#include <vector>
class StateSequence;

class Transitions {
public:
	std::vector<Categorical*> A;

	Transitions() {
	}

	Transitions( std::vector<Categorical*> A ) : A( A ) {
		clearStats();
	}

	void addObs( StateSequence& seq );
	void clearStats();
	void sampleParameters();
	void setNeighbors( StateSequence* seq );
	void addHistory();
	int maximumBlockSize(int maxBlockSize);
	double minimumLogSelfTransition();
	double entry(int row, int col);
	double logEntry(int row, int col);
};

#endif
