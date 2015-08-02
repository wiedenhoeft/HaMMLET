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


#include "Initial.hpp"
#include "StateSequence.hpp"

void Initial::clearStats() {
	Pi->clearObs();
}

void Initial::addObs( StateSequence& seq ) {
	if ( !seq.compress ) {
		for ( size_t i = 0; i < seq.stateSequence.size(); i++ ) {
			Pi->addObs( seq.stateSequence[i], 1 );
		}
	} else {
		for ( size_t i = 0; i < seq.blocks.size(); ++i ) {
			Pi->addObs( seq.stateSequence[i], seq.blocks[i].size );
		}
	}
}


void Initial::sampleParameters() {
	Pi->sampleParameters();
}



void Initial::addHistory() {
	Pi->addHistory();
}



double Initial::pdf( int i ) {
	return Pi->pdf( i );
}
