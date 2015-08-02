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


#include "Transitions.hpp"
#include <fstream>
#include "StateSequence.hpp"


void Transitions::clearStats() {
	for ( size_t i = 0; i < A.size(); i++ ) {
		A[i]->clearObs();
	}
}

double Transitions::entry( int row, int col ) {	
	return A[row]->pdf( col );
}

double Transitions::logEntry( int row, int col ) {	
	return A[row]->logPdf( col );
}


void Transitions::addObs( StateSequence& seq ) {
	if ( !seq.compress ) {
		for ( size_t i = 0; i < seq.stateSequence.size() - 1; i++ ) {
			A[seq.stateSequence[i]]->addObs( seq.stateSequence[i + 1], 1 );
		}
	} else {
		for ( size_t i = 0; i < seq.stateSequence.size() - 1; i++ ) {
			//self transitions
            A[seq.stateSequence[i]]->addObs( seq.stateSequence[i], seq.blocks[i].size-1 );
            //out of this block            
			A[seq.stateSequence[i]]->addObs( seq.stateSequence[i + 1], 1 );
		}

        //self transitions
        A[seq.stateSequence[seq.stateSequence.size() - 1]]->addObs( seq.stateSequence[seq.stateSequence.size() - 1], seq.blocks[seq.stateSequence.size()-1].size);
	}
}

void Transitions::sampleParameters() {
	for ( size_t i = 0; i < A.size(); i++ ) {
		A[i]->sampleParameters();
	}
}

void Transitions::addHistory() {
	for ( size_t i = 0; i < A.size(); i++ ) {
		A[i]->addHistory();
	}
}


