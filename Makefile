###     This file is part of HaMMLET.
### 
###     HaMMLET is free software: you can redistribute it and/or modify
###     it under the terms of the GNU General Public License as published by
###     the Free Software Foundation, either version 3 of the License, or
###     (at your option) any later version.
### 
###     HaMMLET is distributed in the hope that it will be useful,
###     but WITHOUT ANY WARRANTY; without even the implied warranty of
###     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
###     GNU General Public License for more details.
### 
###     You should have received a copy of the GNU General Public License
###     along with HaMMLET.  If not, see <http://www.gnu.org/licenses/>.


CC=g++
CFLAGS= -c  -O3  --std=c++11 
SOURCES= rnglib.cpp pdflib.cpp asa266.cpp Normal.cpp Uniform.cpp Theta.cpp StateSequence.cpp Categorical.cpp Transitions.cpp Initial.cpp Mixture.cpp NormalTransformed.cpp Product.cpp Model.cpp WaveletTree.cpp main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=hammlet

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o  $@

.cpp.o: 
	$(CC) $(CFLAGS) $< -o $@

mode: 
	g++ -o mode mode.cpp

clean:
	rm -rf *o $(EXECUTABLE)
