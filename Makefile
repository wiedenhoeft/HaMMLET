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

# $@ name of target
# $< first dependence
# $+ all dependencies

EXECUTABLE=hammlet

COMPILER=g++
CFLAGS=-c -Werror  --std=c++11   -fmax-errors=1 
# -Wall -Wuninitialized
SOURCES=main.cpp
OBJECTS=$(SOURCES:%.cpp=%.o)

all: CFLAGS +=  -O3
all: hammlet-manpage.hpp $(SOURCES) $(EXECUTABLE)  

debug: CFLAGS += -g
debug:  hammlet-manpage.hpp $(SOURCES) $(EXECUTABLE) 


%.o: %.cpp
	$(COMPILER) $(CFLAGS) $< -o $@


$(EXECUTABLE): $(OBJECTS) 
	$(COMPILER) $(OBJECTS) -o  $@

clean: 
	rm -rf *.o $(EXECUTABLE)


	
	
	
# Create manuals in different formats from doc/manpage.md. It also produces a header file to hard-code the manpage into the executable for platform-independence. 
# This would typically not be called by the end-user, since it requires that the system has pandoc, xxd, pdflatex and man installed. It is provided here for convenience of the developer.
man: manclean hammlet-manpage.hpp doc/hammlet.man doc/hammlet-manpage.txt
	pandoc doc/hammlet-manpage.md -V lang=en -H doc/pandoc.css  -s -t html > doc/hammlet-manpage.html	
	pandoc doc/hammlet-manpage.md -V lang=en -V papersize=a4 -H doc/man-preamble.tex -s -t latex  --variable classoption=landscape,twocolumn --variable geometry={margin=1in}  -o doc/hammlet-manpage-a4.pdf
	pandoc doc/hammlet-manpage.md -V lang=en -V papersize=letter -H doc/man-preamble.tex -s -t latex  --variable classoption=landscape,twocolumn --variable geometry={margin=1in}  -o doc/hammlet-manpage-letter.pdf

manclean:
	rm -f doc/hammlet-manpage.pdf
	rm -f doc/hammlet-manpage.html
	rm -f doc/hammlet-manpage.txt
	rm -f doc/hammlet.man
	rm -f hammlet-manpage.hpp
	
doc/hammlet.man:
	pandoc doc/hammlet-manpage.md -s -t man --variable adjusting=l > doc/hammlet.man


doc/hammlet-manpage.txt: doc/hammlet.man
	MANWIDTH=80 man  doc/hammlet.man > doc/hammlet-manpage.txt

	
# Create .txt version of manpage, as well as a hexdump with C++ declarations, so we can #include the txt version of the manpage into the source at compile time:
hammlet-manpage.hpp:	doc/hammlet-manpage.txt
	xxd -i doc/hammlet-manpage.txt > hammlet-manpage.hpp
	

