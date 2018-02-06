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
# $< first dependency
# $+ all dependencies



COMPILER=g++
CFLAGS=-Werror  --std=c++11   -fmax-errors=1 -Wreturn-type
SRC=./src
TLS=$(SRC)/tools
BIN=./bin
LIB=./lib
DOC=./doc
LOGO=./logo

TOOLS=mapLinesToGenome combineCounts avg maxSegmentation
TOOLSLIST=$(addprefix $(BIN)/, $(TOOLS))

all: CFLAGS +=  -O3
all: $(SRC)/hammlet-manpage.hpp tools hammlet
	chmod ug+x $(BIN)/*

debug: CFLAGS += -g
debug:  $(SRC)/hammlet-manpage.hpp  tools hammlet
	chmod ug+x $(BIN)/*

hammlet: $(SRC)/hammlet-manpage.hpp
	$(COMPILER) $(CFLAGS) $(SRC)/main.cpp  -o $(BIN)/hammlet

	
tools: $(TOOLSLIST)
	

$(BIN)/%:  $(TLS)/%.cpp $(LIB)/gzstream/libgzstream.a
	$(COMPILER) $(CFLAGS)  $< -I$(LIB)/gzstream -L$(LIB)/gzstream -lgzstream -lz -o $@ 

	
# make gzip stream library	
$(LIB)/gzstream/libgzstream.a:
	make -C $(LIB)/gzstream
	
clean: 
	rm -vf $(BIN)/hammlet
	rm -vf $(TOOLSLIST)
	rm -vf $(BIN)/pyhammlet/*.pyc
	rm -vf $(LIB)/gzstream/libgzstream.a
	rm -vf $(LIB)/gzstream/gzstream.o


	
	
# Create manuals in different formats from $(DOC)/manpage.md. It also produces a header file to hard-code the manpage into the executable for platform-independence. 
# This would typically not be called by the end-user, since it requires that the system has pandoc, xxd, pdflatex, awk, base64 and man installed. It is provided here for convenience of the developer.
man: manclean  $(SRC)/hammlet-manpage.hpp $(DOC)/hammlet.man $(DOC)/hammlet-manpage.txt $(DOC)/hammlet-manpage.html $(DOC)/hammlet-manpage-a4.pdf $(DOC)/hammlet-manpage-letter.pdf $(LOGO)/logo-round.base64

	
manclean:
	rm -f $(DOC)/hammlet-manpage-a4.pdf
	rm -f $(DOC)/hammlet-manpage-letter.pdf
	rm -f $(DOC)/hammlet-manpage.html
	rm -f $(DOC)/hammlet-manpage.txt
	rm -f $(DOC)/hammlet.man
	rm -f $(SRC)/hammlet-manpage.hpp

$(LOGO)/logo-round.base64:
	echo -n "data:image/png;base64," > $(LOGO)/logo-round.base64
	base64 -w 0 $(LOGO)/logo-round.png >> $(LOGO)/logo-round.base64
	
$(DOC)/hammlet-manpage-a4.pdf:
	pandoc $(DOC)/hammlet-manpage.md -V lang=english -V papersize=a4 -H $(DOC)/man-preamble.tex -s -t latex  -V classoption=landscape,twocolumn -V geometry={margin=1in}  -o $(DOC)/hammlet-manpage-a4.pdf

$(DOC)/hammlet-manpage-letter.pdf:
	pandoc $(DOC)/hammlet-manpage.md -V lang=english -V papersize=letter -H $(DOC)/man-preamble.tex -s -t latex  -V classoption=landscape,twocolumn -V geometry={margin=1in}  -o $(DOC)/hammlet-manpage-letter.pdf

$(DOC)/hammlet-manpage.html: $(LOGO)/logo-round.base64
	pandoc $(DOC)/hammlet-manpage.md -V lang=english -H $(DOC)/pandoc.css  -s -t html | awk 'BEGIN{getline l < "logo/logo-round.base64"}/\"..\/logo\/logo-round.png\"/{gsub("../logo/logo-round.png",l)}1' > $(DOC)/hammlet-manpage.html
	
$(DOC)/hammlet.man:
	pandoc $(DOC)/hammlet-manpage.md -V adjusting=l -s -t man > $(DOC)/hm
	cat $(DOC)/hm $(LOGO)/logo-boxdrawing-centered.groff > $(DOC)/hammlet.man
	rm $(DOC)/hm
	

$(DOC)/hammlet-manpage.txt: $(DOC)/hammlet.man
	MANWIDTH=80 man  $(DOC)/hammlet.man > $(DOC)/hammlet-manpage.txt
	
# Create .txt version of manpage, as well as a hexdump with C++ declarations, so we can #include the txt version of the manpage into the source at compile time:
$(SRC)/hammlet-manpage.hpp:	$(DOC)/hammlet-manpage.txt
	xxd -i $(DOC)/hammlet-manpage.txt > $(SRC)/hammlet-manpage.hpp
	

