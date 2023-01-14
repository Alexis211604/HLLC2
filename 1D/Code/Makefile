# Compilateur utilisé
CC=g++



# On choisit comment on compile
CXX_FLAGS = $(DEBUG_FLAG)

# Le nom de l'exécutable
PROG = main

# Les fichiers source à compiler
SRC = HLLC.cc timescheme.cpp flux.cpp operation.cpp

# La commande complète : compile seulement si un fichier a été modifié
$(PROG) : $(SRC)
	g++ -O0 -o main -g HLLC.cc timescheme.cpp flux.cpp operation.cpp
# Évite de devoir connaitre le nom de l'exécutable
all : $(PROG)

# Supprime l'exécutable, les fichiers binaires (.o) et les fichiers
# temporaires de sauvegarde (~)
clean_1 :
	rm -f *.o *~ Resultats/Test_1/*.dat $(PROG)

clean_2 :
	rm -f *.o *~ Resultats/Test_2/*.dat $(PROG)

clean_3 :
	rm -f *.o *~ Resultats/Test_3/*.dat $(PROG)

clean_4 :
	rm -f *.o *~ Resultats/Test_4/*.dat $(PROG)
clean_5 :
	rm -f *.o *~ Resultats/Test_5/*.dat $(PROG)

clean_6 :
	rm -f *.o *~ Resultats/Test_6/*.dat $(PROG)

clean_7 :
	rm -f *.o *~ Resultats/Test_7/*.dat $(PROG)


