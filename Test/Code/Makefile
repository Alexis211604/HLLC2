# Compilateur utilisé
CC=g++



# On choisit comment on compile
CXX_FLAGS = $(DEBUG_FLAG)

# Le nom de l'exécutable
PROG = main

# Les fichiers source à compiler
SRC = 2D.cc timescheme.cpp flux.cpp operation.cpp

# La commande complète : compile seulement si un fichier a été modifié
$(PROG) : $(SRC)
	g++ -O0 -o main -g 2D.cc timescheme.cpp flux.cpp operation.cpp
# Évite de devoir connaitre le nom de l'exécutable
all : $(PROG)

# Supprime l'exécutable, les fichiers binaires (.o) et les fichiers
# temporaires de sauvegarde (~)
clean :
	rm -f *.o *~ Affichage/*.dat $(PROG)
