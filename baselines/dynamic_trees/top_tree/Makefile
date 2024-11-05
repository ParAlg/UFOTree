main: main.c tree.h top_tree.h kruskal.h debug.h tree.o top_tree.o kruskal.o debug.o
	gcc -Wall -g main.c tree.o top_tree.o kruskal.o debug.o -o main

tree.o: tree.c tree.h
	gcc -Wall -c -g tree.c

top_tree.o: top_tree.c top_tree.h tree.h
	gcc -Wall -c -g top_tree.c

kruskal.o: kruskal.c kruskal.h
	gcc -Wall -c -g kruskal.c

debug.o: debug.c top_tree.h tree.h
	gcc -Wall -c -g debug.c

clean:
	rm -f *.o main

.PHONY: clean
