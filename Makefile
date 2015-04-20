SuperMaxRepeat: Source.o
	g++ Source.o -o SuperMaxRepeat

Source.o: Source.cpp
	g++ -c Source.cpp
		
clean:
	rm *.o SuperMaxRepeat
