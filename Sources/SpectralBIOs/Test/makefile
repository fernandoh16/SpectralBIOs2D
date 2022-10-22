SOURCE_DIR = /Users/fernandohenriquez/Documents/ETH/UQ_BIOs_Codes/Sources
BIN_DIR    = /Users/FernandoH/Documents/ETH/Codes/SpectralCode/Tests/Bin

FL = -I/Users/fernandohenriquez/Documents/ETH/UQ_BIOs_Codes/Sources -I /Users/fernandohenriquez/Downloads/boost_1_80_0
CPPFLAGS = -I/Users/FernandoH/Documents/ETH/Codes/SpectralCode -I /Users/fernandohenriquez/Downloads/boost_1_80_0 -lfftw3

Point.o: ${SOURCE_DIR}/Point/Point.cpp ${SOURCE_DIR}/Point/Point.hpp
	g++ -std=c++14 -c -O $(FL) ${SOURCE_DIR}/Point/Point.cpp 

Quadrature.o: ${SOURCE_DIR}/Quadrature/Quadrature.cpp  ${SOURCE_DIR}/Quadrature/Quadrature.hpp
	g++ -std=c++14 -c -O ${SOURCE_DIR}/Quadrature/Quadrature.cpp

ParametricBoundaryRepresentation.o: ${SOURCE_DIR}/ParametricBoundaryRepresentation/ParametricBoundaryRepresentation.cpp ${SOURCE_DIR}/ParametricBoundaryRepresentation/ParametricBoundaryRepresentation.hpp ${SOURCE_DIR}/Point/Point.hpp
	g++ -std=c++14 -c -O $(FL) ${SOURCE_DIR}/ParametricBoundaryRepresentation/ParametricBoundaryRepresentation.cpp
	
SpectralQuantity.o: ${SOURCE_DIR}/SpectralQuantity/SpectralQuantity.cpp ${SOURCE_DIR}/SpectralQuantity/SpectralQuantity.hpp
	g++ -std=c++14 -c -O $(FL) ${SOURCE_DIR}/SpectralQuantity/SpectralQuantity.cpp

SpectralBIOs.o: ${SOURCE_DIR}/SpectralBIOs/SpectralBIOs.cpp ${SOURCE_DIR}/SpectralBIOs/SpectralBIOs.hpp ${SOURCE_DIR}/ParametricBoundaryRepresentation/ParametricBoundaryRepresentation.hpp ${SOURCE_DIR}/Point/Point.hpp
	g++ -std=c++14 -c -O  $(FL) ${SOURCE_DIR}/SpectralBIOs/SpectralBIOs.cpp

TestSpectralBIOs.o: TestSpectralBIOs.cpp ${SOURCE_DIR}/ParametricBoundaryRepresentation/ParametricBoundaryRepresentation.hpp ${SOURCE_DIR}/Quadrature/Quadrature.cpp ${SOURCE_DIR}/Point/Point.hpp ${SOURCE_DIR}/SpectralQuantity/SpectralQuantity.cpp
	g++ -std=c++14 -c -O $(FL) TestSpectralBIOs.cpp -o ${SOURCE_DIR}/SpectralBIOs/Test/TestSpectralBIOs.o

TestSpectralBIOs: TestSpectralBIOs.o SpectralBIOs.o ParametricBoundaryRepresentation.o Quadrature.o Point.o SpectralQuantity.o
	g++ -std=c++14 -O -lfftw3 -o TestSpectralBIOs TestSpectralBIOs.o SpectralBIOs.o ParametricBoundaryRepresentation.o Quadrature.o Point.o SpectralQuantity.o

clean:
	rm *.o