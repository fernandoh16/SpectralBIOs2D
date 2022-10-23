include_directories(${Boost_INCLUDE_DIRS})

set(SOURCE_FILES
    ${CMAKE_CURRENT_LIST_DIR}/Point/Point.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ParametricBoundaryRepresentation/ParametricBoundaryRepresentation.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Quadrature/Quadrature.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SpectralBIOs/SpectralBIOs.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SpectralQuantity/SpectralQuantity.cpp
    ${CMAKE_CURRENT_LIST_DIR}/BIOsProblems/DirichletProblem.cpp
    # ${CMAKE_CURRENT_LIST_DIR}/QoIBoundary/QoIBoundary.cpp
    # ${CMAKE_CURRENT_LIST_DIR}/QoISolution/QoISolution.cpp
    ${CMAKE_CURRENT_LIST_DIR}/BoundaryPotentials/BoundaryPotentials.cpp
	#${CMAKE_CURRENT_LIST_DIR}/Halton/Halton.cpp
    PARENT_SCOPE
)

set(HEADER_FILES
    ${CMAKE_CURRENT_LIST_DIR}/Point/Point.hpp
    ${CMAKE_CURRENT_LIST_DIR}/ParametricBoundaryRepresentation/ParametricBoundaryRepresentation.hpp
    ${CMAKE_CURRENT_LIST_DIR}/Quadrature/Quadrature.hpp
    ${CMAKE_CURRENT_LIST_DIR}/SpectralBIOs/SpectralBIOs.hpp
    ${CMAKE_CURRENT_LIST_DIR}/BoundaryRepresentations/Fourier.hpp
    ${CMAKE_CURRENT_LIST_DIR}/BoundaryRepresentations/Kidney.hpp
    ${CMAKE_CURRENT_LIST_DIR}/BoundaryRepresentations/Boomerang.hpp
    ${CMAKE_CURRENT_LIST_DIR}/SpectralQuantity/SpectralQuantity.hpp
    ${CMAKE_CURRENT_LIST_DIR}/BIOsProblems/DirichletProblem.hpp
    # ${CMAKE_CURRENT_LIST_DIR}/QoIBoundary/QoIBoundary.hpp
    # ${CMAKE_CURRENT_LIST_DIR}/QoISolution/QoISolution.hpp
    ${CMAKE_CURRENT_LIST_DIR}/BayesianInversion/BayesianInversion.hpp
    ${CMAKE_CURRENT_LIST_DIR}/BoundaryPotentials/BoundaryPotentials.hpp
	${CMAKE_CURRENT_LIST_DIR}/BayesianInversion/InverseScatteringBI/ForwardModel.hpp
	${CMAKE_CURRENT_LIST_DIR}/BayesianInversion/InverseScatteringBI/PlaneWave.hpp
	${CMAKE_CURRENT_LIST_DIR}/BayesianInversion/InverseScatteringBI/Potential.hpp
	${CMAKE_CURRENT_LIST_DIR}/BayesianInversion/InverseScatteringBI/QuantityOfInterest.hpp
	#${CMAKE_CURRENT_LIST_DIR}/Halton/Halton.hpp
    PARENT_SCOPE
)
