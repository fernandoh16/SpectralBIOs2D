set(UQ_BIOs_Tests
	GenerateMeasurements
    TestBI_IPL_Serial
    #TestBIOperations
    TestDirichletProblem
    #TestForwardUQ_IPL
    #TestForwardUQ_Uniform
    #TestPBR
    #TestQoI
    #TestSpectralBIOs   
    TestSpectralQuantity 
    PARENT_SCOPE
)

set(UQ_BIOs_Tests_Parallel 
    TestBI_IPL_Full
    TestBI_Halton_Full
    PARENT_SCOPE
)

