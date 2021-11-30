add_executable(latticemodelsimulationlib_tests tests.cpp LatticeModelSimulationLib_t.cpp)
target_link_libraries(latticemodelsimulationlib_tests PUBLIC latticemodelsimulationlib Catch2::Catch2)

# allow user to run tests with `make test` or `ctest`
catch_discover_tests(latticemodelsimulationlib_tests)
