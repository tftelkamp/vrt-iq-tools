include(CTest)
include(Catch)

add_executable(tests test_rtlsdr_to_soapy.cpp)
target_link_libraries(tests PRIVATE Catch2::Catch2)

catch_discover_tests(tests ADD_TAGS_AS_LABELS)