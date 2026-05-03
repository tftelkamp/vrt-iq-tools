//
// SPDX-License-Identifier: MIT
//

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_session.hpp>

TEST_CASE( "Assert that something is true (pass)", "[demo]" ) {
    REQUIRE( "1" == "1" );
}

int main(int argc, char** argv)
{
    Catch::Session session;

    int returnCode = session.applyCommandLine( argc, argv );
    if( returnCode != 0 )
        return returnCode;

    return session.run();
}