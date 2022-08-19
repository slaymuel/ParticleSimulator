#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include "../src/particles.h"
/*#include <catch2/reporters/catch_reporter_console.hpp>
#include <catch2/reporters/catch_reporter_registrars.hpp>
namespace Catch
{
	CATCH_REGISTER_REPORTER("console", ConsoleReporter)
}*/

using namespace Simulator;

TEST_CASE( "Particles can be correctly added and removed from Particles object", "[vector]" ) {
    /*int pNum = 1000;
    Particles particles;
    REQUIRE( particles.tot == 0 );

    SECTION( "Randomly adding and removing particles" ) {
        for(int i = 0; i < pNum; i++)
            particles.add_random({10.0, 10.0, 10.0});

        REQUIRE( particles.tot == 1000 );

        
        for(int i = 0; i < 729; i++)
            particles.remove_random();

        REQUIRE( particles.tot == 1000 - 729 );
    }

    auto j = GENERATE(100, 250, 500);
    SECTION( "Adding and removing specific particles" ) {
        std::vector< std::vector<double> > com, pos;
        std::vector<double> charges, r, rf, b_min, b_max;
        std::vector<std::string> names;
        Eigen::Vector3d ae, be, qDisp;

        for(int i = 0; i < j; i++){
            pos.push_back(Random::get_vector());
            com.push_back(Random::get_vector());
            charges.push_back(Random::get_random() - 0.5);
            r.push_back(Random::get_random());
            rf.push_back(Random::get_random());
            b_min.push_back(Random::get_random());
            b_max.push_back(Random::get_random());
            names.push_back("test_particle");
        }

        for(unsigned int i = 0; i < pos.size(); i++){
            ae << pos[i][0], pos[i][1], pos[i][2];
            be << com[i][0], com[i][1], com[i][2];
            qDisp = (ae - be);
            particles.add(com[i], pos[i], qDisp, r[i], rf[i], charges[i], qDisp.norm(), b_min[i], b_max[i], names[i]);
        }

        REQUIRE( particles.tot == pos.size() );

        for(unsigned int i = 0; i < particles.tot; i++){
            REQUIRE(particles[i]->pos[0] == pos[i][0]);
            REQUIRE(particles[i]->pos[1] == pos[i][1]);
            REQUIRE(particles[i]->pos[2] == pos[i][2]);

            REQUIRE(particles[i]->com[0] == com[i][0]);
            REQUIRE(particles[i]->com[1] == com[i][1]);
            REQUIRE(particles[i]->com[2] == com[i][2]);

            REQUIRE(particles[i]->q == charges[i]);
            REQUIRE(particles[i]->q == charges[i]);
            REQUIRE(particles[i]->q == charges[i]);

            REQUIRE(particles[i]->r == r[i]);
            REQUIRE(particles[i]->r == r[i]);
            REQUIRE(particles[i]->r == r[i]);
        }

        REQUIRE(particles.is_valid());
    }*/
}