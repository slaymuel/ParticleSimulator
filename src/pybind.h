#ifdef PY11
PYBIND11_MODULE(particlesimulator, m) {
    py::enum_<Moves::MoveTypes>(m, "MoveTypes")
        .value("TRANSLATE", Moves::MoveTypes::TRANSLATE)
        .value("GCSINGLEADD", Moves::MoveTypes::GCSINGLEADD)
        .value("GCSINGLEREMOVE", Moves::MoveTypes::GCSINGLEREMOVE)
        .value("ROTATE", Moves::MoveTypes::ROTATE)
        .value("SWAP", Moves::MoveTypes::SWAP)
        .value("SINGLESWAP", Moves::MoveTypes::SINGLESWAP)
        .value("VOLUMEMOVE", Moves::MoveTypes::VOLUMEMOVE)
        .value("CHARGETRANS", Moves::MoveTypes::CHARGETRANS)
        .value("CHARGETRANSRAND", Moves::MoveTypes::CHARGETRANSRAND)
        .value("CLUSTER", Moves::MoveTypes::CLUSTER)
        .value("WIDOMINSERTION", Moves::MoveTypes::WIDOMINSERTION)
        .value("WIDOMDELETION", Moves::MoveTypes::WIDOMDELETION)
        .value("GCADD", Moves::MoveTypes::GCADD)
        .value("GCREMOVE", Moves::MoveTypes::GCREMOVE)
        .value("CHARGETRANSLATE", Moves::MoveTypes::CHARGETRANSLATE);

    py::enum_<Samplers::SamplerTypes>(m, "SamplerTypes")
        .value("DENSITY_X", Samplers::SamplerTypes::DENSITY_X)
        .value("DENSITY_Y", Samplers::SamplerTypes::DENSITY_Y)
        .value("DENSITY_Z", Samplers::SamplerTypes::DENSITY_Z)
        .value("ENERGY", Samplers::SamplerTypes::ENERGY)
        .value("WIDOMHS", Samplers::SamplerTypes::WIDOMHS)
        .value("QDIST", Samplers::SamplerTypes::QDIST)
        .value("XDR", Samplers::SamplerTypes::XDR)
        .value("NUMIONS", Samplers::SamplerTypes::NUMIONS)
        .value("PRESSURE", Samplers::SamplerTypes::PRESSURE)
        .value("PRESSUREV", Samplers::SamplerTypes::PRESSUREV)
        .value("FORCEPRESSURE", Samplers::SamplerTypes::FORCEPRESSURE)
        .value("FORCE", Samplers::SamplerTypes::FORCE)
        .value("CLIFFPRESSURE", Samplers::SamplerTypes::CLIFFPRESSURE)
        .value("MODIFIEDWIDOM", Samplers::SamplerTypes::MODIFIEDWIDOM)
        .value("MODIFIEDWIDOMCOULOMB", Samplers::SamplerTypes::MODIFIEDWIDOMCOULOMB);
    
    
    py::class_<Simulator>(m, "Simulator")
        .def(py::init<double, double, std::string>())
        .def("run", &Simulator::run)
        //.def("add_move", &Simulator::add_move)
        .def("add_move", static_cast<void (Simulator::*)(Moves::MoveTypes, std::vector<double>)>(&Simulator::add_move))
        .def("add_sampler", static_cast<void (Simulator::*)(Samplers::SamplerTypes, 
                                                            std::vector<double>)>(&Simulator::add_sampler))
        .def("set_temperature", &Simulator::set_temperature)
        .def("finalize", &Simulator::finalize)
        .def_readwrite("state", &Simulator::state);


    py::class_<State>(m, "State")
        .def("set_geometry", &State::set_geometry)
        .def("load_cp", &State::load_cp)
        .def("set_energy", &State::set_energy)
        .def("equilibrate", &State::equilibrate)
        .def("load_spline", &State::load_spline)
        .def("reset_energy", &State::reset_energy)
        .def_readwrite("particles", &State::particles)
        .def_readonly("energy", &State::energy)
        .def_readonly("cummulativeEnergy", &State::cummulativeEnergy);


    py::class_<Particles>(m, "Particles")
        .def_readonly("pModel", &Particles::pModel)
        .def_readonly("nModel", &Particles::nModel)
        .def("create", &Particles::create, py::arg("pNum"), py::arg("nNum"), py::arg("params"));

    py::class_<Particle>(m, "Particle")
        .def_readonly("com", &Particle::com)
        .def_readonly("pos", &Particle::pos)
        .def_readonly("qDisp", &Particle::qDisp)
        .def_readonly("q", &Particle::q)
        .def_readonly("b", &Particle::b)
        .def_readonly("r", &Particle::r)
        .def_readonly("rf", &Particle::rf);
}
#endif