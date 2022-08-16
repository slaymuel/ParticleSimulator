#ifdef PY11
PYBIND11_MODULE(particlesimulator, m) {
    py::enum_<MoveTypes>(m, "MoveTypes")
        .value("TRANSLATE", MoveTypes::TRANSLATE)
        .value("GCSINGLEADD", MoveTypes::GCSINGLEADD)
        .value("GCSINGLEREMOVE", MoveTypes::GCSINGLEREMOVE)
        .value("ROTATE", MoveTypes::ROTATE)
        .value("SWAP", MoveTypes::SWAP)
        .value("SINGLESWAP", MoveTypes::SINGLESWAP)
        .value("VOLUMEMOVE", MoveTypes::VOLUMEMOVE)
        .value("CHARGETRANS", MoveTypes::CHARGETRANS)
        .value("CHARGETRANSRAND", MoveTypes::CHARGETRANSRAND)
        .value("CLUSTER", MoveTypes::CLUSTER)
        .value("WIDOMINSERTION", MoveTypes::WIDOMINSERTION)
        .value("WIDOMDELETION", MoveTypes::WIDOMDELETION)
        .value("GCADD", MoveTypes::GCADD)
        .value("GCREMOVE", MoveTypes::GCREMOVE)
        .value("CHARGETRANSLATE", MoveTypes::CHARGETRANSLATE);

    

    py::class_<Simulator>(m, "Simulator")
        .def(py::init<double, double, std::string>())
        .def("run", &Simulator::run)
        //.def("add_move", &Simulator::add_move)
        .def("add_move", static_cast<void (Simulator::*)(MoveTypes, std::vector<double>)>(&Simulator::add_move))
        //.def("add_move",  static_cast<void (Simulator::*)(int, std::vector<double>)>(&Simulator::add_move))
        .def("add_sampler", &Simulator::add_sampler, py::arg("i"), py::arg("interval"), py::arg("ds") = 0.05)
        .def("set_temperature", &Simulator::set_temperature)
        .def("set_cp", &Simulator::set_cp)
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
        .def_readonly("particles", &Particles::particles)
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