#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "boltzmann_wealth.cpp"

namespace py = pybind11;

PYBIND11_MODULE(boltzmann_wealth, m) {
    // Agent class bindings
    py::class_<Agent>(m, "Agent")
        .def(py::init<double, double, double>(),
             py::arg("x"), py::arg("y"), py::arg("initial_wealth"))
        .def_readwrite("x", &Agent::x)
        .def_readwrite("y", &Agent::y)
        .def_readwrite("wealth", &Agent::wealth)
        .def_readwrite("wealth_class", &Agent::wealth_class)
        .def("determine_wealth_class", &Agent::determine_wealth_class)
        .def("get_level", &Agent::get_level, py::arg("class_name"))
        .def("transact_risk", &Agent::transact_risk, py::arg("other"), py::arg("risk_factor"))
        .def("transact_class", &Agent::transact_class, py::arg("other"), py::arg("class_factor"))
        .def("movement", &Agent::movement, py::arg("new_x"), py::arg("new_y"));

    // BoltzmannWealthModel class bindings
    py::class_<BoltzmannWealthModel>(m, "BoltzmannWealthModel")
        .def(py::init<double, double, int, double, const std::vector<std::tuple<double, double, double>>&>(),
             py::arg("width"), py::arg("height"), py::arg("num_agents"),
             py::arg("initial_wealth"),
             py::arg("additional_agents") = std::vector<std::tuple<double, double, double>>())
        .def("initialize_agents", &BoltzmannWealthModel::initialize_agents)
        .def_static("set_wealth_thresholds", &BoltzmannWealthModel::set_wealth_thresholds,
                          py::arg("thresholds"))
        .def("step_risk_transaction", &BoltzmannWealthModel::step_risk_transaction,
             py::arg("neighbourhood"), py::arg("risk_factor"))
        .def("step_class_transaction", &BoltzmannWealthModel::step_class_transaction,
             py::arg("neighbourhood"), py::arg("class_factor"))
        .def("step_moving_house", &BoltzmannWealthModel::step_moving_house,
             py::arg("neighbourhood"), py::arg("moving_cost"))
        .def("collect_wealth", &BoltzmannWealthModel::collect_wealth)
        .def("collect_x_coordinate", &BoltzmannWealthModel::collect_x_coordinate)
        .def("collect_y_coordinate", &BoltzmannWealthModel::collect_y_coordinate);
}
