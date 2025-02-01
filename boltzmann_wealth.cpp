#include <unordered_map>
#include <vector>
#include <string>
#include <tuple>
#include <cmath>
#include <algorithm>
#include <random>
#include <cassert>
using namespace std;

// This code is authored by: Peijing Li 
// Contact: pl2932@columbia.edu
// Initial Version: January 25, 2025
//
// Description:
// This source code implements an agent-based model for wealth distribution,
// inspired by the Boltzmann distribution. The simulation involves agents
// that interact, transact, and evolve their wealth over time.
//
// Classes:
// - Agent: Represents an individual agent with its physical coordinates, wealth, and a wealth class.
// - BoltzmannWealthModel: Simulates the collective behaviour of agents, including transactions
//   and movement within a defined spatial region.

class Agent {
public:
    double x, y;                 // Spatial coordinates
    double wealth;               // Wealth
    std::string wealth_class;    // Wealth class 

    static std::unordered_map<std::string, int> wealth_class_map; 
    static std::vector<double> wealth_thresholds;

    // Constructor
    // Initializes an agent with a given position and wealth level.
    Agent(double x, double y, double initial_wealth)
        : x(x), y(y), wealth(initial_wealth){
        wealth_class = determine_wealth_class();
    }

    // Determines the wealth class based on the agent's wealth.
    std::string determine_wealth_class() {
        if (wealth <= wealth_thresholds[0]) {
            return "poverty";
        } else if (wealth <= wealth_thresholds[1]) {
            return "lower middle";
        } else if (wealth <= wealth_thresholds[2]) {
            return "upper middle";
        } else {
            return "upper";
        }
    }


    int get_level(const std::string& class_name) const {
        auto it = wealth_class_map.find(class_name);
        if (it != wealth_class_map.end()) {
            return it->second;
        }
        return -1;
    }



    // Simulates a transaction between the current agent and another agent.
    // The amount of transaction is determined by the risk factor.
    void transact_risk(Agent& other, double risk_factor) {
        double wealth_gap = std::abs(wealth - other.wealth + 0.1);
        double exchange_amt = risk_factor * wealth_gap;
        exchange_amt = std:: min(exchange_amt,std::min(wealth,other.wealth));
        wealth -= exchange_amt;
        other.wealth += exchange_amt;
        wealth_class = determine_wealth_class();
        other.wealth_class = other.determine_wealth_class();
    }


    // Simulates a transaction between the current agent and another agent 
    // according to their gap in wealth level
    void transact_class(Agent& other, double class_factor){
        double overall_amt = wealth;
        int level_diff = abs(get_level(wealth_class) - get_level(other.wealth_class));
        double exchange_amt = overall_amt*class_factor 
                        * (3 - min(3,level_diff))/3;
        exchange_amt = std::max(0.0, std::min(exchange_amt, wealth));
        wealth -= exchange_amt;
        other.wealth += exchange_amt;
        wealth_class = determine_wealth_class();
        other.wealth_class = other.determine_wealth_class();
    }

    // function determine wealth level

    // Updates the spatial coordinates of the agent.
    void movement(double new_x, double new_y){
        x = new_x;
        y = new_y;
    }
};

std::unordered_map<std::string, int> Agent::wealth_class_map = {
    {"poverty", 0},
    {"lower middle", 1},
    {"upper middle", 2},
    {"upper", 3}
};

std::vector<double> Agent::wealth_thresholds = {100.0,100.0,100.0};

class BoltzmannWealthModel{
public:
    // Constructor
    // Initializes the model with a spatial region, agent count, and wealth parameters.
    BoltzmannWealthModel(double width, double height, int num_agents, double initial_wealth,
                         const std::vector<std::tuple<double, double, double>>& additional_agents = {})
        : width(width), height(height), num_agents(num_agents),
         initial_wealth(initial_wealth){
        initialize_agents(); // Generate default agents.

        // Add user-provided agents.
        for (const auto& data : additional_agents) {
            double x, y, wealth;
            std::tie(x, y, wealth) = data;
            agents.emplace_back(x, y, wealth);
        }
    }
    // Initializes agents with random positions and wealth.
    void initialize_agents() {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dist_x(0, width);
        std::uniform_real_distribution<> dist_y(0, height);

        for (int i = 0; i < num_agents; i++) {
            agents.emplace_back(dist_x(gen), dist_y(gen), initial_wealth);
        }
    }

    static void set_wealth_thresholds(const std::vector<double>& thresholds) {
        Agent::wealth_thresholds = thresholds;
    }

    // simulation step under different transaction method
    void step_risk_transaction(double neighbourhood, double risk_factor) {
        transact_with_method(neighbourhood, "risk", risk_factor);
    }
    void step_class_transaction(double neighbourhood, double class_factor) {
        transact_with_method(neighbourhood, "level", class_factor);
    }
    void step_moving_house(double neighbourhood, double moving_cost){
        moving_house(neighbourhood,moving_cost);
    }

    // Collects the wealth of all agents for analysis.
    std::vector<double> collect_wealth() const {
        std::vector<double> wealths;
        wealths.reserve(agents.size());
        for (const auto& agent : agents) {
            wealths.push_back(agent.wealth);
        }
        return wealths;
    }
    std::vector<double> collect_x_coordinate() const {
        std::vector<double> x_coordiantes;
        x_coordiantes.reserve(agents.size());
        for (const auto& agent : agents) {
            x_coordiantes.push_back(agent.x);
        }
        return x_coordiantes;
    }
    std::vector<double> collect_y_coordinate() const {
        std::vector<double> y_coordiantes;
        y_coordiantes.reserve(agents.size());
        for (const auto& agent : agents) {
            y_coordiantes.push_back(agent.y);
        }
        return y_coordiantes;
    }

private:
    double width, height;             // Dimensions of the spatial region.
    int num_agents;                   // Total number of agents.
    double initial_wealth;            // Initial wealth for all agents.
    std::vector<Agent> agents;        // Collection of agents.
    std::vector<double> wealth_thresholds; // wealth thresholds

    // Determine the neighbourhood of the agent, 
    // randomly pick one neighbour to transact and 
    // determine the transaction method
    void transact_with_method(double neighbourhood, 
            const std::string& method, double param){
        std::random_device rd;
        std::mt19937 gen(rd());

        for(auto& agent: agents){
             std::vector<Agent*> neighbours = get_neighbours(agent,neighbourhood);
             Agent* other = nullptr;
            if (!neighbours.empty()){
                std::uniform_int_distribution<> dist(0, neighbours.size() - 1);
                other = neighbours[dist(gen)];
            }  else {
                // Robustness for sparse spatial distribution
                 std::vector<Agent*> all_agents;
                    for (auto& a : agents) {
                        if (&a != &agent) {
                            all_agents.push_back(&a);
                        }
                    }
                    if (!all_agents.empty()) {
                    std::uniform_int_distribution<> dist(0, all_agents.size() - 1);
                    other = all_agents[dist(gen)];
                    }
                }
            if (method == "risk"){
                agent.transact_risk(*other, param);
            } else if (method == "level"){
                agent.transact_class(*other, param);
            } else {
                throw std::invalid_argument("unknown transaction method:" + method);
            }
            }
        }
    
    void moving_house(double neighbourhood, 
            double moving_cost){
        std::random_device rd;
        std::mt19937 gen(rd());

        for(auto& agent: agents){
            if(agent.get_level(agent.wealth_class) > 1){
            std::vector<Agent*> neighbours = get_neighbours(agent,neighbourhood);

            Agent* other = nullptr;
            double total_wealth = 0.0;
            double average_wealth = 0.0;
            double x_weighted = 0.0;
            double y_weighted = 0.0;

            for(const auto& neighbour:neighbours){
                if(neighbour != nullptr){
                double weight = neighbour->wealth;
                total_wealth += weight;
                x_weighted += weight * neighbour->x;
                y_weighted += weight * neighbour->y;
                }
            }

            if(total_wealth == 0){
                continue;
            }

            double new_x = x_weighted / total_wealth;
            double new_y = y_weighted / total_wealth;

            double dx = agent.x - new_x;
            double dy = agent.y - new_y;
            double dr = std::sqrt(dx * dx + dy * dy);
            double moving_cost_total = dr*moving_cost;
            if(agent.wealth > 5 * moving_cost_total){
                agent.movement(new_x,new_y);
                agent.wealth -= moving_cost_total; 
            }
          }
        }
    }


    // Finds neighbors of a given agent within a specified distance.
    std::vector<Agent*> get_neighbours(Agent& agent, double neighbour_dist) {
        std::vector<Agent*> neighbours;
        for (auto& other : agents) {
            if (&other != &agent && distance(agent, other) <= neighbour_dist) {
                    neighbours.push_back(&other);
            }
        }
        return neighbours;
    }
    // Computes the Euclidean distance between two agents.
    double distance(const Agent& a, const Agent& b) const {
        return std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
    }
};
