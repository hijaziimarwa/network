#pragma once

class ThreadControl;
class Settings;
class MDTimer;
class Random;

#include <complex>
#include <fstream>
#include <vector>
#define EMPTY -1

using namespace std;

//#######################################################################

class Cluster
{
private:
public:
    Cluster()
    {
    }
    int head_vertex;
    int num_vertices;
};

//#######################################################################

class Bead
{
private:
public:
    Bead()
    {
    }
    Random *rnd;
    int *neighbours;
    double *eq_dist;
    double mass;
    int id;
    int num_neighbours;
    int num_active_neighbours;
    int max_neighbours;
    int num_crosslinks;
    bool is_mobile;
    bool is_active;

    double position[3], velocity[3];
    double previous_position[3], previous_velocity[3];

    double new_position[3], new_velocity[3];

    double force[3];
    int cluster;

    bool checked;
    int crosslink_node;
    int fiber;
    int next_in_cluster, previous_in_cluster;
    void generate_bead(Settings &settings, double x, double y, double z);
    void deactivate();
};

//#######################################################################

class Fiber
{
private:
public:
    Fiber()
    {
    }

    Random *rnd;

    int id;
    int num_beads, max_beads;
    double r0[3], u[3];
    double length;
    int *neighbour_fibers;
    int num_crosslinks;
    double interbead_eqlength;
    bool is_active;

    void generate_fiber(Settings &settings);
};

class Node
{
private:
public:
    Node()
    {
    }

    int id;
    int *beads;
    int num_beads;
    int *neighbours;
    int num_neighbours;
    double mass;
    double position[3];
    double neighbor_position[3];
    double velocity[3];
};

//##############################################################################

class System
{
private:
public:
    System() {}

    Settings *settings;
    Random *rnd;

    //##############################################################################

    void setup(Settings *settings_);
    void compute_bending_forces();

    void initialize();
    void generate_fibers();
    void generate_beads();
    void activate_fibers();
    void check();
    void initialize_crosslinks();
    void get_distance(double &smin, double &tmin, double &dr2, int n, int m);
    void define_nodes();
    void connect_nodes();
    void distance_between_nodes();
    void reset_forces();
    void compute_forces();
    void compute_stretching_forces();
    void compute_langevin_forces();
    void update_springs();

    void langevin_noise();
    void output_data();


    void deactivate_beads();

    void relabel_beads();
    void find_clusters();

    void clear_clusters();
    void move_between_clusters(int i, int cto);
    void insert_to_cluster(int c, int i);
    void count_clusters();

    void compute_rdf();

    void half_kick();
    void full_kick();
    void move();
    void step();
    void integrate();
    void md_verlet();
    void velocity_initialization();

    Fiber **fibers;
    Bead **beads;
    Cluster **clusters;
    Node **nodes;
    int num_fibers, num_active_fibers;
    int num_beads, num_active_beads;
    int num_crosslinks, num_active_crosslinks;
    int checked_sides,average_of_sides;
    double distance_of_pairs;

    int num_nodes;
    double density;
    //double kc;
    double t, tf, dt, dt_half;
    double fiber_diameter;
    double fs;
    int steps;
    double final_time, time_steps;
    int max_beads_per_fiber;
    int largest_cluster;
    double potential_energy, kinetic_energy;
    double temperature;
    double potential_energy_stretching, potential_energy_bending;
    double total_potential_energy, total_energy;
    double cross_thresh;
    double average_position[3], average_velocity[3], average_velocities[3];

};
