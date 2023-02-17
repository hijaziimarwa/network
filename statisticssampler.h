#pragma once

#include <stdio.h>
#include "system.h"
#include <fstream>

class System;
class Settings;

class StatisticsSampler
{
private:
    System *system;
    Settings *settings;

public:
    StatisticsSampler(System *system);
    void sample();
    void write_data();
    void reset_data_array();
    void write_network();
    void write_network_properties();
    void blue_print();
    void lets_try();
    void sample_average_fiber_length();
    void sample_distance_between_nodes();
    void velocities_of_beads();
    void velocities_of_nodes();
    void sample_total_force();
    void compute_rdf();
    void initialize_timeseriesfile();
    void distance_between_node_and_its_neighbours();
 
    double neighbour_distance;
   

    void write_nodes();

    double Fx, Fy, Fz;
    double v_x, v_y, v_z;
    double global_average_fiber_length;
    double delta;
    double distance_between_nodes;
    double *distance_node;
    int num_sampled_nodes;
    int num_sampled_neighbour_nodes;
    double *velocity_node;


    int id_sampled_bead;

    double r_max = 10;
    double dr = 0.1;
    int nbin = 100;
    double *rbin;
    double *gbin;
    double *hbin;

    int num_data;
    int curr_index;
    double **data_array;
};
