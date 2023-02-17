#include <string>

#pragma once

using namespace std;

class Settings
{
public:
    Settings(char parameter_file[200]);

    int seed;
    int nsteps;
    int num_fibers;
    int num_beads_per_fiber;
    double tf, dt;
    double interbead_distance;
    double cross_thresh;
    double crosslink_threshold;
    double kf, kc,ko;
    double Lx,Ly,Lz;
    double Df;
    double Lf;
    int write_interval, sample_interval;
    double temperature;
    
    double kb;
    double bead_mass;
    double gamma;
    double fiber_diameter;

    string get_suffix();

    string directory;
    string suffix;
    string networkfile;
    string nodesfile;
    string networkpropertiesfile;
    string timeseriesfile;
};

