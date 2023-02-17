#include <math.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <cstring>
#include <sstream>
#include "settings.h"

using namespace std;

template <typename T>
string tostr(const T &t)
{
  ostringstream os;
  os << t;
  return os.str();
}

Settings::Settings(char parameter_file[200])
{

  ifstream infile;
  int check;
  char line[200];
  char curr_param[200];
  char curr_string[200];
  float curr_value;

  // infile.open("parameters");
  infile.open(parameter_file);

  while (infile.peek() != EOF)
  {
    infile.getline(line, sizeof(line));
    check = sscanf(line, "%s%f", curr_param, &curr_value);

    if (strcmp(curr_param, "seed") == 0)
      seed = int(curr_value);
    else if (strcmp(curr_param, "num_fibers") == 0)
      num_fibers = int(curr_value);
    else if (strcmp(curr_param, "num_beads_per_fiber") == 0)
      num_beads_per_fiber = double(curr_value);
    else if (strcmp(curr_param, "interbead_distance") == 0)
      interbead_distance = double(curr_value);
    else if (strcmp(curr_param, "bead_mass") == 0)
      bead_mass = double(curr_value);
    else if (strcmp(curr_param, "fiber_diameter") == 0)
      fiber_diameter = double(curr_value);
    else if (strcmp(curr_param, "num_steps") == 0)
      nsteps = int(curr_value);
    else if (strcmp(curr_param, "final_time") == 0)
      tf = double(curr_value);
    else if (strcmp(curr_param, "timestep") == 0)
      dt = double(curr_value);
    else if (strcmp(curr_param, "temperature") == 0)
      temperature = double(curr_value);
    else if (strcmp(curr_param, "friction_coefficient") == 0)
      gamma = double(curr_value);
    else if (strcmp(curr_param, "cross_thresh") == 0)
      cross_thresh = double(curr_value);
    else if (strcmp(curr_param, "Lx") == 0)
      Lx = double(curr_value);
    else if (strcmp(curr_param, "Ly") == 0)
      Ly = double(curr_value);
    else if (strcmp(curr_param, "Lz") == 0)
      Lz = double(curr_value);
    else if (strcmp(curr_param, "fiber_stiffness") == 0)
      kf = double(curr_value);
    else if (strcmp(curr_param, "bending_stiffness") == 0)
      ko = double(curr_value);
    else if (strcmp(curr_param, "crosslink_stiffness") == 0)
      kc = double(curr_value);
    else if (strcmp(curr_param, "write_interval") == 0)
      write_interval = int(curr_value);
    else if (strcmp(curr_param, "sample_interval") == 0)
      sample_interval = int(curr_value);
  }

  infile.close();

  directory += "../raw_data/";
  suffix += get_suffix();

  networkfile += directory + "network_" + suffix;
  nodesfile += directory + "nodes_" + suffix;
  networkpropertiesfile += directory + "networkproperties_" + suffix + ".dat";
  timeseriesfile += directory + "timeseries_" + suffix + ".dat";

  //kb is 1,38×10^−14
  kb=1.38*(pow(10,-14));
  Lf = interbead_distance * (num_beads_per_fiber - 1); // fix the number of beads and l0 and you calculate Lf

}

string Settings::get_suffix()
{

  string num_fibers_str = tostr(num_fibers);
  string nb_str = tostr(num_beads_per_fiber);
  string l0_str = tostr(interbead_distance);
  string ko_str = tostr(ko);
  string seed_str = tostr(seed);

  string suffix;

  suffix += "Nf" + num_fibers_str+ "_nb" + nb_str + "_l0" + l0_str +"_ko"+ko_str+ "_seed" + seed_str;
  cout << suffix << endl;

  return suffix;
}

