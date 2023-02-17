#include "statisticssampler.h"
#include "settings.h"
#include "system.h"
#include <iomanip>
#include <string.h>
#include <iostream>
#include <numeric>
#include <vector>
using namespace std;

StatisticsSampler::StatisticsSampler(System *system_)
{
  system = system_;
  settings = system->settings;

  num_sampled_nodes = 3;
  num_sampled_neighbour_nodes = 0;
  for (int i = 0; i < num_sampled_nodes; i++)
  {
    num_sampled_neighbour_nodes += system->beads[i]->num_neighbours;
    // cout<<num_sampled_neighbour_nodes<<endl;
  }
  distance_node = new double[num_sampled_neighbour_nodes];

  for (int i = 0; i < num_sampled_nodes; i++)
  {
    distance_node[i] = 0;
  }

  num_data = 4;
  data_array = new double *[settings->write_interval];
  for (int i = 0; i < settings->write_interval; i++)
  {
    data_array[i] = new double[num_data];
  }
  reset_data_array();

  rbin = new double[nbin];
  gbin = new double[nbin];
  hbin = new double[nbin];
  for (int i = 0; i < nbin; i++)
  {
    rbin[i] = dr * (i + 0.5);
    // cout << rbin[i] << endl;
    gbin[i] = 0;
    hbin[i] = 0;
  }

  initialize_timeseriesfile();
}

void StatisticsSampler::reset_data_array()
{
  curr_index = 0;
  for (int i = 0; i < settings->write_interval; i++)
  {
    for (int j = 0; j < num_data; j++)
    {
      data_array[i][j] = 0;
    }
  }
}

void StatisticsSampler::initialize_timeseriesfile()
{
  remove(settings->timeseriesfile.c_str());

  ofstream datafile;
  datafile.open(settings->timeseriesfile, std::ios_base::app);

  datafile << "# TIMESERIES DATA" << endl;
  datafile << "#Column 1 : Time " << endl;
  datafile << "#Column 2 : Kinetic energy " << endl;
  
}

void StatisticsSampler::write_network()
{
  int j;

  // Create file name

  string curr_networkfile;
  stringstream step_string;
  step_string << system->steps;
  string step_str = step_string.str();

  curr_networkfile += settings->networkfile + "_" + step_str + ".dat";
  // cout << curr_networkfile << endl;

  remove(curr_networkfile.c_str());

  ofstream datafile;
  datafile.open(curr_networkfile.c_str(), std::ios_base::app);

  // Write coordinate of beads
  for (int i = 0; i < system->num_beads; i++)
  {
    if (system->beads[i]->is_active)
    {
      datafile << system->beads[i]->fiber << " " << i << " " << system->beads[i]->position[0] << " " << system->beads[i]->position[1] << " " << system->beads[i]->position[2] << endl;

    }
  }

  datafile << endl;

  // Delimiter between fibers and crosslinks

  for (int i = 0; i < system->num_beads; i++)
  {
    if (system->beads[i]->is_active)
    {
      datafile << i << " ";
      for (int k = 0; k < system->beads[i]->num_neighbours; k++)
      {
        if (system->beads[system->beads[i]->neighbours[k]]->is_active)
        {
          datafile << system->beads[i]->neighbours[k] << " ";
        }
      }
      datafile << endl;
    }
  }

  datafile.close();
}

void StatisticsSampler::write_nodes()
{
  int j;
  // Create file name

  string curr_nodesfile;
  stringstream step_string;
  step_string << system->steps;
  string step_str = step_string.str();

  curr_nodesfile += settings->nodesfile + "_" + step_str + ".dat";

  remove(curr_nodesfile.c_str());

  ofstream datafile;
  datafile.open(curr_nodesfile.c_str(), std::ios_base::app);

  // Write coordinate of beads
  for (int i = 0; i < system->num_nodes; i++)
  { // cout<<system->num_nodes<<endl;
    system->nodes[i]->position[0] = 0;
    system->nodes[i]->position[1] = 0;
    system->nodes[i]->position[2] = 0;

    for (int j = 0; j < system->nodes[i]->num_beads; j++)
    { // cout << system->nodes[i]->num_beads << endl;
      system->nodes[i]->position[0] += system->beads[system->nodes[i]->beads[j]]->position[0] / double(system->nodes[i]->num_beads);
      system->nodes[i]->position[1] += system->beads[system->nodes[i]->beads[j]]->position[1] / double(system->nodes[i]->num_beads);
      system->nodes[i]->position[2] += system->beads[system->nodes[i]->beads[j]]->position[2] / double(system->nodes[i]->num_beads);

      // cout<<"   b   "<<system->beads[system->nodes[i]->beads[j]]->position[0]<< " " <<system->beads[system->nodes[i]->beads[j]]->position[1]<< " "<<system->beads[system->nodes[i]->beads[j]]->position[2]<<endl;
      // cout<<system->nodes[i]->position[0]<<system->nodes[i]->position[1]<<system->nodes[i]->position[2]<<endl;
      // cout<<system->beads[system->nodes[i]->beads[j]]->velocity[0]<<system->beads[system->nodes[i]->beads[j]]->velocity[1]<<system->beads[system->nodes[i]->beads[j]]->velocity[2]<<endl;
    }

    datafile << system->nodes[i]->id << " " << system->nodes[i]->position[0] << " " << system->nodes[i]->position[1] << " " << system->nodes[i]->position[2] << endl;
  }
  datafile << endl;

  // Delimiter between fibers and crosslinks

  for (int i = 0; i < system->num_nodes; i++)
  {
    // cout << i << " " << system->nodes[i]->num_neighbours << endl;
    datafile << system->nodes[i]->id << " ";
    for (int k = 0; k < system->nodes[i]->num_neighbours; k++)
    {
      datafile << system->nodes[i]->neighbours[k] << " ";
    }
    datafile << endl;
  }
  datafile << endl;
  // system->distance_between_nodes();

  datafile.close();
}

void StatisticsSampler::write_network_properties()
{

  ofstream datafile;
  datafile.open(settings->networkpropertiesfile, std::ios_base::app);
  datafile << system->num_active_crosslinks << " " << system->num_active_fibers << " " << system->checked_sides << endl;
  datafile.close();
}

void StatisticsSampler::write_data()
{
  ofstream datafile;
  datafile.open(settings->timeseriesfile, std::ios_base::app);
  for (int i = 0; i < curr_index; i++)
  {
    for (int j = 0; j < num_data; j++)
    {
      datafile << data_array[i][j] << " ";
    }

    datafile << endl;
  }
  datafile.close();

  reset_data_array();
}

void StatisticsSampler::sample_total_force()
{
  Fx = Fy = Fz = 0;
  for (int i = 0; i < system->num_beads; i++)
  {
    if (system->beads[i]->is_active && system->beads[i]->is_mobile)
    {
      Fx += system->beads[i]->force[0];
      Fy += system->beads[i]->force[1];
      Fz += system->beads[i]->force[2];
    }
  }
}

void StatisticsSampler::velocities_of_beads()
{
  for (int i = 0; i < system->num_beads; i++)
  {
    if (system->beads[i]->is_active)
    {
      v_x += system->beads[i]->velocity[0];
      v_y += system->beads[i]->velocity[1];
      v_z += system->beads[i]->velocity[2];
    }
  }
}

void StatisticsSampler::sample()
{

  data_array[curr_index][0] = system->t;

  data_array[curr_index][1] = system->kinetic_energy;
  data_array[curr_index][2] = system->total_potential_energy;
  data_array[curr_index][3] = system->total_energy;
  
 
  curr_index++;
}

void StatisticsSampler::sample_distance_between_nodes()
{

  for (int n = 0; n < num_sampled_nodes; n++)
  {
    n = 0;

    
    double delta_x = system->nodes[system->nodes[n]->neighbours[0]]->position[0] - system->nodes[n]->position[0];
    double delta_y = system->nodes[system->nodes[n]->neighbours[0]]->position[1] - system->nodes[n]->position[1];
    double delta_z = system->nodes[system->nodes[n]->neighbours[0]]->position[2] - system->nodes[n]->position[2];

    double delta_r = sqrt(delta_x * delta_x + delta_y * delta_y + delta_z * delta_z);

 
    distance_node[n] = delta_r;
  }
}



void StatisticsSampler::velocities_of_nodes()
{

  for (int i = 0; i < system->num_nodes; i++)
  {
    system->nodes[i]->velocity[0] = 0;
    system->nodes[i]->velocity[1] = 0;
    system->nodes[i]->velocity[2] = 0;

    for (int j = 0; j < system->nodes[i]->num_beads; j++)
    {
      system->nodes[i]->velocity[0] += system->beads[system->nodes[i]->beads[j]]->velocity[0] / double(system->nodes[i]->num_beads);
      system->nodes[i]->velocity[1] += system->beads[system->nodes[i]->beads[j]]->velocity[1] / double(system->nodes[i]->num_beads);
      system->nodes[i]->velocity[2] += system->beads[system->nodes[i]->beads[j]]->velocity[2] / double(system->nodes[i]->num_beads);
    }

    double v_x = system->nodes[i]->velocity[0];
    double v_y = system->nodes[i]->velocity[1];
    double v_z = system->nodes[i]->velocity[2];
    double v = sqrt(v_x * v_x + v_y * v_y + v_z * v_z);
  }
}


void StatisticsSampler::distance_between_node_and_its_neighbours()
{

  for (int i = 0; i < system->num_nodes; i++)
  {
    system->nodes[i]->position[0] = 0;
    system->nodes[i]->position[1] = 0;
    system->nodes[i]->position[2] = 0;

    for (int j = 0; j < system->nodes[i]->num_beads; j++)
    {
      system->nodes[i]->position[0] += system->beads[system->nodes[i]->beads[j]]->position[0] / double(system->nodes[i]->num_beads);
      system->nodes[i]->position[1] += system->beads[system->nodes[i]->beads[j]]->position[1] / double(system->nodes[i]->num_beads);
      system->nodes[i]->position[2] += system->beads[system->nodes[i]->beads[j]]->position[2] / double(system->nodes[i]->num_beads);
    }
  }
  int m = 0;
  for (int n = 0; n <num_sampled_nodes; n++)
  {
    for (int i = 0; i < system->nodes[n]->num_neighbours; i++)
    {
      double neigh_x = system->nodes[system->nodes[n]->neighbours[i]]->position[0];
      double neigh_y = system->nodes[system->nodes[n]->neighbours[i]]->position[1];
      double neigh_z = system->nodes[system->nodes[n]->neighbours[i]]->position[2];

      double node_x = system->nodes[n]->position[0];
      double node_y = system->nodes[n]->position[1];
      double node_z = system->nodes[n]->position[2];

      neighbour_distance = sqrt((neigh_x - node_x) * (neigh_x - node_x) + (neigh_y - node_y) * (neigh_y - node_y) + (neigh_z - node_z) * (neigh_z - node_z));

     
      // cout << system->nodes[n]->neighbours[i] << endl;
      distance_node[m] = neighbour_distance;

      m++;
      
    }
  }
}
