#include <iostream>
#include <math.h>
#include <time.h>
#include <fstream>
#include <iomanip>
#include "system.h"
#include "statisticssampler.h"
#include "settings.h"

using namespace std;

int main(int args, char *argv[])
{
  char *parameter_file;
  parameter_file = argv[1];

  cout << "### LOADING PARAMETERs ###" << endl;
  Settings *settings = new Settings(parameter_file); // reset all parameters and loads them from parameter file
  cout << endl;

  System *system = new System(); // create new System
  system->setup(settings);
  StatisticsSampler *sampler = new StatisticsSampler(system); // creates the samlping class

  system->velocity_initialization();
  cout << "### RUNNING SIMULATION ###" << endl;
  //create a file 
  ofstream energyfile;
  energyfile.open("/home/hijazi/Documents/Collagen_Fiber_Network/code_static/energy.txt");

  system->t = 0;
  //system->kc = 0;
  int t_counter = 0;
  while (system->t <= settings->tf)
  {
    system->t += system->dt;
    system->steps += 1;
   // cout<<"hello "<<endl;
    system->compute_forces();
    system->integrate();
    system->total_potential_energy = system->potential_energy_stretching + system->potential_energy_bending;

    system->total_energy = system->kinetic_energy + system->total_potential_energy;

    system->temperature = 2 * system->kinetic_energy / (3 * system->num_active_beads * settings->kb);
    energyfile<<system->t<<" " << system->temperature<<endl;
    t_counter++;

    if (t_counter == settings->sample_interval)
    {

     sampler->sample();


      t_counter = 0;
      if (sampler->curr_index == settings->write_interval)
      { 
        sampler->write_data();
        
      }
    }
  }
  

  cout << "==> Simulation successfully performed." << endl;

  return 0;
}
