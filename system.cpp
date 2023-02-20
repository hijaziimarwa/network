#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fstream>
#include <string.h>
#include <math.h>
#include "system.h"
#include "settings.h"
#include "random.h"
#include "statisticssampler.h"
#include "geometry.h"
#include <iomanip>

#define EMPTY -1

using namespace std;

void System::setup(Settings *settings_)
{
  settings = settings_;              // the System class is defined from a set of settings.
  rnd = new Random(-settings->seed); // random seed

  initialize(); // intialize the system
  // compute_forces();
}

void System::initialize()
{
  // Temporal quantitites
  dt = settings->dt;

  dt_half = 0.5 * dt;
  t = 0;
  steps = 0;

  // Structural quantitites
  num_crosslinks = 0;                                  // total number of crosslinks in the system
  num_fibers = settings->num_fibers;                   // total number of fibers in the system
  max_beads_per_fiber = settings->num_beads_per_fiber; // the number of beads per fiber is set constant here
                                                       // the mass of each bead is set constant here
  num_beads = 0;
  density;

  // Initializing the fibers
  cout << "### GENERATING BEADS ###" << endl;
  cout << "Number of fibers = " << num_fibers << endl;
  cout << "Number of beads per fiber = " << settings->num_beads_per_fiber << endl;
  cout << "inter_bead_distance = " << settings->interbead_distance << endl;
  generate_fibers();
  generate_beads();

  cout << "Total number of beads = " << num_beads << endl;
  cout << "===> Beads generated." << endl
       << endl;

  // Initializing the crosslinks
  cout << "### ESTABLISHING CROSSLINKS###" << endl;
  initialize_crosslinks();
  for (int i = 0; i < num_beads; i++) // Loop over all beads of the system
  {
    /*if (beads[i]->is_active)
    {
      cout << beads[i]->position[0] << " " << beads[i]->position[1] << " " << beads[i]->position[2] << endl;
    }*/
    beads[i]->num_active_neighbours = beads[i]->num_neighbours;
    beads[i]->mass = settings->bead_mass;
  }
  cout << "Number of crosslinks = " << num_crosslinks << endl;
  /**
  for (int i=0; i<num_fibers;i++){
    total_crosslinks += fibers[i]->num_crosslinks;
  }

  cout << "Number of crosslinks (second method) = " << total_crosslinks*0.5 << endl;
  **/
  cout << "===> Crosslinks established." << endl
       << endl;

  // Cleaning up the network
  cout << "### CLEANING UP THE NETWORK###" << endl;
  num_active_beads = num_beads;
  int old_num_active_beads = 0;
  while (old_num_active_beads != num_active_beads)
  {                                          // Check if the number of beads after the iteration is the same as the number of beads before
    old_num_active_beads = num_active_beads; // Save the number of beads before the iteration
    find_clusters();                         // Find the connectivity clusters
    // cout << num_active_beads << endl;
    deactivate_beads();
    activate_fibers();
  }

  cout << num_active_beads << " remaining interconnected beads" << endl;
  cout << num_active_fibers << " remaining interconnected fibers" << endl;
  cout << num_active_crosslinks << " remaining interconnected crosslinks" << endl;
  cout << "===> Network finalized." << endl
       << endl;

 
  define_nodes();
  

}

void System::generate_fibers()
{
  fibers = new Fiber *[num_fibers]();
  for (int n = 0; n < num_fibers; n++)
  {
    fibers[n] = new Fiber();
    fibers[n]->id = n;
    fibers[n]->generate_fiber(*settings);
    num_beads += fibers[n]->num_beads;
  }
}

void System::generate_beads()
{
  double x, y, z;
  int curr_id;

  beads = new Bead *[num_beads]();

  for (int n = 0; n < num_fibers; n++)
  {
    x = fibers[n]->r0[0];
    y = fibers[n]->r0[1];
    z = fibers[n]->r0[2];
    for (int i = 0; i < fibers[n]->num_beads; i++)
    {
      curr_id = fibers[n]->id * fibers[n]->num_beads + i;
      beads[curr_id] = new Bead();
      beads[curr_id]->id = curr_id;
      beads[curr_id]->max_neighbours = settings->num_fibers + 2;
      beads[curr_id]->fiber = n;
      beads[curr_id]->generate_bead(*settings, x, y, z); // initialize bead
      // Shift new bead along the predefined direction
      x += fibers[n]->u[0];
      y += fibers[n]->u[1];
      z += fibers[n]->u[2];
    }

    for (int i = 0; i < fibers[n]->num_beads - 1; i++)
    {
      curr_id = fibers[n]->id * fibers[n]->num_beads + i;
      beads[curr_id]->neighbours[beads[curr_id]->num_neighbours] = curr_id + 1;
      beads[curr_id + 1]->neighbours[beads[curr_id + 1]->num_neighbours] = curr_id;

      beads[curr_id]->eq_dist[beads[curr_id]->num_neighbours] = settings->interbead_distance;
      beads[curr_id + 1]->eq_dist[beads[curr_id + 1]->num_neighbours] = settings->interbead_distance;

      beads[curr_id]->num_neighbours++;
      beads[curr_id + 1]->num_neighbours++;
    }
  }
}

void System::get_distance(double &t1min, double &t2min, double &dr2, int n, int m)
{
  double un[3], um[3], Delta[3];
  double un2 = 0;
  double um2 = 0;
  double und, umd;
  int num_beads_n, num_beads_m;

  num_beads_n = fibers[n]->num_beads;
  num_beads_m = fibers[m]->num_beads;

  for (int d = 0; d < 3; d++)
  {
    un[d] = fibers[n]->u[d] * (fibers[n]->num_beads - 1);
    um[d] = fibers[m]->u[d] * (fibers[m]->num_beads - 1);
    Delta[d] = fibers[n]->r0[d] - fibers[m]->r0[d];
  }

  t1min = 0;
  t2min = 0;
  dr2 = 0;
  distance_between_segments(un, um, Delta, t1min, t2min, dr2);
}

void System::initialize_crosslinks()
{
  double t1min, t2min, dr2, dr;
  double x0n, x1n, y0n, y1n, z0n, z1n;
  double x0m, x1m, y0m, y1m, z0m, z1m;
  int num_beads_n, num_beads_m;
  int bead_n1, bead_n2;

  cross_thresh = settings->cross_thresh;

  for (int n1 = 0; n1 < num_fibers; n1++)
  {
    for (int n2 = 0; n2 < n1; n2++)
    {

      t1min = 0;
      t2min = 0;
      get_distance(t1min, t2min, dr2, n1, n2);

      if (dr2 < cross_thresh)
      {

        //  Identifies closest beads.
        // bead_n1 = n1 * settings->num_beads_per_fiber + int(t1min * settings->interbead_distance * (fibers[n1]->num_beads - 1) + 0.5);
        // bead_n2 = n2 * settings->num_beads_per_fiber + int(t2min * settings->interbead_distance * (fibers[n2]->num_beads - 1) + 0.5);

        bead_n1 = n1 * settings->num_beads_per_fiber + int(t1min * (fibers[n1]->num_beads - 1) + 0.5);
        bead_n2 = n2 * settings->num_beads_per_fiber + int(t2min * (fibers[n2]->num_beads - 1) + 0.5);

        double dx = beads[bead_n1]->position[0] - beads[bead_n2]->position[0];
        double dy = beads[bead_n1]->position[1] - beads[bead_n2]->position[1];
        double dz = beads[bead_n1]->position[2] - beads[bead_n2]->position[2];
        double real_dr2 = dx * dx + dy * dy + dz * dz;

        if (real_dr2 < cross_thresh)
        {
          beads[bead_n1]->neighbours[beads[bead_n1]->num_neighbours] = bead_n2;
          beads[bead_n2]->neighbours[beads[bead_n2]->num_neighbours] = bead_n1;

          beads[bead_n1]->eq_dist[beads[bead_n1]->num_neighbours] = sqrt(real_dr2);
          beads[bead_n2]->eq_dist[beads[bead_n2]->num_neighbours] = sqrt(real_dr2);

          beads[bead_n1]->num_neighbours++;
          beads[bead_n2]->num_neighbours++;

          beads[bead_n1]->num_crosslinks++;
          beads[bead_n2]->num_crosslinks++;

          fibers[n1]->num_crosslinks++;
          fibers[n2]->num_crosslinks++;

          num_crosslinks++;
        }
      }
    }
  }
}

void System::deactivate_beads()
{
  bool deactivate;

  for (int i = 0; i < num_beads; i++) // Loop over all beads of the system
  {
    if (beads[i]->is_active)
    {
      if (beads[i]->cluster != largest_cluster)
      {
        for (int k = 0; k < beads[i]->num_neighbours; k++)
        {
          beads[beads[i]->neighbours[k]]->num_active_neighbours--;
        }
        beads[i]->is_active = false;
      }
    }
  }

  int num_disconnected_beads = 1;
  while (num_disconnected_beads > 0)
  {                                     // Check if the number of disconnected beads in the previous iteration is greater than 0
    num_disconnected_beads = 0;         // Initialize the number of disconnected beads of the current iteration
    for (int i = 0; i < num_beads; i++) // Loop over all beads of the system
    {
      if (beads[i]->is_active)
      {
        deactivate = false; // By default, assume it should not be disconnected

        if (beads[i]->is_mobile) // If the bead is mobile (i.e. in the simulation box)
        {
          if (beads[i]->num_active_neighbours == 1) // Disconnect if it has only 1 neighbour
          {
            deactivate = true;
          }
        }
        else
        {                    // If the bead is not mobile
          deactivate = true; // First, assume it should be disconnected.
          for (int k = 0; k < beads[i]->num_neighbours; k++)
          { // However, if it has one mobile and active neighbour, it should NOT be disconnected
            if (beads[beads[i]->neighbours[k]]->is_mobile & beads[beads[i]->neighbours[k]]->is_active)
            {
              deactivate = false;
            }
          }
        }
        if (deactivate)
        {
          for (int k = 0; k < beads[i]->num_neighbours; k++)
          {
            beads[beads[i]->neighbours[k]]->num_active_neighbours--;
          }
          beads[i]->is_active = false;
          num_disconnected_beads++;
        }
      }
    }
  }

  num_active_beads = 0;
  num_active_crosslinks = 0;
  for (int i = 0; i < num_beads; i++) // Loop over all beads of the system
  {
    if (beads[i]->is_active)
    {
      num_active_beads++;
      num_active_crosslinks += beads[i]->num_crosslinks;
    }
  }
  num_active_crosslinks /= 2;
}

void System ::activate_fibers()
{
  num_active_fibers = 0;
  for (int i = 0; i < num_beads; i++)
  {
    if (beads[i]->is_active)
    {
      fibers[beads[i]->fiber]->is_active = true;
    }
  }

  for (int i = 0; i < num_fibers; i++)
  {
    if (fibers[i]->is_active)
    {
      num_active_fibers++;
    }
  }
}

void System::clear_clusters()
{

  clusters = new Cluster *[num_beads]();

  for (int i = 0; i < num_beads; i++)
  {
    clusters[i] = new Cluster();
    clusters[i]->head_vertex = EMPTY;
    clusters[i]->num_vertices = 0;
  }

  for (int i = 0; i < num_beads; i++)
  {
    beads[i]->cluster = -1;
    beads[i]->next_in_cluster = EMPTY;
    beads[i]->previous_in_cluster = EMPTY;
  }
}

void System::count_clusters()
{
  unsigned int nmax = 0;
  for (int i = 0; i < num_beads; i++)
  {
    if (clusters[i]->num_vertices > nmax)
    {
      nmax = clusters[i]->num_vertices;
      largest_cluster = i;
    }
  }
}

void System::move_between_clusters(int i, int cto)
{

  long cfrom = beads[i]->cluster;

  if (beads[i]->previous_in_cluster != EMPTY)
  {
    beads[beads[i]->previous_in_cluster]->next_in_cluster = beads[i]->next_in_cluster;
  }
  else
  {
    clusters[cfrom]->head_vertex = beads[i]->next_in_cluster;
  }

  if (beads[i]->next_in_cluster != EMPTY)
  {
    beads[beads[i]->next_in_cluster]->previous_in_cluster = beads[i]->previous_in_cluster;
  }

  beads[i]->previous_in_cluster = EMPTY;
  beads[i]->next_in_cluster = clusters[cto]->head_vertex;
  if (beads[i]->next_in_cluster != EMPTY)
  {
    beads[beads[i]->next_in_cluster]->previous_in_cluster = i;
  }
  clusters[cto]->head_vertex = i;
  beads[i]->cluster = cto;

  clusters[cfrom]->num_vertices--;
  clusters[cto]->num_vertices++;
}

void System::insert_to_cluster(int c, int i)
{

  beads[i]->previous_in_cluster = EMPTY;
  beads[i]->next_in_cluster = clusters[c]->head_vertex;
  if (beads[i]->next_in_cluster != EMPTY)
  {
    beads[beads[i]->next_in_cluster]->previous_in_cluster = i;
  }
  clusters[c]->head_vertex = i;
  beads[i]->cluster = c;
  clusters[c]->num_vertices++;
}

void System::find_clusters()
{

  int curr_cluster;
  int curr_bead;
  int c = 0;

  clear_clusters();

  for (int i = 0; i < num_beads; i++)
  {
    if (beads[i]->is_active)
    {
      insert_to_cluster(c, i);
      c++;
    }
  }

  for (int i = 0; i < num_beads; i++)
  {
    if (beads[i]->is_active)
    {
      for (int k = 0; k < beads[i]->num_neighbours; k++)
      {
        if (beads[beads[i]->neighbours[k]]->is_active)
        {
          curr_cluster = beads[beads[i]->neighbours[k]]->cluster; // Identify the cluster C of the current neighbour N of B
          if (beads[i]->cluster != curr_cluster)
          {
            curr_bead = clusters[curr_cluster]->head_vertex;
            while (curr_bead != EMPTY)
            {
              move_between_clusters(curr_bead, beads[i]->cluster);
              curr_bead = clusters[curr_cluster]->head_vertex;
            }
          }
        }
      }
    }
  }

  count_clusters();
}

void System::relabel_beads()
{

  int index;

  int *new_id, *old_id;
  int *new_neighbours;

  new_id = new int[num_beads]();
  old_id = new int[num_active_beads]();

  // First, copy active beads into a new array, change their id accordingly and create an arraw that mapy the original id to the new one.

  for (int i = 0; i < num_beads; i++)
  {
    new_id[i] = -1;
  }

  index = 0;
  for (int i = 0; i < num_active_beads; i++)
  {
    while (!beads[index]->is_active)
    {
      index++;
    }
    beads[i] = beads[index];
    new_id[index] = i;
    old_id[i] = index;
    index++;
  }

  // Then, change the id of neighbours to the new ones.
  for (int i = 0; i < num_active_beads; i++)
  {
    index = 0;
    for (int k = 0; k < beads[old_id[i]]->num_active_neighbours; k++)
    {
      while (new_id[beads[old_id[i]]->neighbours[index]] == -1)
      {
        index++;
      }
      beads[i]->neighbours[k] = new_id[beads[old_id[i]]->neighbours[index]];
      index++;
    }
    beads[i]->num_neighbours = beads[old_id[i]]->num_active_neighbours;
    beads[i]->id = i;
  }

  num_beads = num_active_beads;

  delete[] new_id;
  delete[] old_id;
}

void System::define_nodes()
{

  int j, num_beads_in_node, jj, m, new_beads;

  int *list_of_beads;
  list_of_beads = new int[num_beads];

  for (int i = 0; i < num_beads; i++)
  {
    beads[i]->checked = false;
    beads[i]->crosslink_node = -1;
  }

  nodes = new Node *[num_crosslinks]();
  num_nodes = 0;
  for (int i = 0; i < num_beads; i++)
  {
    if (beads[i]->num_crosslinks > 0 & !beads[i]->checked & beads[i]->is_active) // beads[i]->nc >o means the bead has a certain number
    // of crosslinks that means I am encountring a new node
    {
      // Find first bead of node

      num_beads_in_node = 0;
      for (int i = 0; i < num_beads; i++)
      {
        list_of_beads[i] = -1;
      }
      list_of_beads[0] = i;
      beads[i]->checked = true;
      new_beads = 1;

      // Add new beads to the node

      while (new_beads > 0)
      {
        m = num_beads_in_node + new_beads;
        new_beads = 0;
        for (int k = num_beads_in_node; k < m; k++)
        {
          j = list_of_beads[k];
          if (beads[j]->is_active)
          {
            for (int k = 0; k < beads[j]->num_neighbours; k++)
            {
              jj = beads[j]->neighbours[k];
              if (beads[j]->fiber != beads[jj]->fiber & !beads[jj]->checked & beads[jj]->is_active)
              {
                list_of_beads[m + new_beads] = jj;
                beads[jj]->checked = true;
                new_beads++;
              }
            }
          }
        }
        num_beads_in_node = m;
      }

      // Create the node and define its properties

      nodes[num_nodes] = new Node();
      nodes[num_nodes]->id = num_nodes;
      nodes[num_nodes]->num_beads = num_beads_in_node;

      nodes[num_nodes]->beads = new int[num_beads_in_node];
      nodes[num_nodes]->mass = 0;

      for (int k = 0; k < num_beads_in_node; k++)
      {
        j = list_of_beads[k];
        nodes[num_nodes]->beads[k] = j;
        beads[j]->crosslink_node = num_nodes;
        nodes[num_nodes]->mass += beads[j]->mass;
      }
      num_nodes++;
    }
  }

  delete[] list_of_beads;

  connect_nodes();
  // list of neighbours of a specific node
}

void System::connect_nodes()
{

  int current_bead, new_current_bead;
  int neighbouring_node, num_neighbours;
  bool broken_loop;
  int i, j, jj;

  int *list_of_neighbours;
  list_of_neighbours = new int[num_nodes];

  for (int n = 0; n < num_nodes; n++)
  {

    for (int i = 0; i < num_nodes; i++)
    {
      list_of_neighbours[i] = -1;
    }
    num_neighbours = 0;

    for (int l = 0; l < nodes[n]->num_beads; l++)
    {
      i = nodes[n]->beads[l];

      for (int k = 0; k < beads[i]->num_neighbours; k++)
      {
        j = beads[i]->neighbours[k];
        if (beads[j]->is_active & beads[i]->fiber == beads[j]->fiber)
        {

          for (int i = 0; i < num_beads; i++)
          {
            beads[i]->checked = false;
          }

          current_bead = j;
          beads[i]->checked = true;
          beads[j]->checked = true;

          broken_loop = false;
          while (beads[current_bead]->num_crosslinks == 0)
          { /////////////hhhhhhhhh///////////////
            // cout << n << " " << i << " " << j << " " << current_bead << " " << beads[current_bead]->num_crosslinks << " " << beads[current_bead]->num_neighbours << endl;
            if (!beads[current_bead]->is_mobile)
            {
              broken_loop = true;
              break;
            }
            for (int kk = 0; kk < beads[current_bead]->num_neighbours; kk++)
            {
              jj = beads[current_bead]->neighbours[kk];
              // cout << jj << endl;
              if (beads[jj]->is_active & beads[current_bead]->fiber == beads[jj]->fiber & !beads[jj]->checked)
              {
                new_current_bead = jj;
                beads[jj]->checked = true;
              }
            }
            current_bead = new_current_bead;
          }

          if (broken_loop)
          {
            neighbouring_node = -1;
          }
          else
          {
            list_of_neighbours[num_neighbours] = beads[current_bead]->crosslink_node;
            num_neighbours++;
          }
        }
      }
    }

    nodes[n]->neighbours = new int[num_neighbours];
    nodes[n]->num_neighbours = num_neighbours;
    // position nieghbours in ascending order
    for (int i = 0; i < num_neighbours; i++)
    {
      nodes[n]->neighbours[i] = list_of_neighbours[i];

      // cout << n << " " << i << " " << list_of_neighbours[i] << endl;
      // cout << " Node " << n << " has " << num_neighbours << " neighbours " << list_of_neighbours[i] << endl;
    }
    // cout << nodes[n]->position[0] << " " << nodes[n]->position[1] << " " << nodes[n]->position[2] << endl;
  }

  delete[] list_of_neighbours;
}

void System::reset_forces()
{
  potential_energy_stretching = 0;
  potential_energy_bending = 0;

  for (int i = 0; i < num_beads; i++)
  {

    beads[i]->force[0] = 0;
    beads[i]->force[1] = 0;
    beads[i]->force[2] = 0;
  }
}

void System::compute_forces()
{
  reset_forces();
  compute_stretching_forces();
  compute_bending_forces();
  compute_langevin_forces();


}

void System::compute_stretching_forces()
{
  int j;
  double dx, dy, dz, dr2, dr, delta;
  double radial_force;
  ofstream myfile;
  ofstream myfile1;
  myfile.open("/home/hijazi/Documents/Collagen_Fiber_Network/code_static/id_same.txt");
  myfile1.open("/home/hijazi/Documents/Collagen_Fiber_Network/code_static/id_cross.txt");

  for (int i = 0; i < num_beads; i++)
  {

    if (beads[i]->is_active && beads[i]->is_mobile)
    {

      for (int k = 0; k < beads[i]->num_neighbours; k++)
      {

        j = beads[i]->neighbours[k];
        if (beads[j]->is_active && beads[j]->is_mobile)
        {
          if (beads[i]->id < beads[j]->id)
          {
            dx = beads[j]->position[0] - beads[i]->position[0];
            dy = beads[j]->position[1] - beads[i]->position[1];
            dz = beads[j]->position[2] - beads[i]->position[2];
            dr2 = dx * dx + dy * dy + dz * dz;
            dr = sqrt(dr2);

            delta = dr - beads[i]->eq_dist[k];

            if (beads[i]->fiber == beads[j]->fiber)
            {
             // myfile << beads[i]->id << " " << beads[j]->id << endl;
              radial_force = settings->kf * delta / dr;
              potential_energy_stretching += 0.5 * settings->kf * delta * delta;
              myfile<<delta<<" "<<dr<<" "<<beads[i]->eq_dist[k]<<" " <<endl;
            }
            else
            {
             // myfile1 << beads[i]->id << " " << beads[j]->id << endl;
              radial_force = settings->kc * delta / dr;
              if (0.5 * settings->kc * delta * delta > 0.1)
              {
                // cout << "CROSSLINK " << steps << " " << i << " " << j << " " << 0.5 * settings->kc * delta * delta << " " << dr << " " << beads[i]->eq_dist[k] << " " << radial_force << endl;
              }
              potential_energy_stretching += 0.5 * settings->kc * delta * delta;
              myfile1<<delta<<" "<<dr<<" "<<beads[i]->eq_dist[k]<<" " <<endl;
            }

            beads[i]->force[0] += radial_force * dx;
            beads[i]->force[1] += radial_force * dy;
            beads[i]->force[2] += radial_force * dz;

            beads[j]->force[0] -= radial_force * dx;
            beads[j]->force[1] -= radial_force * dy;
            beads[j]->force[2] -= radial_force * dz;
          }
        }
      }
    }
  }

  myfile.close();
  myfile1.close();
}

void System::compute_bending_forces()
{
  int j1, j2;
  double t1x, t1y, t1z, t2x, t2y, t2z;
  double norm1, norm2;
  double scalar_product;
  double t, cos_theta, theta, delta_tetha;
  double bending_force;
  double F1x, F1y, F1z, F2x, F2y, F2z, F3x, F3y, F3z;

  for (int i; i < num_beads; i++)
  {
    if (beads[i]->is_active && beads[i]->is_mobile)
    {
      for (int k1 = 0; k1 < beads[i]->num_neighbours; k1++)
      {
        j1 = beads[i]->neighbours[k1];
        if (beads[j1]->is_active && beads[j1]->is_mobile)
        {
          // cout<<beads[i]->fiber<<beads[j1]->fiber<<endl;
          if (beads[i]->fiber == beads[j1]->fiber)
          {

            t1y = beads[j1]->position[1] - beads[i]->position[1];
            t1z = beads[j1]->position[2] - beads[i]->position[2];

            t1x = beads[i]->position[0] - beads[j1]->position[0];
            t1y = beads[i]->position[1] - beads[j1]->position[1];
            t1z = beads[i]->position[2] - beads[j1]->position[2];

            norm1 = sqrt(t1x * t1x + t1y * t1y + t1z * t1z);

            for (int k2 = 0; k2 < k1; k2++)
            {
              j2 = beads[i]->neighbours[k2];
              if (beads[j2]->is_active && beads[j2]->is_mobile)
              {
                if (beads[i]->fiber == beads[j2]->fiber)
                {
                  t2x = beads[j2]->position[0] - beads[i]->position[0];
                  t2y = beads[j2]->position[1] - beads[i]->position[1];
                  t2z = beads[j2]->position[2] - beads[i]->position[2];

                  norm2 = sqrt(t2x * t2x + t2y * t2y + t2z * t2z);
                  // scalar_product=inner_product(t1x, t1y, t1z, t2x, t2y, t2z);
                  t = t1x * t2x + t1y * t2y + t1z * t2z;
                  // cout<<"t2x"<<t2x<<endl;
                  cos_theta = t / (norm1 * norm2);

                  theta = acos(cos_theta);

                  // delta_tetha = (cos_theta - cos_theta_0) * (cos_theta - cos_theta_0)

                  potential_energy_bending += settings->ko * (1 - cos_theta);

                  F1x = (settings->ko / (norm1)) * (-t2x / (norm2) + t1x / (norm1) * (t / (norm1 * norm2)));
                  F3x = (settings->ko / (norm2)) * (t1x / (norm1)-t2x / (norm2) * (t / (norm1 * norm2)));
                  F2x = -F1x - F3x;

                  F1y = (settings->ko / (norm1)) * (-t2y / (norm2) + t1y / (norm1) * (t / (norm1 * norm2)));
                  F3y = (settings->ko / (norm2)) * (t1y / (norm1)-t2y / (norm2) * (t / (norm1 * norm2)));
                  F2y = -F1y - F3y;

                  F1z = (settings->ko / (norm1)) * (-t2z / (norm2) + t1z / (norm1) * (t / (norm1 * norm2)));
                  F3z = (settings->ko / (norm2)) * (t1z / (norm1)-t2z / (norm2) * (t / (norm1 * norm2)));
                  F2z = -F1z - F3z;

                  beads[j1]->force[0] += F1x;
                  beads[j1]->force[1] += F1y;
                  beads[j1]->force[2] += F1z;

                  beads[j2]->force[0] += F3x;
                  beads[j2]->force[1] += F3y;
                  beads[j2]->force[2] += F3z;

                  beads[i]->force[0] += F2x;
                  beads[i]->force[1] += F2y;
                  beads[i]->force[2] += F2z;
                }
              }
            }
          }
        }
      }
    }
  }
}

void System::compute_langevin_forces()
{

  double C;

  C = sqrt(2 * settings->temperature * settings->kb * settings->gamma / settings->dt);
  for (int i = 0; i < num_beads; i++)
  {
    if (beads[i]->is_active and beads[i]->is_mobile)
    {

      beads[i]->force[0] += -settings->gamma * beads[i]->velocity[0] + C * rnd->nextGauss();
      beads[i]->force[1] += -settings->gamma * beads[i]->velocity[1] + C * rnd->nextGauss();
      beads[i]->force[2] += -settings->gamma * beads[i]->velocity[2] + C * rnd->nextGauss();
    }
  }
}

void System::update_springs()
{
  for (int i = 0; i < num_beads; i++)
  {
    if (beads[i]->is_active)
    {
      for (int k = 0; k < beads[i]->num_neighbours; k++)
      {
        int j = beads[i]->neighbours[k];
        if (beads[j]->is_active)
        {
          if (beads[i]->fiber != beads[j]->fiber)
          {
            if (beads[i]->eq_dist[k] > 0)
            {
              beads[i]->eq_dist[k] -= 0.00001;
            }
            else if (beads[i]->eq_dist[k] < 0)
            {
              beads[i]->eq_dist[k] = 0;
            }
          }
        }
      }
    }
  }
}

void Fiber::generate_fiber(Settings &settings)
{
  double theta, phi;

  rnd = new Random(-settings.seed - id); // Maps the random number generator to the global one.

  num_beads = settings.num_beads_per_fiber;         // constant number of beads per fiber
  interbead_eqlength = settings.interbead_distance; // Equilibrium distance between beads in a fiber.
  num_crosslinks = 0;
  is_active = false;
  // Define a random orientation
  theta = M_PI * rnd->nextDouble();
  phi = 2 * M_PI * rnd->nextDouble();
  u[0] = interbead_eqlength * cos(phi) * sin(theta);
  u[1] = interbead_eqlength * sin(phi) * sin(theta);
  u[2] = interbead_eqlength * cos(theta);

  // Define the position of the first bead of the fiber
  r0[0] = settings.Lx * rnd->nextDouble();
  r0[1] = settings.Ly * rnd->nextDouble();
  r0[2] = settings.Lz * rnd->nextDouble();

  length = (num_beads - 1) * interbead_eqlength; // total length of the fiber
}

void Bead::generate_bead(Settings &settings, double x, double y, double z)
{

  rnd = new Random(-settings.seed - settings.num_fibers - id); // Maps the random number generator to the global one.

  position[0] = x;
  position[1] = y;
  position[2] = z;

  // Boundary conditions : if the bead lies out of the box it is kept immobile.
  is_mobile = true;
  if (x < 0 || x > settings.Lx || y < 0 || y > settings.Ly || z < 0 || z > settings.Lz)
    is_mobile = false;

  is_active = true;

  // Initialize crosslink structure
  num_neighbours = 0;
  num_crosslinks = 0;
  neighbours = new int[max_neighbours](); // list of beads from other fibers connected to the current bead with crosslinks
  eq_dist = new double[max_neighbours](); // list of beads from other fibers connected to the current bead with crosslinks
}

void System ::velocity_initialization()
{

  kinetic_energy = 0;

  // Velocity is drawn from a Boltzmann distribution

  for (int i = 0; i < num_beads; i++)
  {
    if (beads[i]->is_active)
    {
      beads[i]->velocity[0] = sqrt(settings->kb * settings->temperature / beads[i]->mass) * rnd->nextGauss();
      beads[i]->velocity[1] = sqrt(settings->kb * settings->temperature / beads[i]->mass) * rnd->nextGauss();
      beads[i]->velocity[2] = sqrt(settings->kb * settings->temperature / beads[i]->mass) * rnd->nextGauss();

      beads[i]->previous_position[0] = beads[i]->position[0] - beads[i]->velocity[0] * settings->dt;
      beads[i]->previous_position[1] = beads[i]->position[1] - beads[i]->velocity[1] * settings->dt;
      beads[i]->previous_position[2] = beads[i]->position[2] - beads[i]->velocity[2] * settings->dt;

      kinetic_energy = kinetic_energy + 0.5 * beads[i]->mass * (beads[i]->velocity[0] * beads[i]->velocity[0] + beads[i]->velocity[1] * beads[i]->velocity[1] + beads[i]->velocity[2] * beads[i]->velocity[2]);

    }
    
  }
  cout<<kinetic_energy<<endl;
  temperature = 2 * kinetic_energy / (3 * num_active_beads * settings->kb);
  cout<<temperature<<endl;
 
}

void System ::integrate()
{
  kinetic_energy = 0;
  double dt2 = dt * dt;

  for (int i; i < num_beads; i++)
  {
    if (beads[i]->is_active & beads[i]->is_mobile)
    {
      beads[i]->new_position[0] = 2 * beads[i]->position[0] - beads[i]->previous_position[0] + dt2 * beads[i]->force[0] / beads[i]->mass;
      beads[i]->new_position[1] = 2 * beads[i]->position[1] - beads[i]->previous_position[1] + dt2 * beads[i]->force[1] / beads[i]->mass;
      beads[i]->new_position[2] = 2 * beads[i]->position[2] - beads[i]->previous_position[2] + dt2 * beads[i]->force[2] / beads[i]->mass;

      //    velocity at time t #####
      beads[i]->velocity[0] = (beads[i]->new_position[0] - beads[i]->previous_position[0]) / (2 * dt);
      beads[i]->velocity[1] = (beads[i]->new_position[1] - beads[i]->previous_position[1]) / (2 * dt);
      beads[i]->velocity[2] = (beads[i]->new_position[2] - beads[i]->previous_position[2]) / (2 * dt);

      kinetic_energy = kinetic_energy + 0.5 * beads[i]->mass * (beads[i]->velocity[0] * beads[i]->velocity[0] + beads[i]->velocity[1] * beads[i]->velocity[1] + beads[i]->velocity[2] * beads[i]->velocity[2]);

      //   update old position

      beads[i]->previous_position[0] = beads[i]->position[0];
      beads[i]->previous_position[1] = beads[i]->position[1];
      beads[i]->previous_position[2] = beads[i]->position[2];

      // update current position #####

      beads[i]->position[0] = beads[i]->new_position[0];
      beads[i]->position[1] = beads[i]->new_position[1];
      beads[i]->position[2] = beads[i]->new_position[2];
    }
  }
  temperature = 2 * kinetic_energy / (3 * num_active_beads * settings->kb);
}
