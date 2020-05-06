/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

#define EPS 0.00001

using namespace std;
using std::string;
using std::vector;
using std::default_random_engine;
using std::normal_distribution;
static default_random_engine gen;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */

  num_particles = 1000;  // TODO: Set the number of particles
  // This line creates a normal (Gaussian) distribution for x
  normal_distribution<double> dist_x(0, std[0]);
  normal_distribution<double> dist_y(0, std[1]);
  normal_distribution<double> dist_theta(0, std[2]);

  for (int i = 0; i < num_particles; ++i) {
    //   double sample_x, sample_y, sample_theta;
    //   sample_x = dist_x(gen);
    //   sample_y = dist_y(gen);
    //   sample_theta = dist_theta(gen);   
    // }
    Particle p;

    //Initialize all particles to 
    //*   first position (based on estimates of x, y, theta and their uncertainties
    //*   from GPS) and all weights to 1. 
    p.id = i;
    p.x = x;
    p.y = y;
    p.theta = theta;
    p.weight = 1.0;

    // Add random Gaussian noise to each particle
    p.x += dist_x(gen);
    p.y += dist_y(gen);
    p.theta += dist_theta(gen);

    // Set of current particles
    particles.push_back(p);
    
	  }

  // is_initialized
  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */


  // This line creates a normal (Gaussian) distribution for x
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);

  // Apply bicycle model calculation here. 
  for (int i = 0; i < num_particles; ++i) {

    // If angle is zero
    if ( fabs(yaw_rate) == 0.0 ) { 
      particles[i].x += velocity * delta_t * cos( particles[i].theta );
      particles[i].y += velocity * delta_t * sin( particles[i].theta );
      
    } 
    // If angle is available
    else {
      particles[i].x += velocity / yaw_rate * ( sin( particles[i].theta + yaw_rate * delta_t ) - sin( particles[i].theta ) );
      particles[i].y += velocity / yaw_rate * ( cos( particles[i].theta ) - cos( particles[i].theta + yaw_rate * delta_t ) );
      particles[i].theta += yaw_rate * delta_t;
    }
  
    //add random Gaussian noise
    particles[i].x += dist_x(gen);
    particles[i].y += dist_y(gen);
    particles[i].theta += dist_theta(gen);

  }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

  //apply Nearest Neighbour
  for (unsigned int i = 0; i < observations.size(); i++) {

    // Get current observation that already transformed 
    LandmarkObs transformed_obs_coor = observations[i];

    //Must put a high value, I'll use the sensor_range which is 50
    double min_dist = 50;
    int map_id = -1;
    
    for (unsigned int j = 0; j < predicted.size(); j++) {
      // Get current prediction
      LandmarkObs prediction_coor = predicted[j];
      
      // get distance between observation and predicted landmarks
      double cur_dist = dist(transformed_obs_coor.x, transformed_obs_coor.y, prediction_coor.x, prediction_coor.y);

      // Search for closest predicted landmark
      if (cur_dist < min_dist) {
        min_dist = cur_dist;
        map_id = prediction_coor.id;
      }
    }

    // Assign observation's id to the nearest predicted landmark's id for each particles
    observations[i].id = map_id;
  }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  //Update each particle, for every particle in p, do this:
  for (int i = 0; i < num_particles; i++) {

    //get back particles gps and theta
    double x = particles[i].x;
    double y = particles[i].y;
    double theta = particles[i].theta;

    //Transform the observation coordinates from sensors for each particle.
    //Then stored in transformed_observations (observation now in MAP'S coordinate system)
    vector<LandmarkObs> transformed_observations;
    for(unsigned int j = 0; j < observations.size(); j++) {
      double xx = cos(theta)*observations[j].x - sin(theta)*observations[j].y + x;
      double yy = sin(theta)*observations[j].x + cos(theta)*observations[j].y + y;
      transformed_observations.push_back(LandmarkObs{ observations[j].id, xx, yy });
    }

    //Search for landmarks in sensor's range. (sensor is at each particle coordinate)
    //The only landmarks in search region are then stored in predictions
    vector<LandmarkObs> predictions;
    for(unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {
      float landmark_x = map_landmarks.landmark_list[j].x_f;
      float landmark_y = map_landmarks.landmark_list[j].y_f;
      int landmark_id = map_landmarks.landmark_list[j].id_i;

      if (dist(x, y, landmark_x, landmark_y) <= sensor_range) {
        predictions.push_back(LandmarkObs{ landmark_id, landmark_x, landmark_y });
      }
    }
    // Call dataAssociate
    dataAssociation(predictions, transformed_observations);

    // Reset weight.
    particles[i].weight = 1.0;

    //Search for associate obs with landmark
    for(unsigned int j = 0; j < transformed_observations.size(); j++) {
      double observation_x = transformed_observations[j].x;
      double observation_y = transformed_observations[j].y;
      int associate_id = transformed_observations[j].id;
      double associate_x, associate_y;

      for (unsigned int k = 0; k < predictions.size(); k++) {
        if (predictions[k].id == associate_id) {
          associate_x = predictions[k].x;
          associate_y = predictions[k].y;
        }
      }

      // Calculate weight using multivariate Gaussian
      double dX = observation_x - associate_x;
      double dY = observation_y - associate_y;

      double weight = ( 1/(2*M_PI*std_landmark[0]*std_landmark[1])) * exp( -( dX*dX/(2*std_landmark[1]*std_landmark[0]) + (dY*dY/(2*std_landmark[1]*std_landmark[1])) ) );
      if (weight == 0) {
        particles[i].weight *= EPS;
      } else {
        particles[i].weight *= weight;
      }
    }
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

  vector<double> weights;

  for(int i = 0; i < num_particles; i++) {
    weights.push_back(particles[i].weight);
    }

  // Normalize weights
  double sum = 0;
  for (int i = 0; i < weights.size(); i++) sum += weights[i];
  for (int i = 0; i < weights.size(); i++) weights[i] /= sum;

  // particles = resampled_particles;
  std::discrete_distribution<int> dist(weights.begin(), weights.end());
  std::vector<Particle> new_particles;
	new_particles.reserve(particles.size());

	for (int i = 0; i < particles.size(); i++) {
		int sampled_index = dist(gen);
		new_particles.push_back(particles[sampled_index]);
	}

	particles = new_particles;

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}