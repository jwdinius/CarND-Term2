/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <math.h> 

// additional include for updateWeights
#include <map>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	default_random_engine gen;

	// this was the minimum number I found gave decent performance and wasn't terribly slow
	num_particles = 60;

	// Pre-allocate now that we know the right size
	particles.resize(num_particles);
	weights.resize(num_particles);

	// setup normal distributions
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_th(theta, std[2]);

	// initialize the particles and weights
	for (size_t i = 0; i < num_particles; ++i) {		
		particles[i].id     = i;
		particles[i].x      = dist_x(gen);
		particles[i].y      = dist_y(gen);
		particles[i].theta  = dist_th(gen);
		particles[i].weight = 1.;

		weights[i] = 1.;

	}
	
	is_initialized = true;

	return;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;

	// setup normal distributions (for adding noise of zero mean to the motion process below)
	normal_distribution<double> dist_x( 0., std_pos[0]);
	normal_distribution<double> dist_y( 0., std_pos[1]);
	normal_distribution<double> dist_th(0., std_pos[2]);

	// loop over particles...
	for (size_t i = 0; i < particles.size(); i++) {
		// if yaw rate is not too small, use the full EOM
		if (fabs(yaw_rate) > .001) {
			particles[i].x += velocity / yaw_rate * ( sin(particles[i].theta + yaw_rate * delta_t)
				                                    - sin(particles[i].theta) );
			particles[i].y += velocity / yaw_rate * ( cos(particles[i].theta) 
				                                    - cos(particles[i].theta + yaw_rate * delta_t));
			particles[i].theta += yaw_rate * delta_t;
		}
		// otherwise use the simplification
		else {
			particles[i].x += velocity * cos(particles[i].theta) * delta_t;
			particles[i].y += velocity * sin(particles[i].theta) * delta_t;
		}

		// add noise to the particle based on process model
		particles[i].x     += dist_x(gen);
		particles[i].y     += dist_y(gen);
		particles[i].theta += dist_th(gen);

	}

	return;

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	double d_min, d;
	int id_min;
	// loop over all observations
	for (size_t i = 0; i < observations.size(); i++) {
		// (re)initialize distance to something large and index to nonsense (it'll get updated)
		d_min  = 1.e5;
		id_min = -1;

		// loop over all predictions to find the closest (nearest neighbor) for association
		for (size_t j = 0; j < predicted.size(); j++) {
			d = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
			// if distance is smaller than our current smallest, update
			if (d < d_min) {
				d_min  = d;
				id_min = predicted[j].id;
			}
		}

		// now assign the identity of the best prediction for that observation
		observations[i].id = id_min;
	}

	return;
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	
	// Create a std::map to easily reference Map::map_landmarks by its id_i field (notation is a bit
	// confusing here)
	std::map<int, Map::single_landmark_s> landmarks_id_map;

	for (unsigned int i = 0; i < map_landmarks.landmark_list.size(); i++) {

		landmarks_id_map.insert(std::make_pair(map_landmarks.landmark_list[i].id_i,
		                                       map_landmarks.landmark_list[i]));
	}

	double ct, st, xm, ym, d;
	double x, y, mux, muy, md, sx, sy, p;
	
	// loop over particles...
	for (unsigned int i = 0; i < num_particles; i++) {

		// temp variables to clean up local(vehicle)-to-global(map) transformation
		ct = cos(particles[i].theta);
		st = sin(particles[i].theta);

		// local-to-global coordinate transformation (loop over all observations)
		std::vector<LandmarkObs> obs_global;
		obs_global.resize(observations.size());
		for (unsigned j = 0; j < observations.size(); j++) {

			xm = particles[i].x + observations[j].x * ct - observations[j].y * st; 
			ym = particles[i].y + observations[j].x * st + observations[j].y * ct;
			obs_global[j] = LandmarkObs{(int)j, xm, ym};	

		}

		// Find landmarks that are within sensor_range
		std::vector<LandmarkObs> marks_in_range;
		for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {

			d  = dist(particles[i].x, particles[i].y,
			          map_landmarks.landmark_list[j].x_f,
					  map_landmarks.landmark_list[j].y_f);

			if (d < sensor_range) {
				marks_in_range.push_back(LandmarkObs{map_landmarks.landmark_list[j].id_i,
				                                     map_landmarks.landmark_list[j].x_f,
												     map_landmarks.landmark_list[j].y_f });
			}

		}

		// Calculate weights based upon marks within range
		if (marks_in_range.size() > 0) {

			// Associate each measurement with the nearest neighbor landmark that is in range
			dataAssociation(marks_in_range, obs_global);

			// reset the particle weight
			particles[i].weight = 1.;

			// this is just a trial for c++ niceness, usually I would do something like:
			// for (size_t i = 0; i < obs_global.size(); i++)
			for (const auto obs:obs_global) {
				x = landmarks_id_map[obs.id].x_f;
				y = landmarks_id_map[obs.id].y_f;
				mux = obs.x;
				muy = obs.y;
				sx = std_landmark[0];
				sy = std_landmark[1];
				// compute mahalonobis distance
				md = dist(x/sx, y/sy, mux/sx, muy/sy);
				// compute likelihood
				p  = exp(-md/2.);
				p /= 2. * M_PI * sx * sy;
				
				particles[i].weight *= p;

			}

			// assign weight based upon likelihood calculated above
			weights[i] = particles[i].weight;

		} 
		else {
			// if no landmarks were observed for the current particle, set its weight to 0
			weights[i] = 0.0;

		}

	}

	return;
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	// define parameters needed for random draws
	std::default_random_engine gen;
	std::discrete_distribution<int> dist(weights.begin(), weights.end());

	// create new container for resampled particles (and reserve a size for it)
	std::vector<Particle> new_particles;
	new_particles.reserve(particles.size());

	// resample baesd on weights (see definition of dist above)
	for (int i = 0; i < particles.size(); i++) {
		int sampled_index = dist(gen);
		new_particles.push_back(particles[sampled_index]);
	}

	// reset particle container to be the resampled particle container
	particles = new_particles;

	return;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
