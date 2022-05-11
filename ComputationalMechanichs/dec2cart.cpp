#include "dec2cart.h"
#include <cmath>

using namespace C;
const double MU = 4e6;
const double PI = 3.14159265358979323846;

inline void dec2cart::orbital_momentum(){
	this->h1 = this->y * this->v_z - this->z * this->v_y;
	this->h2 = this->z * this->v_x - this->x * this->v_z;
	this->h3 = this->x * this->v_y - this->y * this->v_x;
}

inline void dec2cart::eccentricity_vector(){
	double r = std::sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
	this->e1 = (this->v_y * this->h3 - this->v_z * this->h2) / MU - this->x / r;
	this->e2 = (this->v_z * this->h1 - this->v_x * this->h3) / MU - this->y / r;
	this->e3 = (this->v_x * this->h2 - this->v_y * this->h1) / MU - this->z / r;
	this->e = std::sqrt(this->e1 * this->e1 + this->e2 * this->e2 + this->e3 * this->e3);
}

inline void dec2cart::ascending_node(){
	this->n1 = -h2;
	this->n2 = h1;
}

inline void dec2cart::true_anomaly(){
	double rdr = this->x * this->v_x + this->y * this->v_y + this->z * this->v_z;
	double r = std::sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
	if (rdr >= 0) {
		this->nu = std::acos((this->e1 * this->x + this->e2 * this->y + this->e3 * this->z) / (r * this->e));
	}
	else {
		this->nu = 2 * PI - std::acos((this->e1 * this->x + this->e2 * this->y + this->e3 * this->z) / (r * this->e));
	}
}

inline void dec2cart::eccentric_anomaly(){
	double k = std::sqrt((1 + this->e) / (1 - this->e));
	this->E = 2 * std::atan(std::tan(nu / 2) / k);
}

void dec2cart::Omega_calc() {
	if (this->n2 >= 0) {
		this->Omega = std::acos(this->n1 / std::sqrt(this->n1 * this->n1 + this->n2 * this->n2));
	}
	else {
		this->Omega = 2 * PI - std::acos(this->n1 / std::sqrt(this->n1 * this->n1 + this->n2 * this->n2));
	}
}

void dec2cart::omega_calc() {
	if (this->e3 >= 0) {
		this->omega = std::acos((this->n1 * this->e1 + this->n2 * this->e2) / (this->e * std::sqrt(this->n1 * this->n1 + this->n2 * this->n2)));
	}
	else {
		this->omega = 2 * PI - std::acos((this->n1 * this->e1 + this->n2 * this->e2) / (this->e * std::sqrt(this->n1 * this->n1 + this->n2 * this->n2)));
	}
}

dec2cart::dec2cart(double x, double y, double z, double Vx, double Vy, double Vz) {
	this->x = x;
	this->y = y;
	this->z = z;
	this->v_x = Vx;
	this->v_y = Vy;
	this->v_z = Vz;
	orbital_momentum();
	eccentricity_vector();
	ascending_node();
	true_anomaly();
	eccentric_anomaly();
	omega_calc();
	Omega_calc();
}

std::vector<double> dec2cart::keplerian_elements(double time) {
	double i = std::acos(this->h3 / std::sqrt(this->h1 * this->h1 + this->h2 * this->h2 + this->h3 * this->h3));

	double r = std::sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
	double v_sq = this->v_x * this->v_x + this->v_y * this->v_y + this->v_z * this->v_z;
	double a = 1 / (2 / r - v_sq / MU);

	double tau = time - (this->E - this->e * std::sin(this->E)) / std::sqrt(MU / (a * a * a));

	std::vector<double> cont = { a, this->e, i, this->Omega, this->omega, tau };
	return cont;
}
