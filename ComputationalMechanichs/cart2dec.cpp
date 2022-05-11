#include "cart2dec.h"
#include <cmath>
#include <vector>
#include <iostream>

const double MU = 4e6;
const double PI = 3.14159265358979323846;
using namespace N;

inline void cart2dec::mean_movement() {
	this->n = std::sqrt(MU / (this->a * this->a * this->a));
	//std::cout << this->n << std::endl;
}

inline void cart2dec::mean_anomaly(double time) {
	this->M = this->n * (time - this->tau);
	//std::cout << this->M << std::endl;
}

inline void cart2dec::eccentric_anomaly(double eps) {
	double tmpE_0 = this->M;
	double tmpE_1 = this->M + this->e * std::sin(tmpE_0);
	//int count = 0;
	while (std::abs(tmpE_1 - tmpE_0) > eps) {// and count < 2500) {
		tmpE_0 = tmpE_1;
		tmpE_1 = this->M + this->e * std::sin(tmpE_0);
		//count++;
		//std::cout << tmpE_1 << std::endl;
	}
	this->E = tmpE_1;
	//std::cout << tmpE_1 << std::endl;
}

inline void cart2dec::true_anomaly() {
	//std::cout << 0 << std::endl;
	double k = std::sqrt((1 + this->e) / (1 - this->e));
	//std::cout << 0 << std::endl;
	this->nu = 2 * std::atan(k * std::tan(this->M / 2));
	//std::cout << 0 << std::endl;
}

inline void cart2dec::radius() {
	this->r = this->a * (1 - this->e * std::cos(this->E));
}

cart2dec::cart2dec(double a, double e, double i, double Omega, double omega, double tau) {
	this->a = a;
	this->e = e;
	this->i = i;
	this->Omega = Omega;
	this->omega = omega;
	this->tau = tau;
	cart2dec::mean_movement();
}

std::vector<double> cart2dec::coordinates_velocities(double time) {
	cart2dec::mean_anomaly(time);
	//std::cout << 4 << std::endl;
	cart2dec::eccentric_anomaly(1e-7);
	//std::cout << 1 << std::endl;
	cart2dec::true_anomaly();
	//std::cout << 2 << std::endl;
	cart2dec::radius();
	//std::cout << 3 << std::endl;
	this->u = this->omega + this->nu;

	double x = this->r * (std::cos(this->u) * std::cos(this->Omega) - std::sin(this->u) * std::sin(this->Omega) * std::cos(this->i));
	double y = this->r * (std::cos(this->u) * std::sin(this->Omega) + std::sin(this->u) * std::cos(this->Omega) * std::cos(this->i));
	double z = this->r * (std::sin(this->u) * std::sin(this->i));

	double p = this->a * (1 - this->e * this->e);
	double V_r = std::sqrt(MU / p) * this->e * std::sin(this->nu);
	double V_n = std::sqrt(MU / p) * (1 + this->e * std::cos(this->nu));

	double V_1 = V_r * std::cos(this->u) - V_n * std::sin(this->u);
	double V_2 = V_r * std::sin(this->u) + V_n * std::cos(this->u);

	double Vx = (V_1 * std::cos(this->Omega) - V_2 * std::sin(this->Omega) * std::cos(this->i));
	double Vy = (V_1 * std::sin(this->Omega) + V_2 * std::cos(this->Omega) * std::cos(this->i));
	double Vz = (V_2 * std::sin(this->i));

	std::vector<double> coordinates_velocities = { x, y, z, Vx, Vy, Vz };
	return coordinates_velocities;
}
