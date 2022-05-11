#pragma once
#include <vector>
namespace C {
	class dec2cart
	{
	private:
		double x, y, z, v_x, v_y, v_z;
		double h1, h2, h3;
		double e1, e2, e3;
		double n1, n2;
		double nu, e, E;
		double omega, Omega;

		inline void orbital_momentum();

		inline void eccentricity_vector();

		inline void ascending_node();

		inline void true_anomaly();

		inline void eccentric_anomaly();

		inline void Omega_calc();

		inline void omega_calc();

	public:

		dec2cart(double x, double y, double z, double Vx, double Vy, double Vz);

		std::vector<double> keplerian_elements(double time);
	};
}