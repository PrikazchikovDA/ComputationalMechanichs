#pragma once
#include <vector>

namespace N 
{
	class cart2dec
	{
	private:
		double a, e, i, Omega, omega, tau;
		double n, M, E, nu, r, u;

		inline void mean_movement();

		inline void mean_anomaly(double time = 0);

		inline void eccentric_anomaly(double eps = 1e-4);

		inline void true_anomaly();

		inline void radius();
	public:
		cart2dec(double a, double e, double i, double Omega, double omega, double tau);

		std::vector<double> coordinates_velocities(double time);

	};
}

