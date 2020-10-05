#include<iostream>
#include<fstream>
#include<vector>

constexpr size_t T_before = 1000;
constexpr size_t M = 1000000;
constexpr size_t T = 1;
constexpr long double dt = 0.001;

std::vector<long double> f(std::vector<long double> coor)
{
	long double x = coor[0];
	long double y = coor[1];
	long double z = coor[2];
	std::vector<long double> result = {
		10 * (y - x),
		28 * x - y - x * z,
		-8 / 3 * z + x * y
	};
	return result;
}

std::vector<long double> f_variation(std::vector<long double> coor, std::vector<long double> variation)
{
	long double x = coor[0];
	long double y = coor[1];
	long double z = coor[2];
	long double x_ = variation[0];
	long double y_ = variation[1];
	long double z_ = variation[2];
	std::vector<long double> result = {
		-10 * x_ + 10 * y_,
		(28 - z) * x_ - y_ - x * z_,
		-8 / 3 * z_ + x * y_ + x_ * y
	};
	return result;
}

void CountNextCoor(std::vector<long double>& coor)
{
	std::vector<long double> k1, k2, k3, k4, current_coor;
	current_coor = coor;
	k1 = f(current_coor);
	for (size_t i = 0; i < current_coor.size(); i++)
	{
		current_coor[i] = coor[i] + k1[i] * dt / 2;
	}
	k2 = f(current_coor);
	for (size_t i = 0; i < current_coor.size(); i++)
	{
		current_coor[i] = coor[i] + k2[i] * dt / 2;
	}
	k3 = f(current_coor);
	for (size_t i = 0; i < current_coor.size(); i++)
	{
		current_coor[i] = coor[i] + k3[i] * dt;
	}
	k4 = f(current_coor);
	for (size_t i = 0; i < coor.size(); i++)
	{
		coor[i]+= dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
	}
}

void CountNextVariation(std::vector<long double> coor, std::vector<long double>& variation)
{
	std::vector<long double> k1, k2, k3, k4, current_variation;
	current_variation = variation;
	k1 = f_variation(coor, current_variation);
	for (size_t i = 0; i < current_variation.size(); i++)
	{
		current_variation[i] = variation[i] + k1[i] * dt / 2;
	}
	k2 = f_variation(coor, current_variation);
	for (size_t i = 0; i < current_variation.size(); i++)
	{
		current_variation[i] = variation[i] + k2[i] * dt / 2;
	}
	k3 = f_variation(coor, current_variation);
	for (size_t i = 0; i < current_variation.size(); i++)
	{
		current_variation[i] = variation[i] + k3[i] * dt;
	}
	k4 = f_variation(coor, current_variation);
	for (size_t i = 0; i < variation.size(); i++)
	{
		variation[i] += dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
	}
}


int main()
{
	//std::ofstream fout;
	//fout.open("../wwwroot/output/test/result.csv");//Введите свой путь
	std::vector<long double> coor = { 0.1, 0.1, 0.1 };
	for (size_t i = 0; i < T_before; i++)
		CountNextCoor(coor);

	std::vector<long double> variation_1 = { 1, 0, 0 };
	std::vector<long double> variation_2 = { 0, 1, 0 };
	std::vector<long double> variation_3 = { 0, 0, 1 };
	long double log_sum_1 = 0;
	long double log_sum_2 = 0;
	long double log_sum_3 = 0;

	for (size_t i = 0; i < M; i++)
	{
		for (size_t j = 0; j < T; j++)
		{
			CountNextVariation(coor, variation_1);
			CountNextVariation(coor, variation_2);
			CountNextVariation(coor, variation_3);
			CountNextCoor(coor);
		}
		long double norm_1 = sqrtl(variation_1[0] * variation_1[0] + variation_1[1] * variation_1[1] + variation_1[2] * variation_1[2]);
		log_sum_1 += logl(norm_1);
		//Нормировка по вектору скорости
		/*variation_1 = { 10*(coor[1]-coor[0]), 28*coor[0] - coor[1] - coor[0]*coor[2], -8/3*coor[2]+coor[0]*coor[1] };
		norm_1 = sqrtl(variation_1[0] * variation_1[0] + variation_1[1] * variation_1[1] + variation_1[2] * variation_1[2]);*/
		for (auto& el : variation_1)
			el /= norm_1;
		/*
		long double dot_21 = variation_1[0] * variation_2[0] + variation_1[1] * variation_2[1] + variation_1[2] * variation_2[2];
		for (size_t i = 0; i < 3; i++)
			variation_2[i] -= variation_1[i] * dot_21;
		long double norm_2 = sqrtl(variation_2[0] * variation_2[0] + variation_2[1] * variation_2[1] + variation_2[2] * variation_2[2]);
		log_sum_2 += logl(norm_2);
		for (auto& el : variation_2)
			el /= norm_2;
		long double dot_31 = variation_3[0] * variation_1[0] + variation_3[1] * variation_1[1] + variation_3[2] * variation_1[2];
		long double dot_32 = variation_3[0] * variation_2[0] + variation_3[1] * variation_2[1] + variation_3[2] * variation_2[2];
		for (size_t i = 0; i < 3; i++)
			variation_3[i] -= variation_1[i] * dot_31 + variation_2[i] * dot_32;
		long double norm_3 = sqrtl(variation_3[0] * variation_3[0] + variation_3[1] * variation_3[1] + variation_3[2] * variation_3[2]);
		log_sum_3 += logl(norm_3);
		for (auto& el : variation_3)
			el /= norm_3;
		*/
	}
	std::cout << "Lyapunov exponents:" << log_sum_1 / M / T / dt << " | " << log_sum_2 / M / T / dt << " | " << log_sum_3 / M / T / dt;
}