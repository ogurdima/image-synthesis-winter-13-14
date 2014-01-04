#pragma once
#include <vector>
using std::vector;

class Chi2Inv
{
public:
	static double chi2inv(double p, int v);

private:
	static vector<vector<double>> _chi2inv;
	Chi2Inv();
};

