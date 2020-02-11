#include "jkr_lookup_table.h"

void LUT::init(double d_delta_delta_c){
	delta_delta_c[0] = -1.;
	a_a0[0] = calculate_a_a0(delta_delta_c[0]);
	fn_fc[0] = calculate_fn_fc(a_a0[0]);
	double alpha = 1.;
	for(int i=0;i<points-1;++i) {
		//If next delta_delta_c is over the high resolution region increase d_delta_delta_c by factor alpha_l
		if(delta_delta_c[i]+alpha*d_delta_delta_c>crit_lim){
			alpha = alpha_l;
		}
		
		delta_delta_c[i+1] = delta_delta_c[i]+alpha*d_delta_delta_c;

		if ((delta_delta_c[i+1] > 1.) & (delta_delta_c[i] < 1.))
		{
			delta_delta_c[i+1] = 1.;
		}

		a_a0[i+1] = calculate_a_a0(delta_delta_c[i+1]);
		fn_fc[i+1] = calculate_fn_fc(a_a0[i+1]);
		//Adding extra points for neg delta_delta_c values due to large changes in fn_fc
		if (delta_delta_c[i+1] < 0.)
		{
			while(fabs(fn_fc[i+1]-fn_fc[i]) > d_delta_delta_c*fn_fc_tol_factor)
			{
				delta_delta_c[i+1] = 0.5*(delta_delta_c[i+1]+delta_delta_c[i]);
				a_a0[i+1] = calculate_a_a0(delta_delta_c[i+1]);
				fn_fc[i+1] = calculate_fn_fc(a_a0[i+1]);
			}
		}
	}
	delta_delta_c_max = delta_delta_c[points-1];
}

std::pair<int,int> LUT::get_closest_indexes(double delta_delta_c_val)
{
		std::pair<int,int> ind = std::make_pair(0,points-1);
		int mid_old = 0;
		int mid;
		while(true)
		{
				mid = (ind.first+ind.second)/2;
				if (mid_old == mid)
				{
						break;
				}
				if (delta_delta_c_val < delta_delta_c[mid])
				{
						ind.second = mid;
				}
				else
				{
						ind.first = mid;
				}
				mid_old = mid;
		}
		return ind;
}

double LUT::get_a_a0(double delta_delta_c_val)
{
	double a_a0_val;

	//If point is larger than delta_max then calculate value (Slower performance)
	if ( delta_delta_c_val > delta_delta_c_max )
	{         
		//add point to delta
		a_a0_val = calculate_a_a0(delta_delta_c_val);
		if (outside_LUT == false)
		{
			std::cout << "Simulation outside LUT table, increase the resolution for speed! delta_delta_c = " 
								<< delta_delta_c_val << std::endl;
			outside_LUT = true;
		}
	}
	else if (delta_delta_c_val < -1.)
	{
		a_a0_val = 0.;
	}
	else
	{
		std::pair<int,int> ind = get_closest_indexes(delta_delta_c_val);
		//Linear interpolation
		double a,b,f;
		a = a_a0[ind.first];
		b = a_a0[ind.second];
		f = (delta_delta_c_val - delta_delta_c[ind.first])/(delta_delta_c[ind.second]-delta_delta_c[ind.first]);
		a_a0_val = (a * (1.0f - f)) + (b * f);
	}
	return a_a0_val;
}
double LUT::get_fn_fc(double delta_delta_c_val)
{
	double fn_fc_val;
	double a_a0_val;

	//If point is larger than delta_max then calculate value (Slower performance)
	if ( delta_delta_c_max < delta_delta_c_val )
	{         
		//add point to delta
		a_a0_val = calculate_a_a0(delta_delta_c_val);
		fn_fc_val  = calculate_fn_fc(a_a0_val);
		
		if (outside_LUT == false)
		{
			std::cout << "Simulation outside LUT table, increase resolution for speed! delta_delta_c = " 
								<< delta_delta_c_val << std::endl;
			outside_LUT = true;
		}
	}
	else if (delta_delta_c_val < -1.)
	{
		fn_fc_val = 0.;
	}
	else
	{
		std::pair<int,int> ind = get_closest_indexes(delta_delta_c_val);
		//Linear interpolation
		double a,b,f;
		a = fn_fc[ind.first];
		b = fn_fc[ind.second];
		f = (delta_delta_c_val - delta_delta_c[ind.first])/(delta_delta_c[ind.second]-delta_delta_c[ind.first]);
		fn_fc_val = (a * (1.0f - f)) + (b * f);
	}
	return fn_fc_val;
}
//Returns Fh/Fc where Fh is apperent hertz force according to thornton 2015 eq 3.76
double LUT::get_fh_fc(double delta_delta_c_val)
{
	double fn_fc = get_fn_fc(delta_delta_c_val);
	return calculate_fh_fc(fn_fc);
}
double LUT::calculate_fh_fc(double fn_fc)
{
	return (fn_fc+2.+pow(4.*fn_fc+4.,0.5));
}
double LUT::calculate_a0(double Wij,double reff,double Yeff)
{
	return pow(9./2.*Wij*M_PI*pow(reff,2.)/(Yeff),1./3.);
}
double LUT::calculate_delta_c(double a_0,double reff)
{
	return pow(a_0,2.)/(2.*pow(6,1./3.)*reff);
}
double LUT::calculate_fc(double Wij,double reff)
{
	return 3./2.*M_PI*Wij*reff;
}

//writing out to csv file delta_delta_c,a_a0,fn_fc
void LUT::write_csv_table()
{
	std::ofstream outfile;
	outfile.open("LUT_TABLE.csv");

	for(int i=0;i<points;++i)
	{
		outfile << std::setprecision(24) << std::fixed << delta_delta_c[i] 
						<< "," << a_a0[i] << "," << fn_fc[i] << std::endl;
	}
}

//prints out table info (points) and max delta_delta_c
void LUT::print_table_info()
{
	std::cout << "LUT table has " << points << " points, with max(delta_delta_c) = " 
						<< delta_delta_c_max << std::endl;
}

//prints out full table 
void LUT::print_table()
{
	std::cout << "delta_delta_c\t\ta_a0\t\tfn_fc" << std::endl;
	for(int i=0;i<points;++i) {
		std::cout << delta_delta_c[i] << "\t\t" << a_a0[i] << "\t\t" << fn_fc[i] << std::endl;
	}
	std::cout << std::endl;
}

double LUT::get_max_delta_delta_c()
{
	return delta_delta_c_max;
}
int LUT::get_nr_of_points()
{
	return points;
}

//If input is a vector loop over vector and calculate for each value a_a0
void calculate_a_a0(std::vector<double> &delta_delta_c,std::vector<double> &a_a0)
{
	std::vector<double> a_a0_array(delta_delta_c.size());
	for (unsigned int i = 0;i < delta_delta_c.size();i = i + 1)
	{
		a_a0[i] = calculate_a_a0(delta_delta_c[i]);
		std::cout << "Calculating a_a0 not using lookup!" << std::endl;
	}
}

double calculate_a_a0(double delta_delta_c)
{
	int Degree = 4;
	double a_a0 = 0.0;
	double op[Degree+1], zeroi[Degree], zeror[Degree]; // Coefficient vectors

	//Want to solve 0 = A*x^2+Bx^0.5+C but substitute y=x^0.5 and get 0 = A*y^4+B*y+C
	op[0] = 2.*pow(6., 1./3.);      // A*y^4
	op[1] = 0.;                                       // 0*y^3     
	op[2] = 0.;                                       // 0*y^2
	op[3] = -4./3.*pow(6., 1./3.);  // B*y^1
	op[4] = -1*delta_delta_c;                                 // C*y^0

	rpoly_ak1(op, &Degree, zeror, zeroi);

	for (int rt = 0; rt < 4; rt++) 
	{
		//Root must be larger than zero and for negeative values of deltan it should be the last root.
		if (zeror[rt] >= 0 && zeroi[rt] == 0 )
		{
		  //Only allow solution y>0 and imag(y) == 0 since y=x^0.5
		  a_a0 = pow(zeror[rt],2.0);
		}
	}
	return a_a0;
}
double calculate_fn_fc(double a_a0)
{
	return 4.*(pow(a_a0,3.)-pow(a_a0,(3./2.)));
}


