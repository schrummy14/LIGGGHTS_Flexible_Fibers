#ifndef JKR_LOOKUP_TABLE_H
#define JKR_LOOKUP_TABLE_H 1

#include <iostream>
#include "PolynomialRootFinder.h"       //Added for solving JKR polynomial
#include <cmath>
#include <vector>
#include <algorithm>
#include <utility>
#include <fstream>
#include <iomanip>

double calculate_a_a0(double delta_delta_c);
double calculate_fn_fc(double a_a0);
void calculate_a_a0(std::vector<double> &delta_delta_c,std::vector<double> &a_a0);

class LUT {
    private:
        const static int points = int(3e5);
        double crit_lim = 1.;
        double alpha_l = 2.;
        double fn_fc_tol_factor = 1.;
        double delta_delta_c_max;
        double delta_delta_c[points];
        double a_a0[points];
        double fn_fc[points];
        bool outside_LUT = false;

    public:
        LUT() {};
        //~LUT(){delete delta_delta_c;delete a_a0;delete fn_fc;};
        void init(double d_delta_delta_c);
        std::pair<int,int> get_closest_indexes(double delta_delta_c_val);
        double get_a_a0(double delta_delta_c_val);
        double get_fn_fc(double delta_delta_c_val);
        double get_fh_fc(double delta_delta_c_val);
        double calculate_fh_fc(double fn_fc);
        
        double calculate_a0(double Wij,double reff,double Yeff);
        double calculate_delta_c(double a_0,double reff);
        double calculate_fc(double Wij,double reff);
        
        double get_max_delta_delta_c();
        int get_nr_of_points();

        //writing out to csv file delta_delta_c,a_a0,fn_fc
        void write_csv_table();
        //prints out table info (points) and max delta_delta_c
        void print_table_info();
        //prints out full table 
        void print_table();
};

#endif
