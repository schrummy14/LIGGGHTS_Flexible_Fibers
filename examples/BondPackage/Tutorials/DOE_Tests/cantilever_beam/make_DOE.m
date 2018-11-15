close all
clear
clc

N_bond_damp_vals = 2;
N_bond_Youngs_modulus = 2;

bond_damp_vals = linspace(10,50,N_bond_damp_vals);
bond_Youngs_modulus_vals = linspace(1e9,1e10,N_bond_Youngs_modulus);

variable_names = [
    "run_num"
    "bond_damp_val"
    "bond_youngs_modulus"
    ];

[bvs,bYms] = ndgrid(bond_damp_vals,bond_Youngs_modulus_vals);

flag = make_DOE_txt(variable_names,bvs(:),bYms(:));