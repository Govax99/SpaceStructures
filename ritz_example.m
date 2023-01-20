clear;
close all;
clc;

%% MSS_2021_a1

%define variables and data
syms x
syms a2 a3
syms da2 da3
ritz_coeff = [a2 a3];
syms l EJ k kt F
var = [l EJ k kt F];
data = [500 1e10 80 2e7 2500];

%build relevant displacements
w1 = a2*(x/l)^2 + a3*(x/l)^3;
w2 = (a2 + a3)*(x/l)^2;
dw1 = da2*(x/l)^2 + da3*(x/l)^3;
dw2 = (da2 + da3)*(x/l)^2;

%write principle of virtual work
pvw = int(diff(dw1,x,2)*EJ*diff(w1,x,2) + diff(dw2,x,2)*2*EJ*diff(w2,x,2),x,0,l) + ...
    subs(dw2*k*w2,x,l) + subs(diff(dw2,x)*kt*diff(w2,x),x,l) == subs(dw2,x,l)*F;

% compute and show solution
[sol_sym, sol] = ritz_sym(pvw,var,data,ritz_coeff);
sol_sym
sol