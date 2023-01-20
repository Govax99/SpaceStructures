function [symsol,sol] = ritz_sym(pvw,var,data,coeff,dcoeff)
%RITZ_SYM solve pvw approximated with ritz symbolically
%
% PROTOTYPE:
%     a = ritz_sym(pvw,var,data,coeff,dcoeff)
%
% INPUT:
%    pvw[-]             one line principle of virtual work [J]
%    var[1+]            data parameter variables[-]
%    data[1+]           data numeric values [-]
%    coeff[1+]          ritz coefficent names [-]
%    dcoeff[1+]         OPTIONAL - ritz coefficent delta names [-]
%
% OUTPUT:
%    symsol[dim]        contains symbolic expression of the ritz coefficients [-]
%    sol[dim]           contains numerical values of ritz coefficients
%
% CONTRIBUTORS:
%   Davide Zamblera - Politecnico di Milano - Space eng.
%
% VERSIONS
% 2022-08-20: First version
%

if nargin < 5
    %if dcoeff is not specify assume a -> da
    coeff_str = arrayfun(@string, coeff);
    dcoeff_str = arrayfun(@(x) "d"+x,coeff_str);
    dcoeff = str2sym(dcoeff_str);
end
sys = lhs(pvw) - rhs(pvw);
sys = children(collect(sys,dcoeff));
ch = [sys{:}];
ch = ch(ch ~= dcoeff); %case where there is only 1 variable has wrong childrens


sys_eq = subs(ch == 0,dcoeff,ones(size(dcoeff)));


[A,B] = equationsToMatrix(sys_eq, coeff);
symsol = linsolve(A,B);
sol = double(subs(symsol,var,data));

end