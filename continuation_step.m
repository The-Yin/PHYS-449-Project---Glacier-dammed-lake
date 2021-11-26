function [v_out, Dv_out, faux] = continuation_step(v_0,pars,functions_structure)
%Function whose zeros define arc length continuation step for finding
%       v_0:            (2*ndegf+1)-by-1 input,[initial guess of closed orbit] + [ntial
%                       guess of fixed point] + [updated parameter value]
%       parameters:     parameter structure to be passed to be
%                       parallel_shoot_par
%       v_prev:         previous solution v_0
%       dv_norm:        norm of required step length
%       scalefactor:    metric for computing step length
%       scale:          step size computed from change in parameter v_0(4)
%                       ('linear', default) or from change in
%                       log(parameter) ('log')
%Output is
%       v_out:          function whose zeros define arc length continuation
%                       solution
%       Dv_out:         Jacobian matrix of that function
%       faux:           Optional structure with additional information on closed
%                       orbit, with fields
%           DvP:        When v_out is a fixed point, Jacobian of Poincare map
%                       (extended to three dimensions with null space in
%                       direction of flow at fixed point)
%           t_orbit:    period of orbit
%       
%WCode is identical to v3 except that it invokes parallel_shoot_par_v5

%initialize
if nargin > 1 && isfield(pars,'parameters')
    parameters = pars.parameters;
end
if nargin > 1 && isfield(pars,'dv_norm')
    dv_norm = pars.dv_norm;
else
    dv_norm = 1e-2;
end
if nargin > 1 && isfield(pars,'scalefactor')
    scalefactor = pars.scalefactor;
else
    scalefactor = ones(size(v_0));
end
if nargin > 1 && isfield(pars,'scale')
    scale = pars.scale;
else
    scale = 'linear';
end
if nargin > 1 && isfield(pars,'vprev')
    vprev = pars.vprev;
else
    vprev = ones(size(v_0));
end

v_out = zeros(size(v_0));
Dv_out = zeros(length(v_0));
[v_out1 Dv_out1 faux] = shoot_method_para(v_0,parameters,functions_structure);
v_out(1:end-1) = v_out1;
Dv_out(1:end-1,:) = Dv_out1;
switch scale
    case 'linear'
        v_out(end) = norm((v_0-vprev).*(scalefactor.^(1/2))) - dv_norm;
        Dv_out(end,:) = (v_0-vprev).*scalefactor./norm((v_0-vprev).*(scalefactor.^(1/2)));
    case 'log'
        v_out(end) = norm((log(v_0)-log(vprev)).*(scalefactor.^(1/2))) - dv_norm;
        Dv_out(end,:) = (log(v_0)-log(vprev)).*scalefactor./v_0./norm((log(v_0)-log(vprev)).*(scalefactor.^(1/2)));
end