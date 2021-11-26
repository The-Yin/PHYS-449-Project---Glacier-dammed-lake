function fout = arclength_continuation_method(v_0,functions_structure, lmin,ds,scalefactor,scale,direction,npts,i_crossing,parameters,srchparams)
%fout = parallel_periodic_arclength_nullcline(v_0,functions_structure, lmin,ds,scalefactor,scale,direction,npts,i_crossing,parameters,srchparams)
%Arc length continuation computation of closed orbits using
% @shoot_method, @continuation_step, @shoot_method_para, and @Newton_single
%Input variables (all optional) are:
%       v_0:    2*ndegf-by-one vector consisting of initial guess for initial
%               point on a closed  orbit, concatenated with
%               an initial guess for fixed point of dynamical system
%       functions_structure: structure of dynamic systems function and event function
%       lmin:   first parameter value to be used
%       ds:     Arc length size
%       scalefactor: (diagonal) metric for computation of arc length, takes
%               the form of a (2*ndegf+1)-by-1 vector so that arc
%               length is computed as \sum_i a_i dv_i^2 if scale is
%               'linear', and as  \sum_i a_i d(log(v_i))^2 if scale is
%               'log'
%       scale:  step size ('linear') or ('log')  log(parameter)
%       direction:  initial direction for changes in l, choices +1 and -1
%       npts:   Number of points at which solution should be computed
%       i_crossing: sets index of variable that will return to its original
%               value to compute Poincare map. 
%       parameters: parameter structure
%       srchparams: search parameters to be passed to Newton solver
%Output variables
%       fout:   output structure containing fields
%           v:  ndegf-by-npts array of initial points lying on closed orbits
%           v_steady:   ndegf-by-npts array of fixed points
%           l:  1-by-npts vector of corresponding parameter values
%           error:  1-by-npts vector of booleans, indicating true if
%               computation did not converge, false otherwise
%           Dv: ndegf-by-ndegf-by-npts array of Jacobians of Poincare map v_in ->
%               v_f, to compute stability of closed orbit
%           eigmax: maximum absolute value of eigenvalues of Dv, should
%                determine stability of orbit
%           t_orbit: orbital period
%           parameters: the parameters structure


%set i_crossing
parameters.i_crossing = i_crossing

%set up higher level parameters structure to pass to 
pars.dv_norm = ds;
pars.scalefactor = scalefactor;
pars.scale = scale;
pars.parameters = parameters;

%initialize output
fout.v = zeros(parameters.ndegf,npts);
fout.v_steady = zeros(parameters.ndegf,npts);
fout.l = zeros(1,npts);
fout.error = zeros(1,npts);
fout.Dv = zeros(parameters.ndegf,parameters.ndegf,npts);
fout.eigmax = zeros(1,npts);

fout.t_orbit = zeros(1,npts);
fout.parameters = parameters;



%compute first point: find steady state for imposed parameter value ---
%improve this to make more general using anonymous functions
%vaux = ttt_steady(v_0(parameters.ndegf+1:2*parameters.ndegf),parameters,srchparams);

%amend initial guess
%v_0(parameters.ndegf+1:2*parameters.ndegf) = vaux;

% I let steady solution to be zero always, the orginal code is above
vaux = [0;0];
v_0(parameters.ndegf+1:2*parameters.ndegf) = vaux;



%now compute corresponding closed orbit --- ditto for making more general
[vtemp, errorflag, faux] = Newton_single(@shoot_method,v_0(1:parameters.ndegf),parameters,srchparams,functions_structure);



fout.v(:,1) = vtemp;
fout.v_steady(:,1) = vaux;
fout.l(1) = lmin;
fout.error(1) = errorflag;

fout.Dv(:,:,1) = faux.DvP;
fout.eigmax(1) = max(real(eig(faux.DvP)));

fout.t_orbit(1) = faux.t_orbit;

%amend initial guess again

v_0(1:parameters.ndegf) = vtemp;

%loop
vtemp = [v_0; lmin];   %change variable vector to include parameter guess



dv = zeros(2*parameters.ndegf+1,1);
switch scale
    case 'linear'
        dv(end) = ds*direction/scalefactor(end)^(1/2);
    case 'log'
        dv(1:2*parameters.ndegf) = 1;
        dv(end) = exp(direction*ds/scalefactor(end)^(1/2));
end

for ii=2:npts
    ii
    pars.parameters.ii = ii
    vprev = vtemp;
    switch scale
        case 'linear'
            vtemp = vprev+dv;
        case 'log'
            vtemp = vprev.*dv;
    end
    pars.vprev = vprev;
    %can probably make this more general by putting anonymous funciton
    %handle in pars structure and having parallel_shoot_arclength_v4 the
    %same code (not problem-specific)
    [vtemp, errorflag, faux] =  Newton_single(@continuation_step,vtemp,pars,srchparams,functions_structure);
    fout.v(:,ii) = vtemp(1:parameters.ndegf);
    fout.v_steady(:,ii) = vtemp(parameters.ndegf+1:2*parameters.ndegf);
    fout.l(ii) = vtemp(2*parameters.ndegf+1);
    fout.error(ii) = errorflag;
    fout.Dv(:,:,ii) = faux.DvP;
    fout.eigmax(ii) = max(abs(eig(faux.DvP)));
    fout.t_orbit(ii) = faux.t_orbit;
    switch scale
        case 'linear'
            dv = vtemp-vprev;
        case 'log'
            dv = vtemp./vprev;
    end
end

end
