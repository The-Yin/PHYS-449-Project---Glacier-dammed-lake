function [v_0,function_structure,lmin,ds,scalefactor,scale,direction,npts,i_crossing,parameters,srchparams] = argument_generator()
% Output the arguments needed to start arc length continuation method by @arclength_continuation_method
% change the code of this function directly if need to modify the argument
% function_structure is main part of this code:
%             evolve: n = 0 - right hand side of equations of a specific system
%                     n = 1 - Jacobian matrix of the system
%                     n = 2 - Hessian matrix of the system
%             events: n = 0 - event function
%                     n = 1 - derivative of event function with respect to x_i
%             steppar: n = 0 - derivative of right hand side of system's ODEs with resprect to parameter
%                      n = 1 - derivative of Jacobian with resprect to parameter
%             mass: mass function
% Just change the function_structure in this code if we want to analysis different system

% initiate fragmented  arguments
v_0 = [1/2;(sqrt(7)-2)/2;0;0]
lmin = 1  
ds = 1e-2
scalefactor = [1 1 1 1 1]'
scale = 'linear'
direction = +1
npts = 10
i_crossing = 1

% initiate function sturucture
function_structure.evolve = @evolvefunction
function_structure.steppar = @stepparfunction
function_structure.mass = @massfunction
function_structure.events = @eventsfunction

% initiate parameter sturucture
parameters.ndegf = 2
parameters.type = 'c'
parameters.c = 1
parameters.solver.method = '15s'
parameters.solver.RelTol = 1e-12
parameters.solver.AbsTol = 1e-12
parameters.solver.delta_v_0 =  100*eps*norm(v_0)
parameters.t_span = [0 40]
parameters.ii =1

% initiate srchparams sturucture
srchparams.itmax = 25;
srchparams.tolF = sqrt(eps);
srchparams.toldelta = sqrt(eps);
srchparams.verbose = 0;
end


% the functions below are determined by a specific system
function result = evolvefunction(t,v,n,params)
x = v(1)
y = v(2)
c = params
switch n
    case 0
        result = [x-y-c*(x^2+y^2)*x;y+x-c*(x^2+y^2)*y]
    case 1
        result = [1-c*(3*x^2+y^2),-1-c*2*x*y;1-c*2*x*y,1-c*(x^2+3*y^2)]
    case 2
        result = [-1*c*6*x,-1*c*2*y;-1*c*2*y,-1*c*2*x]
        result(:,:,2) = [-1*c*2*y,-1*c*2*x;-1*c*2*x,-1*c*6*y]
end
end

function result = massfunction(t,params)
result = [1,0;0,1]
end

function result = eventsfunction(t,v,n,params)
x = v(1)
y = v(2)
c = params
if n == 0
    event_aux = [x-y-c*(x^2+y^2)*x;y+x-c*(x^2+y^2)*y];
    event_out = event_aux(1);
    event_aux(1) = [];
    if [1] == sign(event_aux)
        isterminal = true;
    else
        isterminal = false;
    end
    direction = 0;
    result = [event_out,isterminal,direction]
else
    result = [1-c*(3*x^2+y^2),-1-c*2*x*y]
end
end

function result = stepparfunction(t,v,n)
x = v(1)
y = v(2)
if n == 0
    result = [-(x^2+y^2)*x;-(x^2+y^2)*y]
else
    result = [-(3*x^2+y^2),-2*x*y;-2*x*y,-(x^2+3*y^2)]
end
end

