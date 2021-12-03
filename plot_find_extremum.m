function [f_min,f_max] = plot_find_extremum(x_set,y_set,period_set,c_set,function_structure)
% Find extremum value of x_i in the closed orbit with given parameter
n = length(x_set)
f_min = zeros(1,n)
f_max = zeros(1,n)
for i = 1:n
    x = x_set(i)
    y = y_set(i)
    c = c_set(i)
    period = period_set(i)
    tspan = [0 period]
    w0 = [x,y,c]
    [t,w] = ode15s(@evolve_system,tspan,w0)

    f_min(i) = min(w(:,1))
    f_max(i) = max(w(:,1))
end

function v_out = evolve_system(t,v_in)
%change function_structure.evolve(n=0) to the form that suitable for ode15s
%initialize
v_out = zeros(3,1);
%unpack parameters
c = v_in(3);
%unpack input
x = v_in(1);
y = v_in(2);
%compute output by evolve function
result = function_structure.evolve(t,[x,y],0,c)
v_out(1) = result(1)
v_out(2) = result(2)
v_out(3) = 0
end
end