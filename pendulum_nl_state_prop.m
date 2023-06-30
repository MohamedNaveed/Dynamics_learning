function [state_n] = pendulum_nl_state_prop(t, state, U, model)

t_span = [(t-1)*model.dt, t*model.dt];
options = odeset('RelTol',1e-8,'AbsTol',1e-10);

[temp, X_out] = ode45(@(t,y) pendulum_nl_ode(t,y,U,model), t_span, state, options);
%X_out(end,1) = atan2(sin(X_out(end,1)), cos(X_out(end,1)));

state_n = X_out(end,:);

end