function xdot = objectEOM(t,x,rho,Cd,A,m,g,wind_vel)
    % Inputs: t: time 
    % unpacking position and velocity
    pE_E = x(1:3);
    vE_E = x(4:6);

    v_E = vE_E - wind_vel;
    v_a = norm(v_E);

    g_matrix = [0, 0, g]' ;

    a_drag = (-1/2 * rho * v_a * A * Cd / m .* v_E );
    a_E = a_drag + g_matrix;
    xdot = [v_E; a_E];
end
