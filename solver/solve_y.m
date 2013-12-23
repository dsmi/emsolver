function Y = solve_y(mesh, contacts, opts)
% Y = solve_y(mesh, contacts, opts)
%
% Given a system of conductors with N terminals, this function calculates
% the admittance matrix Y of size N-by-N by applying the unity voltage to
% each of the terminals one by one and calling solve_c to find the resulting
% currents. The opts parameters is the solver options and is passed
% to solve_c.
%

[ M, T ] = mklhsmat(mesh, contacts, opts);

nports = length(contacts);
Y = [];
for ip=1:nports,
    potentials = zeros(nports,1);
    potentials(ip) = 1;
    % Find currents for the given contact potentials
    currents = solve_c(mesh, contacts, opts, M, T, potentials);
    % Minus is here because we want to treat the contact as a positive
    % terminal of the port and the ground as the negative one; for the
    % positive terminal the positive current direction is into the
    % terminal; while the solver returns the currents flowing out of
    % the contacts.
    Y = [ Y -currents ];
end
