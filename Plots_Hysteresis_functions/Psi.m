function [ psi ] = Psi( xi, epsilon )

N = length(xi);

psi = zeros(N,1);

psi = epsilon*tan(0.5*pi*xi);

end

