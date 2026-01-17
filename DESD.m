% This function is called by Dynamics.m
function  [Psi_dot] = DESD(Psi,sigma,Deqn,Omega)

Psi_dot = Omega - sigma*transpose(Deqn)*sin(Deqn*Psi);

end