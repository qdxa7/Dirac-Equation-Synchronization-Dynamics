function  [Psi_dot] = Dirac(Psi,sigma, Deqn,Omega)


Psi_dot = Omega - sigma*transpose(Deqn)*sin(Deqn*Psi);

end