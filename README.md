# Oscullations
Codes that implement the evolution of a binary orbit at leading order in the Post-Newtonian expansion, via the method of oscullations. The time evolution of the orbital elements is computed at linear order by $Element(t+\delta t)=Element(t)+(dElement/dt) \delta t$, where $(dElement/dt)$ is computed at linear order in [1].


**orbit(omega, Omega, l, e, a, m1, m2, f=0)**:
Encodes the orbit of a binary pair with masses m1 and m2. The orbit is parametrized by the masses of the bodies (in solar masses),
            the eccentricity e, and the semi-major axis a. The orbital plane is described by three angles Omega, omega and l.
            The evolution of the orbit is measured with respect to the true anomaly f. By default, the starting anomaly is set to f=0. 
            Distances are given in m. Time in seconds, Masses in solar masses.



            
[1] Gravity, Newtonian, post-Newtonian, Relativistic. E. Poisson and C. M. Will. Cambridge University. Chapter 3.
    
