import numpy as np

class orbit:
    """ 
        orbit(omega, Omega, l, e, a, m1, m2, f=0):
    
            Encodes the orbit of a binary pair with masses m1 and m2. The orbit is parametrized by the masses of the bodies (in solar masses),
            the eccentricity e, and the semi-major axis a. The orbital plane is described by three angles Omega, omega and l.
            The evolution of the orbit is measured with respect to the time t. By default, the starting anomaly is set to f=0. 
            Distances are given in m. Time in seconds, Masses in solar masses.
            
            see [Gravity, Newtonian, post-Newtonian, Relativistic. E. Poisson and C. M. Will. Cambridge University. Chapter 3]
    
        ************************************************************************************************************************************
            
            Methods:
            
                distance():
                    Returns the lenght of the radial vector.
                    output: float 32
                
                position_unit():
                    Returns the unit vector along the radial direction.
                    output: float 32, shape (3,)
                
                position():
                    Returns the position vector.
                    output: float 32, shape (3,)
            
                periastrum_unit():
                    Returns the unit vector along the direction to the periastrum.
                    output: float 32, shape (3,)
                    
                lambdavector():
                    Returns the unit vector lambda orthogonal to the position vector in the orbital plane.
                    output: float 32, shape (3,)
                    
                zvector():
                    Returns the unit vector orthogonal to the orbital plane.
                    output: float 32, shape (3,)
            
                time_step(delta, steps=1):
                    input: float, int
                    Performs *steps* time steps in orbital evolution, each with time lenght delta. This is done by using the oscullating formalism,
                    with the orbital parameters evolving according to param(t+delta)=(dparam/dt)*delta.
                    
                set_force(R,S,W):
                    input: function, function, function
                    Sets the force perturbation to a specific function shape. The force vector is decomposed according to the criteria in
                    [Gravity, Newtonian, post-Newtonian, Relativistic. E. Poisson and C. M. Will. Cambridge University. Chapter 3]
                
                force():
                    Returns the force vector.
                    output: float 32, shape (3,)
            
     
    """
    
    def __init__(self, omega, Omega, l, e, a, m1, m2, f=0):
        self.initial=[omega, Omega, l, e, a, m1, m2, f]
        self.l=l
        self.Omega=Omega
        self.omega=omega
        self.f=f
        self.e=e
        self.a=a
        self.p=a*(1-e**2)        
        self.t=0
        self.Ms=1.98847*10**30
        self.GM=6.674*10**(-11)*self.Ms*(m1+m2)
        
        
        
        self.R=lambda : 0
        self.S=lambda : 0
        self.W=lambda : 0
        
    def reset(self):
        [omega, Omega, l, e, a, m1, m2, f]=self.initial
        self.t=0
        self.l=l
        self.Omega=Omega
        self.omega=omega
        self.f=f
        self.e=e
        self.a=a
        self.p=a*(1-e**2)        
        self.t=0
        
    def distance(self):
        return self.p/(1+self.e*np.cos(self.f))
    
    def position(self):
        return self.distance()*self.position_unit()
    
    def position_unit(self):
        return np.array([self.__xcoord(),self.__ycoord(),self.__zcoord()],dtype='float32')
    
    def __xcoord(self):
        return np.cos(self.Omega)*np.cos(self.omega+self.f)-np.cos(self.l)*np.sin(self.Omega)*np.sin(self.omega+self.f)
    
    def __ycoord(self):
        return np.sin(self.Omega)*np.cos(self.omega+self.f)+np.cos(self.l)*np.cos(self.Omega)*np.sin(self.omega+self.f)
    
    def __zcoord(self):
        return np.sin(self.l)*np.sin(self.omega+self.f)        
    
    def lambdavector(self):
        return np.array([-np.cos(self.Omega)*np.sin(self.omega+self.f)-np.cos(self.l)*np.sin(self.Omega)*np.cos(self.omega+self.f),
                 -np.sin(self.Omega)*np.sin(self.omega+self.f)+np.cos(self.l)*np.cos(self.Omega)*np.cos(self.omega+self.f),
                  np.sin(self.l)*np.cos(self.omega+self.f)
                 ],dtype='float32')
    
    def zvector(self):
        return np.array([np.sin(self.l)*np.sin(self.Omega),
                        -np.sin(self.l)*np.cos(self.Omega),
                         np.cos(self.l)
                        ],dtype='float32')
    
    def periastrum_unit(self):
        return np.array([np.cos(self.Omega)*np.cos(self.omega)-np.cos(self.l)*np.sin(self.Omega)*np.sin(self.omega),
                        np.sin(self.Omega)*np.cos(self.omega)+np.cos(self.l)*np.cos(self.Omega)*np.sin(self.omega),
                        np.sin(self.l)*np.sin(self.omega)],dtype='float32')
    
    
    
    
    def time_step(self, delta, steps=1, deltas=0):
        for step in range(steps):
            
            deltaf=np.sqrt(self.GM/self.p**3)*(1+self.e*np.cos(self.f))**2 
            +(1/self.e)*np.sqrt(self.p/self.GM)*(np.cos(self.f)*self.R() - self.S()*np.sin(self.f)*(2+self.e*np.cos(self.f))/(1+self.e*np.cos(self.f)))
            
            deltap=2*np.sqrt(self.p**3/(self.GM))*1/(1+self.e*np.cos(self.f))*self.S()
            
            deltae=np.sqrt(self.p/(self.GM))*(np.sin(self.f)*self.R() +(2*np.cos(self.f)+self.e*(1+np.cos(self.f)**2))/(1+self.e*np.cos(self.f))*self.S())
            
            deltal=np.sqrt(self.p/(self.GM))*(np.cos(self.omega+self.f)*self.W()/(1+self.e*np.cos(self.f)))
            
            deltaomega=np.sqrt(self.p/(self.GM))*(1/self.e)*(-np.cos(self.f)*self.R() + (2+self.e*np.cos(self.f))*np.sin(self.f)*self.S()/(1+self.e*np.cos(self.f)) )
            
            if self.W()!=0:
                deltaomega+=-np.sqrt(self.p/(self.GM))*(1/self.e)*((self.e/np.tan(self.l))*self.W()*np.sin(self.omega+self.f)/(1+self.e*np.cos(self.f)) )
                                                                  
            
            if np.sin(self.l)==0:
                deltaOmega=0
            else:
                deltaOmega=(np.sqrt(self.p/(self.GM))*self.W()*np.sin(self.omega+self.f)/(1+self.e*np.cos(self.f)))/np.sin(self.l)
            
            self.t+=delta
            self.f+=deltaf*delta
            self.p+=deltap*delta
            self.e+=deltae*delta
            self.l+=deltal*delta
            self.Omega=deltaOmega*delta
            self.omega=deltaomega*delta
            self.a=self.p/(1-self.e**2)
            
            if deltas==1:
                return {'deltaf': deltaf,'deltap':deltap, 'deltae': deltae, 'deltal': deltal,'deltaOmega': deltaOmega,'deltaomega': deltaomega}
    

    
    def set_force(self,R,S,W):
        if R is not None:
            self.R=R
        if S is not None:
            self.S=S
        if W is not None:
            self.W=W
            
    def force(self):
        return self.R()*self.position_unit()+self.S()*self.lambdavector()+self.W()*self.zvector()
