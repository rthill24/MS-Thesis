#Paik and Thayamballi's steel regression equation for uniaxial panel strength
# (c) 2009 University of Michigan, Dr. Matthew Collette NA&ME Dept.


class PaikCReg:
    """
    Class for calculating the compressive strength of steel stiffened panels
    by Paik and Thayamballi's regression equation, limited by the Euler
    buckling strength of the material, as given in "Ultimate strength of
    aluminum plates and stiffened panels for marine applications" Paik and
    Duran, Marine Technology, Vol 41, No. 3, July 2004. 
    """
    
    def ucs(self, panel):
        '''
        Returns the ultimate strength of the panel by the equation
        '''
        b = panel.getBeta()
        l = panel.getLambda()
        ys = panel.getYsavg()
        

        #Formula should be limited by Euler buckling stress
        limit = ys/(l**2)

        #Caclulate regression formula
        bottom = 0.995 + 0.936*l**2 + 0.170*b**2 + 0.188*l**2*b**2 - 0.067*l**4
        bottom = bottom**(0.5)

        #Limit yield strength to proof strength
        if bottom > 1:
            us = ys/bottom
        else:
            us = ys

        #Limit to Euler buckling load
        if limit < us:
            us = limit

        return us


    
        
