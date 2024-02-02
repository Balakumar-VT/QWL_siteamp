import numpy as np

##Density and shear wave velocity at source
vr = 3.5
pr = 2.75

class Amplification():

    def __init__(self,vs,depth,freq = np.logspace(-1,2,100)):
        self.vs = vs
        self.depth = depth
        self.freq = freq

        self.thickness = np.insert(np.diff(self.depth),0,self.depth[0])
        self.site_amp = self.calc_QWL_amp()

    def calc_tt(self):
        self.tt = 1/(4*self.freq)

    def calc_density(self,v):
        return 1.742 + 0.2875*v


    def calc_vs_eql(self):
        
        z_eql = []
        for i in self.tt:
            tt_prof = 0
            z = 0
            for t,v in zip(self.thickness,self.vs):
                if tt_prof + t/v <= i:
                    tt_prof = tt_prof + t/v
                    z = z+t
                else:
                    diff_tt = i - tt_prof
                    z = z+ diff_tt*v
                    break
                
            z_eql.append(z)

        self.z_eql = np.array(z_eql)
        self.vs_eql = self.z_eql/self.tt


    def calc_density_eql(self):
        density_eql = []
        for i in self.z_eql:
            z = 0
            d = 0
            j = 0
            while z+self.thickness[j] < i:
                
                d = d+self.calc_density(self.vs[j])*self.thickness[j]
                z = z+self.thickness[j]
        
                j = j+1
                
            diff_z = i-z
            d = d+(self.calc_density(self.vs[j]))*diff_z
        
            density_eql.append(d/i)

        self.density_eql = np.array(density_eql)

    def calc_QWL_amp(self):

        self.calc_tt()
        self.calc_vs_eql()
        self.calc_density_eql()

        amp = np.sqrt((pr*vr)/(self.density_eql*self.vs_eql))

        return amp
                