import numpy as np
import pandas as pd

##Density and shear wave velocity at source
vr = 3.5
pr = 2.75

class get_Vs():

    def __init__(self,freqs,site_amp):

        self.freqs = freqs
        self.site_amp = site_amp

        self.get_shallow_layer()
        self.get_vs_profile()


    def solve_eqn(self,a,b,c):
        discriminant = np.sqrt(b**2 - 4*a*c)
        
        sol1 = (-b+discriminant)/(2*a)
        sol2 = (-b-discriminant)/(2*a)
        
        solution = max(sol1,sol2)
        return solution

    def calc_density(self,v):
        return 1.742 + 0.2875*v

    def get_shallow_layer(self):

        site_amp_high = self.site_amp[0]
        freq_high = self.freqs[0]
        self.vs_shallow = self.solve_eqn(0.2875,1.742,-(pr*vr)/site_amp_high**2)
        self.depth_shallow = self.vs_shallow/(4*freq_high)


    def get_vs_profile(self):
        
        tt_old = 1/(4*self.freqs[0])
        depth_old = self.depth_shallow
        density_old = self.calc_density(self.vs_shallow)
        velocity = [self.vs_shallow]
        depth = [self.depth_shallow]
        
        for f,s in zip(self.freqs[1:],self.site_amp[1:]):
            tt_new = 1/(4*f)
            del_tt = tt_new - tt_old
        
            a = 0.2875*del_tt
            b = 1.742*del_tt
            c = depth_old*density_old - (pr*vr*tt_new)/(s**2)
        
            vs_current = self.solve_eqn(a,b,c)
            velocity.append(vs_current)
        
            del_depth = vs_current*del_tt
            density_new = ((depth_old*density_old)+(del_depth*(1.742+0.285*vs_current)))/(depth_old+del_depth)
        
            depth_old = depth_old + del_depth
            density_old = density_new
            tt_old = tt_new
            depth.append(depth_old)

        self.vs_profile = pd.DataFrame([depth,velocity]).T
        self.vs_profile.columns = ['Depth (km)','Vs (km/s)']
