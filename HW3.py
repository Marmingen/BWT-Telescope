#####################################################################
# M. Marminge
#####################################################################
## IMPORTS

import math

#####################################################################
## TELESCOPE CLASS

class Telescope():
    """initialized with UU-BWT information"""
    def __init__(self, theta_blaze=66.5, f_col=0.8, L_grate=0.2, 
                 gpm=72e3, diameter=0.9, CCD_pixelsize=15e-6,
                 focal_ratio=0.25):
        
        #### input data ####
        self.theta_blaze = theta_blaze*math.pi/180 # from degrees to radians
        self.f_col = f_col # m
        self.L_grate = L_grate # m
        self.gpm = gpm # grooves per meter
        self.diameter = diameter # m
        self.CCD_pixelsize = CCD_pixelsize # m
        self.focal_ratio = focal_ratio
        self.wavelengths = []
        
        #### calculated data ####
        self.N = L_grate*gpm # amount of grooves
        self.delta = 1/gpm # m
        self.W = L_grate # m
    
    def set_seeing(self, seeing):
        """input seeing as arcsec"""
        self.seeing = seeing/60**2/180*math.pi # from arcsec to radians
        self.image = self.seeing*self.diameter/self.focal_ratio
        
    def set_sampling(self, sampling):
        self.target_image = self.CCD_pixelsize*sampling
        self.sampling = sampling
        
    def add_wavelengths(self, lambdas):
        self.wavelengths += [Wavelength(lam, self) for lam in lambdas]
        
    def calc_true_res(self):
        if self.image and self.target_image:
            for wave in self.wavelengths:
                wave.calc_f_cam(self)
                wave.calc_true_resolution(self)
        else:
            print("set local seeing and desired sampling rate first")
            
    def print_wavelength_data(self):
        if self.image and self.target_image:
            print("data for the wavelengths")
            print("|","-"*85, "|")
            print(f"| {'wavelength [Å]':<15} | {'optimal m':<10} | {'optimal s [microns]':<20} | {'f_cam [m]':<14} | {'true res':<14} |")
            for wave in self.wavelengths:
                print(f"| {wave.wavelength/1e-10:<15.3f} | {wave.m:<10} | {wave.optimal_s/1e-6:<20.3f} | {wave.f_cam:<14.3f} | {wave.true_resolution:<14.2f} |")
            print("|", "-"*85, "|")
        else:
            print("set local seeing and desired sampling rate first")
            
    def print_telescope_data(self):
        print("data for the telescope")
        print("|", "-"*85, "|")
        print(f"| {'diameter [m]':<20} | {'f_tel [m]':<20} | {'f_col [m]':<20} | {'blaze angle [deg]':<17}|")
        print(f"| {self.diameter:<20} | {self.diameter/self.focal_ratio:<20.2f} | {self.f_col:<20} | {math.degrees(self.theta_blaze):<17.2f}|")
        print("|", "-"*85, "|")
        print(f"| {'L_grate [m]':<20} | {'gpm [1/m]':<20} | {'pixelsize [microns]':<20} | {'N':<17}|")
        print(f"| {self.L_grate:<20} | {self.gpm:<20} | {self.CCD_pixelsize/1e-6:<20.2f} | {self.N:<17}|")
        print("|", "-"*85, "|")
        
#####################################################################
## WAVELENGTH CLASS

class Wavelength():
    
    def __init__(self, wavelength, telescope):
        
        #### input data ####
        self.telescope = telescope
        self.wavelength = wavelength
        
        #### calculated data ####
        self.__calc_m(self.telescope)
        self.optimal_resolution = self.m*telescope.N
        self.__calc_optimal_s(self.telescope)
        self.__calc_beta(self.telescope)
        
    def __calc_m(self, tel):
        """calculates the optical order"""
        self.m = int(2*tel.delta/self.wavelength*math.sin(tel.theta_blaze)+0.5)

    def __calc_optimal_s(self, tel):
        """calculates the optimal slit width"""
        self.optimal_s = self.wavelength*tel.f_col/(tel.W*math.cos(tel.theta_blaze))
        
    def __calc_beta(self, tel):
        """calculates the reflection angle"""
        self.beta = math.asin(self.m*self.wavelength/tel.delta - math.sin(tel.theta_blaze))
        
    def calc_f_cam(self, tel):
        """calculates the optimal focal length of the camera"""
        self.f_cam = tel.f_col*tel.target_image/tel.image*math.cos(self.beta)/math.cos(tel.theta_blaze)
    
    def calc_true_resolution(self, tel):
        """calculates the true resolution based on sampling and seeing"""
        self.true_resolution = self.wavelength*self.m*self.f_cam/(tel.target_image*tel.delta*math.cos(self.beta))

#####################################################################
## MAIN

def main():
    
    lam1 = 4000e-10
    lam2 = 6000e-10
    lam3 = 8000e-10
    lam4 = 256000e-10
    lam5 = 1e-7
    
    lambs = [lam1, lam2, lam3]
    
    BWT = Telescope()
    
    BWT.diameter = 39

    BWT.set_seeing(2)
    BWT.set_sampling(3)
    
    
    BWT.add_wavelengths(lambs)
    
    BWT.calc_true_res()
    
    
    BWT.print_wavelength_data()
    
    BWT.print_telescope_data()

#####################################################################
## RUN CODE

if __name__ == "__main__":
    main()