import numpy as np
import subprocess
import pandas as pd
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def miller_ind(x):
        return '_'.join(map(str, tuple(x)))

def modify_unit_vector(x):
    if len(x[x % 2 ==1]) == 3: # if all 3 odd for BCC/// If 2 are even for FCC
        return x/2
    else:
        return x
    
    
def read_lmp_data(filepath):
    
    with open(filepath, 'r') as file:
        lines = file.readlines()
        atoms = lines[15:]

    with open('.core', 'w') as file:
        file.write('id type x y z\n')
        file.writelines(atoms)
    return pd.read_csv(filepath_or_buffer='core', sep='\s+')

    
class GenDipole:
    
    def __init__(self, angle, line_direction, lattice_const, elastic_const={'c11':243.0, 'c12': 145.0, 'c44':119.00}, atomsk_bin_path=f'{os.getenv("CONDA_PREFIX")}/bin/atomsk',propfile='Fe_prop.txt', z_basis = np.array([1, 0, 1]), Lx_max = 950, Ly_max = 100, Lz_max = 200):
        self.theta = angle
        self.y_basis = line_direction
        self.alat = lattice_const
        self.propfile = propfile
        self.Lx_max = Lx_max
        self.Ly_max = Ly_max
        self.Lz_max = Lz_max
        self.z_basis = z_basis
        temp_x_basis = np.cross(self.y_basis, self.z_basis)
        self.x_basis = (temp_x_basis/np.gcd.reduce(list(temp_x_basis))).astype(int)
        self.bin_exec = atomsk_bin_path

        docstring = f"""#Properties 
#Elastic constants (GPa)
elastic Voigt
{elastic_const['c11']} {elastic_const['c11']} {elastic_const['c11']}  #C11 C22 C33
{elastic_const['c12']} {elastic_const['c12']} {elastic_const['c12']}   #C23 C31 C12
{elastic_const['c44']} {elastic_const['c44']} {elastic_const['c44']}  #C44 C55 C66
"""
    
        with open(self.propfile, 'w') as file:
            file.write(docstring)
            file.write('#Crystallographic orientation of the block (this will rotate the elastic tensor accordingly)\n')
            file.write('orientation\n')
            file.write(miller_ind(self.x_basis)+'\n')
            file.write(miller_ind(self.y_basis)+'\n')
            file.write(miller_ind(self.z_basis)+'\n')
    
        self.b = np.sqrt(3)/2*self.alat
        self.bx = self.b*np.cos((90-self.theta)*np.pi/180)   ####
        self.by = self.b*np.sin((90-self.theta)*np.pi/180)   ####
        self.bz = 0

        ux_basis = modify_unit_vector(self.x_basis)
        uy_basis = modify_unit_vector(self.y_basis)
        uz_basis = modify_unit_vector(self.z_basis)
        
        ux = np.linalg.norm(ux_basis)*self.alat
        uy = np.linalg.norm(uy_basis)*self.alat
        uz = np.linalg.norm(uz_basis)*self.alat
        
        self.Nx = int(self.Lx_max/ux)+1
        self.Ny = int(self.Ly_max/uy)+1
        self.Nz = int(self.Lz_max/uz)+1
        
        self.Lx = self.Nx*ux
        self.Ly = self.Ny*uy
        self.Lz = self.Nz*uz
        
        
    
    def gen_same_glideplane_struc(self): #, b, bx, by, Nx, Ny, Nz, Lx):
        
        miller_x = miller_ind(self.x_basis) 
        miller_y = miller_ind(self.y_basis) 
        miller_z = miller_ind(self.z_basis) 

        outfile = f"dipole_same_glideplane_Fe_{self.theta}.lmp"

        if os.path.isfile(outfile):
            os.remove(outfile)

        atomsk_command = f'{self.bin_exec} --create bcc {self.alat} Fe orient [{miller_x}] [{miller_y}] [{miller_z}] \
                    -dup {self.Nx} {self.Ny} {self.Nz} -prop {self.propfile} \
                    -disloc loop 0.501*box 0.501*box 0.501*box Z {-0.251*self.Lx} {self.bx} {self.by} 0 0.3 \
                    {outfile}'
        
        args = atomsk_command.split()
        self.result = subprocess.run(args)
        return self.result
    
    def gen_diff_glideplane_struc(self, n_image):
        """ will use 'n_image' no. of images to sum periodic image contributions in the dipole displacement """
        
        miller_x = miller_ind(self.x_basis) 
        miller_y = miller_ind(self.y_basis) 
        miller_z = miller_ind(self.z_basis)
        
        
        disloc_str = ''

        loc_z1 = 0.75
        loc_z2 = 0.25
        loc_x1 = 0.50
        
    
        if (n_image==0) or (self.bx<np.finfo(float).eps):
            disloc_str+=(f'-disloc {loc_z1*self.Lz} {loc_x1*self.Lx} mixed y z {self.bx} {self.by} 0 ')
            disloc_str+=(f'-disloc {loc_z2*self.Lz} {loc_x1*self.Lx} mixed y z {-self.bx} {-self.by} 0 ')
            print (self.bx)

        else:
            for im_z in range(-n_image , n_image):
                p_z1 = im_z+loc_z1
                p_z2 = im_z+loc_z2
            
                for im_x in range(-n_image, n_image):
                                 
                    p_x = im_x+loc_x1
    
                    disloc_str+=(f'-disloc {p_z1*self.Lz} {p_x*self.Lx} mixed y z {self.bx} {self.by} 0 ')
                    disloc_str+=(f'-disloc {p_z2*self.Lz} {p_x*self.Lx} mixed y z {-self.bx} {-self.by} 0 -center com ')
        
        #bin_exec = '/users/sbagchi/atomsk/atomsk'
        outfile = f"dipole_diff_glideplane_Fe_{self.theta}.lmp"

        if os.path.isfile(outfile):
            os.remove(outfile)

        atomsk_command = f'{self.bin_exec} --create bcc {self.alat} Fe orient [{miller_x}] [{miller_y}] [{miller_z}] \
                    -dup {self.Nx} {self.Ny} {self.Nz} -prop {self.propfile} ' + disloc_str + outfile

            
        args = atomsk_command.split()
        self.result = subprocess.run(args)
            
        return #self.result
    
    
    
    def create_prsitine_struc(self, outfile='simple_pristine.lmp'):
        
        miller_x = miller_ind(self.x_basis) 
        miller_y = miller_ind(self.y_basis) 
        miller_z = miller_ind(self.z_basis)
        
        #bin_exec = '/users/sbagchi/atomsk/atomsk' 

        if os.path.isfile(outfile):
            os.remove(outfile)

        atomsk_command = f'{self.bin_exec} --create bcc {self.alat} Fe orient [{miller_x}] [{miller_y}] [{miller_z}] \
                        -dup {self.Nx} {self.Ny} {self.Nz} {outfile}'

        args = atomsk_command.split()
        result = subprocess.run(args)

        return result.stdout
    
   

    def periodic_corr(self, outfile, reset_ids=True):
        
        self.create_prsitine_struc()
        
        dipole_file=f'dipole_diff_glideplane_Fe_{self.theta}.lmp'
        pristine_file='simple_pristine.lmp' 
        # outfile='im_corr_struc_{:.6f}.lmp'.format(self.theta)
        
        #---load dipole and pristine dataframes----
        X = read_lmp_data(dipole_file)
        X0 = read_lmp_data(pristine_file)
        

        #-----Search for the corners------
        corner_1 = int(X0[(X0.x==0) & (X0.y==0) & (X0.z==0)].id)  # (0, 0, 0)
        
        temp = round(self.Lx  - self.Lx/self.Nx,5)
        corner_2 = int(X0[(abs(X0.x-temp)<1e-4) & (X0.y==0) & (X0.z==0)].id) # (Lx, 0, 0)
        
        temp = round(self.Ly  - self.Ly/self.Ny,5)
        corner_3 = int(X0[(abs(X0.y-temp)<1e-4) & (X0.x==0) & (X0.z==0)].id) # (0, Ly, 0)
        
        temp = round(self.Lz  - self.Lz/self.Nz,5)
        corner_4 = int(X0[(abs(X0.z-temp)<1e-4) & (X0.y==0) & (X0.x==0)].id) # (0, 0, Lz)

        #-------- Obtain image-summed Displacements ----------------

        X['ux'] = X.x - X0.x
        X['uy'] = X.y - X0.y
        X['uz'] = X.z - X0.z

        #--------- evaulate tensor [d] (u_sum = u_pbc + [d]*x)

        d = np.zeros((3,3))

        d[0][0] = (float(X[X.id==corner_1].ux) - float(X[X.id==corner_2].ux))/(float(X0[X0.id==corner_1].x) - float(X0[X0.id==corner_2].x))
        d[1][0] = (float(X[X.id==corner_1].uy) - float(X[X.id==corner_2].uy))/(float(X0[X0.id==corner_1].x) - float(X0[X0.id==corner_2].x))
        d[2][0] = (float(X[X.id==corner_1].uz) - float(X[X.id==corner_2].uz))/(float(X0[X0.id==corner_1].x) - float(X0[X0.id==corner_2].x))

        d[0][1] = (float(X[X.id==corner_3].ux) - float(X[X.id==corner_1].ux))/(float(X0[X0.id==corner_3].y) - float(X0[X0.id==corner_1].y))
        d[1][1] = (float(X[X.id==corner_3].uy) - float(X[X.id==corner_1].uy))/(float(X0[X0.id==corner_3].y) - float(X0[X0.id==corner_1].y))
        d[2][1] = (float(X[X.id==corner_3].uz) - float(X[X.id==corner_1].uz))/(float(X0[X0.id==corner_3].y) - float(X0[X0.id==corner_1].y))


        d[0][2] = (float(X[X.id==corner_4].ux) - float(X[X.id==corner_1].ux))/(float(X0[X0.id==corner_4].z) - float(X0[X0.id==corner_1].z))
        d[1][2] = (float(X[X.id==corner_4].uy) - float(X[X.id==corner_1].uy))/(float(X0[X0.id==corner_4].z) - float(X0[X0.id==corner_1].z))
        d[2][2] = (float(X[X.id==corner_4].uz) - float(X[X.id==corner_1].uz))/(float(X0[X0.id==corner_4].z) - float(X0[X0.id==corner_1].z))
    
    

    
        # Obtain PBC displacements by subtracting the error term from u_sum (u_pbc = u_sum - d*r)

        ux_pbc = X.ux - (d[0,0]*X0.x + d[0,1]*X0.y + d[0,2]*X0.z)
        uy_pbc = X.uy - (d[1,0]*X0.x + d[1,1]*X0.y + d[1,2]*X0.z)
        uz_pbc = X.uz - (d[2,0]*X0.x + d[2,1]*X0.y + d[2,2]*X0.z)

        # Add corrected PBC displacement to atom positions

        Xn = X0
        Xn.x = X0.x + ux_pbc
        Xn.y = X0.y + uy_pbc
        Xn.z = X0.z + uz_pbc
    
        # Remove the extra half-plane near left (right) edge of the box
        if self.bx>np.finfo(float).eps:
            Xnn = Xn[(Xn.x>=min(Xn.x) + self.bx)] # Xn.x <=max(Xn.x) - self.bx)
        else:
            Xnn = Xn

        # Resetting ids in case of deleted atoms
        if reset_ids:
            Xnn.id = np.arange(1,len(Xnn.id)+1, dtype=int)

        # Write lmp output

        with open(pristine_file, 'r') as file:
            lines = file.readlines()

        with open(outfile,'w') as file:
            file.write('# PBC corrected BCC Fe\n')
            file.write('\n')
            file.write(f'      {len(Xnn.id)} atoms\n')
            file.writelines(lines[3:15])
            Xnn.to_csv(file, mode='a', header=None, index=None, sep=' ')
            
        return
    
    
    def with_prism_loop(self):
        
        miller_x = miller_ind(self.x_basis) 
        miller_y = miller_ind(self.y_basis) 
        miller_z = miller_ind(self.z_basis) 
        
        #bin_exec = '/users/sbagchi/atomsk/atomsk' 

        outfile = "prismatic_loop.lmp"

        if os.path.isfile(outfile):
            os.remove(outfile)

        atomsk_command = f'{self.bin_exec} --create bcc {self.alat} Fe orient [{miller_x}] [{miller_y}] [{miller_z}] \
                    -dup {self.Nx} {self.Ny} {self.Nz} -prop {self.propfile} \
                    -disloc loop 0.501*box 0.501*box 0.501*box X {-0.251*self.Lx} {self.bx} {self.by} 0 0.3 \
                    {outfile}'
        
        args = atomsk_command.split()
        self.result = subprocess.run(args)
        return self.result


    def apply_prestress(self, filename, sig11=0.,sig22=0.,sig33=0.,sig12=0.,sig13=0.,sig23=0.):
        
        stress_label = round(np.linalg.norm(np.array([sig11, sig22, sig33, sig12, sig13, sig23])), 2)

        #bin_exec = '/users/sbagchi/atomsk/atomsk'
        atomsk_command = f'{self.bin_exec} {filename} -prop {self.propfile} -stress xx {sig11} -stress yy {sig22} -stress zz {sig33} -stress xy {sig12} \
                -stress xz {sig13} -stress yz {sig23} \
                prestress_{stress_label}.lmp'
        
        args = atomsk_command.split()
        self.result = subprocess.run(args)
        return self.result

    def shear_loop(self, burgers_vec, glide_plane_normal='Y', radius=150, res_shear=0.0, shear_comp='xy'):
        
        miller_x = miller_ind(self.x_basis)
        miller_y = miller_ind(self.y_basis)
        miller_z = miller_ind(self.z_basis)

        #bin_exec = '/users/sbagchi/atomsk/atomsk'

        outfile = "Shear_Loop_{res_shear}.lmp"

        if os.path.isfile(outfile):
            os.remove(outfile)

        atomsk_command = f'{self.bin_exec} --create bcc {self.alat} Fe orient [{miller_x}] [{miller_y}] [{miller_z}] \
                    -dup {self.Nx} {self.Ny} {self.Nz} -prop {self.propfile} -disloc loop 0.501*box 0.501*box 0.501*box {glide_plane_normal} {radius}\
                    {burgers_vec[0]} {burgers_vec[1]} {burgers_vec[2]} 0.33\
                    -stress {shear_comp} {res_shear}\
                    {outfile}'

        args = atomsk_command.split()
        self.result = subprocess.run(args)
        return self.result
