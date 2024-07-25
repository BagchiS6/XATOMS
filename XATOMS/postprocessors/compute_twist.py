from ovito.data import CutoffNeighborFinder
from ovito.modifiers import ExpressionSelectionModifier, DeleteSelectedModifier
import numpy as np

# Taken from Soumendu's code
def cart2pol(X):
    x = X[0]
    y = X[1]
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)


def extract_local_angle(finder, index, central_atom_pos):

        neigh_index = []
        neigh_pos_polar = []
        neigh_pos_cart = []

        for neigh in finder.find(index):
                neigh_index.append(neigh.index)
                neigh_pos_cart.append(np.array([neigh.delta[0]+central_atom_pos[0], neigh.delta[1]+central_atom_pos[1]]))
                neigh_pos_polar.append(cart2pol(np.array([neigh.delta[0], neigh.delta[1]])))


        neigh_pos_polar = np.array(neigh_pos_polar)
        neigh_pos_cart = np.array(neigh_pos_cart)


        neigh_pos_polar_sorted = neigh_pos_polar[abs(neigh_pos_polar[:,1]).argsort()]
        neigh_pos_cart_sorted = neigh_pos_cart[abs(neigh_pos_polar[:,1]).argsort()]

        neigh_index = np.array(neigh_index)
        neigh_index_sorted = neigh_index[abs(neigh_pos_polar[:,1]).argsort()]

        is_zero = True # In case interested to implement local twists in future

        return {'index' :neigh_index_sorted[0], 'polar_angle': neigh_pos_polar_sorted[0,1], \
                'rel_cart_pos': neigh_pos_cart_sorted[0,:], 'zero_axis': is_zero}


def extract_layer_alignment(finder, initial_index, num_iter=10):
        ic = 0
        lattice_axis_info = {'index': [], 'polar_angle':[], \
                             'rel_cart_pos': [], \
                             'zero_axis': []}
        index = initial_index
        central_atom_pos= np.array([0,0])
        while ic < num_iter:
                local_angle_info = extract_local_angle(finder, index, central_atom_pos)

                # might be smarter way to do the following, for now bootstrap!
                lattice_axis_info['index'].append(local_angle_info['index'])
                lattice_axis_info['polar_angle'].append(local_angle_info['polar_angle'])
                lattice_axis_info['rel_cart_pos'].append(local_angle_info['rel_cart_pos'])
                lattice_axis_info['zero_axis'].append(local_angle_info['zero_axis'])

                index = local_angle_info['index']
                central_atom_pos = local_angle_info['rel_cart_pos']
                print (f'Done with iteration number: {ic}')
                ic+=1

        return lattice_axis_info



def compute_interlayer_angle(layer_1_info, layer_2_info):

        from scipy import stats
        from scipy.stats import t

        x1 = x2 = y1 = y2 = np.array([0.])

        x1=np.append(x1, np.array(layer_1_info['rel_cart_pos'])[:,0][layer_1_info['zero_axis']])
        y1=np.append(y1, np.array(layer_1_info['rel_cart_pos'])[:,1][layer_1_info['zero_axis']])

        x2=np.append(x2, np.array(layer_2_info['rel_cart_pos'])[:,0][layer_2_info['zero_axis']])
        y2=np.append(y2, np.array(layer_2_info['rel_cart_pos'])[:,1][layer_2_info['zero_axis']])

        res1=stats.linregress(x1,y1)
        res2=stats.linregress(x2,y2)

        # plt.plot(x1,y1, 'bo')
        # plt.plot(x1, res1.intercept + res1.slope*x1, 'b-')
        # plt.plot(x2,y2, 'ro')
        # plt.plot(x2, res2.intercept + res2.slope*x2, 'r-')

        tinv = lambda p, df: abs(stats.t.ppf(p/2, df))
        ts1 = tinv(0.05, len(x1)-2)
        err_1 = ts1*res1.stderr

        ts2 = tinv(0.05, len(x2)-2)
        err_2 = ts2*res2.stderr

        num_x = abs(res1.slope -res2.slope)
        denom_y = abs(1 + res1.slope*res2.slope)

        twist_angle = np.arctan2(num_x, denom_y)*180/np.pi

        return twist_angle, err_1, err_2

def get_interlayer_twist(data, cutoff, id_1, id_2, reference_particle_type=2, num_iter=10):

        modifier_select = ExpressionSelectionModifier(expression=f"ParticleType=={reference_particle_type}")
        modifier_delete = DeleteSelectedModifier()
        data.data.apply(modifier_select)
        data.data.apply(modifier_delete)

        index_1 = np.argwhere(data.data.particles["Particle Identifier"]==id_1)

        assert data.data.particles["Particle Identifier"][index_1] == id_1

        index_2 = np.argwhere(data.data.particles["Particle Identifier"]==id_2)
        assert data.data.particles["Particle Identifier"][index_1] == id_1

        finder = CutoffNeighborFinder(cutoff, data.data)
        layer_1_info = extract_layer_alignment(finder, index_1, num_iter)
        layer_2_info = extract_layer_alignment(finder, index_2, num_iter)
        
        return compute_interlayer_angle(layer_1_info, layer_2_info)

