import os, sys, math, json


import numpy as np
from sklearn.neighbors import KDTree



atomic_radii = {
	"default": {
	        "H": 1.200,
	        "HE": 1.400,
	        "C": 1.700,
	        "N": 1.550,
	        "O": 1.520,
	        "F": 1.470,
	        "NA": 2.270,
	        "MG": 1.730,
	        "P": 1.800,
	        "S": 1.800,
	        "CL": 1.750,
	        "K": 2.750,
	        "CA": 2.310,
	        "NI": 1.630,
	        "CU": 1.400,
	        "ZN": 1.390,
	        "SE": 1.900,
	        "BR": 1.850,
	        "CD": 1.580,
	        "I": 1.980,
	        "HG": 1.550,
	        "other": 2.000,
	}
}




class SASA(object):
	def __init__(self, ball_radius=1.40, n_points=100, radii_dict="default", **kwargs):

		assert  ball_radius > 0
		assert n_points > 1

		self.ball_radius = ball_radius
		self.n_points = n_points
		global atomic_radii
		self.radii_dict = atomic_radii[radii_dict]

		self._sphere = self._compute_sphere()


        



	def _compute_sphere(self):
		n = self.n_points

		dl = np.pi * (3 - 5**0.5)
		dz = 2.0 / n

		longitude = 0
		z = 1 - dz / 2

		coords = np.zeros((n, 3), dtype=np.float32)
		for k in range(n):
			r = (1 - z * z) ** 0.5
			coords[k, 0] = math.cos(longitude) * r
			coords[k, 1] = math.sin(longitude) * r
			coords[k, 2] = z
			z -= dz
			longitude += dl

		return coords


	def compute(self, entity, **kwargs):
		print("Computing ASA...")

		
		kdt = KDT(entity, **kwargs)

		radii_list = []
		n_atoms = 0
		for a in kdt.atoms:
			n_atoms += 1
			if a in self.radii_dict:
				radii_list.append(self.radii_dict[a.element])
			else:
				radii_list.append(self.radii_dict["other"])




		radii = np.array(radii_list, dtype=np.float64)

		radii += self.ball_radius
		twice_maxradii = np.max(radii) * 2


		asa_array = np.zeros((n_atoms, 1), dtype=np.int64)
		#ptset = set(range(self.n_points))

		#print(ptset)

		for i in range(n_atoms):
			#exposed_points = ptset.copy()

			i_radii = radii[i]
			s_on_i = (np.array(self._sphere, copy=True) * i_radii) + kdt.coords[i]
			
			sphere_kdt = KDT(s_on_i, leaf_size=10)

			i_neighbours = kdt.of(coords=s_on_i, radius=twice_maxradii, distances=False, unique=True)


			#print(i, len(i_neighbours))
			neighbour_radii = np.array([radii_list[j] for j in i_neighbours])
			neighbour_coords = np.array([kdt.coords[j] for j in i_neighbours])

			overlap_indexes = sphere_kdt.of(neighbour_coords, radius=neighbour_radii, unique=True)

			#print(i, len(overlap_indexes))
			
			asa_array[i] = self.n_points - len(overlap_indexes)


		f = radii * radii * (4 * np.pi / self.n_points)
		asa_array = asa_array * f[:, np.newaxis]
		#print(asa_array)


		for atom, asa in zip(entity.atoms(), asa_array):
			atom.set_misc("SASA", float(asa[0]))
		return asa_array





class KDT(object):
    def __init__(self, coords_or_entity, leaf_size=10, **kwargs):
    	from .base import BiopythonOverlayClass
    	if isinstance(coords_or_entity, BiopythonOverlayClass):
    		atoms = coords_or_entity.atoms()
    		coords = np.array([a.coord for a in atoms], dtype=np.float64)
    	else:
    		atoms = None
    		coords = np.array(coords_or_entity)

    	self.atoms = atoms
    	self.coords = coords

    	self.tree = tree = KDTree(self.coords, leaf_size=leaf_size)




    def of(self, coords, radius=10, distances=False, unique=False):

        if np.isscalar(coords[0]):
            coords = [coords]
        coords = np.array(coords)
        neigh_indexes = []
        out = self(coords, radius=radius, distances=distances)
        if distances:
            neigh_distances = []
            [neigh_indexes.extend(n) for n in out[0]]
            [neigh_indexes.extend(n) for n in out[1]]
            return neigh_indexes, neigh_distances
        else:
            if unique:
                [neigh_indexes.extend(n) for n in out]
                neigh_indexes = [int(i) for i in set(neigh_indexes)]
                return neigh_indexes

            else:
                return out



    def atom_of(self, item):
        assert self.atoms is not None
        return self.atoms[item]

    def coord_of(self, item):
        return self.coords[item]

    def __call__(self, item, radius, distances=False):
        return self.tree.query_radius(item, r=radius, return_distance=distances)









