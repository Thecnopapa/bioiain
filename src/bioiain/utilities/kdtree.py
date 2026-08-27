import os, sys, math, json
import numpy as np

from sklearn.neighbors import KDTree



class KDT(object):
    def __init__(self, coords_or_entity, leaf_size=10, quiet=False, params=None, symops=None, com=None, mode="atoms", auto_parse_symmetry=True, **kwargs):
        from ..base import PseudoAtom
        self.entity = None
        self.mode = None
        self.has_symmetries = False
        self.quiet = quiet

        if not quiet:
            log(3, "Building KDT...")
        from ..base import Entity
        if isinstance(coords_or_entity, Entity):
            self.entity = coords_or_entity
            self.mode = mode
            if not quiet:
                log(4, f"Entity: {self.entity}")
                log(4, f"Mode: {self.mode}")
            if mode == "atoms":
                atoms = coords_or_entity.atoms(hetatm=True)
                coords_or_entity._kdtrees["atom"] = self
            elif mode == "ca":
                atoms = [res.ca for res in coords_or_entity.residues()]
                coords_or_entity._kdtrees["ca"] = self
            coords = np.array([a.coord for a in atoms], dtype=np.float64)
            if auto_parse_symmetry:
                if params is None:
                    params = self.entity.params()
                if symops is None:
                    symops = self.entity.symops()
        else:
            atoms = [PseudoAtom(c) if not isinstance(c, PseudoAtom) else c for c in coords_or_entity]
            coords = np.array([a.coord if isinstance(a, PseudoAtom) else a for a in atoms])
        if not quiet:
            log(4, f"N={len(coords)}")

        operations = [1]*len(coords)
        positions = [None]*len(coords)

        if params is not None and symops is not None:
            self.has_symmetries = True
            asu_atoms = atoms
            atoms = []
            coords = []
            operations = []
            positions = []
            for atom in asu_atoms:
                for op, (coord, pos) in atom.all(centre=com, params=params, symops=symops).items():
                    atoms.append(atom)
                    coords.append(coord)
                    operations.append(op)
                    positions.append(pos)
            if not quiet:
                log(4, f"N+symm={len(coords)}")

        self.atoms = atoms
        self.coords = coords
        self.operations = operations
        self.positions = positions

        self.tree = KDTree(self.coords, leaf_size=leaf_size)

    def __repr__(self):
        if self.entity is not None:
            return f"<bi.{self.__class__.__name__}: {self.entity} ({self.mode}) N={len(self.coords)} symmetries={self.has_symmetries}>"
        else:
            return f"<bi.{self.__class__.__name__}: N={len(self.coords)} symmetries={self.has_symmetries}>>"


    def neighbours(self, coords, n_neighbours=2, distances=False, unique=False):
        from ..base import PseudoAtom

        if isinstance(coords, PseudoAtom) or np.isscalar(coords[0]):
            coords = [coords]
        coords = np.array([a.coord if isinstance(a, PseudoAtom) else a for a in coords])
        neigh_indexes = []
        out = self._nearest(coords, n_neighbours=n_neighbours, distances=distances)
        if distances:
            neigh_distances = []
            [neigh_indexes.extend(n) for n in out[1]]
            [neigh_distances.extend(n) for n in out[0]]
            return neigh_indexes, neigh_distances
        else:
            if unique:
                [neigh_indexes.extend(n) for n in out]
                neigh_indexes = [int(i) for i in set(neigh_indexes)]
                return neigh_indexes
            else:
                return out

    def of(self, *args, **kwargs):
        return self.radius(*args, **kwargs)

    def radius(self, coords, radius=10, distances=False, unique=False):
        from ..base import PseudoAtom
        if isinstance(coords, PseudoAtom) or np.isscalar(coords[0]):
            coords = [coords]
        coords = np.array([a.coord if isinstance(a, PseudoAtom) else a for a in coords])
        neigh_indexes = []
        out = self._radius(coords, radius=radius, distances=distances)
        if distances:
            neigh_distances = []
            [neigh_indexes.extend(n) for n in out[0]]
            [neigh_distances.extend(n) for n in out[1]]
            return neigh_indexes, neigh_distances
        else:
            if unique:
                [neigh_indexes.extend(n) for n in out]
                neigh_indexes = [int(i) for i in set(neigh_indexes)]
                return neigh_indexes

            else:
                return out

    def __len__(self):
        return len(self.coords)

    def __iter__(self):
        self.i = 0
        return self

    def __next__(self):
        if self.i >= len(self):
            raise StopIteration
        r = self[self.i]
        self.i += 1
        return r
    def __getitem__(self, item):
        return {"atom": self.atoms[item], "coord": self.coords[item], "op": self.operations[item], "pos": self.positions[item]}

    def atom_of(self, item):
        return self.atoms[item]

    def coord_of(self, item):
        return self.coords[item]

    def pos_of(self, item):
        return self.positions[item]

    def op_of(self, item):
        return self.operations[item]

    def _radius(self, item, radius, distances=False):
        return self.tree.query_radius(item, r=radius, return_distance=distances)

    def _nearest(self, item, n_neighbours,  distances=False):
        return self.tree.query(item, k=n_neighbours, return_distance=distances)





