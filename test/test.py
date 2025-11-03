import sys
import Bio.PDB as bp

sys.path.append('.')
import src.bioiain as bi
print(bi.bp)



bi.log("warning", "works")


pdb = bi.biopython.loadPDB("./test/1M2Z.pdb")
cif = bi.biopython.loadPDB("./test/1M2Z.cif")
print(bi.biopython.Structure)




print(pdb.child_class)
print(cif.__dict__.keys())




