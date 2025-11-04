import sys
sys.path.append('.')
import Bio.PDB as bp



from src.bioiain.biopython import downloadPDB
import src.bioiain as bi



downloadPDB("./test/data", "test_list", ["5JJM", "6nwl"], file_path="test/pdb_list.txt", file_format="pdb")





