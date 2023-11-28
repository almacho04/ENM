from prody import *
import argparse
import os

# ! Script to visualize the modes of a PDB file using VMD

def extract_filename(file_path):
    filename = os.path.basename(file_path)
    filename_without_extension, _ = os.path.splitext(filename)
    return filename_without_extension
def main(args):
    file_name = extract_filename(args)
    print("Starting ANM analysis on", file_name)
    anm = ANM(file_name)
    prot = parsePDB(args)
    calphas = prot.select("calpha")
    print("building Hessian...")
    anm.buildHessian(calphas, cutoff=11.)
    print("Hessian built.")
    anm.calcModes()
    print("Modes calculated.")
    writeNMD(file_name + ".nmd", anm[:3], calphas)
    print(f"Use \033[1;31m vmd -e {file_name}.nmd\033[0m to visualize the modes")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform ANM analysis on a PDB file")
    parser.add_argument("--pdb", required=True, help="Path to the PDB file for ANM analysis")
    args = parser.parse_args()
    main(args.pdb)
