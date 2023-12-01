from prody import *
from enm import *
from pylab import *
import argparse
import os

# ! Script to visualize the modes of a PDB file using VMD

def extract_filename(file_path):
    filename = os.path.basename(file_path)
    filename_without_extension, _ = os.path.splitext(filename)
    return filename_without_extension
# def main(args):
def main():
    file_name = "./DESRES-Trajectory_sarscov2-10880334-no-water-no-ion-glueCA.pdb"
    just_name = extract_filename(file_name)
    print("Starting ANM analysis on", file_name)
    enm = ENM(file_name)
    enm.filter('CA')
    print(format_colored(len(enm.atoms)))
    b = enm.getCoords()
    H = enm.getHessian()
    print("Hessian built.")
    anm = ANM("ANM analsysis on " + file_name)
    anmENM = ANM("anm with hessian from enm")

    anmENM.setHessian(H)
    prot = parsePDB(file_name)
    calphas = prot.select("protein and name CA")
    print("building Hessian...")
    anm.buildHessian(calphas, cutoff=11.)
    print("Hessian built.")
    anm.calcModes()
    print("Modes calculated.")
    anmENM.calcModes()
    print(f"Comparison of shapes: ANM {anm.getHessian().shape}, ENM {H.shape}")
    writeNMD(just_name + "FULL.nmd", anm[:3], calphas)
    writeNMD(just_name + ".nmd", anmENM[:3], enm)
    print(f"Use \033[1;31m vmd -e {just_name}FULL.nmd\033[0m to visualize the modes")
    print(f"Use \033[1;31m vmd -e {just_name}.nmd\033[0m to visualize the modes")

if __name__ == "__main__":
    # parser = argparse.ArgumentParser(description="Perform ANM analysis on a PDB file")
    # parser.add_argument("--pdb", required=True, help="Path to the PDB file for ANM analysis")
    # args = parser.parse_args()
    # main(args.pdb)
    main()
