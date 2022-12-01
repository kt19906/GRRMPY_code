import grrmpy
from pathlib import Path
from ase import Atoms
import ase
import pickle

def read(filename:str):
    return read_grrmpy_obj(filename)

def write(filename:str,data, format:str=None):
    if format is None:
        p = Path(filename)
        format = p.suffix[1:]
    if format == "xyz":
        write_xyz(filename,data)
    elif format == "cif":
        write_cif(filename,data)
    elif format == "vasp":
        write_vasp(filename,data)
    elif format == "vasp-xdatcar":
        write_vasp_xdatcar(filename,data)
    elif format == "traj":
        write_traj(filename,data)
    elif format == "pickle" or format == "pkl":
        write_grrmpy_obj(filename,data)
    else:
        raise Exception("formatを指定してください")

def _to_atoms(data):
    if isinstance(data,grrmpy.structure.structures.Structures):
        atoms = data.atoms_list
    elif isinstance(data,grrmpy.structure.structures.Structure):
        atoms = data.atoms
    elif type(data)==Atoms:
        atoms = data
    elif type(data) == list:
        atoms = [atoms if atoms else Atoms() for atoms in data]
    else:
        raise TypeError("Structures, EQList, TSList, Aromsリスト, PTList, Structure, EQ, TS, PT, Atomsのいずれかです")
    return atoms

def _to_1d_atoms(data):
    if isinstance(data,grrmpy.structure.structures.Structure):
        atoms = data.atoms
    elif type(data)==Atoms:
        atoms = data
    else:
        raise TypeError("Structure, EQ, TS, PT, Atomsのいずれかです")
    return atoms

def _to_2d_atoms(data):
    if isinstance(data,grrmpy.structure.structures.Structures):
        atoms = data.atoms_list
    elif type(data) == list:
        atoms = [atoms if atoms else Atoms() for atoms in data]
    else:
        raise TypeError("Structures, EQList, TSList, PTList, Structure, Aromsリストのいずれかです")
    return atoms

def write_xyz(filename,data):
    atoms = _to_atoms(data)
    ase.io.write(filename,atoms,format="xyz")
    
def write_cif(filename,data):
    atoms = _to_atoms(data)
    ase.io.write(filename,atoms,format="cif")
    
def write_vasp(filename,data):
    atoms = _to_1d_atoms(data)
    ase.io.write(filename,atoms,format="vasp")
    
def write_vasp_xdatcar(filename,data):
    atoms = _to_2d_atoms(data)
    ase.io.write(filename,atoms,format="vasp-xdatcar")

def write_traj(filename,data):
    atoms = _to_atoms(data)
    ase.io.write(filename,atoms,format="traj")
    
def write_grrmpy_obj(filename,data):
    dct = data.todict()
    with open(filename,"wb") as f:
        pickle.dump(type(data).__name__, f)
        pickle.dump(dct, f)
        
def read_grrmpy_obj(filename):
    with open(filename, "rb") as f:
        obj_type = pickle.load(f)
        dct = pickle.load(f)
    if obj_type == "GrrmData":
        return grrmpy.GrrmData.fromdict(dct)
    elif obj_type == "EQList":
        return grrmpy.structure.EQList.fromdict(dct)
    elif obj_type == "TSList":
        return grrmpy.structure.TSList.fromdict(dct)
    elif obj_type == "PTList":
        return grrmpy.structure.PTList.fromdict(dct)
    elif obj_type == "EQ":
        return grrmpy.structure.EQ.fromdict(dct)
    elif obj_type == "TS":
        return grrmpy.structure.TS.fromdict(dct)
    elif obj_type == "PT":
        return grrmpy.structure.PT.fromdict(dct)
    elif obj_type == "COM":
        return grrmpy.structure.COM.fromdict(dct)
    elif obj_type == "Geometry":
        return grrmpy.geometry.Geometry.fromdict(dct)
    elif obj_type == "Geometries":
        return grrmpy.geometry.Geometries.fromdict(dct)