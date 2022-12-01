import numpy as np
from ase import Atoms
from pathlib import Path
from grrmpy.io.read_poscar import get_cell
from grrmpy.geometry import Geometry
from grrmpy.structure.comfile import COM
from ase.units import kJ,Hartree,mol

class Structure():
    def __init__(self,atoms:Atoms=None):
        self.set_atoms(atoms)
        self._geometry = None

    @property
    def atoms(self):
        if self._atoms == Atoms():
            return None
        else:
            return self._atoms
        
    def get_atoms(self):
        return self.atoms
    
    @atoms.setter
    def atoms(self,atoms):
        if atoms is None:
            self._atoms = Atoms()
        elif type(atoms) == Atoms:
            self._atoms = atoms
        else:
            raise TypeError("atomsに設定できるのはAtomsオブジェクトです")
                
    def set_atoms(self,atoms):
        self.atoms = atoms
            
    @property
    def positions(self):
        return self._atoms.get_positions()
    
    @property
    def geometry(self):
        return self._geometry
    
    @geometry.setter
    def geometry(self,geometry):
        self._geometry = geometry
        
    def set_geometry(self,target0=None,target1=[],target2=[],mult=1.0,**kwargs):
        geometry = Geometry(self.atoms,target0,target1,target2,mult,**kwargs) # self._atomsにしない
        self.geometry = geometry
    
    @geometry.deleter
    def geometry(self):
        self._geometry = None
    
    @property  
    def cell(self):
        return self._atoms.get_cell()
    
    def get_cell(self):
        return self._atoms.get_cell()
    
    @cell.setter
    def cell(self,cell):
        self._atoms.cell = cell
        
    def set_cell(self,cell,scale_atoms:bool=False,apply_constraint:bool=True):
        if type(cell) == str or type(cell) == Path:
            cell = get_cell(cell)
            self._atoms.cell = cell
        else:
            self._atoms.set_cell(cell,scale_atoms,apply_constraint)
        
    @property
    def pbc(self):
        return self._atoms.get_pbc()
    
    def get_pbc(self):
        return self._atoms.get_pbc()
    
    @pbc.setter
    def pbc(self,pbc):
        self._atoms.pbc = pbc
        
    def set_pbc(self,pbc):
        self.pbc = pbc
    
    @property
    def mol(self):
        if self._geometry is None:
            raise RuntimeError(f"molを呼び出すには{__class__.__name__}にgeometryを設定してください")
        else:
            return self.geometry.mol
        
    @property
    def smiles(self):
        if self._geometry is None:
            raise RuntimeError(f"smilesを呼び出すには{__class__.__name__}にgeometryを設定してください")
        else:
            return self.geometry.smiles
        
    def __copy__(self):
        new_obj = self.__class__()
        new_obj._atoms = self._atoms.copy() if not self._atoms is None else None
        return new_obj
    
    def copy(self):
        return self.__copy__()
        
    def __len__(self):
        return len(self.get_atoms())
    
    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(atoms={self.atoms})"
    
    def __bool__(self):
        if self.atoms is None:
            return False
        return True
    
    def __eq__(self, other):
        return self.atoms == other.atoms
        
    def __ne__(self, other: object) -> bool:
        return not self == other
    
    @classmethod
    def fromdict(cls,dct):
        atoms = dct["atoms"]
        if atoms:
            atoms = Atoms.fromdict(atoms)
        geometry = dct["geometry"]
        if geometry:
            geometry = Geometry.fromdict(geometry)
        new_obj = cls(atoms=atoms)
        new_obj._geometry = geometry
        return new_obj
    
    def todict(self):
        atoms = self.atoms.todict() if self.atoms else None
        geometry = self._geometry.todict() if self._geometry else None
        return {"atoms":atoms,"geometry":geometry}
        
        
class EQ(Structure):
    def __init__(self,energy:float=None,atoms:Atoms=None,frozen_atoms=None):
        super().__init__(atoms)
        self._energy = energy
        self.set_frozen_atoms(frozen_atoms)
        
    @property
    def energy(self):
        return self._energy
    
    @property
    def frozen_atoms(self):
        if self._frozen_atoms == Atoms():
            return None
        else:
            return self._frozen_atoms
        
    def get_frozen_atoms(self):
        return self.frozen_atoms # self._frozen_atomsにしない
    
    @frozen_atoms.setter
    def frozen_atoms(self, atoms:Atoms):
        if type(atoms) == Atoms:
            self._frozen_atoms = atoms
        elif atoms is None:
            self._frozen_atoms = Atoms()
        else:
            raise TypeError("Atomsオブジェクトを設定してください")
        
    def set_frozen_atoms(self,data):
        """FrozenAtomsを設定する
        
        | 2種類の方法で設定できる.
        | 1つはFrozenAtomsのAtomsオブジェクトを直接与える方法
        | 1つはcomファイルのパスを与える方法

        Parameters:
        
        data: Atoms or Path object or str
            Atomsオブジェクトまたは, comファイルのパス
        """
        if type(data) == Atoms:
            self.frozen_atoms = data
        elif type(data) == str or type(data) == Path:
            com = COM(data)
            self.frozen_atoms = com._atoms_dict["Frozen Atoms"]
        elif data is None:
            self.frozen_atoms = Atoms()
        else:
            raise TypeError("Atoms,COMファイルパス,Noneのいずれかでず")
    
    def get_atoms(self,frozen_atoms=True):
        """親クラスを上書きしている"""
        if frozen_atoms and self.frozen_atoms:
            return self.atoms + self.frozen_atoms
        else:
            return self.atoms
    
    def __copy__(self):
        new_obj = super().__copy__()
        new_obj._energy = self.energy if not self.energy is None else None
        return new_obj
      
    def __str__(self):
        if self._atoms == Atoms():
            return f"{self.__class__.__name__}()"
        else:
            tokens = []
            N = len(self)
            if N <= 60:
                symbols = self._atoms.get_chemical_formula('reduce')
            else:
                symbols = self._atoms.get_chemical_formula('hill')
            tokens.append("symbols='{0}'".format(symbols))

            for name in sorted(self._atoms.arrays):
                if name in ['numbers', 'positions']:
                    continue
                tokens.append('{0}=...'.format(name))
                
            tokens.append("energy={0}".format(self.energy))

            return '{0}({1})'.format(self.__class__.__name__, ', '.join(tokens))
    
    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(energy={self.energy},atoms={self.atoms})"
    
    @classmethod
    def fromdict(cls,dct):
        energy = dct["energy"]
        frozen_atoms = dct["frozen_atoms"]
        if frozen_atoms:
            frozen_atoms = Atoms.fromdict(frozen_atoms)
        new_obj = super().fromdict(dct)
        new_obj._energy = energy
        new_obj.frozen_atoms = frozen_atoms
        return new_obj
    
    def todict(self):
        dct = super().todict()
        energy = self.energy
        frozen_atoms = self.frozen_atoms.todict() if self.frozen_atoms else None
        dct.update({"energy":energy,"frozen_atoms":frozen_atoms})
        return dct

class TS(EQ):
    def __init__(self, energy:float=None, atoms:Atoms=None, connection:list=None,frozen_atoms=None):
        super().__init__(energy, atoms)
        self.connection =connection
        self._ini_eq = None
        self._fin_eq = None
        self.set_frozen_atoms(frozen_atoms)
        
    @property
    def connection(self):
        return self._connection
    
    @connection.setter
    def connection(self, connection):
        if connection is None:
            self._connection = None
        elif type(connection) == list:
            self._connection = np.array(connection)
        elif type(connection) == np.ndarray:
            self._connection = connection
        else:
            raise TypeError("connectionsはlistまたはnp.arrayです")
     
    @property
    def ini_eq(self):
        return self._ini_eq
    
    @ini_eq.setter
    def ini_eq(self,ini_eq):
        if type(ini_eq) == EQ:
            self._ini_eq = ini_eq
        elif ini_eq is None:
            self._ini_eq = None
        else:
            raise TypeError("EQオブジェクトを設定してください")
    
    @property
    def fin_eq(self):
        return self._fin_eq
    
    @fin_eq.setter
    def fin_eq(self,fin_eq):
        if type(fin_eq) == EQ:
            self._fin_eq = fin_eq
        elif fin_eq is None:
            self._fin_eq = None
        else:
            raise TypeError("EQオブジェクトを設定してください")
    
    @property
    def forward_energy(self):
        return self.get_forward_energy()
    
    def get_forward_energy(self,unit="Hartree"):
        if self.ini_eq and self.fin_eq:
            energy = self.energy - self.fin_eq.energy
            if unit=="Hartree":
                return energy
            elif unit == "kJ/mol":
                return energy*Hartree*mol/kJ
        return None
        
    @property
    def reverse_energy(self):
        return self.get_reverse_energy()
    
    def get_reverse_energy(self,unit="Hartree"):
        if self.ini_eq and self.fin_eq:
            energy = self.energy - self.ini_eq.energy
            if unit=="Hartree":
                return energy
            elif unit == "kJ/mol":
                return energy*Hartree*mol/kJ
        return None
    
    def __str__(self):
        text = super().__str__()[:-1] #最後の")"を外す
        text += f", CONNECTION : {self.connection[0]} - {self.connection[1]})"
        return text
    
    def __call__(self): #geometryで同じ反応か比較できるようにするため
        return self.ini_eq, self.fin_eq
        
    def __copy__(self):
        energy = self.energy if not self.energy is None else None
        atoms = self.atoms.copy() if not self.atoms is None else None
        connection = self.connection.copy() if not self.connection is None else None
        ini_eq = self.ini_eq.copy() if not self.ini_eq is None else None
        fin_eq = self.fin_eq.copy() if not self.fin_eq is None else None
        new_obj = __class__(energy,atoms,connection)
        new_obj.ini_eq = ini_eq
        new_obj.fin_eq = fin_eq
        return new_obj
    
    @classmethod
    def fromdict(cls,dct):
        connection = dct["connection"]
        ini_eq = dct["ini_eq"]
        if ini_eq:
            ini_eq = EQ.fromdict(ini_eq)
        fin_eq = dct["fin_eq"]
        if fin_eq:
            fin_eq = EQ.fromdict(fin_eq)
        new_obj = super().fromdict(dct)
        new_obj.connection = connection
        new_obj.ini_eq = ini_eq
        new_obj.fin_eq = fin_eq
        return new_obj
    
    def todict(self):
        dct = super().todict()
        connection = self.connection
        ini_eq = self.ini_eq.todict() if self.ini_eq else None
        fin_eq = self.fin_eq.todict() if self.fin_eq else None
        dct.update({"connection":connection,"ini_eq":ini_eq,"fin_eq":fin_eq})
        return dct
    
    
class PT(TS):
    def __init__(self, energy=None, atoms=None, connection=None,frozen_atoms=None):
        super().__init__(energy, atoms, connection,frozen_atoms)
        
