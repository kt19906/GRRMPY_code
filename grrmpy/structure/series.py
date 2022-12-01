import numpy as np
import pandas as pd
from ase import Atoms

class A:
    def a(self):
        print("ddd")

class Series(np.ndarray):
    def __new__(cls, data):
        new_data = []
        for atoms in data:
            if atoms is not None and type(atoms)!=Atoms:
                raise TypeError("AtomsオブジェクトまたはNoneを要素とするリストしか受け付けません")
            elif atoms is None:
                atoms = Atoms()
            else:
                atoms = atoms
            new_data.append(atoms)
        obj = np.asarray(new_data,dtype=object).view(cls)
        return obj
    
    @property
    def calc(self):
        return self._calc
    
    @calc.setter
    def calc(self,calc):
        self._calc = calc
        for atoms in self:
            atoms.clac = self.calc()
            
    def set_calculator(self, calc):
        return self.calc(calc)
    
    def get_chemical_symbols(self):
        return [atoms.get_chemical_symbols() for atoms in self]
    
    