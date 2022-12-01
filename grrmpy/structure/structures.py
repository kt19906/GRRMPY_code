
from ..io.read_listlog import _open_file,_read_mark_of_list,_logtext2energies,_read_connections
from ..io.read_poscar import get_cell,get_cell_and_pbc
from ..conv.log2atoms import _log2atoms
from .structure import EQ,TS,PT,Structure,COM
import grrmpy.geometry.geometries as gg
import numpy as np
import pandas as pd
from pathlib import Path
import copy
from ase import Atoms
from ase.units import kJ,Hartree,mol

"""
EQListのプロパティはは実質的にself._strcturesとself._geometries,comの3つのみ
その他のプロパティはこの2つのプロパティから生成される

マージした際,COMオブジェクトは1個目のものが採用される
マージの際, Atomsの原子の数や種類などの確認は行なわない
"""

class Structures():
    _element = Structure
    def __init__(self,atoms_list=None):
        self._geometries = None
        self._strctures = None
        if atoms_list is not None:
            self._check_same_len(atoms_list)
            self._strctures = [self._element(atoms) for atoms in atoms_list]
    
    def _check_same_len(self,atoms_list):
        """atoms_listの原子が全て同じ数か確認"""
        len_atoms = np.array([len(atoms) for atoms in atoms_list if atoms is not None and atoms!= Atoms()])
        return np.all(len_atoms > len_atoms[0]) # 全部同じ原子数の場合True
        
    @property
    def atoms_list(self):
        if self._strctures is None:
            return None
        else:
            return [eq.atoms for eq in self._strctures]
    
    def get_atoms_list(self):
        return self.atoms_list
    
    @property
    def positions_list(self):
        if self._strctures is None:
            return None
        else:
            return np.array([eq.positions for eq in self._strctures])
    
    @property
    def geometries(self):
        return self._geometries
    
    @geometries.setter
    def geometies(self,geometries):
        if geometries is None:
            self._geometries = None
        elif type(geometries) == gg.Geometries:
            self._geometries = geometries
        else:
            raise TypeError("Geometriesオブジェクトを設定してください.")
        
    @geometries.deleter
    def geometries(self):
        self._geometries = None
        
    def set_geometries(self,target0=None,target1=[],target2=[],mult=1.0,**kwargs):
        geometries = gg.Geometries(self.atoms_list,target0,target1,target2,mult,**kwargs)
        self.geometies = geometries
        
    @property
    def cells(self):
        return self._strctures[0].get_cell()
    
    def get_cells(self):
        return self._strctures[0].get_cell()
    
    @cells.setter
    def cells(self,cell):
        for structure in [structure for eq in self._strctures]:
            structure.cell = cell
    
    def set_cells(self,cell, scale_atoms=False, apply_constraint=True):
        if type(cell) == str or type(cell) == Path:
            """cell=POSCARファイルを指定した場合"""
            cell = get_cell(cell)
            for structure in [structure for structure in self._strctures]:# in self.atoms_listではダメ(Atoms()でなくNoneのため)
                structure.cell = cell
        else:
            for structure in [structure for structure in self._strctures]: 
                structure.set_cell(cell,scale_atoms,apply_constraint)
                
    def _set_cells_and_pbcs(self,poscar=None):
        if poscar:
            cell,pbc = get_cell_and_pbc(poscar)
            for structure in [structure for structure in self._strctures]: 
                structure.cell = cell
                structure.pbc = pbc
                
    @property
    def pbcs(self):
        return self._strctures[0].get_pbc()
    
    def get_pbcs(self):
        return self._strctures[0].get_pbc()
    
    @pbcs.setter
    def pbcs(self,pbc):
        for structure in [structure for structure in self._strctures]: 
            structure.pbc = pbc
            
    def set_pbcs(self,pbc):
        self.pbcs = pbc
            
    @property
    def mols(self):
        if self.geometries is None:
            raise RuntimeError(f"molsを呼び出すには{self.__class__.__name__}にgeometriesを設定して下さい")
        else:
            return [geo.mols for geo in self.geometries]
        
    @property
    def smileses(self):
        if self.geometries is None:
            raise RuntimeError(f"smilesesを呼び出すには{self.__class__.__name__}にgeometriesを設定して下さい")
        else:
            return [geo.smiles for geo in self.geometries]
        
    @property
    def group(self):
        if self.geometries is None:
            raise RuntimeError(f"groupを呼び出すには{self.__class__.__name__}にgeometriesを設定して下さい")
        else:
            return self.geometries.group
        
    @property
    def cluster(self):
        if self.geometries is None:
            raise RuntimeError(f"clusterを呼び出すには{self.__class__.__name__}にgeometriesを設定して下さい")
        else:
            return self.geometries.cluster
        
    def build_from_structure(self,structure_list):
        """EQオブジェクトのリストからEQListオブジェクトを作成する
        
        Copyされない(ミュータブル)

        Parameters:
        
        structure_list: list of EQ
            EQオブジェクトを格納したリスト

        Returns:
            EQList: EQListオブジェクト
        """
        new_obj = self.__class__()
        new_obj._strctures = structure_list
        return new_obj
      
    def append(self,strcture):
        """EQオブジェクトを追加する

        Parameters:
        
        strcture: EQ
            EQオブジェクト

        Note:
            geometriesプロパティはリセットされる(Noneになる)
        """
        if type(strcture) == self._element:
            self._strctures.append(copy.copy(strcture))
            self.del_geometries()
        else:
            raise TypeError(f"{self._element}オブジェクトを指定してください")
    
    def insert(self,index:int,strcture):
        """指定した位置にEQオブジェクトを挿入する

        Parameters:
        
        strcture: EQ
            EQオブジェクト
        index: int
            挿入する位置. Index番号

        Note:
            geometriesプロパティはリセットされる(Noneになる)
        """
        if type(strcture) == self._element:
            self._strctures.insert(index,copy.copy(strcture))
            self.del_geometries()
        else:
            raise TypeError(f"{self._element}オブジェクトを指定してください")
        
    def __getitem__(self, item):
        if isinstance(item, slice):
            return self.build_from_structure(self._strctures[item])
        item = np.array(item)
        if item.ndim == 0 and item.dtype == int:
            """index番号で指定"""
            return self._strctures[item]
        elif item.ndim == 1 and item.dtype == int:
            """ファンシーインデックスで指定"""
            return self.build_from_structure([self._strctures[i] for i in item])
        elif item.ndim == 1 and item.dtype == bool:
            """ブーリアンインデックスで指定"""
            if len(item) == len(self):
                return self.build_from_structure([self._strctures[i] for i,b in enumerate(item) if b])
            else:
                raise IndexError("ブーリアンインデックスの要素数が一致しません")
        else:
            raise IndexError("引数はint, list(ブーリアンorファンシー), sliceのいずれかです")
        
    def __len__(self):
        """構造の数(EQの数)を返す"""
        return len(self._strctures)
    
    def __iter__(self):
        return iter(self._strctures)
    
    def __copy__(self):
        new_obj = self.__class__()
        new_obj.log = self.log.copy()
        new_obj._strctures = [copy.copy(strcut) for strcut in self._strctures]
        return new_obj
    
    def copy(self):
        """複製する"""
        return self.__copy__()
    
    def __add__(self,other):
        """データをマージする"""
        if type(other)==type(self):
            structures1 = [copy.copy(strcuture) for strcuture in self._strctures]
            structures2 = [copy.copy(strcuture) for strcuture in other._strctures]
            new_obj = self.__class__()
            new_obj.log = self.log.copy() + other.log.copy()
            new_obj._strctures = structures1 + structures2
            return new_obj
        else:
            raise TypeError(f"{self.__class__.__name__}と{other.__class__.__name__}でマージはできません")
        
    def __iadd__(self,other):
        return self + other
    
    def __bool__(self):
        if self._strctures:
            return True
        return False
    
    @classmethod
    def fromdict(cls,dct):
        strctures = dct["strctures"]
        if strctures:
            strctures = [cls._element.fromdict(structure) for structure in strctures]
        geometries = dct["geometries"]
        if geometries:
            geometries = gg.Geometries.fromdit(geometries)
        new_obj = cls()
        new_obj._strctures = strctures
        new_obj._geometries = geometries
        return new_obj
    
    def todict(self):
        strctures = [strcture.todict() for strcture in self._strctures] if self._strctures else None
        geometries = self._geometries.todict() if self._geometries else None
        return {"strctures":strctures,"geometries":geometries}
    
   
class EQList(Structures):
    """\*EQ_list.logの情報をまとめたクラス
    
    Parameters:
    
    log: str
        \*EQ_list.logのパス
        
        
    Properties:
    
    energies: list of float
        logファイルから読み取ったエネルギー
    atoms_list: list of Atoms
        Atomsオブジェクトのリスト
    positions_list:
        座標のリスト
    geometries: Geometries
        Geometriesオブジェクト
    mols:
        Molオブジェクトのリスト
    smileses: 
        SMILESのリスト
    group: list of integers
        | group名(数字)のリスト,同じ数字の構造は同じ構造をしている
        | 例えば [0,0,1,2,2,3] の場合, EQ0とEQ1, EQ3とEQ4がそれぞれ同じ構造である
    cluster: 2D list of integers
        | group毎にeq番号をまとめる. 各要素内の構造は同じ構造をしている
        | 例えば [[0,1],[2],[3,4],[5]] の場合,EQ0とEQ1, EQ3とEQ4がそれぞれ同じ構造である
        
    Note:
        mols, smileses, group, clusterのプロパティーはgeometriesにGeometriesを設定してからでないと呼び出ない
    """
    _element = EQ
    def __init__(self,log=None,com=None,poscar=None):
        self.log = [log]
        self._strctures = None
        self._geometries = None
        self.com = COM(com)
        if log:
            logtext = _open_file(log) # list.logファイルを読み込む
            self.file_check(logtext[0]) # logファイルの1行目からEQ_list.logであることを確認
            hash_idx, energy_idx = _read_mark_of_list(logtext) # logファイル中の"#",や"Engergy="の位置を取得
            atoms_list_gen = _log2atoms(logtext,hash_idx,energy_idx) # エネルギーを返すジェネレーター
            energies_gen = _logtext2energies(logtext, energy_idx) # 座標を返すジェネレーター
            self._strctures = self.build_structure_obj(energies_gen, atoms_list_gen, logtext) # EQクラスのリスト
            self._set_cells_and_pbcs(poscar)

    def file_check(self,text):
        token = ""
        if text != "List of Equilibrium Structures\n": # EQファイルであるかチェック    
            if text == "List of Path Top (Approximate TS) Structures\n":
                token = "PT_list.logを指定している??"
            elif text == "List of Transition Structures\n":
                token = "TS_list.logを指定している??"
            raise Exception(f"EQ_list.logファイルを読み込めません\n{token}")
        
    def get_atoms_list(self, frozen_atoms=True):
        """親クラスを上書きしている"""
        return [structure.get_atoms(frozen_atoms) for structure in self._strctures]
        
    @property
    def frozen_atoms_list(self):
        return self[0].frozen_atoms
    
    def get_frozen_atoms_list(self):
        return self.frozen_atoms_list
    
    @frozen_atoms_list.setter
    def frozen_atoms_list(self, atoms:Atoms):
        if type(atoms) != Atoms:
            raise TypeError("Atomsオブジェクトを設定して下さい")
        for structure in self._strctures:
            structure.frozen_atoms = atoms
        self.com = COM()
            
    def set_frozen_atoms_list(self,data):
        if type(data) == str or type(data) == Path:
            self.com = COM(data)
            fozen_atoms = self.com._atoms_dict["Frozen Atoms"]
            for structure in self._strctures:
                structure.frozen_atoms = fozen_atoms
        else:
            self.frozen_atoms_list = data
            
    @property
    def energies(self):
        return np.array([eq.energy for eq in self._strctures])
        
    def build_structure_obj(self, energies_gen, atoms_list_gen, _): # _にはlogtextが入るがこれはTSList,PTList用(子クラス)
        return [self._element(energy,atoms,self.com.frozen_atoms) for energy,atoms in zip(energies_gen,atoms_list_gen)]

    @property
    def summary(self):
        if self.geometies is None:
            group = [None for _ in range(len(self))]
        else:
            group = self.group
        summary = pd.DataFrame(
            data = {"node":[i for i in range(len(self))],
                    "name":[f"EQ{i}" for i in range(len(self))],
                    "group":group,
                    "E/Hartree":self.energies,
                    "E/kJmol-1":[i*Hartree*mol/kJ for i in self.energies]
                    }
        )
        return summary
    
    def __str__(self):
        tokens = []
        if self._strctures is not None:
            N = len(self._strctures[0])
            if N <= 60:
                symbols = self._strctures[0].atoms.get_chemical_formula('reduce')
            else:
                symbols = self._strctures[0].atoms.get_chemical_formula('hill')
            tokens.append("symbols='{0}'".format(symbols))

            for name in sorted(self._strctures[0].atoms.arrays):
                if name in ['numbers', 'positions']:
                    continue
                tokens.append('{0}=...'.format(name))

            if len(self) > 5:
                tokens.append(f"energies={str(self.energies[:5])[:-1]}...]")
            else:
                tokens.append(f"energies={self.energies}")

        return '{0}({1})'.format(self.__class__.__name__, ', '.join(tokens))
    
    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(log={self.log},com={self.com.file})"
    
    def __add__(self,other):
        new_obj = super().__add__(other)
        new_obj.log = self.log.copy() + other.log.copy()
        new_obj.com = self.com.copy()
        return new_obj
    
    def __copy__(self):
        new_obj = super().__copy__()
        new_obj.com = self.com.copy()
        return new_obj
    
    @classmethod
    def fromdict(cls,dct):
        new_obj = super().fromdict(dct)
        new_obj.com = COM.fromdict(dct["com"])
        return new_obj
    
    def todict(self):
        dct = super().todict()
        dct.update({"com":self.com.todict()})
        return dct
          
                      
class TSList(EQList):
    _element = TS
    def __init__(self,log=None,com=None,poscar=None):
        super().__init__(log,com,poscar)
        
    @property
    def connections(self):
        return np.array([strcture.connection for strcture in self._strctures])
    
    @connections.setter
    def connections(self,connections):
        connections = np.array(connections)
        for strcture,connection in zip(self._strctures,connections):
            strcture.connection = connection
        
    def file_check(self,text):
        token = ""
        if text != "List of Transition Structures\n": # TSファイルであるかチェック
            if text == "List of Equilibrium Structures\n":
                token = "EQ_list.logを指定している??"
            elif text == "List of Path Top (Approximate TS) Structures\n":
                token = "PT_list.logを指定している??"
            raise Exception(f"TS_list.logファイルを読み込めません\n{token}")
        
    def build_structure_obj(self, energies_gen, atoms_list_gen, logtext):
        connections = _read_connections(logtext)
        return [self._element(energy,atoms,c,self.com.frozen_atoms) for energy,atoms,c in zip(energies_gen,atoms_list_gen,connections)]
    
    def _set_ini_eq(self,ini_eq_list:list):
        for structure,ini_eq in zip(self._strctures,ini_eq_list):
            structure.ini_eq = ini_eq
            
    def _set_fin_eq(self,fin_eq_list:list):
        for structure,fin_eq in zip(self._strctures,fin_eq_list):
            structure.fin_eq = fin_eq
    
    @property
    def summary(self):
        summary = pd.DataFrame(
            data = {"edge":[i for i in range(len(self))],
                    "name":[f"{self._element.__name__}{i}" for i in range(len(self))],
                    "source":self.connections[:,0],
                    "target":self.connections[:,1],
                    "E/Hartree":self.energies,
                    "E/kJmol-1":[i*Hartree*mol/kJ for i in self.energies],
                    "forward/kJmol-1":[strcture.get_forward_energy("kJ/mol") for strcture in self._strctures],
                    "reverse/kJmol-1":[strcture.get_reverse_energy("kJ/mol") for strcture in self._strctures],
                    }
        )
        return summary
            
    def __str__(self):
        text = super().__str__()[:-1]
        if len(self) > 5:
            text += f", connections = {str(self.connections[:5])[:-1]},...])"
        else:
            text += f", connections = {self.connection})"
        return text
    
    def insert(self, index:int, strcture):
        self.connections[:index] += 1
        return super().insert(index, strcture)
            
        
class PTList(TSList):
    _element = PT
    def __init__(self,log=None,com=None,poscar=None):
        super().__init__(log,com,poscar)     
        
    def file_check(self,text):
        token = ""
        if text != "List of Path Top (Approximate TS) Structures\n": # PTファイルであるかチェック
            if text == "List of Equilibrium Structures\n":
                token = "EQ_list.logを指定している??"
            elif text == "List of Transition Structures\n":
                token = "TS_list.logを指定している??"
            raise Exception(f"TS_list.logファイルを読み込めません\n{token}")
        
