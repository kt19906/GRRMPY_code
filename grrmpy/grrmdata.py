from typing import List
import numpy as np
from copy import deepcopy
from pathlib import Path

from .structure import EQList, TSList, PTList, COM
from .geometry.geometries import Geometries

class GrrmData():
    """GRRMデータの情報をまとめたクラス
    
    | 作成方法は主に2つある
    | - \*list.logファイルを使う方法
    | - EQList,TSList,PTListを使う方法
    | - EQクラス,TSクラス,PTクラスのリストを使う方法(未実装)
    
    Parameters:
    
    eq: str or list of atoms or EQList
        \*EQ_list.logのパス または EQList
    ts: str or list of atoms or TSList
        \*TS_list.logのパス または TSList
    pt: str or list of atoms or EQList
        \*PT_list.logのパス または PTList
    comfile: str
        \*comファイルのパス
    poscar: str
        POSCARのパス
    """
    def __init__(self,eq=None,ts=None,pt=None,comfile=None,poscar=None):
        if type(eq) == EQList:
            self.__eq = eq
        else:
            self.__eq = EQList(eq,comfile,poscar)
        if type(ts) == TSList:
            self.__ts = ts
        else:
            self.__ts = TSList(ts,comfile,poscar)
            if self.ts:
                self.ts._set_ini_eq([self.eq[i] for i in self.ts.connections[:,0]])
                self.ts._set_fin_eq([self.eq[i] for i in self.ts.connections[:,1]])
        if type(pt) == PTList:
            self.__pt = pt
        else:
            self.__pt = PTList(pt,comfile,poscar)
            if self.pt:
                self.pt._set_ini_eq([self.eq[i] for i in self.pt.connections[:,0]])
                self.pt._set_fin_eq([self.eq[i] for i in self.pt.connections[:,1]])
        
    @property
    def eq(self):
        return self.__eq
    
    @eq.setter
    def eq(self,eq):
        self.set_eq(eq)
    
    def set_eq(self,eq):
        if type(eq) == EQList:
            self.__eq = eq
        else:
            raise TypeError("EQListを指定して下さい")
            
    @property
    def ts(self):
        return self.__ts
    
    @ts.setter
    def ts(self,ts):
        self.set_ts(ts)
    
    def set_ts(self,ts):
        if type(ts) == TSList:
            self.__ts = ts
        else:
            raise TypeError("TSListを指定して下さい")
    
    @property
    def pt(self):
        return self.__pt
    
    @pt.setter
    def pt(self,pt):
        self.set_pt(pt)
    
    def set_pt(self,pt):
        if type(pt) == PTList:
            self.__pt = pt
        else:
            raise TypeError("PTListを指定して下さい")
    
    @property
    def com(self):
        if self.eq:
            return self.eq.com
        if self.ts:
            return self.ts.com
        if self.pt:
            return self.pt.com
        
    @property
    def eqlog(self):
        if self.eq:
            return self.eq.log
        return None
    
    @property
    def tslog(self):
        if self.ts:
            return self.ts.log
        return None
    
    @property
    def ptlog(self):
        if self.pt:
            return self.pt.log
        return None
       
    @property
    def cells(self):
        return self.eq.get_cells()
    
    def get_cells(self):
        return self.eq.get_cells()
    
    @cells.setter
    def cells(self, cells):
        if not self.eq.atoms_list is None:
            self.eq.cells = cells
        if not self.ts.atoms_list is None:
            self.ts.cells = cells
        if not self.pt.atoms_list is None:
            self.pt.cells = cells
            
    def set_cell(self, cell, scale_atoms=False, apply_constraint=True):
        if not self.eq.atoms_list is None:
            self.eq.set_cells(cell,scale_atoms,apply_constraint)
        if not self.ts.atoms_list is None:
            self.ts.set_cells(cell,scale_atoms,apply_constraint)
        if not self.pt.atoms_list is None:
            self.pt.set_cells(cell,scale_atoms,apply_constraint)
            
    @property
    def pbcs(self):
        return self.eq.get_pbcs()
    
    def get_pbcs(self):
        return self.eq.get_pbcs()
    
    @pbcs.setter
    def pbcs(self,pbc):
        if not self.eq.atoms_list is None:
            self.eq.pbcs(pbc)
        if not self.ts.atoms_list is None:
            self.ts.pbcs(pbc)
        if not self.pt.atoms_list is None:
            self.pt.pbcs(pbc)
    
    def set_pbc(self,pbc):
        self.pbc = pbc
            
    def search_path(ini:int,fin:int=None,group:bool=True,pt:int=1,priority:int=0,pseudo_energy:bool=True):
        """反応を検索する
        
        Parameters:
        
        ini: int
            group=Trueの場合,始状態のグループ番号, group=Falseの場合,始状態のEQ番号
        fin: int
            | group=Trueの場合,終わり状態のグループ番号, group=Falseの場合,始状態のEQ番号
            | Noneの場合全ての経路を検索する
        group: bool
            | Trueの場合,Group化を考慮して経路検索を行なう.
            | 同じgroup内のEQ間は自由に行き来できる
        pt: int
            | - pt = 0
            |     TSのみで検索する
            | - pt = 1
            |     PTも含め検索するが,同じEQ間でTSもPTも存在する場合にはTSを優先する
            | - pt = 2
            |     PTも含め検索するが,同じEQ間でTSもPTも存在する場合にはよりエネルギーの低い方を優先する
        pseudo_energy: bool
            | Trueの場合,
        """
        pass
    
    def __add__(self,other):
        new_obj = self.__class__()
        new_obj.set_eq(self.eq + other.eq)
        ts_connections = np.concatenate([self.ts.connections, other.ts.connections+len(self.eq)], 0) 
        pt_connections = np.concatenate([self.pt.connections, other.pt.connections+len(self.eq)], 0) 
        new_obj.set_ts(self.ts + other.ts)
        new_obj.ts.connections = ts_connections
        new_obj.set_pt(self.pt + other.pt)
        new_obj.pt.connections = pt_connections
        return new_obj
    
    def __iadd__(self,other):
        return self + other
        
    def __copy__(self):
        new_obj = self.__class__() 
        new_obj.eq = self.eq.copy()
        new_obj.ts = self.ts.copy()
        new_obj.pt = self.pt.copy()
        return new_obj
    
    def copy(self):
        return self.__copy__()
    
    def __getitem__(self,item):
        if item == "EQ" or item == "eq":
            return self.eq
        elif item == "TS" or item == "ts":
            return self.ts
        elif item == "PT" or item == "pt":
            return self.pt
        else:
            raise KeyError("'EQ','TS','PT'のいずれかです")
        
    def __repr__(self) -> str:
        tokens = []
        for eq,ts,pt in zip(self.eqlog,self.tslog,self.ptlog):
            eq = f"'{eq}'" if eq is not None else None
            ts = f"'{ts}'" if ts is not None else None
            pt = f"'{pt}'" if pt is not None else None
            com = f"'{self.com}'" if self.com is not None else None
            tokens.append(f"{self.__class__.__name__}(eqlog={eq},tslog={ts},ptlog={pt},com={com})")
        return " + ".join(tokens)
    
    def __bool__(self):
        if self.eq:
            return True
        return False
    
    @classmethod
    def fromdict(cls,dct):
        eq = dct["eq"]
        ts = dct["ts"]
        pt = dct["pt"]
        new_obj = cls()
        new_obj.eq = EQList.fromdict(eq)
        new_obj.ts = TSList.fromdict(ts)
        new_obj.pt = PTList.fromdict(pt)
        return new_obj
        
    def todict(self):
        eq = self.eq.todict()
        ts = self.ts.todict()
        pt = self.pt.todict()
        return {"eq":eq,"ts":ts,"pt":pt}