import pandas as pd
from ase.units import eV, kJ, mol, Hartree
from ase import Atoms
from ase.constraints import FixAtoms
import pickle

#user
from grrmpy.path.functions import to_excell,to_fig, to_html, to_plotly
from grrmpy.calculator import pfp_calculator

class ReactPath():
    """
    
        Parameters:
        
        data: dict or DataFrame
            energyとnameまたはatomsとnameの2つのキーを持つ辞書.

            - 'energy': list of float
                | エネルギーをリストで与える
                | 単位はkJ/mol, eV単位で与える場合,unit='eV'とする
            - 'atoms': list of Atoms
                | Atomsのリスト
            - 'name' : list of str
                | Ex) ['EQ0','TS1','EQ3','TS3','EQ2','PT2','EQ5']
                | 'TS','PT'の文字列を含んでいる場合,それぞれグラフ描写時に赤色,緑色で表示される
                
        title: srt
            | Pathを識別するためにタイトルを付ける事が可能.
            | デフォルトはNone
        positions: list of bool
            | solidラインの位置をカスタマイズする
            | 使用方法は   を参照
            | 指定しない場合はすべてTrueリストが作成される
        unit: string
            dataの設定の際に'energy'で設定した場合,与えたエネルギーの単位を記入する
            'kJ/mol','Hartree','eV'のいずれか
        clac_funct: function object
            | dataの設定の際に'atoms'で設定した場合,calculatorを与える関数が必要.
            | デフォルトではPFP.
            
        Note:
            constraintsはFixAtomsだけfromdict,frompkl変換できる
    """
    def __init__(self,data:dict,positions:list=None,title=None,unit="eV",calc_func=pfp_calculator):
        """

        Parameters:
        
        data: dict or DataFrame
            energyとnameまたはatomsとnameの2つのキーを持つ辞書.
            
            - 'energy': list of float
                | エネルギーをリストで与える
                | 単位はkJ/mol, eV単位で与える場合,unit='eV'とする
            - 'atoms': list of Atoms
                Atomsのリスト
            - 'name' : list of str
                | Ex) ['EQ0','TS1','EQ3','TS3','EQ2','PT2','EQ5']
                | 'TS','PT'の文字列を含んでいる場合,それぞれグラフ描写時に赤色,緑色で表示される
                
        title: srt
            | Pathを識別するためにタイトルを付ける事が可能.
            | デフォルトはNone
        positions: list of bool
            | solidラインの位置をカスタマイズする
            | 使用方法は   を参照
            | 指定しない場合はすべてTrueリストが作成される
        unit: string
            dataの設定の際に'energy'で設定した場合,与えたエネルギーの単位を記入する
            'kJ/mol','Hartree','eV'のいずれか
        clac_funct: function object
            | dataの設定の際に'atoms'で設定した場合,calculatorを与える関数が必要.
            | デフォルトではPFP.
        """
        self.data = data
        self.title = title
        self.calc_func = calc_func

        if positions is None:
            positions = [True for _ in range(len(self.data))]
        self._positions = positions
        
        if 'atoms' in self.data.columns:
            for atoms in self.data["atoms"]:
                atoms.calc = self.calc_func()
            self.data["energy"] = [atoms.get_potential_energy()*mol/kJ for atoms in self.data["atoms"]]
        else:
            if unit == "eV":
                self.data["energy"] = [i*mol/kJ for i in self.data["energy"]]
            elif unit == "Hartree":
                self.data["energy"] = [i*mol/(kJ*Hartree) for i in self.data["energy"]]
        
    @property
    def data(self):
        return self._data
    
    @data.setter
    def data(self,data):
        if type(data) == dict:
            if not "name" in data.keys():
                raise KeyError("'name'キーがありません")
            self._data = pd.DataFrame(data)
        elif type(data) == pd.DataFrame:
            if not "name" in data.columns:
                raise KeyError("'name'キーがありません")
            self._data = data
        self._positions = [True for _ in range(len(self.data))]
        self._data = self._data.reset_index(drop=True)
    
    def get_energy(self):
        return self.data["energy"].to_list()
    
    def get_name(self):
        return self.data["name"].to_list()
    
    @property
    def positions(self):
        return self._positions
    
    @positions.setter
    def positions(self, positions):
        n_true = len([i for i in positions if i]) # Trueの個数
        if not n_true == len(self.data):
            raise ValueError("Trueの数がデータの数と一致しません")
        self._positions = positions
        
    def get_solid_df(self):
        solid_xini = [i*2+1 for i,b in enumerate(self.positions) if b]
        solid_xfin = [i+1 for i in solid_xini]
        std_e = self.get_energy()[0]
        solid_yini = [i-std_e for i in self.get_energy()]
        df = pd.DataFrame({"solid_xini":solid_xini,
                           "solid_xfin":solid_xfin,
                           "name":self.get_name(),
                           "solid_yini":solid_yini,
                           "solid_yfin":solid_yini
                           })
        return df
    
    def get_dot_df(self):
        solid_df = self.get_solid_df()
        dot_xini = solid_df["solid_xfin"][:-1].to_list()
        dot_xfin = solid_df["solid_xini"][1:].to_list()
        dot_yini = solid_df["solid_yini"][:-1].to_list()
        dot_yfin = solid_df["solid_yini"][1:].to_list()
        dot_name = [f"{self.data['name'][i]}-{self.data['name'][i+1]}" 
                    for i in range(len(self.data)-1)]
        df = pd.DataFrame({"dot_xini":dot_xini,
                           "dot_xfin":dot_xfin,
                           "dot_name":dot_name,
                           "dot_yini":dot_yini,
                           "dot_yfin":dot_yfin})
        return df
    
    def get_ea(self):
        """律速段階の活性化障壁と反応を出力"""
        ts = self.data['name'].str.contains('TS')
        pt = self.data['name'].str.contains('PT')
        i_ts_pt = [i for i,(t,p) in enumerate(zip(ts,pt)) if any([t,p])] # TS,PTの位置
        e = [self.data['energy'][i]-self.data['energy'][i-1] for i in i_ts_pt]
        ea = max(e)
        ea_idx = i_ts_pt[e.index(ea)]
        reac_name = [self.data['name'][ea_idx-1],self.data['name'][ea_idx],self.data['name'][ea_idx+1]]
        return max(e),reac_name
          
    def write_excel(self,outfile:str,**kwargs):
        """Excellにグラフを作成する
        
        その他の引数はreaction_path.functions.to_excell()を参照

        Parameters:
        
        outfile: str
            Excelファイル名

        """
        table_data = self.data
        solid_data = self.get_solid_df()
        dot_data = self.get_dot_df()
        to_excell(outfile,table_data,solid_data,dot_data,**kwargs)
        
    def write_html(self,outfile=None,**kwargs):
        """反応座標グラフをhtmlに保存する

        Parameters:
        
        outfile: str or Path
            | 出力先のhtmlファイル
            | Noneの場合は標準出力する
        """
        solid_data = self.get_solid_df()
        dot_data = self.get_dot_df()
        if outfile is None:
            return to_plotly(solid_data,dot_data)
        html_text = to_html(solid_data,dot_data,**kwargs)
        with open(outfile,"w") as f:
            f.write(html_text)
        
    def preview(self):
        solid_data = self.get_solid_df()
        dot_data = self.get_dot_df()
        to_fig(solid_data,dot_data,xlabel="Energy [kJ/mol]")
        
    def get_fig(self):
        solid_data = self.get_solid_df()
        dot_data = self.get_dot_df()
        return to_fig(solid_data,dot_data)
    
    def __str__(self):
        p = self.get_name()
        e = self.get_energy()
        ea = self.get_ea()
        return f"Title:{self.title}\nPath{p},Energy:{e},Ea:{ea}"
    
    def todict(self):
        energy = [e for e in self.data["energy"]]
        name = [n for n in self.data["name"]]
        data_dict = {"name":name,"energy":energy,"positions":self.positions,"title":self.title}
        if 'atoms' in self.data.columns:
            atoms = [atoms.todict() for atoms in self.data["atoms"]]
            data_dict["atoms"] = atoms
        return data_dict
    
    @classmethod
    def fromdict(cls,data):
        energy = [i for i in data["energy"]]
        name = [i for i in data["name"]]
        data_dict = {"name":name,"energy":energy}
        if "atoms" in data:
            constraints_list = []
            for atoms_dict in data["atoms"]:
                constraints_list.append(atoms_dict.pop("constraints", []))
            atoms = [Atoms.fromdict(i) for i in data["atoms"]]
            for a,c in zip(atoms,constraints_list):
                a.set_constraint(c) # ASEのバグ？によってfromdictが使えないのでconstraintsは手動で設定
            data_dict["atoms"] = atoms
        positions = data["positions"]
        title = data["title"]
        unit="kJ/mol"
        return cls(data_dict,positions,title,unit)
    
    def topkl(self,outfile):
        with open(outfile, "wb") as f:
            pickle.dump(self.todict(),f)
        
    @classmethod
    def frompkl(self,openfile):
        with open(openfile, "rb") as f:
            data_dict = pickle.load(f)
        return self.fromdict(data_dict)