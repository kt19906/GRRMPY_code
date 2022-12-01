from ase import Atoms
from itertools import combinations
# USER
from grrmpy.calculator import pfp_calculator

###作成途中
def atoms2log(name,*data,indices=None,calc_func=None):
    """AtomsオブジェクトからXXX_EQ_list.lig, XXX_TS_list.logファイルを作成する

    Parameters:
    
    name: str
        XXX_EQ_list.lig, XXX_TS_list.logのXXX部分.lig
    data:
        方法1
            **eq_list**: list of Atoms
                | EQのAtomsのリスト
            **ts_list**: list of Atoms
                | TSのAtomsのリスト
            **connections**: list of list of int
                | Ex) [[0,1],[1,2],[4,'??']]
            
            >>> atoms2log("Sample", eq_list, ts_list,connections)
        方法2
            **ini_list**: list of Atoms
                | iniのAtomsのリスト
            **fin_list**: list of Atoms
                | finのAtomsのリスト
            **ts_list**: list of Atoms
                | TSのAtomsのリスト
            
            >>> atoms2log("Sample", ini_list, fin_list, ts_list)
    indices: list of integers
        | 全ての原子を書き込むと容量が大きくなってしまうのでEQ,TS_listに保存する原子のindexを指定する.
        | Noneの場合は全ての原子が_list.logファイルに書き込まれる
    calc_func: functions object
        | Calculatorを返す関数
        | デフォルトはNone. Noneの場合は予めAtomsにclaculatorを設定しておく必要がある.
    """
    if type(data[2][0]) == Atoms:
        """方法2"""
        eq_list = data[0]
        ts_list = data[1]
        connections = data[2]
    else:
        """方法1"""
        ini_list = data[0]
        fin_list = data[1]
        ts_list = data[2]
        all_eq_list = ini_list + fin_list
        
        ## なんか
    _atoms2log(name,eq_list,ts_list,connections,indices,calc_func)

def _atoms2log(name,eq_list,ts_list,connections,indices,calc_func):
    eq_energy_list = get_energy(eq_list,calc_func)
    _write_eq_list(name,eq_list,eq_energy_list,indices)
    ts_energy_list = get_energy(ts_list,calc_func)
    _write_ts_list(name,ts_list,connections,ts_energy_list,indices)
    
def get_energy(atomslist,calc_func):
    if calc_func:
        for atoms in atomslist:
            atoms.calc = calc_func()
    return [atoms.get_potentail_enegy() for atoms in atomslist]

def _write_eq_list(name,eq_list,eq_energy_list,indices):         
    with open(f"{name}_EQ_list.log","w") as f:
        f.write("List of Equilibrium Structures\n\n")  
        for atoms,energy in zip(eq_list,eq_energy_list):
            atoms = atoms[indices]
            chemical_symbols = atoms.get_chemical_symbols()
            positions = atoms.get_positions()
            ## なんか
            
def _write_ts_list(name,ts_list,connections,ts_energy_list,indices):
    with open(f"{name}_EQ_list.log","w") as f:
        f.write("List of Transition Structures\n\n")  
        for atoms,energy,connection in zip(ts_list,ts_energy_list,connections):
            atoms = atoms[indices]
            chemical_symbols = atoms.get_chemical_symbols()
            positions = atoms.get_positions()
            ## なんか

