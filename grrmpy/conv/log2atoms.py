from ase import Atoms
from grrmpy.io.read_listlog import _open_file,_read_mark_of_list,_get_chemical_symbols,_int2float

def log2atoms(logfile):
    """logファイルをAtomsオブジェクトのリストに変換する

    Parameters:
    
    logfile: str
        List.logファイルのパス

    Returns:
        list of Atoms: Atomsオブジェクトのリスト
    """
    logtext = _open_file(logfile)
    hash_idx, energy_idx = _read_mark_of_list(logtext)
    atoms_list = [atoms for atoms in _log2atoms(logtext,hash_idx,energy_idx)]
    return atoms_list

def _log2atoms(logtext:str,hash_idx:list,energy_idx:list):
    chemical_symbols = _get_chemical_symbols(logtext, energy_idx)
    natoms = len(chemical_symbols)
    positions_list = (_int2float(logtext[i+1:i+natoms+1]) for i in hash_idx) # 座標を取得
    atoms_list = (Atoms(chemical_symbols,p) for p in positions_list)
    return atoms_list # ジェネレーター