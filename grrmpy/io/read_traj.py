from pathlib import Path
from natsort import natsorted
from ase.io import read,iread

def read_traj(folder):
    """
    
    | フォルダ内のtrajファイルを検索し,Atomsのリストする.
    | trajファイルは0.traj 1.traj 2.traj...のようになっている必要がある.
    | また例えば,0.traj 2.traj のように1.trajが抜けていた場合はNoneの要素になる.
    | trajの中身がimagesの場合は別の関数である,iread_trajを用いる


    Parameters:
    
    folder: str or Path
        フォルダ名

    Returns:
        list: Atomsのリスト
    """
    p = Path(folder)
    traj_files = natsorted([p for p in p.glob("*traj")],lambda x:x.name)
    max_i = int(traj_files[-1].stem)
    atoms_list = [read(p.joinpath(f"{i}.traj")) 
                  if p.joinpath(f"{i}.traj").exists() 
                  else None 
                  for i in range(max_i)]
    return atoms_list

def iread_traj(folder):
    """read_trajのimages版"""

    p = Path(folder)
    traj_files = natsorted([p for p in p.glob("*traj")],lambda x:x.name)
    max_i = int(traj_files[-1].stem)
    atoms_list = [[atoms for atoms in iread(p.joinpath(f"{i}.traj"))]
                  if p.joinpath(f"{i}.traj").exists() 
                  else None 
                  for i in range(max_i)]
    return atoms_list