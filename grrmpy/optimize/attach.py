from ase.neb import  NEB
from ase.io import Trajectory
from math import sqrt
from functools import partial 

# User
from grrmpy.calculator import pfp_calculator
from grrmpy.functions import get_diff,get_fmax

###他に参照しているところが無ければ消したい##########
def get_fmax(atoms):
    return sqrt((atoms.get_forces()**2).sum(axis=1).max())
########################################################

def optimize_eq(neb,calc_func=pfp_calculator):
    """iniとfinを最適化しながらNEB計算を行なう

    Parameters:
    
    neb: neb object
        NEBオブジェクト
    calc_func: function object
        calculatorを返す関数
    """    
    if get_fmax(neb.images[0]) > get_fmax(neb.images[1]):
        ini = neb.images[1].copy()
        ini.calc = calc_func()
        neb.images[0] = ini
    if get_fmax(neb.images[-1]) > get_fmax(neb.images[-2]):
        fin = neb.images[-2].copy()
        fin.calc = calc_func()
        neb.images[-1] = fin

     
neb_maxstep_climb_false = {
    10:0.01,
    5:0.05,
    2:0.07,
    1:0.1,
    0.1:0.2,
    0.07:0.3,
    0.06:0.4,
    0:0.5,
}

neb_maxstep_climb_true = {
    1:0.05,
    0:0.01,
}

opt_maxstep = {
    10:0.1,
    5:0.2,
    2:0.3,
    0:0.35,
}

def automate_maxstep(opt,maxstep=None):
    """maxstepの値を変更しながら最適化を行なう

    Parameters:
    
    opt: optimizer object
        Optimizerオブジェクト
    maxstep (_type_, optional)
        | {fmax:maxstep}の辞書で与える.
        | Noneの場合はNEB(climb=True,False),opt等に合わせて自動で設定する
    """
    if maxstep is None:
        if type(opt.atoms) == NEB:
            maxstep = neb_maxstep_climb_true if opt.atoms.climb else neb_maxstep_climb_false
        else:
            maxstep = opt_maxstep
            
    maxstep = sorted(maxstep.items(),key=lambda x:x[0],reverse=True)
    now_fmax = get_fmax(opt.atoms)
    for fmax,step in maxstep:
        if now_fmax > fmax:
            opt.maxstep = step
            break


def write_traj(outfile,opt,dist=0.25,mic=None,method=0):
    """最適化中の構造をある条件に従って保存する.
    
    | IRC Pathのようなものを作成することができる
    | 
    | method=0 : TS構造からの距離がdistずつ離れるたびに記録
    | method=1 : 緩和中の移動距離がdistずつ移動するたびに記録

    Parameters:
    
    outfile: str
        保存先のtrajファイルパス.既に存在するファイルを設定するとエラー
    opt: Optimizer object
        Optimizer
    dist: float
        dist Å構造が変化する度に構造を保存する.
    mic:
        | 周期境界条件で最小移動規則を適用し距離を算出する場合,True
        | Noneの場合,計算する構造が周期境界を持っている場合は自動でTrueにする.
    method: Int
        メソッド
        
    Note:
        | 他のattachとは異なりlambda文でのattachを行なわないことに注意する.
        | Exampleを参照する.
        
    Example:
    
        >>> opt = LBFGS(atoms)
        >>> opt.attach(write_traj('forward_path.traj',opt))
    """
    return partial(_write_traj, outfile, HoldAtoms(opt,dist,mic), dist, method)

def _write_traj(outfile,holdatoms,dist,method):
    """
    maxstepを変えて再計算した時にも追記できるような発展的な使用にはこれを使った方がいいかも
    """
    if holdatoms.first_write:
        mode = "w"
        holdatoms.first_write = False
    else:
        mode = "a"
    atoms1 = holdatoms.atoms
    atoms2 = holdatoms.opt.atoms
    if method == 0:
        method0(atoms1,atoms2,outfile,holdatoms,dist,mode)
    elif method == 1:
        method1(atoms1,atoms2,outfile,holdatoms,dist,mode)
    
def method0(atoms1,atoms2,outfile,holdatoms,dist,mode):
    if get_diff(atoms1,atoms2,mic=holdatoms.mic) > holdatoms.dist:
        Trajectory(outfile, mode=mode, atoms=atoms2).write()
        holdatoms.dist += dist
        
def method1(atoms1,atoms2,outfile,holdatoms,dist,mode):
    if get_diff(atoms1,atoms2,mic=holdatoms.mic) > dist:
        Trajectory(outfile, mode=mode, atoms=atoms2).write()
        holdatoms.atoms = atoms2.copy()

class HoldAtoms():
    """最適化時に構造を保持するためのクラス
    """
    def __init__(self,opt,dist,mic):
        self.first_write = True
        self.opt = opt
        self.atoms = opt.atoms.copy()
        self.dist = dist
        if mic is None:
            self.mic= True if any(self.atoms.get_pbc()) else False
        else:
            self.mic = mic