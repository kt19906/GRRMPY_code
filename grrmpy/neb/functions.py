from ase.geometry import find_mic
import numpy as np

# User
from grrmpy.calculator import pfp_calculator
from grrmpy.functions import to_html_energy_diagram, calc_curvature_circle_fitting

def to_html_nebgraph(neb_obj,calc_func=pfp_calculator,
                     full_html=False,
                     unit="kJ/mol",
                     title="NEB Energy Diagram",
                     xaxis_title="Reaction Coordinate",
                     yaxis_title=None,
                     **kwargs):
    """NEBのエネルギーダイアグラムのhtmlテキストを作成する
    
    Parameters:
    
    unit: string
        'eV', 'kJ/mol', 'Hartree', 'kcal/mol'のいずれか
    yaxis_title: string
        Noneの場合,Energy({unit})
    kwargs:
        plotpyのto_htmlの引数を参照
    """
    return to_html_energy_diagram(neb_obj.images,
                                  calc_func,
                                  full_html,
                                  unit,title,
                                  xaxis_title,
                                  yaxis_title,
                                  **kwargs)
    
def insert_image(neb, mic=False, i=None, a=0.01, clac_func=pfp_calculator):
    """TS周りに新たにimageを挿入する.
    
    TS構造の両側に新たなImageを挿入する.
    
    Parameters:
    
    neb: NEB object
        NEBオブジェクト
    mic: bool
        周期境界条件で最小画像規則を使用する場合は,True
    i: int or None
        | i番目のimageの両側に新たなイメージを作成する.
        | Noneの場合,imax(TS)に挿入する.
    a: float or list of float
        | 挿入する構造のTS構造との変位
        | i番目の構造からa[Å]離れた構造を挿入する
        | 2要素のリストで与えた場合,i-1番目に1番目要素,i+1番目に2番目の要素を適用する
    clac_func: fnction object
        claculatorを返す関数. デフォルトはpfpのcalculator
    """
    if i is None:
        i = neb.imax
    if type(a) == int or type(a) == float:
        a = [a,a]
    elif type(a) == list or type(a) == tuple:
        a = a
    else:
        raise TypeError("aはfloatまたは要素数2のリストです")
    images = neb.images
    ts = neb.images[i].copy()
    def insert(idx,ins_idx,a):
        """
        idx: tsの隣の構造のidx番号
        ins_idx: 挿入する位置
        """
        pos1 = images[idx].get_positions()
        pos2 = ts.get_positions()
        d = pos2 - pos1
        if mic:
            d = find_mic(d, ts.get_cell(), ts.pbc)[0]
        n = np.linalg.norm(d)/a
        d /= (n - 1.0)
        new_pos = pos1 + (n-2) * d # n-2(ts付近)
        unconstrained_image = ts.copy()
        unconstrained_image.set_positions(new_pos,apply_constraint=True)
        unconstrained_image.calc = clac_func()
        images.insert(ins_idx,unconstrained_image)
    if i != len(images)-1:
        insert(i+1,i+1,a[1])
    if i != 0:
        insert(i-1,i,a[0])


            