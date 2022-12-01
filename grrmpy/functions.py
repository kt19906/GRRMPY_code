import math
from math import sqrt
import numpy as np
import plotly.graph_objects as go
import plotly.io as pyi
import networkx as nx
from ase.units import kJ,Hartree,mol,kcal
from ase.geometry import find_mic
from ase.geometry.analysis import Analysis
from ase.build.rotate import rotation_matrix_from_points
from ase.neighborlist import build_neighbor_list,natural_cutoffs
# User
from grrmpy.calculator import pfp_calculator

def get_fmax(atoms):
    """fmaxを返す"""
    return sqrt((atoms.get_forces()**2).sum(axis=1).max())

def draw_graph(images,
               calc_func=pfp_calculator,
               unit="kJ/mol",
               title="Energy Diagram",
               xaxis_title="",
               yaxis_title=None):
    """エネルギーダイアグラムのグラフを作成する
    
    Parameters:
    
    images: list of Atoms
        Atomsオブジェクトのリスト
    full_html: bool
        | <html>タグたか始まる完全なhtmlテキストを出力する場合,True
        | Falseの場合<div>タグのみ
    unit: string
        'eV', 'kJ/mol', 'Hartree', 'kcal/mol'のいずれか
    xaxis_title: str
        x軸のタイトル  
    yaxis_title: string
        Noneの場合,Energy({unit})
    
    """        
    x = [i for i in range(len(images))]
    try:
        y = [atoms.get_potential_energy() for atoms in images]
    except:
        for image in images:
            image.calc = calc_func()
        y = [atoms.get_potential_energy() for atoms in images]
           
    y = [i-y[0] for i in y] # iniのエネルギーを0スタートで表記
    if unit == "kJ/mol":
        y = [i*mol/kJ for i in y]
    elif unit == "Hartree":
        y = [i/Hartree for i in y]
    elif unit == "kcal/mol":
        y = [i*mol/kcal for i in y]
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(x=x, y=y,
                   line=dict(width=2,color="black"))
    )
    if yaxis_title is None:
        yaxis_title = f"Energy({unit})"
    fig.update_layout(
        title=title,
        xaxis_title=xaxis_title,
        yaxis_title=unit
    )
    return fig
    
def to_html_energy_diagram(images,calc_func=pfp_calculator,
                           full_html=False,
                           unit="kJ/mol",
                           title="Energy Diagram",
                           xaxis_title="",
                           yaxis_title=None,
                           **kwargs):
    """エネルギーダイアグラムのhtmlテキストを作成する
    
    Parameters:
    
    images: list of Atoms
        Atomsオブジェクトのリスト
    full_html: bool
        | <html>タグたか始まる完全なhtmlテキストを出力する場合,True
        | Falseの場合<div>タグのみ
    unit: string
        'eV', 'kJ/mol', 'Hartree', 'kcal/mol'のいずれか
    xaxis_title: str
        x軸のタイトル  
    yaxis_title: string
        Noneの場合,Energy({unit})
    
    """
    fig = draw_graph(
        images,
        calc_func=calc_func,
        unit=unit,
        title=title,
        xaxis_title=xaxis_title,
        yaxis_title=yaxis_title)
    return pyi.to_html(fig,full_html=full_html,**kwargs)


def circle_fitting(x, y):
    """Circle Fitting with least squared
        input: point x-y positions  
        output  cxe x center position
                cye y center position
                re  radius of circle 
    """

    sumx = sum(x)
    sumy = sum(y)
    sumx2 = sum([ix ** 2 for ix in x])
    sumy2 = sum([iy ** 2 for iy in y])
    sumxy = sum([ix * iy for (ix, iy) in zip(x, y)])

    F = np.array([[sumx2, sumxy, sumx],
                  [sumxy, sumy2, sumy],
                  [sumx, sumy, len(x)]])

    G = np.array([[-sum([ix ** 3 + ix * iy ** 2 for (ix, iy) in zip(x, y)])],
                  [-sum([ix ** 2 * iy + iy ** 3 for (ix, iy) in zip(x, y)])],
                  [-sum([ix ** 2 + iy ** 2 for (ix, iy) in zip(x, y)])]])

    try:
        T = np.linalg.inv(F).dot(G)
    except np.linalg.LinAlgError:
        return 0, 0, float("inf")

    cxe = float(T[0] / -2)
    cye = float(T[1] / -2)

    try:
        re = math.sqrt(cxe ** 2 + cye ** 2 - T[2])
    except np.linalg.LinAlgError:
        return cxe, cye, float("inf")
    return cxe, cye, re


def calc_curvature_circle_fitting(x, y, npo=1):
    """各点での曲率を求める
    Calc curvature
    x,y: x-y position list
    npo: the number of points using Calculation curvature
    ex) npo=1: using 3 point
        npo=2: using 5 point
        npo=3: using 7 point
    """
    cv = []
    n_data = len(x)

    for i in range(n_data):
        lind = i - npo
        hind = i + npo + 1

        if lind < 0:
            lind = 0
        if hind >= n_data:
            hind = n_data

        xs = x[lind:hind]
        ys = y[lind:hind]
        (cxe, cye, re) = circle_fitting(xs, ys)

        if len(xs) >= 3:
            # sign evaluation
            c_index = int((len(xs) - 1) / 2.0)
            sign = (xs[0] - xs[c_index]) * (ys[-1] - ys[c_index]) - (
                    ys[0] - ys[c_index]) * (xs[-1] - xs[c_index])

            # check straight line
            a = np.array([xs[0] - xs[c_index], ys[0] - ys[c_index]])
            b = np.array([xs[-1] - xs[c_index], ys[-1] - ys[c_index]])
            theta = math.degrees(math.acos(
                np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))))

            if theta == 180.0:
                cv.append(0.0)  # straight line
            elif sign > 0:
                cv.append(1.0 / -re)
            else:
                cv.append(1.0 / re)
        else:
            cv.append(0.0)

    return cv

def n_sampling(lst, n, end=False):
    """listを与えたとき,いい感じに均等になるようにn個サンプリングする
    
    Parameters:
        lst: list
            リスト
        n: int
            抽出する数
        end: bool
            始めと終わりの構造を必ず含める場合True
    """
    division = len(lst) / n
    a = [lst[round(division * i):round(division * (i + 1))][0] for i in range(n)]
    if end:
        original_idx = [i for i in range(len(lst))]
        idx = [original_idx[round(division * i):round(division * (i + 1))][0] for i in range(n)]
        if original_idx[-1] != idx[-1]:
            a[-1] = lst[-1]
    return a

def get_diff(atoms1,atoms2,mic=False):
    """atoms1とatoms2の距離の差分を算出する

    Parameters:
    
    atoms1: Atoms
        Atomsオブジェクト
    atoms2: Atoms
        Atomsオブジェクト
    mic: bool
        周期境界条件で最小移動規則を適用する場合.True

    Returns:
        float: 距離
    """
    pos1 = atoms1.get_positions()
    pos2 = atoms2.get_positions()
    d = pos2 - pos1
    if mic:
        d = find_mic(d, atoms2.get_cell(), atoms2.pbc)[0]
    d = np.linalg.norm(d)
    return d


def minimize_rotation_and_translation_for_specified_indices_only(target, atoms, indices=None):
    """atomsの原子の位置をtargetの位置と近くなるように(最小二乗法)配置する.
    
    ase.build.rotate.minimize_rotation_and_translationでindexを指定できるようにした関数
    indices=Noneの時はASEのminimize_rotation_and_translationと同じで全ての原子を動かす.
    """
    if indices is None:
        indices = [i for i in range(len(target))]
    p = atoms.get_positions()[indices]
    p0 = target.get_positions()[indices]

    # centeroids to origin
    c = np.mean(p, axis=0)
    p -= c
    c0 = np.mean(p0, axis=0)
    p0 -= c0

    # Compute rotation matrix
    R = rotation_matrix_from_points(p.T, p0.T)

    atoms.positions[indices] = np.dot(p, R.T) + c0
    
def get_bond(atoms, unique=True, mult=1.0, **kwargs):
    """結合している原子のindex番号のリストを返す
    
    Parameters:
    
    atoms: Atoms
        対象のAtoms
    unique: bool
        Trueの場合,A-B,B-Aのいずれかのみを返す
    mult: float
        大きい程,離れていても結合していると判定される.
    kwargs:
        | 元素毎で解離の共有結合半径を指定できる
        | ex) H=0.5 で水素の共有結合半径を0.5Åに変更できる.
    """
    cutoff = natural_cutoffs(atoms, mult=mult,**kwargs)
    nl = build_neighbor_list(atoms,cutoff)
    ana = Analysis(atoms, nl=nl)
    if unique:
        bonds = ana.unique_bonds[0]
    else:
        bonds = ana.all_bonds[0]
    return bonds


def connected_components(atoms,indices=None,mult=1.0,**kwargs):
    """atoms中に存在する分子をindex番号毎にまとめたジェネレーターを返す
    
    Exampleを参照
    
    Parameters:
    
    atoms: Atoms
        対象のAtoms
    indices: list of int
        特定の原子のみcomponentを調べる時に指定する.
        Noneの場合,全ての原子を対象とする.
    mult: float
        大きい程,離れていても結合していると判定される.
    kwargs:
        | 元素毎で解離の共有結合半径を指定できる
        | ex) H=0.5 で水素の共有結合半径を0.5Åに変更できる.
        
        
    Examples:
    
        例えばatomsの[5,6,7,8,9]がメタン分子,[10,11]がNO分子だった場合,
        分子ごとに分離することができる
        
        >>> gen = connected_components(atoms, indices=[5,6,7,8,9,10,11])
        >>> for i in gen:
        >>>     print(i)
        >>> #--> {5,6,7,8,9}
        >>> #--> {10,11}
    """
    if indices is None:
        indices = [i for i in range(len(atoms))]
    G = nx.Graph()
    bonds = get_bond(atoms,mult=mult,**kwargs)
    for i in indices:
        for j in bonds[i]:
            if j in indices:
                G.add_edge(i, j)
    return nx.connected_components(G)