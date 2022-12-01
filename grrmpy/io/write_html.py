from ase.neb import NEB
from ase import Atoms
import ase.units as units
from ase.vibrations import Vibrations
from pathlib import Path

#User
from grrmpy.calculator import pfp_calculator
from grrmpy.vibrations.functions import to_html_table_and_imode,to_html_graph
from grrmpy.neb.functions import to_html_nebgraph
from grrmpy.functions import to_html_energy_diagram

def write_html(html,obj,*args,**kwargs):
    """htmlファイルとして保存する
    
    Parameters:
    
    html: string
        保存ファイル名
    obj: object
        | objの引数には次のオブジェクトを指定できる
        | - NEB
        |     NEBのエネルギーダイアグラムを作成
        | - Vibrations
        |     エネルギーダイアグラムと振動数表を作成
        | - trajファイル
        |     エネルギーダイアグラムを作成
    """
    if type(obj) == NEB:
        write_neb_graph(html,obj, calc_func=pfp_calculator,**kwargs)
    elif type(obj) == Vibrations:
        write_vibtb_and_vibgraph(html,obj, calc_func=pfp_calculator,**kwargs)
    elif type(obj) == Path or type(obj) == str:
        write_vib_graph(html,obj,calc_func=pfp_calculator,**kwargs)
    elif type(obj) == list:
        if type(obj[0]) == Atoms:
            write_energy_diagram(html,obj,calc_func=pfp_calculator,**kwargs)
    
def write_vib_table(html:str,vib_obj,full_html:bool=True):
    """振動数の表をhtmlのstrを出力する&虚振動の振動モード番号を出力する

    Parameters:
    
    html: string or None
        | 出力するhtmlファイル名.
        | Noneの場合,htmlテキストを戻り値として出力する
    vib_obj: Viblations object
        Viblationsオブジェクト
    full_html: boolean
        | <html>タグから始まる,完全なhtmlを出力する場合True
        | Falseの場合<div>タグから始まるテキストを出力
    """
    html_text,_ = to_html_table_and_imode(vib_obj,full_html=full_html,include_plotlyjs="cdn")
    if not html:
        return html_text        
    with open(html,"w") as f:
        f.write(html_text)
        
def write_vib_graph(html,vib_obj, n:int, outfile=None, calc_func=pfp_calculator, kT=units.kB * 300, nimages=30,full_html=True):
    """エネルギーダイアグラムを作成する

    Parameters:
    
    html: string or None
        | 出力するhtmlファイル名.
        | Noneの場合,htmlテキストを戻り値として出力する
    vib_obj: Viblations object
        Viblationsオブジェクト
    n: integer
        振動モード番号
    outfile: str
        | 振動の構造をtrajで出力する場合,trajファイル名.
        | Noneの場合出力しない
        | (write_mode()でtrajファイルを出力する)
    calc_func: function object
        Claculatorを返す関数
    nimages: integer
        | イメージの数, デフォルトは30
        | 偶数で設定する事を推奨する
    full_html:
        | <html>タグから始まる,完全なhtmlを出力する場合True
        | Falseの場合<div>タグから始まるテキストを出力
    """
    html_text = to_html_graph(vib_obj,n,outfile,calc_func,kT,nimages,full_html,include_plotlyjs="cdn")
    if not html:
        return html_text  
    with open(html,"w") as f:
        f.write(html_text)
        

def write_vibtb_and_vibgraph(html,vib_obj, outfile=None, calc_func=pfp_calculator, kT=units.kB * 300, nimages=30):
    """振動数の表とエネルギーダイアグラムのhtmlテキストを出力する.

    Parameters:
    
    html: string or None
        | 出力するhtmlファイル名.
        | Noneの場合,htmlテキストを戻り値として出力する
    vib_obj: Viblations object
        Viblationsオブジェクト
    outfile: str
        | 振動の構造をtrajで出力する場合,trajファイル名.
        | Noneの場合出力しない
        | (write_mode()でtrajファイルを出力する)
    calc_func: function object
        Claculatorを返す関数
    nimages: integer
        | イメージの数, デフォルトは30
        | 偶数で設定する事を推奨する
        
    Note:
        | 虚振動がないor複数ある場合は,エネルギーダイアグラムは出力しない
    """
    tb_txt,imode = to_html_table_and_imode(vib_obj,full_html=False,include_plotlyjs="cdn")
    if imode:
        fig_txt = write_vib_graph(None,vib_obj,imode,outfile,calc_func,kT,nimages,full_html=False)
    else:
        fig_txt = ""
    if not html:
        return tb_txt+fig_txt
    with open(html) as f:
        f.write(tb_txt+fig_txt)
        
def write_neb_graph(html,neb_obj,calc_func=pfp_calculator,full_html=True, unit="kJ/mol",**kwargs):
    html_txt = to_html_nebgraph(neb_obj,calc_func,full_html,unit,full_html,include_plotlyjs="cdn",**kwargs)
    if not html:
        return html_txt
    with open(html,"w") as f:
        f.write(html_txt)
        
def write_energy_diagram(html,images,calc_func=pfp_calculator,full_html=True, unit="kJ/mol",**kwargs):
    html_txt = to_html_energy_diagram(images,calc_func,full_html,unit,full_html,include_plotlyjs="cdn",**kwargs)
    if not html:
        return html_txt
    with open(html,"w") as f:
        f.write(html_txt)
        
