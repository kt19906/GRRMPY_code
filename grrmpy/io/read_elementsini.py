import pandas as pd

def read_elementsini(elements_ini_path):
    """構造をViewerで可視化する際の原子の色を規定する辞書を作成する.
    
    Parameters:
    
    elements_ini_path: Union[str, Path]
        elements.ini(VESTAで使われるやつ)のパス

    Returns:
        dict: 元素名,hexカラーコードの辞書.{'H':'0xffffff'}
    """     
    # iniファイルの読み込み   
    df = pd.read_table(elements_ini_path,header=None,engine='python',delim_whitespace=True)
    element = df[1].tolist()
    r = df[5].tolist()
    g = df[6].tolist()
    b = df[7].tolist()
    color_dict = {e:"0x"+format(hex(int(255*R)).replace('0x', ''), '0>2')+ # RGBをHexに変換
                  format(hex(int(255*G)).replace('0x', ''), '0>2')+
                  format(hex(int(255*B)).replace('0x', ''), '0>2') 
                  for e,R,G,B in zip(element,r,g,b)} # {"H":"0xffffff"...}
    return color_dict
    # grrmpy.visualize.functions.generate_js_codeも参照するとよい
    
def read_csv(filename):
    """原子の色を規定したcsvファイルを読み込みcolor_dictを返す

    Args:
        filename (_type_): csvファイル(jmol_color.csv)
    Returns:
        dict: 元素名,hexカラーコードの辞書.{'H':'0xffffff'}
    """
    df = pd.read_csv(filename)
    element = df.iloc[:,0].tolist()
    r = df.iloc[:,1].tolist()
    g = df.iloc[:,2].tolist()
    b = df.iloc[:,3].tolist()
    color_dict = {e:"0x"+format(hex(int(R)).replace('0x', ''), '0>2')+ # RGBをHexに変換
                  format(hex(int(G)).replace('0x', ''), '0>2')+
                  format(hex(int(B)).replace('0x', ''), '0>2') 
                  for e,R,G,B in zip(element,r,g,b)} # {"H":"0xffffff"...}
    return color_dict