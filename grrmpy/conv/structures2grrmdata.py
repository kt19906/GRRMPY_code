from grrmpy.structure import EQList, PTList, TSList
from grrmpy.grrmdata import GrrmData

def structures2grrmdata(eqlist:EQList=None,tslist:TSList=None,ptlist:PTList=None):
    """EQList,TSList,PTListからGrrmDataオブジェクトを作成する
    
    引数にはNoneを指定することも可能

    Parameters:
    
    eqlist: EQList
        EQListオブジェクト
    tslist: TSList
        TSListオブジェクト
    ptlist: PTList
        PTListオブジェクト
        
    Returns:
    
    GrrmData: GrrmDataオブジェクト
    """
    grrmdata = GrrmData()
    if type(eqlist) == EQList:
        grrmdata.eq = eqlist.copy()
        com = eqlist.com.copy()
    else:
        raise TypeError("eqlistにはEQlistオブジェクトを指定します")
    if type(tslist) == TSList:
        grrmdata.ts = tslist.copy()
        com = tslist.com.copy()
    else:
        raise TypeError("tslistにはTSlistオブジェクトを指定します")
    if type(ptlist) == PTList:
        grrmdata.pt = ptlist.copy()
        com = ptlist.com.copy()
    else:
        raise TypeError("ptlistにはPTlistオブジェクトを指定します")
    grrmdata.com = com
    return grrmdata