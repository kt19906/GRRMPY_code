from grrmpy.optimize.attach import optimize_eq,automate_maxstep,write_traj
from grrmpy.optimize.optimizer import CustumOptimizer
import warnings

try:
    from matlantis_features.ase_ext.optimize import FIRELBFGS #今後matlantis_featuresのアップデートの際に場所が変更される恐れがあるため.
    __all__ = ["optimize_eq","automate_maxstep","write_traj","FIRELBFGS","CustumOptimizer"]
    
except:
    warnings.warn('matlantis_featuresのFIRELBFGSのディレクトリの位置が変更されたためインポートできませんでした')
    __all__ = ["optimize_eq","automate_maxstep","write_traj","CustumOptimizer"]