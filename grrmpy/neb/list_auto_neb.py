from pathlib import Path
import warnings

from grrmpy.neb.auto_neb import AutoNEB

class ListAutoNEB():
    def __init__(self,images_list,clac_func,logfolder="NEBLog"):
        self.logfolder = Path(logfolder)
        if self.logfolder.exists():
            warnings.warn("既に{logfolder}フォルダが存在します")