import warnings
from ase.optimize import FIRE
from ase.io import write,iread
from pathlib import PurePath

# USER
from grrmpy.optimize import FIRELBFGS
from grrmpy.optimize.attach import automate_maxstep
from grrmpy.io.read_listlog import log2atoms,read_connections
from grrmpy import pfp_calculator

class ReCalcGRRM():
    def __init__(self,
                 eq_listlog,
                 ts_listlog=None,
                 pt_listlig=None,
                 constraints=[],
                 calc_functions = pfp_calculator,
                 ):
        """

        Parameters:
        
        *data:
            指定の仕方は3通りある.
            
            方法1:
                eq_listとconnectionsを与える.
                >>> eq_list = log2atoms('XXX_EQ_list.log','XXX.com','POSCAR')
                >>> connections = read_connections('XXX_TS_list.log')
                >>> auto_neb = AutoNEB(eq_list,connections,16)
            方法2:
                data1にEQのAtomsのリスト,data2にconnections情報(List[List[int]])を与える.
                connections情報はgrrmpy.io.read_listlog.read_connectionsで取得可能.
            方法3:
                data1にiniの構造のAtomsのリスト,data2にfinの構造のAtomsのリストを与える.
                つまり,data1の0番目の要素とdata2の0番目の要素でNEB計算が行なわれる.
            方法4:
                imagesを与える場合
        
        """
        if isinstance(data1,str) or isinstance(data1,PurePath):
            # data1がEQ_list.logだった場合.
            eq_list = log2atoms(data1,comfile,poscar,constraints)
            if isinstance(data2,str) or isinstance(data2,PurePath):
                # data2 がTS_list.logの場合.
                connections = read_connections(data2)
            else:
                connections = data2
        else:
            # data1がAotomsリストだった場合
            if type(data2) == list:
                # data2がAtomsリストだった場合
                ini_list = data1
                fin_list = data2
            else:
                # data1がeq_listのAtomsリスト,data2がconnectionsの場合.
                eq_list = data1
                connections = data2

    @property 
    def dafault_param(self):
        """run()のデフォルトのparameterを返す.
        CALC_MODE:
            "opt_only": bool
                最適化計算のみを行なう場合True.
        OPT:
            "optimizer": Optimizer
                構造最適化計算時に使用するOptimizer. Noneの場合構造最適化計算を行なわない.
            "fmax": float
                収束条件
            "maxstep_list": list of float
                maxstepsをリスト化したもの."steps_list"と要素数を合わせる必要がある.
                Noneの要素を指定した場合,automate_maxstepに合わせて自動的にmaxstepsを変化させる.
                optimizerがFIRELBFGSの場合は要素をNoneすることはできない.
            "automate_maxstep": dict
                {10:0.1, 5:0.2, 2:0.3, 0:0.35}の場合,fmaxが10以上の時maxstep=0.1にするという意味.
                keyが0の要素を必ず含める必要があるので注意.
            "steps_list":
                stepsのリスト.
            "logfile_folder":
                logファイルを保存するフォルダ名.Noneの場合,logファイルを作成しない.
            "trajectory_folder":
                trajctoryファイルを保存するフォルダ名.Noneの場合作成しない.
            "run_keyword":dict
                その他,最適化計算のrunの際のキーワード引数を指定できる.
            "attach_list":list of dict
                [{"function":lambda:XXX, "interval":XXX},]の形式で与える.
        IRC: 
            基本的に"OPT"と同じ.
            
            "optimizer": Optimizer
                構造最適化計算時に使用するOptimizer.
            "fmax": float
            "maxstep_list": list of float
            "automate_maxstep": dict
            "steps_list":
            "logfile_folder":
            "trajectory_folder":
        """
        return {
            "CALC_MODE":{
                "opt_only":False,
            },
            "OPT":{
                "optimizer":LBFGS,
                "fmax":0.001,
                "maxstep_list":[0.05,None],
                "automate_maxstep":{10:0.1, 5:0.2, 2:0.3, 0:0.35},
                "steps_list":[200,10000],
                "logfile_folder":"log",
                "trajectory_folder":None,
                "save_folder":"OPT",
                "run_keyword":{},
                "attach_list":[],
            },
            "NEB":{
                "nimages_list":[20,20],
                "parallel":True,
                "optimizer":FIRE,
                "run_keyword":{},
                "attach_list":[],
            },
            "SNEB":{
                
            },
            "IRC":{
                "optimizer":FIRELBFGS,
                "fmax":0.001,
                "maxstep_list":[0.05,0.2],
                "automate_maxstep":{10:0.1, 5:0.2, 2:0.3, 0:0.35},
                "steps_list":[200,10000],
                "logfile_folder":"log",
                "trajectory_folder":None,
                "save_folder":"IRC",
                "run_keyword":{},
                "attach_list":[],
            },
        }
        
    def set_param(self,param):
        self.opt = OPT(**param["OPT"])
        self.irc = IRC(**param["IRC"])

    def run(self,param):
        self.set_param(param)

class OPT():
    def __init__(self,
                 optimizer,
                 fmax,
                 maxstep_list,
                 automate_maxstep,
                 steps_list,
                 logfile_folder = None,
                 trajectory_folder=None,
                 save_folder = "OPT",
                 run_keyword = {},
                 attach_list = [],
                 calc_func = pfp_calculator,
                 ):
        self.optimizer = optimizer
        self.fmax = fmax
        self.maxstep_list = maxstep_list
        self.automate_maxstep = automate_maxstep
        self.steps_list = steps_list
        self.logfile_folder = logfile_folder
        self.trajectory_folder = trajectory_folder
        self.save_folder = save_folder
        self.run_keyword = run_keyword
        self.attach_list = attach_list
        self.calc_func = calc_func

    def _run(self,atoms,filename):
        for maxstep,steps in zip(self.maxstep_list,self.steps_list):
            atoms.calc = self.calc_func()
            # logfile引数
            if self.logfile_folder:
                logfile = f"{self.logfile_folder}/{filename}.log"
            else:
                logfile = None
            # trajectory引数
            if self.trajectory_folder:
                trajectory = f"{self.trajectory_folder}/{filename}.traj"
            else:
                trajectory = None
            # その他引数
            if self.optimizer == FIRELBFGS:
                kwargs = {"maxstep_fire":maxstep if maxstep else 0.2,
                          "maxstep_lbfgs":maxstep if maxstep else 0.2,
                          "logfile":logfile,
                          "trajectory":trajectory}
            else:
                kwargs = {"maxstep":maxstep if maxstep else 0.2,
                          "logfile":logfile,
                          "trajectory":trajectory}
            kwargs.update(self.run_keyword)
            
            opt = self.optimizer(atoms,**kwargs)
            # attach
            if maxstep is None:
                opt.attach(lambda:automate_maxstep(opt,self.automate_maxstep))
            for attach in self.attach_list:
                opt.attach(**attach)
            # run
            opt.run(fmax=self.fmax,steps=steps)
        return True if opt.converged else False
    
    def run(self,atoms,filename):
        converged = self._run(atoms,filename)
        if converged:
            write(f"{self.save_folder}/{filename}.traj",atoms)

class IRC(OPT):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
    
    def irun(self,atoms,filename):
        return super().irn(atoms,filename)
    
    def _run(self,vib_images,filename1,filename2):
        ini = vib_images[14].copy() # 正確に何番か！！！！！
        fin = vib_images[16].copy() # 正確に何番か！！！！！
        ini.calc = self.calc_func()
        fin.calc = self.calc_func()
        converged1 = self.irun(ini,filename1)
        converged2 = self.irun(fin,filename1)
        return converged1,converged2
    
