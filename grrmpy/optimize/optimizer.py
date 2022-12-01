from functools import partial

def CustumOptimizer(optimizer,atoms,**kwargs):
    """クラス風の関数"""
    custun_optimizer = type(f"Custum{optimizer.__name__}", (optimizer,),{})
    custun_optimizer.stop = False # stopがTrueになると計算を止める
    custun_optimizer.run = run # runメソッドを上書き
    return custun_optimizer(atoms,**kwargs)

def run(self, fmax=0.05, steps=None):
    """self.stop=Trueになると計算を停止する.(ASEのoptimizerをほぼコピペ)"""  
    def irun(self):
        self.atoms.get_forces()
        yield False

        if self.nsteps == 0:
            self.log()
            self.call_observers()

        while not self.converged() and self.nsteps < self.max_steps and not self.stop: # この部分だけ変更

            self.step()
            self.nsteps += 1

            yield False

            self.log()
            self.call_observers()
            self.stop_func()

        yield self.converged()
        
    self.fmax = fmax
    if steps:
        self.max_steps = steps
        
    for converged in irun(self):
        pass
    return converged