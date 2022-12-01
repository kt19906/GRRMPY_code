import sys
import warnings
import traceback
from pathlib import Path
from typing import Any,Dict,List,Optional,Union
from io import StringIO
import threading
import time
import numpy as np
from ase import Atoms
from ase.constraints import FixAtoms

from ase.io import Trajectory, write
from ase.visualize.ngl import NGLDisplay
from ipywidgets import (Button, Checkbox, Output,
                        Text, BoundedFloatText,RadioButtons,Image,
                        HBox,VBox,Tab,Dropdown,Layout)
from ipyevents import Event  
from traitlets import Bunch

import nglview as nv
from nglview import NGLWidget
from nglview.component import ComponentViewer
from nglview.color import ColormakerRegistry

# USER
from grrmpy.visualize.functions import (update_tooltip_atoms,generate_js_code,
                                        get_struct,add_force_shape,rotate_view,spin_view)
from grrmpy.visualize.color import default,vesta,jmol

class Viewer(NGLDisplay):
    def __init__(
        self,
        atoms: Union[Atoms, Trajectory, List[Atoms]],
        xsize: int = 500,
        ysize: int = 500,
        ):
        super().__init__(atoms, xsize=xsize, ysize=ysize)
        self.v = self.gui.view  # For backward compatibility...
        # del self.gui # デフォルトのGUIを削除
        # self.gui = HBox([self.view, VBox()]) # GUIを再設定
        
        # Make useful shortcuts for the user of the class
        self.gui.view = self.view
        self.gui.control_box = self.gui.children[1]
        self.gui.custom_colors = self.custom_colors
        
        ####### Property #####################################
        self.replace_structure = False
        self._use_struct_cache = True
        self._struct_cache = []
        self._force_components = []
        self.force_color = [1, 0, 0]  # Red vector for force color.
        self.pre_label = False
        if isinstance(atoms, Atoms):
            self._struct_cache = [None]
        else:
            # atoms is Trajectory or List[Atoms]
            self._struct_cache = [None for _ in range(len(atoms))]
        #色の設定
        self.cm = ColormakerRegistry
        default_jscode = generate_js_code(default)
        vesta_jscode = generate_js_code(vesta)
        jmol_jscode = generate_js_code(jmol)
        self.cm.add_scheme_func('default',default_jscode)
        self.cm.add_scheme_func('vesta',vesta_jscode)
        self.cm.add_scheme_func('jmol',jmol_jscode)
        ###################################################

        # ---原子上にマウスを置いたときに,原子のindexと位置を表示する
        update_tooltip_atoms(self.view, self._get_current_atoms())
        
        # GUI作成&表示
        self.build_gui()
        
        # 初期表示
        self.view.camera = "orthographic" if self.camera_style=='平行投影' else "perspective"
        self.view.add_spacefill()
        self.view.add_ball_and_stick()
        self.view.add_label(
            color="black",
            labelType="text",
            labelText=["" for _ in range(len(self._get_current_atoms()))],
            zOffset=2.0,
            attachment="middle_center",
            radius=1,
        )
        self._update_repr()
        
        self.view.unobserve(NGLWidget._on_frame_changed)
        self.view.observe(self._on_frame_changed, names=["frame"])
        
    @property
    def camera_style(self):
        return self.gui.camera_radio_btn.value
    
    @property
    def show_force(self):
        return self.gui.show_force_checkbox.value
        
    @property
    def show_charge(self):
        return self.gui.show_charge_checkbox.value
    
    
    def build_gui(self):
        ###カメラ(RadioButton)###
        self.gui.camera_radio_btn = RadioButtons(
            options=['平行投影','透視投影'],
            value='平行投影',
            description='カメラ')
        self.gui.camera_radio_btn.observe(self.change_camera)
        ###ラベル(RadioButton)###
        self.gui.label_radio_btn = RadioButtons(
            options=['なし','インデックス','電荷','FixAtoms'],
            value='なし',
            description='ラベル')
        self.gui.label_radio_btn.observe(self.change_label)
        ###セル(CheckBox)####
        self.gui.cell_check_box = Checkbox(value=True,description="セルユニット",)
        self.gui.cell_check_box.observe(self.show_unitcell)
        ###カラースキーム###
        self.csel = Dropdown(options=["default","vesta","element","jmol"],
                             value='default', description='色')
        self.csel.observe(self._update_repr)
        ###モデル###
        self.gui.model_radio_btn = RadioButtons(
            options=['球棒モデル','空間充填モデル'],
            value='球棒モデル',
            description='モデル')
        self.gui.model_radio_btn.observe(self._update_repr)
        ###再配置(チェックボックス)###
        self.gui.replace_structure_checkbox = Checkbox(
            value=self.replace_structure,
            description="再配置")
        self.gui.replace_structure_checkbox.observe(self.change_replace_structure)
        ###アウトプット(エラー表示)###
        self.gui.out_widget = Output(layout={"border": "0px solid black"})
        ##IMG##
        self.gui.a = Button(description='A',tooltip='A軸方向から見る',layout = Layout(width='30px'))
        self.gui.a.on_click(lambda e,x=180,y=-90,z=90: self.rotate_view(e,x,y,z))
        self.gui.b = Button(description='B',tooltip='B軸方向から見る',layout = Layout(width='30px'))
        self.gui.b.on_click(lambda e,x=90,y=0,z=-90: self.rotate_view(e,x,y,z))
        self.gui.c = Button(description='C',tooltip='C軸方向から見る',layout = Layout(width='30px'))
        self.gui.c.on_click(lambda e,x=180,y=0,z=0: self.rotate_view(e,x,y,z))
        self.gui.up = Button(description='x',tooltip=' x軸周りで回転(上)',layout = Layout(width='30px'))
        self.gui.up.on_click(lambda e,axis=0, angle=10: self.spin_view(e,axis,angle))
        self.gui.down = Button(description='x*',tooltip='x軸周りで回転(下))',layout = Layout(width='30px'))
        self.gui.left = Button(description='y*',tooltip='y軸周りで回転(左)',layout = Layout(width='30px'))
        self.gui.right = Button(description='y',tooltip='y軸周りで回転(右)',layout = Layout(width='30px'))
        self.gui.cc = Button(description='z',tooltip='z軸周りで回転(時計周り)',layout = Layout(width='30px'))
        self.gui.rc = Button(description='z*',tooltip='z軸周りで回転(反時計周り)',layout = Layout(width='30px'))
        self.gui.step = BoundedFloatText(value=10,min=0,max=360,step=10,
                                         description='Step(*):',
                                         layout = Layout(width='150px'))
        ##ダウンロード##
        self.gui.download = Button(description='PNGをダウンロード',
                                   tooltip='ローカルPCにPNGをダウンロードする')
        ##保存##
        self.gui.save = Button(description='ファイルへ保存',
                               tooltip='.png,.htmlとその他ASEのwriteで指定できる拡張子で保存できる.')
        ###表示####
        # r = list(self.gui.control_box.children)
        img1 = HBox([
            self.gui.a,
            self.gui.b,
            self.gui.c,
            self.gui.step])
        img2 = HBox([
            self.gui.up,
            self.gui.down,
            self.gui.right,
            self.gui.left,
            self.gui.cc,
            self.gui.rc,
            ])
        general = VBox([
            img1,
            img2,
            HBox([self.gui.download,self.gui.save]),
            self.gui.label_radio_btn,])
        other = VBox([
            self.csel,
            self.rad,
            self.gui.model_radio_btn,
            self.gui.camera_radio_btn,
            self.gui.cell_check_box,
        ])
        detail = VBox([
            
        ])
        
        self.tab = Tab([general,other,detail],_titles={0:"基本", 1:"その他",2:"詳細"})
        self.gui.control_box.children = tuple([self.tab,
                                               self.gui.replace_structure_checkbox,
                                               self.gui.out_widget])
        
    def rotate_view(self,e,x,y,z):
        rotate_view(self.view,x=x,y=y,z=z)
        
    def spin_view(self,e, axis, angle):
        spin_view(self.view, axis, angle)
        
    def update_label(self,labelText):
        self.view.update_label(
            color="black",
            labelType="text",
            labelText=labelText,
            zOffset=2.0,
            attachment="middle_center",
            radius=1,
        )
        
    def _change_label(self,atoms,option):
        if option == "インデックス":
            labelText=[str(i) for i in range(len(atoms))]
        elif option == "電荷":
            try:
                labelText = np.round(atoms.get_charges().ravel(), 3).astype("str").tolist()
            except:
                labelText=["" for _ in range(len(atoms))]
                with self.gui.out_widget:
                    raise Exception("Calculatorを設定してください")    
        elif option == "FixAtoms":
            labelText = self._get_fix_atoms_label_text(atoms)
        self.update_label(labelText)
        
    def change_label(self,e=None):
        self.gui.out_widget.clear_output()
        option = self.gui.label_radio_btn.value
        atoms = self._get_current_atoms()
        if option == "なし":
            labelText=["" for _ in range(len(atoms))]
            self.update_label(labelText)
            return 
        else:
            self._change_label(atoms,option)
            return 
           
    def change_camera(self,e=None):
        option = self.gui.camera_radio_btn.value
        if option == "平行投影":
            self.view.camera = 'orthographic'
        elif option == "透視投影":
            self.view.camera = 'perspective'
            
    def _update_repr(self,e=None):
        option = self.gui.model_radio_btn.value
        if option == "球棒モデル":
            self.view.update_spacefill(radiusType='covalent',
                                    radiusScale=self.rad.value,
                                    color_scheme=self.csel.value)#color_scale='rainbow')
            self.view.update_ball_and_stick(color_scheme=self.csel.value)
        elif option == "空間充填モデル":
            self.view.remove_spacefill()
            self.view.add_spacefill()
            self.view.update_spacefill(radiusType="vwf",color_scheme=self.csel.value)
            
    def show_unitcell(self,e=None):
        if self.gui.cell_check_box.value: 
            self.view.add_unitcell() # Cellの表示
        else:
            self.view.remove_unitcell()
                
    def _get_current_atoms(self) -> Atoms:
        if isinstance(self.atoms, Atoms):
            return self.atoms
        else:
            return self.atoms[self.view.frame]
        
    def _get_fix_atoms_label_text(self,atoms):
        indices_list = []
        for constraint in atoms.constraints:
            if isinstance(constraint, FixAtoms):
                indices_list.extend(constraint.index.tolist())
        label_text = []
        for i in range(len(atoms)):
            if i in indices_list:
                label_text.append("FixAtoms")
            else:
                label_text.append("")
        return label_text
        
    def change_replace_structure(self,event: Optional[Bunch] = None):
        if self.gui.replace_structure_checkbox.value:
            self.replace_structure = True
            self._on_frame_changed(None)
        else:
            self.replace_structure = False

    def _on_frame_changed(self, change: Dict[str, Any]):
        """set and send coordinates at current frame"""
        v: NGLWidget = self.view
        atoms: Atoms = self._get_current_atoms()

        if self.replace_structure:
            # set and send coordinates at current frame
            struct = self._struct_cache[v.frame]
            if struct is None:
                struct = get_struct(atoms)
                if self._use_struct_cache:
                    self._struct_cache[v.frame] = struct  # Cache
            v._remote_call("replaceStructure", target="Widget", args=struct)
        else:
            # Only update position info
            v._set_coordinates(v.frame)

        # Tooltip: update `var atoms_pos` inside javascript.
        atoms = self._get_current_atoms()
        if atoms.get_pbc().any():
            _, Q = atoms.cell.standard_form()
        else:
            Q = np.eye(3)
        Q_str = str(Q.tolist())
        var_str = f"this._Q = {Q_str}"
        v._execute_js_code(var_str)
        
        ###ラベルの更新
        option = self.gui.label_radio_btn.value
        if option != "なし":
            atoms = self._get_current_atoms()
            self._change_label(atoms,option)
        
    def _ipython_display_(self, **kwargs):
        """viewプロパティを書かなくてもjupyter上で勝手に表示してくれる"""
        return self.gui._ipython_display_(**kwargs)
        
