class ReactPaths():
    """複数のReactPathをマージしてまとめたクラス
    """
    def __init__(self,*react_path,title=None):
        """
        
        Parameters:
        
        react_path: ReactPath object
            ReactPathオブジェクト. 可変長引数
        title:
            | タイトルを付ける事ができる.デフォルトはNone
            | キーワード引数なので注意(位置引数でない)
        """
        self.title = title
        self.paths = list(react_path)
        
    def get_solid_df(self):
        pass
    
    def get_dot_df(self):
        pass