
import bcat.stage1.chopper_wheel


class stage1_container:
    
    def __init__(self, data):
        self.data = data
        pass

    def chopper_wheel(self):
        return bcat.stage1.chopper_wheel(self.data)

