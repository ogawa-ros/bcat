
import bcat.stage1.chopper_wheel


class stage1_container(object):

    def __init__(self, data):
        self.data = data
        self.cw = bcat.stage1.chopper_wheel.stage1_chopper_wheel(self.data)
        pass

    def get_chopper_wheel_spec(self):
        return self.cw.chopper_wheel()
