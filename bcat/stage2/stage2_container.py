#import bcat.stage1.imaging

class stage2_container:

    def __init__(self, data):
        self.data = data
        pass

    def imaging(self):
        return bcat.stage2.imaging.imaging(self.data)
