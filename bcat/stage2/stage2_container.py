
import pickle
import pathlib

#import bcat.stage1.imaging


class stage2_container:

    def __init__(self, data):
        self.data = data
        pass

    def save(self, basedir='.'):
        filename = self.data.label + '.stage2.pickle'
        saveto = pathlib.Path(basedir) / filename
        pickle.dump(self, open(saveto, 'wb'))
        return
        
    def imaging(self):
        return bcat.stage2.imaging.imaging(self.data)
