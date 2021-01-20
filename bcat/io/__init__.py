
import pickle

def load(filepath):
    return pickle.load(open(filepath, 'rb'))
