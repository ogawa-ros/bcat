import numpy
import pickle
import pathlib


def _write_header(f, dtype, header):
    h = {}
    h.update(header)
    h['dtype'] = dtype
    ph = pickle.dumps(h)
    ph_size = len(ph)
    
    buf = numpy.array(ph_size, dtype='uint32').tobytes()
    buf += ph
    f.write(buf)
    return


def save(path, array, header={}, mode='w'):
    path = pathlib.Path(path)
    
    if (not path.exists()) | (mode=='w'):
        with open(path, 'wb') as f:
            _write_header(f, array.dtype, header)
            return array.tofile(f)
    else:
        with open(path, 'ab') as f:
            return array.tofile(f)
    pass


def _read_header(path):
    headersize_size = 4
    
    with open(path, 'rb') as f:
        header_size = numpy.frombuffer(f.read(headersize_size),
                                       dtype='uint32')[0]
        header = pickle.loads(f.read(header_size))
        pass
    
    start_addr = headersize_size + header_size
    return header, start_addr


def load(path, mode='r+'):
    header, start_addr = _read_header(path)
    ret = {}
    ret.update(header)
    ret['array'] = numpy.memmap(
        path,
        dtype = header['dtype'],
        mode = mode,
        offset = start_addr,
    )
    return ret
