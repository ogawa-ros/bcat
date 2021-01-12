import time
import numpy
import pandas
import necstdb
import astropy
import astroquery.lamda
from astropy.units import GHz
from astropy.units import kHz
from astropy.units import km
from astropy.units import s
from astropy.constants import c

import bcat.structure
import bcat.stage1

class opu1p85(object):

    def __init__(self):
        self.loc1p85 = astropy.coordinates.EarthLocation(
                lon = 138.472153 * astropy.units.deg,
                lat = 35.940874 * astropy.units.deg,
                height = 1386 * astropy.units.m,
                )

        self.param = {
                'az': '1p85m-az',
                'el': '1p85m-el',
                'az_cmd':'1p85m-az_cmd',
                'el_cmd':'1p85m-el_cmd',
                }

        self.boards ={
            'xffts_board01':'board1',
            'xffts_board02':'board2',
            'xffts_board03':'board3',
            'xffts_board04':'board4'
                }

        self.line_board = {'12CO21':'xffts_board02','13CO21':'xffts_board01','C18O21':'xffts_board01',
                      '12CO32':'xffts_board04','13CO32':'xffts_board03','C18O32':'xffts_board03'
                      }


    def open(self,path):
        db_name = path
        db = necstdb.opendb(db_name)
        return db

    def necstdb2pandas(self,db):
        resample_list = []
        _df_resample,table = {},{}
        #いつもの
        for i, j in self.param.items():
            table[i] = db.open_table(j).read(astype='pandas')
            table[i]['timestamp'] = pandas.to_datetime(table[i]['timestamp'], unit='s')
            table[i] = table[i].set_index('timestamp').sort_index().rename(columns={'data': i}).astype(float)
        _df_resample = pandas.concat([table[i] for i in table], axis = 1)
        return _df_resample

    def get_restfreq(self,mol, trans):
        col, trans_, elev = astroquery.lamda.Lamda.query(mol=mol)
        return trans_[trans_['Upper']==trans+1]['Frequency'][0] * GHz

    def get_status(self,db):
        status = db.open_table('otf-status').read(astype='pandas')
        status['timestamp'] = pandas.to_datetime(status['timestamp'], unit='s')
        status = status.set_index('timestamp').sort_index()
        return status

    def read_wcs(self,db,df_resample):
        wcs = db.open_table('necst-telescope-coordinate-wcs').read(astype='array')
        data = wcs['data'][:, 1:3]
        time = pandas.to_datetime(wcs['timestamp'], unit='s')
        wcs_divide = pandas.DataFrame(data=data, columns=['wcs_x','wcs_y'], index=time).sort_index()
        df_resample = pandas.concat([df_resample,wcs_divide], axis=1)
        return df_resample

    def get_linedata(self,db,spec):
        f_spec = numpy.linspace(0,2,2**15)*GHz
        flo_b7 = 325.08*GHz
        flo_b6 = 225*GHz
        slo_lsb = 6.5*GHz
        slo_usb = 4.0*GHz
        slo_4_6 = 4.0*GHz
        slo_19_21 = 19.2*GHz
        rf1 = (f_spec+flo_b6-slo_lsb)
        rf2 = (f_spec+flo_b6+slo_usb)
        rf3 = (f_spec+flo_b7+slo_4_6)
        rf4 = (f_spec+flo_b7+slo_19_21)

        rf12co21 = self.get_restfreq('co',2)
        rf13co21 = self.get_restfreq('13co',2)
        rfc18o21 = self.get_restfreq('c18o',2)
        rf12co32 = self.get_restfreq('co',3)
        rf13co32 = self.get_restfreq('13co',3)
        rfc18o32 = self.get_restfreq('c18o',3)
        df = (2*GHz/2**15)
        dv_12 = df/rf12co21*c
        dv_13 = df/rf13co21*c
        dv_18 = df/rfc18o21*c
        vwidth = 150*km/s
        fwidth_12CO21 = (vwidth*rf12co21/c).to('GHz')
        fwidth_13CO21 = (vwidth*rf13co21/c).to('GHz')
        fwidth_C18O21 = (vwidth*rfc18o21/c).to('GHz')
        fwidth_12CO32 = (vwidth*rf12co32/c).to('GHz')
        fwidth_13CO32 = (vwidth*rf13co32/c).to('GHz')
        fwidth_C18O32 = (vwidth*rfc18o32/c).to('GHz')
        mask_12CO21 = (rf2>rf12co21-fwidth_12CO21)&(rf2<rf12co21+fwidth_12CO21)
        mask_13CO21 = (rf1>rf13co21-fwidth_13CO21)&(rf1<rf13co21+fwidth_13CO21)
        mask_C18O21 = (rf1>rfc18o21-fwidth_C18O21)&(rf1<rfc18o21+fwidth_C18O21)
        mask_12CO32 = (rf4>rf12co32-fwidth_12CO32)&(rf4<rf12co32+fwidth_12CO32)
        mask_13CO32 = (rf3>rf13co32-fwidth_13CO32)&(rf3<rf13co32+fwidth_13CO32)
        mask_C18O32 = (rf3>rfc18o32-fwidth_C18O32)&(rf3<rfc18o32+fwidth_C18O32)

        mask = {'12CO21':mask_12CO21,'13CO21':mask_13CO21,'C18O21':mask_C18O21,
                '12CO32':mask_12CO32,'13CO32':mask_13CO32,'C18O32':mask_C18O32,
                }
        rf = {'12CO21':rf2,'13CO21':rf1,'C18O21':rf1,
                '12CO32':rf4,'13CO32':rf3,'C18O32':rf3,
                }

        xffts= db.open_table(self.line_board[spec]).read(astype='array')
        line_data = pandas.DataFrame(data=(xffts['data'][:,:-1])[:,mask[spec]], index=pandas.to_datetime(xffts['data'][:,-1], unit='s')) # 分光データの最後の timestamp
        f = bcat.structure.freq_axis(2/32768*astropy.units.GHz, rf[spec][mask[spec]][0])
        del xffts
        return line_data,f

    def create_spec(self,db,spec):
        _df_resample = self.necstdb2pandas(db)
        df_resample = self.read_wcs(db,_df_resample)
        line_data,freq = self.get_linedata(db,spec)
        df_spec_2 = pandas.concat([df_resample["wcs_x"][::1],df_resample["wcs_y"][::1], line_data], axis=1)
        df_spec_1 = pandas.concat([df_spec_2["wcs_x"][::1],df_spec_2["wcs_y"][::1]], axis=1).interpolate()
        del df_spec_2
        df_spec = pandas.concat([df_spec_1, line_data], axis=1).dropna(how="any")
        del line_data
        return df_spec,freq

    def create_obsmode(self,df_resample,db):
        status = self.get_status(db)
        hot_s = status[status['data'] == b'hot start\x00'].index
        hot_e = status[status['data'] == b'hot end  \x00'].index
        off_s = status[status['data'] == b'off start\x00'].index
        off_e = status[status['data'] == b'off end  \x00'].index
        on_s = status[status['data'] == b'on start \x00'].index
        on_e = status[status['data'] == b'on finish\x00'].index
        scan_s = status[status['data'] == b'otf line '].index
        obsmode = []
        for d in df_resample.index:
            d_flag = None
            for i in range(len(on_s)):
                if (d > on_s[i]) & (d < on_e[i]):
                    d_flag = 'ON'
                    break
                else:
                    continue
            for i in range(len(off_s)):
                if (d > off_s[i]) & (d < off_e[i]):
                    d_flag = 'OFF'
                    break
                else:
                    continue

            for i in range(len(hot_s)):
                if (d > hot_s[i]) & (d < hot_e[i]):
                    d_flag = 'HOT'
                    break
                else:
                    continue
            if d_flag:
                obsmode.append(d_flag)
            else:
                obsmode.append('TRANS')
        return obsmode

    def create_coord(self,df_spec,frame='fk5'):
        coord = astropy.coordinates.SkyCoord(
            ra = numpy.array(df_spec['wcs_x']) * astropy.units.deg,
            dec = numpy.array(df_sepc['wcs_y']) * astropy.units.deg,
            frame = frame,
            location = self.loc1p85,
            obstime = astropy.time.Time(
                        numpy.array(df_spec.index.astype(numpy.int64)/10**9),
                        format='unix'),
            )
        return coord

    def container(self,label,path,spec='12CO21'):
        db = self.open(path)
        df_spec,freq = self.create_spec(db,spec)
        obsmode = self.create_obsmode(df_spec)
        coord = self.create_coord(df_spec)
        d1_data = bcat.structure.stage1_data(
            label=label,
            obsmode=obsmode,
            coord=coord,
            spectrum=df_spec.drop(['wcs_x', 'wcs_y'], axis=1).values,
            rf=freq)
        d1 = bcat.stage1.cnotainer(d1_data)
        return d1
