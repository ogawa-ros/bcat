import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import tqdm

import astropy
from astropy.units import km, s, deg, m
from astropy import units as u
import astropy.time
import astropy.io.fits as fits
from astropy.wcs import WCS

class speed:
#     def __init__(self):  
#         self.base_data = d2
#         self.header = 
    def smooth(self, spectrum, spectral_axis, rest_freq, target_velocity):
        """
        スペクトルを窓巻数でconvolutionして、なます

        -----------------------------------------------
        引数：
        spectrum：on点のスペクトル
        spectral_axis：周波数（np.interpを使うためは、これは必ず降順でばければならない）
        rest_freq：観測した帯域（物質）の周波数
        target_velocity：regridしたい速度情報
        -----------------------------------------------
        """
        #横軸の間隔を取得
        df = abs(spectral_axis[1] - spectral_axis[0])
        #targetの周波数を計算
        target_freq = target_velocity * rest_freq/ astropy.constants.c + rest_freq
        target_df = target_freq[1] - target_freq[0]
        #targetの周波数の間隔を取得
        #窓巻数でconvolution
        win = signal.windows.boxcar(abs(int(target_df/df + 0.5)))
        win = win/sum(win)
        sp_smoothed = signal.convolve(spectrum, win, mode='same')
        
        return sp_smoothed*u.K
    
    
    def thin(self, flux, spectral_velocity, target_velocity):
        """
        速度のregrid
        target_velocityがないところは、numpyのinterpで線形補間する

        -----------------------------------------------
        引数：
            flux：スペクトルに単位（u*K）つけたもの
            spectral_velocity：周波数を速度に変換して、ドップラーシフトの補正をしたもの
            target_velocity：regridしたい速度
        -----------------------------------------------
        """
    #     print(flux, spectral_velocity, target_velocity)
        sp_thinned = np.interp(
            target_velocity.to('m/s').value, 
            spectral_velocity.to('m/s').value,
            flux.value
        )

        return sp_thinned*flux.unit
    
    def gauss(self, x, kernel_sigma):
        return np.exp(-x**2/(2*(kernel_sigma)**2))
    
    def speed_imaging(self, base_data, target_velocity, v_corre, target_freq):
        """
        ここでは、速度方向のregridを行う。

        -----------------------------------------------
        引数：
            base_data：いわゆるd2(データを取得して、chopper_wheel、containerしたデータ)
            target_velocity：ターゲットの速度
            v_corre：ドップラーシフトの速度補正
            target_freq：観測した周波数帯（物質）、ここでは13CO
        -----------------------------------------------
        """

        thinned_list = []
        velocity_convention = "radio"
        for i in range(len(base_data.data.coord.icrs.ra)):
            spectrum = base_data.data.spectrum[i]
            #numpy interpを使うために、降順にする必要がある。
            rf = base_data.data.get_rf_axis()
            if rf[0] < rf[-1]:
                flux = spectrum[::-1]*u.K
                spectral_axis = rf[::-1]

            else:
                flux = spectrum*u.K
                spectral_axis = rf

            smoothed = self.smooth(flux, spectral_axis, target_freq, target_velocity)

            spectral_velocity = spectral_axis.to('km/s', equivalencies=astropy.units.doppler_radio(target_freq)) + v_corre[i]
            #v_correctionを足す
            #w51は速度50km/sぐらい
            thinned = self.thin(smoothed, spectral_velocity, target_velocity)
            thinned_list.append(thinned)

        return thinned_list


    def calc_v_correction(self, base_data,loc_lon, loc_lat, loc_height, frame='fk4'):
        """
        ドップラーシフトの計算
        -----------------------------------------------
        引数：
            base_data：d2 いろいろな観測値が格納されているもの
            loc_lon：観測地の経度
            loc_lat：観測地の緯度
            loc_height：観測地の高さ
            frame：frameはデフォルトでfk4にする
        -----------------------------------------------
        """
        # 太陽の運動方向を定義
        dir_sun = astropy.coordinates.SkyCoord(
            ra = 18 * 15 * deg,                   # R.A. = 18 h
            dec = 30 * deg,                       # Dec. = 30 deg
            frame = frame,                        
            equinox = astropy.time.Time('B1900'), # 1900 年分点を指定する
        ).galactic

        # 太陽の速さを定義
        v_sun = 20 * km / s

        # 太陽速度を、銀河面に直交した座標系に成分分解し、U, V, W を計算
        U = v_sun * np.cos(dir_sun.b) * np.cos(dir_sun.l)
        V = v_sun * np.cos(dir_sun.b) * np.sin(dir_sun.l)
        W = v_sun * np.sin(dir_sun.b)

        # CartesianDifferential 型にしておく
        v_bary = astropy.coordinates.CartesianDifferential(U, V, W)

        # 観測者の地球上の座標を設定します。例は長野県野辺山です。
        loc1p85 = astropy.coordinates.EarthLocation(
                lon = loc_lon * deg,
                lat = loc_lat * deg,
                height = loc_height* m,
        )

        # 観測時刻を定義します
        tobs = astropy.time.Time(base_data.data.coord.obstime.value, format='unix')

        # 観測方向を定義します。
        target = astropy.coordinates.SkyCoord(
                base_data.data.coord.icrs.ra,
                base_data.data.coord.icrs.dec,
                frame = 'icrs',
                obstime = tobs,     # 時刻を設定する
                location = loc1p85, # 観測者座標も設定する
            )
        # if mode=='radec'
        #     target = astropy.coordinates.SkyCoord(
        #         base_data.data.coord.icrs.ra,
        #         base_data.data.coord.icrs.dec,
        #         frame = 'icrs',
        #         obstime = tobs,     # 時刻を設定する
        #         location = loc1p85, # 観測者座標も設定する
        #     )
        # if mode == 'gal':


        # V_observer (LSR 系からみた観測者の速度) の計算
        v_observer = astropy.coordinates.SkyCoord(
            # 観測者位置を、観測時刻での GCRS 座標系に変換する
            loc1p85.get_gcrs(target.obstime)
        ).transform_to(
            # それを、LSR 座標系に変換する
            astropy.coordinates.LSR(v_bary=v_bary) # 太陽速度を指定する
                                                   # デフォルトでは、Schönrich et al. (2010) が使われてしまう
        ).velocity

        # V_obs のうち、天体方向の成分を取り出す
        v_correction = v_observer.d_x * np.cos(target.icrs.dec) * np.cos(target.icrs.ra) + \
                       v_observer.d_y * np.cos(target.icrs.dec) * np.sin(target.icrs.ra) + \
                       v_observer.d_z * np.sin(target.icrs.dec)

        # tobs, ra, dec には配列も渡せます

        return v_correction
    
    def space_regrid(self, on_x, on_y, regrid_x, regrid_y, thinned_spectrum, rms, header, kernel_sigma):
        """
        空間方向のregridをしている
        regirdしたい点の周りにon点が二つが2個しかないなら、nanを入れる
        -----------------------------------------------
        引数：
           on_x：on点のx座標
           on_y：on点のy座標
           regrid_x：regridしたいx座標
           regrid_y：regridしたいy座標
           thinned_spectrum：間引いたスペクトル
        -----------------------------------------------
        出力：
            list
        """
        result = []
        for i in range(len(regrid_x.ravel())):
            # x, y 座標値
            x = regrid_x.ravel()[i]
            y = regrid_y.ravel()[i]
            #実角
            r = 60*3/3600
            sx = x - r/np.cos(np.deg2rad(y))
            ex = x + r/np.cos(np.deg2rad(y))
            sy = y - r
            ey = y + r

            ind = np.where((sx*u.deg<on_x) & (ex*u.deg>on_x) & (sy*u.deg<on_y) & (ey*u.deg>on_y))[0]
            if len(ind) <= 2:
                result.append(np.ones(header['NAXIS3'])*np.nan)
            else:
                wcs_X = on_x[ind]
                wcs_Y = on_y[ind]
                rms_tempo = rms[ind]
                sp = thinned_spectrum[ind]
                dist= ((wcs_X - x*u.deg)**2 + (wcs_Y - y*u.deg)**2)**(1/2)
                weight = self.gauss(dist, kernel_sigma)/rms_tempo**2
                sp = np.array(sp)
                weight_sum = np.sum(weight)
                weight /= weight_sum
                imaging_array = sp * weight[:,None]
                result.append(np.sum(imaging_array, axis=0))

        return result

  