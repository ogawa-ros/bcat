# bcat
1.85m望遠鏡の観測データを解析するためのソフトウェア

## 仕様概要
bcatで観測データを読み込み、chopper wheel、座標のregrid、baselineを順番に行う。
事前にfitsを作成するためのheaderを用意しておく必要がある。
使用例のnotebookは、〜〜にある。

## 使用例
それぞれのパートの例を示す

"""python
d = bcat.io.opu1p85.opendata(path,vwidth=150,spec='12CO21')
"""
でデータを読み込む。
path : データがあるpath
vwidth : 速度範囲
spec : 指定する輝線