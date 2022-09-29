# bcat
1.85m望遠鏡の観測データを解析するためのソフトウェア

<span style='color; red'> 現在バグが見つかったため、調査中:disappointed: </span>

## 仕様概要
bcatで観測データを読み込み、chopper wheel、座標のregrid、baselineを順番に行う。
事前にfitsを作成するためのheaderを用意しておく必要がある。
使用例は、bcat_sxample.ipynbを参照。

## 使用例
主なパートの例を示す

- データの取得
    ```python
    d = bcat.io.opu1p85.opendata(path,vwidth=150, spec='12CO21')
    ```
    でデータを読み込む。  
    path : データがあるpath  
    vwidth : 速度範囲  
    spec : 指定する輝線  


- chopper wheel
    ```python
    Tas, coord_on = d.get_chopper_wheel_spec()
    ```
    上記の「データの取得」で得たデータでchopper wheelを行う。


- 速度補正やデータの整形
    ```python
    image = bcat.stage2.imaging.speed()
    ```
    インスタンス化して使用する。


- baselinefit
    ```python
    baseline = bcat.stage2.Baseline_Rms.Baseline_Rms(numpy_array型のデータ, 速度)
    ```
    baselinefitを行う関数。  
    https://docs.google.com/presentation/d/1YzC0L4Qq0MzBXo5WYCD1chwpRe4QzZ-UQyatYqaom5g/edit?usp=sharing
    にbaseline関数の詳細が書いている。

## 環境
メモリ：60GBのPCで作成を行なった。（目視では、メモリが30GBもあれば十分）

macでも使用可能 (メモリがある程度消費されため、使用する際はメモリの確保が必要)  
mac環境  
   - プロセッサ：1.2 GHz クアッドコアIntel Core i7
   - メモリ：16 GB 3733 MHz LPDDR4X