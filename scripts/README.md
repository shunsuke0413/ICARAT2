# スクリプト一覧

## delete_all_objects_and_vtk.sh

ICARAT/以下のオブジェクトファイル及び vtk ファイルを全て消します．
**実行する前に，スクリプトの中を見て自分が何をするのかをよく確認してください．**

```
vim delete_all_objects_and_vtk.sh
chmod 777 delete_all_objects_and_vtk.sh
./delete_all_objects_and_vtk.sh
```


## vtk_to_in.cpp

vtk file(ascii)から mpi 並列計算用の in ファイルを生成するスクリプト．
metis が vendor/にインストールされていることが前提(docs/how_to_install_libraries.md)．
適時 main 関数内を編集する必要あり．

```
~ファイル編集~
g++ -std=c++14 vtk_to_in.cpp -o vtk_to_in.out
./vtk_to_in.out
```

## stl_to_vtk_voxel.py

stl ファイルから vtk ファイルのボクセルデータを作ることができる．

```
pip3 install vtk pyevtk
~ファイル編集~
python3 stl_to_vtk_voxel.py
```
