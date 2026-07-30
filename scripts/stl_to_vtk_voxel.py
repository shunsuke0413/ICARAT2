# 参考元: https://qiita.com/torisan_piyopiyo/items/588e607c578b4b00903a
# https://qiita.com/matsxxx/items/aff526ea2e1f09590b34
# https://qiita.com/torisan_piyopiyo/items/f4e54d2c7453598a3a39
# python3.Xのみ対応．"pip install vtk"，"pip install pyevtk"が必要
import vtk
import time
import numpy as np
import vtk.util.numpy_support as vnp

######## インプットデータ ########
# インプットするstlパス
filename_in = "jetEngineDesignDomainFine.stl"
# アウトプットするvtkパス
filename_out = "out.vtk"
# メッシュサイズ
mesh_size = 100
# 許容誤差(10^-7~10^-3程度，小さいほど精度が高い)
tol = 1e-7
# ボクセルの種類（立方体(="cubic")か直方体(="rect")）
cubicORrect = "rect"
##################################

# main routine start
start = time.time()

# stlファイルを読む
reader = vtk.vtkSTLReader()
reader.SetFileName(filename_in)
reader.Update()

closed_poly = reader.GetOutput()

# メッシュグリッドを作成する
# x_min:0 x_max:1, y_min:2,y_max:3,z_min:4,z_max:5
bounds = closed_poly.GetBounds()
max_size = max([bounds[1] - bounds[0], bounds[3] -
                bounds[2], bounds[5] - bounds[4]])
cell_dims = [mesh_size, mesh_size, mesh_size]  # x, y, z

if cubicORrect == "cubic":
    mesh_pitch = [max_size/cell_dims[0],
                  max_size/cell_dims[1],
                  max_size/cell_dims[2]]
else:
    mesh_pitch = [(bounds[1] - bounds[0])/cell_dims[0],
                  (bounds[3] - bounds[2])/cell_dims[1],
                  (bounds[5] - bounds[4])/cell_dims[2]]

mins = [bounds[0], bounds[2], bounds[4]]


px, py, pz = mesh_pitch
mx, my, mz = (cell_dims+np.array([1, 1, 1])) * mesh_pitch  # max
points = vtk.vtkPoints()
coords = np.stack(np.mgrid[:mx:px, :my:py, :mz:pz], -1).reshape(-1, 3) + mins
points.SetData(vnp.numpy_to_vtk(coords))

structured_base_mesh = vtk.vtkStructuredGrid()
structured_base_mesh.SetExtent(
    0, cell_dims[0], 0, cell_dims[1], 0, cell_dims[2])
structured_base_mesh.SetPoints(points)

# 構造格子データを非構造格子データに変換する
append = vtk.vtkAppendFilter()
append.AddInputData(structured_base_mesh)
append.Update()
base_mesh = append.GetOutput()

# Voxelの中心の座標を求める
cell_centers = vtk.vtkCellCenters()
cell_centers.SetInputData(base_mesh)
cell_centers.Update()

poly_points = cell_centers.GetOutput()

# vtkSelectEnclosedPointsで内外判定
select_enclosed = vtk.vtkSelectEnclosedPoints()
select_enclosed.SetInputData(poly_points)
select_enclosed.SetSurfaceData(closed_poly)
select_enclosed.SetTolerance(tol)
select_enclosed.Update()

# フィルタ適用、フィルタが持つ形状にぶら下がっている
# "SelectedPoints"という配列を参照してbaseMeshのCellDataに与える
isInsideOrOutside = select_enclosed.GetOutput(
).GetPointData().GetArray("SelectedPoints")
structured_base_mesh.GetCellData().AddArray(isInsideOrOutside)
# PointDataで配列を取得した後、CellDataにぶらさげるので注意

# 形状のCellDataにぶら下がっている"SelectedPoints"という配列から情報を抽出する
threshold = vtk.vtkThreshold()
threshold.SetInputArrayToProcess(
    0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, "SelectedPoints")
threshold.SetInputData(structured_base_mesh)
threshold.SetLowerThreshold(1)  # (min, max) = (1, 1)
threshold.SetUpperThreshold(1)
threshold.Update()

# VoxelデータをVTKファイルで出力する
writer = vtk.vtkDataSetWriter()
writer.SetFileName(filename_out)
writer.SetInputData(threshold.GetOutput())
writer.Update()

elapsed_time = time.time() - start
print("computation time:{0}".format(elapsed_time) + "[sec]")
