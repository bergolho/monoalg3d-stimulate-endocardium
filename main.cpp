#include <iostream>
#include <string>
#include <map>

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkSmartPointer.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkHexahedron.h>
#include <vtkSphereSource.h>
#include <vtkAppendPolyData.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include <vtkLine.h>
#include <vtkCellLocator.h>
#include <vtkPointLocator.h>
#include <vtkCellCenters.h>

using namespace std;

struct point3d {
    uint32_t id;
    double pos[3];
};

void read_alg_mesh (const char filename[], vector<struct point3d> &mesh_points_alg) {
    printf("[!] Reading ALG mesh ...\n");
    FILE *file = fopen(filename, "r");
    uint32_t k = 0;
    double aux[17];
    while (fscanf(file,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",\
            &aux[0],&aux[1],&aux[2],&aux[3],&aux[4],&aux[5],&aux[6],&aux[7],&aux[8],&aux[9],&aux[10],&aux[11],\
            &aux[12],&aux[13],&aux[14],&aux[15],&aux[16]) != EOF) {
        struct point3d p;
        p.id = k;
        p.pos[0] = aux[0];
        p.pos[1] = aux[1];
        p.pos[2] = aux[2];
        k++;

        mesh_points_alg.push_back(p);
    }
    fclose(file);
}

int main (int argc, char *argv[])
{
    if (argc-1 != 5) {
      printf("-----------------------------------------------------------------------------------------------------------------------------------------\n");
      printf("Usage:> %s <mesh.alg> <mesh_endocardium_LV.vtk> <mesh_endocardium_RV.vtk> <output_endocardium_LVRV.pts> <mode>\n", argv[0]);
      printf("-----------------------------------------------------------------------------------------------------------------------------------------\n");
      printf("<mesh.alg> = Mesh file in ALG format\n");
      printf("<mesh_endocardium_LV.vtk> = LV endocardium surface in VTK format\n");
      printf("<mesh_endocardium_RV.vtk> = RV endocardium surface in VTK format\n");
      printf("<output_endocardium_LVRV.pts> = Output stimulation point in PTS format\n");
      printf("<mode> = 1-> LVRV, 2-> LV, 3-> RV\n");
      printf("-----------------------------------------------------------------------------------------------------------------------------------------\n");
      exit(EXIT_FAILURE);
    }

    string filename_mesh_alg = argv[1];
    string filename_mesh_endocardium_LV = argv[2];
    string filename_mesh_endocardium_RV = argv[3];
    string filename_output_endocardium_LVRV = argv[4];
    int mode = atoi(argv[5]);

    // ------------------------------------------------------------------------------------------------------------------------------
    // 1) Read the mesh
    vector<struct point3d> mesh_points_alg;
    read_alg_mesh(filename_mesh_alg.c_str(), mesh_points_alg);

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (uint32_t i = 0; i < mesh_points_alg.size(); i++) {
        points->InsertNextPoint(mesh_points_alg[i].pos[0],mesh_points_alg[i].pos[1],mesh_points_alg[i].pos[2]);
    }
    vtkSmartPointer<vtkUnstructuredGrid> unstructured_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructured_grid->SetPoints(points);

    // PointLocator (Kd-tree for fast coordinates searching)
    vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
    pointLocator->SetDataSet(unstructured_grid);
    pointLocator->BuildLocator();

    cout << "Mesh number of points = " << mesh_points_alg.size() << endl;

    // 2) Read the LV endocardium mesh
    vtkSmartPointer<vtkGenericDataObjectReader> reader_2 = vtkSmartPointer<vtkGenericDataObjectReader>::New();
    reader_2->SetFileName(filename_mesh_endocardium_LV.c_str());
    reader_2->Update();

    vtkPolyData *polydata_grid_LV = reader_2->GetPolyDataOutput();
    int LV_endo_num_cells = polydata_grid_LV->GetNumberOfCells();
    int LV_endo_num_points = polydata_grid_LV->GetNumberOfPoints();

    cout << "LV endo number of points = " << LV_endo_num_points << endl;
    cout << "LV endo number of cells = " << LV_endo_num_cells << endl;

    // 3) Read the RV endocardium mesh
    vtkSmartPointer<vtkGenericDataObjectReader> reader_3 = vtkSmartPointer<vtkGenericDataObjectReader>::New();
    reader_3->SetFileName(filename_mesh_endocardium_RV.c_str());
    reader_3->Update();

    vtkPolyData *polydata_grid_RV = reader_3->GetPolyDataOutput();
    int RV_endo_num_cells = polydata_grid_RV->GetNumberOfCells();
    int RV_endo_num_points = polydata_grid_RV->GetNumberOfPoints();

    cout << "RV endo number of points = " << RV_endo_num_points << endl;
    cout << "RV endo number of cells = " << RV_endo_num_cells << endl;

    vector<vector<double>> stim_endocardium_points_LVRV;
    if (mode == 1) {
      stim_endocardium_points_LVRV.assign(RV_endo_num_points+LV_endo_num_points,vector<double>());
    }
    else if (mode == 2) {
      stim_endocardium_points_LVRV.assign(LV_endo_num_points,vector<double>());
    }
    else if (mode == 3) {
      stim_endocardium_points_LVRV.assign(RV_endo_num_points,vector<double>());
    }
    else {
      printf("[-] ERROR! Invalid 'mode'!\n");
      exit(EXIT_FAILURE);
    }
    
    if (mode == 1) {
      // Find the closest point to the RV endocardium in the mesh
      cout << "[!] Finding closest point in the RV endocardium ..." << endl;
      for (int i = 0, k = 0; i < RV_endo_num_points; i++, k++) {
        //printf("%d/%d\n", i, RV_endo_num_points);
        double pos[3];
        polydata_grid_RV->GetPoint(i,pos);
        int id = pointLocator->FindClosestPoint(pos);

        // Get the center position from the ALG cell
        double center[3];
        center[0] = mesh_points_alg[id].pos[0];
        center[1] = mesh_points_alg[id].pos[1];
        center[2] = mesh_points_alg[id].pos[2];

        stim_endocardium_points_LVRV[k].push_back(center[0]);
        stim_endocardium_points_LVRV[k].push_back(center[1]);
        stim_endocardium_points_LVRV[k].push_back(center[2]);
      }

      // Find the closest point to the LV endocardium in the mesh
      cout << "[!] Finding closest point in the LV endocardium ..." << endl;
      for (int i = 0, k = RV_endo_num_points; i < LV_endo_num_points; i++, k++) {
        //printf("%d/%d\n", i, LV_endo_num_points);
        double pos[3];
        polydata_grid_LV->GetPoint(i,pos);
        int id = pointLocator->FindClosestPoint(pos);

        // Get the center position from the ALG cell
        double center[3];
        center[0] = mesh_points_alg[id].pos[0];
        center[1] = mesh_points_alg[id].pos[1];
        center[2] = mesh_points_alg[id].pos[2];

        stim_endocardium_points_LVRV[k].push_back(center[0]);
        stim_endocardium_points_LVRV[k].push_back(center[1]);
        stim_endocardium_points_LVRV[k].push_back(center[2]);
      }
    }
    // Only LV
    else if (mode == 2) {
      cout << "[!] Finding closest point in the LV endocardium ..." << endl;
      for (int i = 0, k = 0; i < LV_endo_num_points; i++, k++) {
        double pos[3];
        polydata_grid_LV->GetPoint(i,pos);
        int id = pointLocator->FindClosestPoint(pos);

        // Get the center position from the ALG cell
        double center[3];
        center[0] = mesh_points_alg[id].pos[0];
        center[1] = mesh_points_alg[id].pos[1];
        center[2] = mesh_points_alg[id].pos[2];

        stim_endocardium_points_LVRV[k].push_back(center[0]);
        stim_endocardium_points_LVRV[k].push_back(center[1]);
        stim_endocardium_points_LVRV[k].push_back(center[2]);
      }
    }
    // Only RV
    else if (mode == 3) {
      cout << "[!] Finding closest point in the RV endocardium ..." << endl;
      for (int i = 0, k = 0; i < RV_endo_num_points; i++, k++) {
        double pos[3];
        polydata_grid_RV->GetPoint(i,pos);
        int id = pointLocator->FindClosestPoint(pos);

        // Get the center position from the ALG cell
        double center[3];
        center[0] = mesh_points_alg[id].pos[0];
        center[1] = mesh_points_alg[id].pos[1];
        center[2] = mesh_points_alg[id].pos[2];

        stim_endocardium_points_LVRV[k].push_back(center[0]);
        stim_endocardium_points_LVRV[k].push_back(center[1]);
        stim_endocardium_points_LVRV[k].push_back(center[2]);
      }
    }
  
    FILE *file = fopen(filename_output_endocardium_LVRV.c_str(),"w+");
    if (mode == 1) {
      fprintf(file,"%d\n",RV_endo_num_points+LV_endo_num_points);
    }
    else if (mode == 2) {
      fprintf(file,"%d\n",LV_endo_num_points);
    }
    else if (mode == 2) {
      fprintf(file,"%d\n",RV_endo_num_points);
    }
    for (int i = 0; i < stim_endocardium_points_LVRV.size(); i++) {
      fprintf(file,"%g %g %g\n", stim_endocardium_points_LVRV[i][0], stim_endocardium_points_LVRV[i][1], stim_endocardium_points_LVRV[i][2]);
    }
    fclose(file);
    
    file = fopen("outputs/test.vtk","w+");
    fprintf(file,"# vtk DataFile Version 4.1\n");
    fprintf(file,"vtk output\n");
    fprintf(file,"ASCII\n");
    fprintf(file,"DATASET POLYDATA\n");
    fprintf(file,"POINTS %d double\n",stim_endocardium_points_LVRV.size());
    for (int i = 0; i < stim_endocardium_points_LVRV.size(); i++) {
      fprintf(file,"%g %g %g\n", stim_endocardium_points_LVRV[i][0], stim_endocardium_points_LVRV[i][1], stim_endocardium_points_LVRV[i][2]);
    }
    fprintf(file,"VERTICES %d %d\n",stim_endocardium_points_LVRV.size(),stim_endocardium_points_LVRV.size()*2);
    for (int i = 0; i < stim_endocardium_points_LVRV.size(); i++) {
      fprintf(file,"1 %d\n", i);
    }
    fclose(file);

    return 0;

}
