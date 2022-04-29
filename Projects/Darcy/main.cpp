
#include <cmath>
#include <set>

#include <iostream>
#include <fstream>
#include <string>
#include "pzcmesh.h"
#include "pzgmesh.h"
#include "pzstack.h"
#include "TPZVTKGeoMesh.h"
#include "TPZAnalysis.h"
#include <TPZLinearAnalysis.h> //for TPZLinearAnalysis
#include <TPZSSpStructMatrix.h>
#include "TPZBndCond.h"
#include "DarcyFlow/TPZDarcyFlow.h"
//#include "TPZDarcyMaterial.h"

#include <pzgeoel.h>
#include "pzgeoelbc.h"
#include "pzfmatrix.h"
#include "pzbstrmatrix.h"
#include <TPZGeoElement.h>
#include "TPZVTKGeoMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZGeoLinear.h"
#include "tpzgeoelrefpattern.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "TPZGmshReader.h"
#include <pzgmesh.h> //for TPZGeoMesh
#include <pzcmesh.h> //for TPZCompMesh
#include <TPZGeoMeshTools.h> //for TPZGeoMeshTools::CreateGeoMeshOnGrid
#include <MMeshType.h> //for MMeshType
#include <pzmanvector.h>//for TPZManVector
#include <Poisson/TPZMatPoisson.h> //for TPZMatPoisson
#include <TPZBndCond.h> //for TPZBndCond
#include <TPZLinearAnalysis.h> //for TPZLinearAnalysis
#include <TPZSSpStructMatrix.h> //symmetric sparse matrix storage
#include <pzskylstrmatrix.h> //symmetric skyline matrix storage
#include <pzstepsolver.h> //for TPZStepSolver
#include <TPZSimpleTimer.h>
#include "TPZBndCondT.h"
#include "DarcyFlow/TPZHybridDarcyFlow.h"
class TPZCompMesh;
#include "DarcyFlow/TPZMixedDarcyFlow.h"
class TPZGeoMesh;
using namespace std;
class TPZMaterial;

int bottombc=-1;
int rigthbc = -2;
int topbc = -3;
int leftbc =-4;



int bottombc_slope=-1;
int rigthbc_slope = -2;
int leftbc_slope =-3;
int toprigthbc_slope = -4;
int topleftbc_slope = -5;
int rampbc_slope = -6;

// brief Function to create the geometric mesh
TPZGeoMesh *CreateGMesh (int ref );

TPZGeoMesh *  CreateGMeshSlope ( int ref );

TPZCompMesh * CreateCMeshSlopeFlow ( TPZGeoMesh *gmesh, int pOrder );

TPZCompMesh * CreateCMesh ( TPZGeoMesh *gmesh, int pOrder );

int main()
{
	
	int porder=1;
	
	int ref=0;
	
    //TPZGeoMesh *gmesh = CreateGMesh( ref);
	TPZGeoMesh *gmesh = CreateGMeshSlope( ref);

    std::ofstream meshfile ( "GeoMeshSlope.txt" );
	
    gmesh->Print ( meshfile );

    // Create computational mesh
    //TPZCompMesh *cmesh = CreateCMesh ( gmesh, porder );
	
	TPZCompMesh *cmesh = CreateCMeshSlopeFlow ( gmesh, porder );
	
    std::ofstream comVfile ( "CompMeshSlope.txt" );
	
    cmesh->Print ( comVfile );

    //Resolvendo o Sistema:
    int numthreads = 0;
	
    bool optimizeBandwidth = true; // Prevents of renumbering of the equations (As the same of Oden's result)
    
    TPZLinearAnalysis analysis ( cmesh,optimizeBandwidth ); // Create analysis

    //sets number of threads to be used by the solver

    //defines storage scheme to be used for the FEM matrices
    //in this case, a symmetric skyline matrix is used
    TPZSkylineStructMatrix<STATE> matskl ( cmesh );
	
    matskl.SetNumThreads ( numthreads );
	
    analysis.SetStructuralMatrix ( matskl );

    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect ( ELDLt );
    analysis.SetSolver ( step );
	
	analysis.Assemble();
	
	std::ofstream sout ( "rhs.nb" );
	analysis.Rhs().Print("Rhs = ",sout,EMathematicaInput);

    analysis.Run();

	
	
    ///vtk export
    TPZVec<std::string> scalarVars ( 1 ),vectorVars ( 1 );
    scalarVars[0] = "Pressure";
    vectorVars[0] = "GradU";
    analysis.DefineGraphMesh ( 2,scalarVars,vectorVars,"Darcy.vtk" );
    constexpr int resolution{0};
    analysis.PostProcess ( resolution );

    std::cout << "FINISHED!" << std::endl;


    return 0;
}


TPZGeoMesh *  CreateGMesh ( int ref )
{
    const std::string name ( "Darcy Flow " );

    TPZGeoMesh *gmesh  =  new TPZGeoMesh();

    gmesh->SetName ( name );
	
    gmesh->SetDimension ( 2 );

    TPZVec<REAL> coord ( 2 );

    vector<vector<double>> co= {{0., 0.}, {1., 0.}, {1., 1.}, {0., 1.}};
	
    vector<vector<int>> topol = {{0,  1,  2, 3}};

    gmesh->NodeVec().Resize ( co.size() );

    for ( int inode=0; inode<co.size(); inode++ ) {
        coord[0] = co[inode][0];
        coord[1] = co[inode][1];
        gmesh->NodeVec() [inode] = TPZGeoNode ( inode, coord, *gmesh );
    }
    TPZVec <long> TopoQuad ( 4 );
    for ( int iel=0; iel<topol.size(); iel++ ) {
        TopoQuad[0] = topol[iel][0];
        TopoQuad[1] = topol[iel][1];
        TopoQuad[2] =	topol[iel][2];
        TopoQuad[3] = topol[iel][3];
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( iel, TopoQuad, 1,*gmesh );
    }
    int id = topol.size();
    TPZVec <long> TopoLine ( 2 );
    TopoLine[0] = 0;
    TopoLine[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, bottombc, *gmesh );//bottom

    id++;
    TopoLine[0] = 1;
    TopoLine[1] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, rigthbc, *gmesh );//rigth

    id++;
    TopoLine[0] = 2;
    TopoLine[1] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, topbc, *gmesh );//left

    id++;
    TopoLine[0] = 3;
    TopoLine[1] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, leftbc, *gmesh ); //top left

    gmesh->BuildConnectivity();
    for ( int d = 0; d<ref; d++ ) {
        int nel = gmesh->NElements();
        TPZManVector<TPZGeoEl *> subels;
        for ( int iel = 0; iel<nel; iel++ ) {
            TPZGeoEl *gel = gmesh->ElementVec() [iel];
            gel->Divide ( subels );
        }
    }

    
    return gmesh;
}

TPZCompMesh * CreateCMesh ( TPZGeoMesh *gmesh, int pOrder )
{
    // Creating computational mesh:
    TPZCompMesh * cmesh = new TPZCompMesh ( gmesh );
    
	cmesh->SetDefaultOrder ( pOrder );
    
	int dim = 2 ;
	
	cmesh->SetDimModel ( dim );
    
	cmesh->SetAllCreateFunctionsContinuous();

	int matid=1;
	
    auto *material = new TPZDarcyFlow ( matid,dim );
	
	

	REAL permeability = 1.;
    material->SetConstantPermeability ( permeability );
	
  
  const auto rhs = [](const TPZVec<REAL>&loc, TPZVec<STATE> &u){
        const REAL &x = loc[0];
        const REAL &y = loc[1];
        u[0] = -1;
  };
	material->SetForcingFunction(rhs,pOrder);
	
    cmesh->InsertMaterialObject ( material );

    // Condition of contours
    TPZFMatrix<STATE>  val1 ( 2,2,0. );
	
    TPZManVector<STATE,2> val2 ( 2,0. );


	int dirichlet = 0 ;
	
	auto pressure = [](const TPZVec<REAL>&loc, TPZVec<STATE> &u, TPZFMatrix<REAL> &grad){
		REAL &x = loc[0];
		REAL &y = loc[1];
		REAL &z = loc[2];
        u[0] = x;

  	};
	
    auto * BCond0 = material->CreateBC ( material, bottombc, dirichlet, val1, val2 );

    auto * BCond1 = material->CreateBC ( material, rigthbc, dirichlet, val1, val2 );
	
	auto * BCond2 = material->CreateBC ( material, topbc, dirichlet, val1, val2 );
	//BCond2->SetForcingFunctionBC(pressure);
	
	auto * BCond3 = material->CreateBC ( material, leftbc, dirichlet, val1, val2 );

	cmesh->InsertMaterialObject ( BCond0 );
    cmesh->InsertMaterialObject ( BCond1 );
	cmesh->InsertMaterialObject ( BCond2 );
	cmesh->InsertMaterialObject ( BCond3 );
	
    //Creating computational elements that manage the space of the mesh:
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();

    return cmesh;
}


TPZGeoMesh *  CreateGMeshSlope ( int ref )
{
    const std::string name ( "Darcy Flow Slope" );

    TPZGeoMesh *gmesh  =  new TPZGeoMesh();

    gmesh->SetName ( name );
	
    gmesh->SetDimension ( 2 );

    TPZVec<REAL> coord ( 2 );

vector<vector<double>> co= {{0., 0.}, {75., 0.}, {75., 30.}, {45., 30.}, {35., 40.}, {0.,40.},{35./3., 40.},{2 * 35/3., 40.}, {30., 40.},{30., 30.}, {60.,30.},{2* 35./3.,2* 35/3.}, {45., 2* 35/3.},{35./3., 35/3.}, {60., 35./3.}};

vector<vector<int>> topol = {{0,  1,  14, 13},{1,  2,  10, 14}, {14, 10, 3,  12}, {13, 14, 12, 11},{11, 12, 3,  9}, {9,  3,  4,  8},{11, 9,  8,  7},{13, 11, 7, 6},{0, 13,  6, 5}};

    gmesh->NodeVec().Resize ( co.size() );

    for ( int inode=0; inode<co.size(); inode++ ) {
        coord[0] = co[inode][0];
        coord[1] = co[inode][1];
        gmesh->NodeVec() [inode] = TPZGeoNode ( inode, coord, *gmesh );
    }
    TPZVec <long> TopoQuad ( 4 );
    for ( int iel=0; iel<topol.size(); iel++ ) {
        TopoQuad[0] = topol[iel][0];
        TopoQuad[1] = topol[iel][1];
        TopoQuad[2] =	topol[iel][2];
        TopoQuad[3] = topol[iel][3];
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> ( iel, TopoQuad, 1,*gmesh );
    }
    int id = topol.size();
    TPZVec <long> TopoLine ( 2 );
    TopoLine[0] = 0;
    TopoLine[1] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, bottombc_slope, *gmesh );//bottom

    id++;
    TopoLine[0] = 1;
    TopoLine[1] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, rigthbc_slope, *gmesh );//rigth

    id++;
    TopoLine[0] = 2;
    TopoLine[1] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, toprigthbc_slope, *gmesh );//top-rigth

    id++;
    TopoLine[0] = 3;
    TopoLine[1] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, rampbc_slope, *gmesh ); //ramp
	
	id++;
    TopoLine[0] = 4;
    TopoLine[1] = 5;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, topleftbc_slope, *gmesh ); //top-left
	
	id++;
    TopoLine[0] = 5;
    TopoLine[1] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> ( id, TopoLine, leftbc_slope, *gmesh ); //left

    gmesh->BuildConnectivity();
    for ( int d = 0; d<ref; d++ ) {
        int nel = gmesh->NElements();
        TPZManVector<TPZGeoEl *> subels;
        for ( int iel = 0; iel<nel; iel++ ) {
            TPZGeoEl *gel = gmesh->ElementVec() [iel];
            gel->Divide ( subels );
        }
    }

    
    return gmesh;
}

TPZCompMesh * CreateCMeshSlopeFlow ( TPZGeoMesh *gmesh, int pOrder )
{
	    // Creating computational mesh:
    TPZCompMesh * cmesh = new TPZCompMesh ( gmesh );
    
	cmesh->SetDefaultOrder ( pOrder );
    
	int dim = 2 ;
	
	cmesh->SetDimModel ( dim );
    
	cmesh->SetAllCreateFunctionsContinuous();

	int matid=1;
	
    auto *material = new TPZDarcyFlow ( matid,dim );
	//auto *material = new TPZHybridDarcyFlow ( matid,dim );
	//auto *material = new TPZMixedDarcyFlow ( matid,dim );
	
	material->SetId(matid);
	
	REAL permeability = 1.;
    material->SetConstantPermeability ( permeability );
	
  
  const auto rhs = [](const TPZVec<REAL>&loc, TPZVec<STATE> &u){
        const REAL &x = loc[0];
        const REAL &y = loc[1];
        u[0] = -1;
  };
	material->SetForcingFunction(rhs,pOrder);
	


    // Condition of contours
    TPZFMatrix<STATE>  val1 ( 2,2,0. );
	
    TPZManVector<STATE,1> val2 ( 1,0. );


	int dirichlet = 0 ;
	
	auto pressure = [](const TPZVec<REAL>&loc, TPZVec<STATE> &u, TPZFMatrix<REAL> &grad){
		REAL &x = loc[0];
		REAL &y = loc[1];
		REAL &z = loc[2];
        u[0] = x;

  	};
	
    auto * BCond0 = material->CreateBC ( material, bottombc_slope, dirichlet, val1, val2 );

    auto * BCond1 = material->CreateBC ( material, rigthbc_slope, dirichlet, val1, val2 );
	
	auto * BCond2 = material->CreateBC ( material, toprigthbc_slope, dirichlet, val1, val2 );

	auto * BCond3 = material->CreateBC ( material, rampbc_slope, dirichlet, val1, val2 );
	
	auto * BCond4 = material->CreateBC ( material, topleftbc_slope, dirichlet, val1, val2 );
	
	auto * BCond5 = material->CreateBC ( material, leftbc_slope, dirichlet, val1, val2 );

	cmesh->InsertMaterialObject ( material );
	cmesh->InsertMaterialObject ( BCond0 );
    cmesh->InsertMaterialObject ( BCond1 );
	cmesh->InsertMaterialObject ( BCond2 );
	cmesh->InsertMaterialObject ( BCond3 );
	cmesh->InsertMaterialObject ( BCond4 );
	cmesh->InsertMaterialObject ( BCond5 );
	
    //Creating computational elements that manage the space of the mesh:
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
	cmesh->SetAllCreateFunctionsContinuous();
	
    return cmesh;
}

