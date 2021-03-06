/**
ModalAnalysisWG.cpp

This target demonstrate how to use NeoPZ for the modal analysis of a metal-backed
electromagnetic waveguide. For this example, the WR90 waveguide is used.

A HCurl-conforming approximation space is used for the transverse field components and
a H1-conforming approximation space for the axial component.
***/


#include <Electromagnetics/TPZWaveguideModalAnalysis.h> //for TPZMatWaveguideModalAnalysis
#include <MMeshType.h>                                  //for MMeshType
#include <TPZBndCond.h>                                 //for TPZBndCond
#include <TPZEigenAnalysis.h>                           //for TPZLinearAnalysis
#include <TPZGeoMeshTools.h>      //for TPZGeoMeshTools::CreateGeoMeshOnGrid
#include <TPZKrylovEigenSolver.h> //for TPZKrylovEigenSolver
#include <TPZLapackEigenSolver.h> //for TPZLapackEigenSolver
#include <TPZSpectralTransform.h> //for TPZSTShiftAndInvert
#include <TPZNullMaterial.h>      //for TPZNullMaterial
#include <TPZElectromagneticConstants.h>
#include <TPZSimpleTimer.h>
#include <TPZSkylineNSymStructMatrix.h> //non-symmetric skyline matrix storage
#include <TPZSpStructMatrix.h>          //non-symmetric sparse matrix storage
#include <pzcmesh.h>                    //for TPZCompMesh
#include <pzgmesh.h>                    //for TPZGeoMesh
#include <pzmanvector.h>                //for TPZManVector
#include "TPZVTKGeoMesh.h" //for exporting geomesh to vtk
#include <pzbuildmultiphysicsmesh.h>
#include <pzlog.h>


//! BC to be used if splitting domain in half (width-wise)
enum class ESymType { NONE, PMC, PEC };


/**
   @brief Creates the geometrical mesh associated with a rectangular metallic waveguide. 
The domain can be divided in half (width-wise) for exploring symmetry. */
TPZAutoPointer<TPZGeoMesh>
CreateGMeshRectangularWaveguide(TPZVec<int> &matIdVec, const REAL &scale,
                                const REAL wDomain, const REAL hDomain,
                                TPZVec<int> &nDivs,
                                const MMeshType meshType,
                                const bool usingSymmetry);

/**
   @brief Creates the computational meshes used for approximating the waveguide EVP.
   Three meshes will be created: one for the H1 approximation space, one for the
   HCurl approximation space and one multiphysics mesh combining both spaces.
   @note The material ids should be the same given to CreateGMeshRectangularWaveguide.
*/
TPZVec<TPZAutoPointer<TPZCompMesh>>
CreateCMesh(TPZAutoPointer<TPZGeoMesh> gmesh, int pOrder,
            const TPZVec<int> &matIdVec, const CSTATE ur, const CSTATE er,
            const STATE lambda, const REAL &scale, bool usingSymmetry, ESymType sym);


/** @brief Counts active equations per approximation space*/
void CountActiveEquations(TPZVec<TPZAutoPointer<TPZCompMesh>> meshVec,
                          const std::set<int64_t> &boundConnects,
                          int &neq,
                          int &nH1Equations, int &nHCurlEquations);

/** @brief Gets the indices of the equations associated with dirichlet homogeneous BCs
 so they can be filtered out of the global system*/
void FilterBoundaryEquations(TPZVec<TPZAutoPointer<TPZCompMesh>> meshVec,
                             TPZVec<int64_t> &activeEquations, int &neq,
                             int &neqOriginal, int& nh1, int &nhcurl);


int main(int argc, char *argv[]) {
#ifdef PZ_LOG
  /**if the NeoPZ library was configured with log4cxx,
   * the log should be initialised as:*/
  TPZLogger::InitializePZLOG();
#endif

  /***********************
   * setting the problem *
   ***********************/
  
  // width of the waveguide
  constexpr REAL wDomain{0.02286};
  // height of the waveguide
  constexpr REAL hDomain{0.01016};
  //whether to split the domain in two for exploring symmetries
  constexpr bool usingSymmetry{false};
  //which BC to apply if symmetry is used
  constexpr ESymType sym{ESymType::PEC};
  // n divisions in x direction
  constexpr int nDivX{32};
  // n divisions in y direction
  constexpr int nDivY{16};
  // type of elements
  constexpr MMeshType meshType{MMeshType::ETriangular};
  // TPZManVector<Type,N> is a vector container with static + dynamic storage.
  // one can also use TPZVec<Type> for dynamic storage
  TPZManVector<int, 2> nDivs = {nDivX, nDivY};
  // polynomial order to be used in the approximation
  constexpr int pOrder{1};
  // magnetic permeability
  constexpr CSTATE ur{1};
  // electric permittivity
  constexpr CSTATE er{1};
  // operational frequency
  constexpr STATE freq = 25e9;
  // wavelength
  constexpr STATE lambda{pzeletromag::cZero/freq};
  /*Given the small dimensions of the domain, scaling it can help in 
    achieving good precision. Uing k0 as a scale factor results in 
    the eigenvalues -(propagationConstant/k0)^2 = -effectiveIndex^2*/
  constexpr REAL scale{2*M_PI/lambda};


  /******************
   * solver options *
   ******************/
  
  //number of threads to use
  constexpr int nThreads{4};
  //number of genvalues to be computed
  constexpr int nEigenpairs{10};
  //whether to compute eigenvectors (instead of just eigenvalues)
  constexpr bool computeVectors{true};
  //how to sort the computed eigenvalues
  constexpr TPZEigenSort sortingRule {TPZEigenSort::TargetRealPart};
  /*
   The simulation uses a Krylov-based Arnoldi solver for solving the
   generalised EVP. A shift-and-inverse spectral transform is applied in 
   the system for better convergence of the eigenvalues.
   The target variable should be close to the desired eigenvalue (in this case,
   the effective index neff). We are looking for neff=1, i.e., the modes with
  lowest cutoff frequency*/
  constexpr CSTATE target = -1;
  // Dimension of the krylov space to be used. Suggested to be at least nev * 10
  constexpr int krylovDim{200};


  /*********************
   * exporting options *
   *********************/

  //whether to print the geometric mesh in .txt and .vtk formats
  constexpr bool printGMesh{false};
  //whether to export the solution as a .vtk file
  constexpr bool exportVtk{true};
  //resolution of the .vtk file in which the solution will be exported
  constexpr int vtkRes{1};
  //if true, the real part of the electric fields is exported. otherwise, the magnitude
  constexpr bool printRealPart{true};

  /********************
   * advanced options *
   ********************/
  
  //reorder the equations in order to optimize bandwidth
  constexpr bool optimizeBandwidth{true};
  /*
    The equations corresponding to homogeneous dirichlet boundary condition(PEC)
    can be filtered out of the global system in order to achieve better conditioning.
   */
  constexpr bool filterBoundaryEqs{true};

  /*********
   * begin *
   *********/

  
  //scoped-timer 
  TPZSimpleTimer total("Total");

  /* Vector for storing materials(regions) identifiers for the 
     TPZGeoMeshTools::CreateGeoMeshOnGrid function.
     In this example, we have one material for the interior of the waveguide
     and one for each BC*/
  TPZManVector<int, 5> matIdVec({1,-1,-2,-3,-4});
  //creates geometric mesh
  auto gmesh = CreateGMeshRectangularWaveguide(matIdVec,scale,wDomain,hDomain,nDivs,
                                               meshType,usingSymmetry);

  
  //print gmesh to vtk
  if(printGMesh)
  {
    const std::string gmeshFileNameTxt{"gmesh_" +
    std::to_string(nDivX) + " _" + std::to_string(nDivY) + ".vtk"};
    std::ofstream gmeshFileTxt(gmeshFileNameTxt);
    gmesh->Print(gmeshFileTxt);
    const std::string gmeshFileNameVtk{"gmesh_" +
    std::to_string(nDivX) + " _" + std::to_string(nDivY) + ".vtk"};
    std::ofstream gmeshFileVtk(gmeshFileNameVtk);
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, gmeshFileVtk, true);
  }

  /*
   The problem uses an H1 approximation space for the longitudinal component 
   and a HCurl approximation space for the transversal one. Therefore, three
  computational meshes are generated. One for each space and a multiphysics mesh*/
  auto meshVec = CreateCMesh(gmesh,pOrder,matIdVec,ur,er,lambda,
                             scale,usingSymmetry,sym);
  //gets the multiphysics mesh (main mesh)
  auto cmesh = meshVec[0];

  
  TPZEigenAnalysis an(cmesh, optimizeBandwidth);
  an.SetComputeEigenvectors(computeVectors);

  /**
     When using NeoPZ with MKL, an sparse matrix should be used for better
     performance. Otherwise, the skyline matrix has available inhouse solvers.
  */
  TPZAutoPointer<TPZStructMatrix> strmtrx{nullptr};
#ifdef PZ_USING_MKL
  strmtrx = new TPZSpStructMatrix<CSTATE>(cmesh);
#else
  strmtrx = new TPZSkylineNSymStructMatrix<CSTATE>(cmesh);
#endif
  
  strmtrx->SetNumThreads(nThreads);
  
  TPZVec<int64_t> activeEquations;
  //this value is the total number of dofs including dirichlet bcs
  int neq, neqOriginal, neqH1, neqHCurl;
  
  
  if(filterBoundaryEqs){
    FilterBoundaryEquations(meshVec, activeEquations, neq, neqOriginal, neqH1, neqHCurl);
    std::cout<<"neq(before): "<<neqOriginal
             <<"\tneq(after): "<<neq<<std::endl;
    strmtrx->EquationFilter().SetActiveEquations(activeEquations);
  }else{
    std::set<int64_t> boundConnects;
    CountActiveEquations(meshVec,boundConnects,neqOriginal,neqH1,neqHCurl);
  }
  
  an.SetStructuralMatrix(strmtrx);

  
  
  TPZKrylovEigenSolver<CSTATE> solver;
  TPZSTShiftAndInvert<CSTATE> st;
  solver.SetSpectralTransform(st);
  solver.SetKrylovDim(krylovDim);

  {
    /**this is to ensure that the eigenvector subspace is orthogonal to
       the spurious solutions associated with et = 0 ez != 0*/
    TPZFMatrix<CSTATE> initVec(neq, 1, 0.);
    constexpr auto firstHCurl = TPZWaveguideModalAnalysis::HCurlIndex();
    for (int i = 0; i < neqHCurl; i++) {
      initVec(firstHCurl + i, 0) = 1;
    }
    solver.SetKrylovInitialVector(initVec);
  }

  solver.SetTarget(target);
  solver.SetNEigenpairs(nEigenpairs);
  solver.SetAsGeneralised(true);
  solver.SetEigenSorting(sortingRule);

  an.SetSolver(solver);
  {
    std::cout<<"Assembling..."<<std::flush;
    TPZSimpleTimer assemble("Assemble");
    an.Assemble();
    std::cout<<"\rAssembled!"<<std::endl;
  }
  {
    TPZSimpleTimer solv("Solve");
    std::cout<<"Solving..."<<std::flush;
    an.Solve();
    std::cout<<"\rSolved!"<<std::endl;
  }
  auto ev = an.GetEigenvalues();

  for(auto &w : ev){
    std::cout<<w<<std::endl;
  }
  
  if (!computeVectors && !exportVtk) return 0;
  
  TPZStack<std::string> scalnames, vecnames;
  scalnames.Push("Ez");
  vecnames.Push("Et");
  const std::string plotfile = "fieldPlot.vtk";
  constexpr int dim{2};
  an.DefineGraphMesh(dim, scalnames, vecnames,plotfile);

  auto eigenvectors = an.GetEigenvectors();

  TPZFMatrix<CSTATE> evector(neqOriginal, 1, 0.);
  const auto nev = ev.size();

  TPZManVector<TPZAutoPointer<TPZCompMesh>,2> meshVecPost(2);
  meshVecPost[0] = meshVec[1];
  meshVecPost[1] = meshVec[2];
  
  std::cout<<"Post processing..."<<std::endl;
  
  for (int iSol = 0; iSol < ev.size(); iSol++) {
    const CSTATE currentKz = [&ev,iSol](){
      auto tmp = std::sqrt(-1.0*ev[iSol]);
      constexpr auto epsilon = std::numeric_limits<STATE>::epsilon()/
      (10*std::numeric_limits<STATE>::digits10);
      //let us discard extremely small imag parts
      if (tmp.imag() < epsilon)
        {tmp = tmp.real();}
      return tmp;
    }();
    eigenvectors.GetSub(0, iSol, neqOriginal, 1, evector);
    for(auto id : matIdVec){
      auto matPtr =
        dynamic_cast<TPZWaveguideModalAnalysis *>(cmesh->FindMaterial(id));
      if(!matPtr) continue;
      matPtr->SetKz(currentKz);
      matPtr->SetPrintFieldRealPart(printRealPart);
    }
    std::cout<<"\rPost processing step "<<iSol+1<<" out of "<<ev.size()
             <<"(kz = "<<currentKz<<")"<<std::flush;
    an.LoadSolution(evector);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshVecPost, cmesh);
    an.PostProcess(vtkRes);
  }
  std::cout<<"\rFinished post processing"<<std::endl;
  std::cout<<std::endl;
  return 0;
}






TPZAutoPointer<TPZGeoMesh> CreateGMeshRectangularWaveguide(
    TPZVec<int> &matIdVec, const REAL &scale, const REAL wDomain,
    const REAL hDomain, TPZVec<int> &nDivs, const MMeshType meshType,
    const bool usingSymmetry)
{
  TPZSimpleTimer timer("Create gmesh");
  // dimension of the problem
  constexpr int dim{2};
  // lower left corner of the domain
  TPZManVector<REAL, 3> minX = {0, 0, 0};
  
  // upper right corner of the domain
  TPZManVector<REAL, 3> maxX = {
    usingSymmetry ? wDomain * scale / 2 : wDomain * scale,
    hDomain * scale,
    0};
  
  // whether to create boundary elements
  constexpr bool genBoundEls{true};
  
  return TPZGeoMeshTools::CreateGeoMeshOnGrid(dim,minX,maxX,matIdVec,nDivs,
                                              meshType,genBoundEls);
}

TPZVec<TPZAutoPointer<TPZCompMesh>>
CreateCMesh(TPZAutoPointer<TPZGeoMesh> gmesh, int pOrder,
            const TPZVec<int> &matIdVec, const CSTATE ur, const CSTATE er,
            const STATE lambda, const REAL &scale, bool usingSymmetry, ESymType sym)
{
  TPZSimpleTimer timer ("Create cmesh");
  constexpr int dim = 2;
  constexpr bool isComplex{true};


  /*
   First we create the computational mesh associated with the H1 space
   (ez component)*/
  auto * cmeshH1 =new TPZCompMesh(gmesh,isComplex);
  cmeshH1->SetDefaultOrder(pOrder +1);//for deRham compatibility
  cmeshH1->SetDimModel(dim);

  const int volMatId = matIdVec[0];
  //number of state variables in the problem
  constexpr int nState = 1;

  auto dummyMat = new TPZNullMaterial<CSTATE>(volMatId,dim,nState);
  cmeshH1->InsertMaterialObject(dummyMat);

  
  TPZFNMatrix<1, CSTATE> val1(1, 1, 1);
  TPZManVector<CSTATE,1> val2(1, 0.);
  const int nMats = matIdVec.size();
  TPZBndCond *dummyBC = nullptr;
  for (int i = 1; i <nMats; i++) {
    //0 for dirichlet (PEC) and 1 for neumann (PMC)
    const int bcType = i==1 && usingSymmetry && sym == ESymType::PMC ? 1 : 0;
    dummyBC =
      dummyMat->CreateBC(dummyMat, matIdVec[i], bcType, val1, val2);
    cmeshH1->InsertMaterialObject(dummyBC);
  }

  cmeshH1->SetAllCreateFunctionsContinuous();
  cmeshH1->AutoBuild();
  cmeshH1->CleanUpUnconnectedNodes();

  /*
    Then we create the computational mesh associated with the HCurl space
   */
  auto *cmeshHCurl = new TPZCompMesh(gmesh,isComplex);
  cmeshHCurl->SetDefaultOrder(pOrder);
  cmeshHCurl->SetDimModel(dim);
  
  dummyMat = new TPZNullMaterial<CSTATE>(volMatId,dim,nState);
  cmeshHCurl->InsertMaterialObject(dummyMat);

  
  dummyBC = nullptr;
  for (int i = 1; i <nMats; i++) {
    //0 for dirichlet (PEC) and 1 for neumann (PMC)
    const int bcType = i==1 && usingSymmetry && sym == ESymType::PMC ? 1 : 0;
    dummyBC =
      dummyMat->CreateBC(dummyMat, matIdVec[i], bcType, val1, val2);
    cmeshHCurl->InsertMaterialObject(dummyBC);
  }

  cmeshHCurl->SetAllCreateFunctionsHCurl();
  cmeshHCurl->AutoBuild();
  cmeshHCurl->CleanUpUnconnectedNodes();

  
  auto *cmeshMF =new TPZCompMesh(gmesh,isComplex);
  
  TPZWaveguideModalAnalysis *matWG  = new TPZWaveguideModalAnalysis(
      volMatId, ur, er, lambda, 1. / scale);
  cmeshMF->InsertMaterialObject(matWG);

  TPZBndCond *bcMat = nullptr;
  for (int i = 1; i <nMats; i++) {
    //0 for dirichlet (PEC) and 1 for neumann (PMC)
    const int bcType = i==1 && usingSymmetry && sym == ESymType::PMC ? 1 : 0;
    bcMat =
      matWG->CreateBC(matWG, matIdVec[i], bcType, val1, val2);
    cmeshMF->InsertMaterialObject(bcMat);
  }

  cmeshMF->SetDimModel(dim);
  cmeshMF->SetAllCreateFunctionsMultiphysicElem();

  cmeshMF->AutoBuild();
  cmeshMF->CleanUpUnconnectedNodes();

  TPZManVector<TPZCompMesh*,3> meshVecIn(2);
  meshVecIn[TPZWaveguideModalAnalysis::H1Index()] = cmeshH1;
  meshVecIn[TPZWaveguideModalAnalysis::HCurlIndex()] = cmeshHCurl;

  
  TPZBuildMultiphysicsMesh::AddElements(meshVecIn, cmeshMF);
  TPZBuildMultiphysicsMesh::AddConnects(meshVecIn, cmeshMF);
  TPZBuildMultiphysicsMesh::TransferFromMeshes(meshVecIn, cmeshMF);

  cmeshMF->ExpandSolution();
  cmeshMF->ComputeNodElCon();
  cmeshMF->CleanUpUnconnectedNodes();

  TPZVec<TPZAutoPointer<TPZCompMesh>> meshVec(3,nullptr);

  meshVec[0] = cmeshMF;
  meshVec[1 + TPZWaveguideModalAnalysis::H1Index()] = cmeshH1;
  meshVec[1 + TPZWaveguideModalAnalysis::HCurlIndex()] = cmeshHCurl;
  return meshVec;
}

void CountActiveEquations(TPZVec<TPZAutoPointer<TPZCompMesh>> meshVec,
                          const std::set<int64_t> &boundConnects,
                          int &neq,
                          int &nH1Equations, int &nHCurlEquations)
{
  auto cmesh = meshVec[0];
  neq = nH1Equations = nHCurlEquations = 0;
  auto cmeshHCurl = meshVec[1];
  auto cmeshH1 = meshVec[2];
  
  for (int iCon = 0; iCon < cmesh->NConnects(); iCon++) {
    bool isH1;
    if (boundConnects.find(iCon) == boundConnects.end()) {
      if (cmesh->ConnectVec()[iCon].HasDependency())
        continue;
      int seqnum = cmesh->ConnectVec()[iCon].SequenceNumber();
      int blocksize = cmesh->Block().Size(seqnum);
      if (TPZWaveguideModalAnalysis::H1Index() == 0 && iCon < cmeshH1->NConnects()) {
        isH1 = true;
      } else if (TPZWaveguideModalAnalysis::H1Index() == 1 && iCon >= cmeshHCurl->NConnects()) {
        isH1 = true;
      } else {
        isH1 = false;
      }
      for (int ieq = 0; ieq < blocksize; ieq++) {
        neq++;
        isH1 == true ? nH1Equations++ : nHCurlEquations++;
      }
    }
  }
  std::cout << "------\tactive eqs\t-------" << std::endl;
  std::cout << "# H1 equations: " << nH1Equations << std::endl;
  std::cout << "# HCurl equations: " << nHCurlEquations << std::endl;
  std::cout << "# equations: " << neq << std::endl;
  std::cout << "------\t----------\t-------" << std::endl;
  return;
}

void FilterBoundaryEquations(TPZVec<TPZAutoPointer<TPZCompMesh>> meshVec,
                             TPZVec<int64_t> &activeEquations, int &neq,
                             int &neqOriginal,
                             int &nh1, int &nhcurl)
{
  TPZSimpleTimer timer ("Filter dirichlet eqs");
  auto cmesh = meshVec[0];
  TPZManVector<int64_t, 1000> allConnects;
  std::set<int64_t> boundConnects;

  for (int iel = 0; iel < cmesh->NElements(); iel++) {
    TPZCompEl *cel = cmesh->ElementVec()[iel];
    if (cel == nullptr) {
      continue;
    }
    if (cel->Reference() == nullptr) {//there is no associated geometric el
      continue;
    }
    TPZBndCond *mat = dynamic_cast<TPZBndCond *>(
        cmesh->MaterialVec()[cel->Reference()->MaterialId()]);
    if (mat && mat->Type() == 0) {//check for dirichlet bcs
      std::set<int64_t> boundConnectsEl;
      std::set<int64_t> depBoundConnectsEl;
      std::set<int64_t> indepBoundConnectsEl;
      cel->BuildConnectList(indepBoundConnectsEl, depBoundConnectsEl);
      cel->BuildConnectList(boundConnectsEl);
      for(auto val : boundConnectsEl){
        if (boundConnects.find(val) == boundConnects.end()) {
          boundConnects.insert(val);
        }
      }
    }
  }
  for (int iCon = 0; iCon < cmesh->NConnects(); iCon++) {
    if (boundConnects.find(iCon) == boundConnects.end()) {
      TPZConnect &con = cmesh->ConnectVec()[iCon];
      if (con.HasDependency())
        continue;
      const auto seqnum = con.SequenceNumber();
      const auto pos = cmesh->Block().Position(seqnum);
      const auto blocksize = cmesh->Block().Size(seqnum);
      if (blocksize == 0)
        continue;

      const auto vs = activeEquations.size();
      activeEquations.Resize(vs + blocksize);
      for (auto ieq = 0; ieq < blocksize; ieq++) {
        activeEquations[vs + ieq] = pos + ieq;
      }
    }
  }

  neqOriginal = cmesh->NEquations();
  CountActiveEquations(meshVec, boundConnects, neq, nh1, nhcurl);
}

