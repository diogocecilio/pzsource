/* Generated by Together */// $Id: TPZLadeKim.h,v 1.25 2010-06-11 22:12:14 diogo Exp $
#ifndef TPZLADEKIM_H
#define TPZLADEKIM_H

#include "pzlog.h"
#include "TPZPlasticStep.h"
#include "TPZYCLadeKim.h"
#include "TPZLadeKimThermoForceA.h"
#include "TPZLadeNelsonElasticResponse.h"
#include "pzvec_extras.h"
#include "TPZPlasticStepID.h"

//#ifdef PZ_LOG
//    TPZLogger loggerLadeKim("plasticity.LadeKim");
//#endif

#define LADEKIMPARENT TPZPlasticStep<TPZYCLadeKim, TPZLadeKimThermoForceA, TPZLadeNelsonElasticResponse>

/**
 * This class implements the LADE KIM constitutive model
 * comprised of a LadeKim Plastic flow rule and yield criterium
 * and of a LadeNelson elastic response.
 * The initialization takes care of ensuring that the starting
 * plastic state is consistent with the starting stress state
 * (unstressed)
 * Because of some numerical instabilities presented by this model,
 * the first steps in loading a material may lead to unexpected 
 * problems, generally leading to invalid numerical evaluations
 * and NotANumbers. To avoid this, it is suggested that the material
 * be loaded with low stress state (let's say, 0.001 strain hydrostatic)
 * with error tolerances of the same magnitude. Greater values for
 * these variales may lead to undesired integration deviations
 * and values lower night cause the integration to crash.
 * after this first step any integration tolerance and strain increment
 * sizes might be used. This problem is likely to be related to the
 * fact that the elastic response and the plastic flow rules are 
 * singular at null stresses.
*/
class TPZLadeKim : public LADEKIMPARENT  {

public:

  enum {NYield = TPZYCLadeKim::NYield};
    

public:

    TPZLadeKim():LADEKIMPARENT(), faPa(0.), fInitialEps()
    {
		fMaterialTensionSign  = -1; // internally in this material tension is negative
		fInterfaceTensionSign =  1; // by default
    }
	
    TPZLadeKim(const TPZLadeKim & source):LADEKIMPARENT(source)
    {
		faPa		= source.faPa;
		fInitialEps = source.fInitialEps;
    }

    TPZLadeKim & operator=(const TPZLadeKim & source)
    {
		LADEKIMPARENT::operator=(source);
		faPa		= source.faPa;
		fInitialEps = source.fInitialEps;
		
		return *this;
    }
	
	virtual const char * Name() const override
	{
	   return "TPZLadeKim";	
	}
	
    /**
    SetUp feeds all the parameters necessary to the method, distributing its values
    inside the aggregation hierarchy and computing the correct initial plasticity 
    parameter to ensure the irreversibility effect of plastic deformations.
    Elastic Mudulus:    poisson, M,    lambda
    Failure Criterion:  a,       m,    neta1
    Plastic Potential:  ksi2,    mu
    Hardening Function: C,       p
    Yield Function:     h,       alpha
    Atmospheric pression pa - to dimensionalize/adim. the stresses
    */
    void SetUp(REAL poisson, REAL M, REAL lambda,
                      REAL a, REAL m, REAL neta1,
                      REAL ksi2, REAL mu,
                      REAL C, REAL p,
                      REAL h, REAL alpha,
                      REAL pa)
    {
	   int interfaceCompressionSign = - fInterfaceTensionSign;
       faPa = interfaceCompressionSign * a * fabs(pa);
       REAL ksi1 = 0.00155 * pow(m, -1.27);
       LADEKIMPARENT::fYC.SetUp(ksi1, ksi2, h, alpha, mu, neta1, m, fabs(pa));
       LADEKIMPARENT::fER.SetUp(lambda, M, poisson, fabs(pa));
       LADEKIMPARENT::fTFA.SetUp(ksi1, p, h, C, fabs(pa));
       TPZTensor<REAL> nullSigma, epsA;
       // The function below should apply the plastic loop to the unstressed state
       // nullSigma, which in the scope of this material means application of
       // an isotropic cohesion of 'a'. The internal plastic variable value and 
       // plastic work should automatically be evaluated to meaningful values.
       // Note that the evaluated plastic work MUST equal that one computed above. (alphaN)
       this->ApplyLoad(nullSigma, epsA /* initial total strain */);
	   fInitialEps = LADEKIMPARENT::GetState();
    }
    virtual void SetUp(const TPZTensor<REAL> & epsTotal)  override {
        LADEKIMPARENT::SetUp(epsTotal);
    }
	
	virtual void Print(std::ostream & out) const override
	{
		out << "\n" << this->Name();
		out << "\n Base Class Data:\n";
		LADEKIMPARENT::Print(out);
		out << "\nTPZLadeKim internal members:";
		out << "\n a*Pa = " << faPa;
		out << "\n InitialEps = " << fInitialEps;
		
	}

	public:
int ClassId() const override;

    void Write(TPZStream &buf, int withclassid) const override{
        LADEKIMPARENT::Write(buf, withclassid);

        buf.Write(&faPa, 1);
        buf.Write(&fInitialEps.m_eps_t[0], 6);
        buf.Write(&fInitialEps.m_eps_p[0], 6);
        buf.Write(&fInitialEps.m_hardening, 1);

        buf.Write(&fYC.fKsi1, 1);
        buf.Write(&fYC.fh, 1);
        buf.Write(&fYC.m_hardening, 1);
        buf.Write(&fYC.fKsi2, 1);
        buf.Write(&fYC.fMu, 1);

        fInitialEps.Write(buf, withclassid);
    }

    void Read(TPZStream& buf, void* context) override {
        LADEKIMPARENT::Read(buf, context);

        buf.Read(&faPa, 1);
        buf.Read(&fInitialEps.m_eps_t[0], 6);
        buf.Read(&fInitialEps.m_eps_p[0], 6);
        buf.Read(&fInitialEps.m_hardening, 1);

        buf.Read(&fYC.fKsi1, 1);
        buf.Read(&fYC.fh, 1);
        buf.Read(&fYC.m_hardening, 1);
        buf.Read(&fYC.fKsi2, 1);
        buf.Read(&fYC.fMu, 1);

        fInitialEps.Read(buf, context);
    }
    
    /**
    Set the plastic state variables
    */
	virtual void SetState(const TPZPlasticState<REAL> &state) override
	{
		TPZPlasticState<REAL> temp(state);
		temp += fInitialEps;
		LADEKIMPARENT::SetState(temp);
	}
	
    /**
    Retrieve the plastic state variables
    */
    virtual TPZPlasticState<REAL> GetState () const override
    {
		TPZPlasticState<REAL> temp = LADEKIMPARENT::GetState();
		temp -= fInitialEps;
		return temp;
    }

    /**
    Computes the strain tensor as a function of the stress state.
    This function returns the inverse of function void Sigma(...) using a Newton's scheme.
    It first performs the Stress translation to remove the effect of apparent cohesion,
    making it possible to use the cohesionless Lade Kim model simulating the cohesion effect.
    @param [in] sigma stress tensor
    @param [out] epsTotal deformation tensor
    */
    virtual void ApplyLoad(const TPZTensor<REAL> & sigma, TPZTensor<REAL> &epsTotal) override
    {
       // Deformation translation from the cohesive material to the equivalent cohesionless
       epsTotal.Add(fInitialEps.m_eps_t, +1.);
       TPZTensor<REAL> I, cohesionlessSigma(sigma);
       I.Identity();
       // Stress translation from the cohesive to the equivalent cohesionless material
       cohesionlessSigma.Add(I, faPa);
       LADEKIMPARENT::ApplyLoad(cohesionlessSigma, epsTotal);
       // Deformation translation from the equivalent cohesionless to the cohesive material
       epsTotal.Add(fInitialEps.m_eps_t, -1.);
    }

    /**
    * Load the converged solution, updating the damage variables
    */
    virtual void ApplyStrain(const TPZTensor<REAL> &epsTotal) override
    {
        TPZTensor<REAL> translatedEpsTotal(epsTotal);
        // Deformation translation from the cohesive to the equivalent cohesionless material
        translatedEpsTotal.Add(fInitialEps.m_eps_t, 1.);
        LADEKIMPARENT::ApplyStrain(translatedEpsTotal);
    }

    /**
    * Load the converged solution, updating the damage variables
    */
    virtual void ApplyStrainComputeDep(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma, TPZFMatrix<REAL> &Dep)  override
    {
       TPZTensor<REAL> translatedEpsTotal(epsTotal);
       // Deformation translation from the cohesive to the equivalent cohesionless material
       translatedEpsTotal.Add(fInitialEps.m_eps_t, 1.);
       LADEKIMPARENT::ApplyStrainComputeDep(translatedEpsTotal, sigma, Dep);

       TPZTensor<REAL> I;
       I.Identity();
       // Stress translation from the equivalent cohesionless to the cohesive material
       sigma.Add(I, - faPa);
    }
	
    virtual void ApplyStrainComputeSigma(const TPZTensor<REAL> &epsTotal, TPZTensor<REAL> &sigma,  TPZFMatrix<REAL> * tangent = NULL) override
    {

        bool require_tangent_Q = true;
        if (!tangent) {
            require_tangent_Q = false;
        }
        
#ifdef PZDEBUG
        // Check for required dimensions of tangent
        if (!(tangent->Rows() == 6 && tangent->Cols() == 6)) {
            std::cerr << "Unable to compute the tangent operator. Required tangent array dimensions are 6x6." << std::endl;
            DebugStop();
        }
#endif
        
        if (require_tangent_Q) {
            DebugStop(); // implemented this functionality.
        }
        
       TPZTensor<REAL> translatedEpsTotal(epsTotal);
       // Deformation translation from the cohesive to the equivalent cohesionless material
       translatedEpsTotal.Add(fInitialEps.m_eps_t, 1.);
       LADEKIMPARENT::ApplyStrainComputeSigma(translatedEpsTotal, sigma);

       TPZTensor<REAL> I;
       I.Identity();
       // Stress translation from the equivalent cohesionless to the cohesive material
       sigma.Add(I, - faPa);
    }

    /**
    return the value of the yield functions for the given deformation
     * @param epsTotal [in] deformation tensor (total deformation
     * @param phi [out] vector of yield functions
    */
    virtual void Phi(const TPZTensor<REAL> &epsTotal, TPZVec<REAL> &phi) const override
    {
        TPZTensor<REAL> translatedEpsTotal(epsTotal);
        // Deformation translation from the cohesive to the equivalent cohesionless material
        translatedEpsTotal.Add(fInitialEps.m_eps_t, 1.);
        LADEKIMPARENT::Phi(translatedEpsTotal, phi);
    }

private:

  
    /**
    * variable related to the stress translation to simulate material cohesion
    * Already contains dimensional factor
    */
    REAL faPa;

    /**
    * variable to store the plastic state related to the unstressed state
    * in a framework of a cohesion material.
    * Since LadeKim models the cohesion material with the stress axis translated
    * from a value of a, it is important to keep the related elastic deformation
    * in this variable to subtract it from the epstotal answers. It is done
    * because, intuitively, the initially nonstressed state should relate to an
    * initially undeformed state.
    */
    TPZPlasticState<REAL> fInitialEps;

public:


// The following static members load test data from article 
// Lade, Paul V. Kim, Moon K. Single Hardening Constitutove Model for Soil, Rock and Concrete. 
// Int. Journal of Solid Structures, vol.32, No14. pp 1963-1978. Elsevier Science, 1994

    // Plain Concrete
    static void PlainConcrete(TPZLadeKim & material)
    {
	REAL poisson = 0.18;
	REAL M       = 361800.;
	REAL lambda  = 0.;
	REAL a       = 28.5;
	REAL m       = 1.113;
	REAL neta1   = 159800.;
	REAL ksi2    = -2.92;
	REAL mu      = 5.06;
	REAL C       = 0.712E-12;
	REAL p       = 3.8;
	REAL h       = 1.990;
	REAL alpha   = 0.75;
	REAL pa      = 14.7;
	
	material.fResTol = 1.e-8;
		
	material.SetUp(poisson, M, lambda,
			a, m, neta1,
			ksi2, mu,
			C, p,
			h, alpha,
			pa);
	
    }
    
    static void PlainConcreteMPa(TPZLadeKim & material)
    {
        REAL poisson = 0.18;
        REAL M       = 361800.*0.0068948;
        REAL lambda  = 0.;
        REAL a       = 28.5;//dimensionless
        REAL m       = 1.113;//dimensionless
        REAL neta1   = 159800.;//dimensionless
        REAL ksi2    = -2.92;//dimensionless
        REAL mu      = 5.06;//dimensionless
        REAL C       = 0.712E-12;
        REAL p       = 3.8;
        REAL h       = 1.990;
        REAL alpha   = 0.75;
        REAL pa      = 14.7*0.0068948;
        
        material.fResTol = 1.e-8;
		
        material.SetUp(poisson, M, lambda,
                       a, m, neta1,
                       ksi2, mu,
                       C, p,
                       h, alpha,
                       pa);
        
    }


    // Loose Sacramento River Sand
    static void LooseSacrRiverSand(TPZLadeKim & material)
    {
	REAL poisson = 0.2;
	REAL M       = 500.;
	REAL lambda  = 0.28;
	REAL a       = 0.;
	REAL m       = 0.093;
	REAL neta1   = 28.;
	REAL ksi2    = -3.72;
	REAL mu      = 2.36;
	REAL C       = 0.127E-3;
	REAL p       = 1.65;
	REAL h       = 0.534;
	REAL alpha   = 0.794;
	REAL pa      = 14.7;
	
	material.fResTol = 1.e-8;
		
	material.SetUp(poisson, M, lambda,
			a, m, neta1,
			ksi2, mu,
			C, p,
			h, alpha,
			pa);
		
    }

    // Dense Sacramento River Sand
    static void DenseSacrRiverSand(TPZLadeKim & material)
    {
	REAL poisson = 0.2;
	REAL M       = 900.;
	REAL lambda  = 0.28;
	REAL a       = 0.;
	REAL m       = 0.23;
	REAL neta1   = 80.;
	REAL ksi2    = -3.09;
	REAL mu      = 2.00;
	REAL C       = 0.396E-4;
	REAL p       = 1.82;
	REAL h       = 0.765;
	REAL alpha   = 0.229;
	REAL pa      = 14.7;
		
	material.fResTol = 1.e-8;
	
	material.SetUp(poisson, M, lambda,
			a, m, neta1,
			ksi2, mu,
			C, p,
			h, alpha,
			pa);
		
    }

    // Fine Silica Sand
    static void FineSilicaSand(TPZLadeKim & material)
    {
	REAL poisson = 0.27;
	REAL M       = 440.;
	REAL lambda  = 0.22;
	REAL a       = 0.;
	REAL m       = 0.1;
	REAL neta1   = 24.7;
	REAL ksi2    = -3.69;
	REAL mu      = 2.26;
	REAL C       = 0.324E-3;
	REAL p       = 1.25;
	REAL h       = 0.355;
	REAL alpha   = 0.515;
	REAL pa      = 14.7;
	
	material.fResTol = 1.e-8;
		
	material.SetUp(poisson, M, lambda,
			a, m, neta1,
			ksi2, mu,
			C, p,
			h, alpha,
			pa);
		
    }
    
    // Fine Silica Sand
    static void FineSilicaSandPaperIII(TPZLadeKim & material)
    {
        REAL poisson = 0.2;
        REAL M       = 440.;
        REAL lambda  = 0.22;
        REAL a       = 0.;
        REAL m       = 0.1;
        REAL neta1   = 24.7;
        REAL ksi2    = -3.69;
        REAL mu      = 2.26;
        REAL C       = 0.324E-3;
        REAL p       = 1.25;
        REAL h       = 0.355;
        REAL alpha   = 0.515;
        REAL pa      = 14.7;
        
        material.fResTol = 1.e-8;
		
        material.SetUp(poisson, M, lambda,
                       a, m, neta1,
                       ksi2, mu,
                       C, p,
                       h, alpha,
                       pa);
		
    }
    
    // Loose Santa Monica Beach Sand
    static void LooseSantaMonicaBeachSand(TPZLadeKim & material)
    {
        REAL poisson = 0.26;
        REAL M       = 600.;
        REAL lambda  = 0.27;
        REAL a       = 0.;
        REAL m       = 0.107;
        REAL neta1   = 32.6;
        REAL ksi2    = -3.65;
        REAL mu      = 2.10;
        REAL C       = 0.000204;
        REAL p       = 1.51;
        REAL h       = 0.60;
        REAL alpha   = 0.79;
        REAL pa      = 14.7;
        
        material.fResTol = 1.e-8;
		
        material.SetUp(poisson, M, lambda,
                       a, m, neta1,
                       ksi2, mu,
                       C, p,
                       h, alpha,
                       pa);
		
    }

public:


//////////////////CheckConv related methods/////////////////////

    /**
    number of types of residuals
    */
    inline int NumCases();

    static TPZManVector<REAL,6+6+1+TPZYCLadeKim::NYield> gRefResidual;
    /**
    LoadState will keep a given state as static variable of the class
    */
    inline void LoadState(TPZFMatrix<REAL> &state);

    inline void ComputeTangent(TPZFMatrix<REAL> &tangent, TPZVec<REAL> &, int icase);

    inline void Residual(TPZFMatrix<REAL> &res,int icase);



    static void CheckConv()
    {
	const int nVars = 6+6+1+TPZYCLadeKim::NYield;
        gRefResidual.Resize(nVars);
	
	// Creating the LadeKim obejct
	TPZLadeKim LadeKim;
	// setup 
	
	TPZLadeKim::PlainConcrete(LadeKim);
	
	TPZFNMatrix<nVars> input(nVars,1), Range(nVars,1);
	//plastic strains
	input(_XX_,0) = 0.0017;
	input(_YY_,0) = 0.0013;
	input(_ZZ_,0) = 0.0011;
	input(_XY_,0) = 0.0007 ;
	input(_XZ_,0) = 0.0005 ;
	input(_ZZ_,0) = 0.0003 ;
	
	//alpha
	input(6,0) = 0.00019;
	
	//delgamma
	input(7,0) = 0.000023;
	
	//total strains
	input(8+_XX_,0) = 0.029;
	input(8+_YY_,0) = 0.031;
	input(8+_ZZ_,0) = 0.037;
	input(8+_XY_,0) = 0.0041;
	input(8+_XZ_,0) = 0.0047;
	input(8+_YZ_,0) = 0.0049;
	
	Range = input * (1./19.);
	TPZVec< REAL > Coefs(1,1.);
	CheckConvergence(LadeKim, input, Range, Coefs);
	
    }

    virtual int GetNYield() const {
        return as_integer(NYield);
    }
};

//////////////////CheckConv related methods/////////////////////



inline int TPZLadeKim::NumCases()
{
    return 1;
}

inline void TPZLadeKim::LoadState(TPZFMatrix<REAL> &state)
{
  int i;
  const int nVars = 6+6+1+TPZYCLadeKim::NYield;
  for(i=0; i<nVars; i++) gRefResidual[i] = state(i,0);

//#ifdef PZ_LOG
//  {
//    std::stringstream sout;
//    sout << "Tension " << state;
//    LOGPZ_DEBUG(loggerLadeKim,sout.str().c_str());
//  }
//#endif
}

inline void TPZLadeKim::ComputeTangent(TPZFMatrix<REAL> &tangent, TPZVec<REAL> &, int icase)
{
  const int nyield = TPZYCLadeKim::NYield;
  const int nVars = 7+nyield+6; // also includes epsTotal
  typedef TFad<nVars,REAL> TFAD;
  int i, j;

  TPZPlasticState<TFAD> Np1_FAD;
  TPZManVector<TFAD, nVars-6> epsRes_FAD(nVars-6);
  TPZManVector<TFAD, nyield > delGamma_FAD(nyield);
  REAL normEpsPErr;

  switch(icase)
  {
    case 0:
      {
	for(i=0;i<6;i++)
	{
		Np1_FAD.m_eps_p[i] = gRefResidual[i];
		Np1_FAD.m_eps_p[i].diff(i,nVars);
	}
	Np1_FAD.m_hardening = gRefResidual[i];
	Np1_FAD.m_hardening.diff(i++,nVars);
	for(j=0;j<nyield;j++)
	{
		delGamma_FAD[j] = TFAD(gRefResidual[i]);
		delGamma_FAD[j].diff(i++,nVars);
	}
	for(j=0;j<6;j++)
	{
		Np1_FAD.m_eps_t[j] = gRefResidual[i];
		Np1_FAD.m_eps_t[j].diff(i++,nVars);
	}

	
	//Compute Residual
    PlasticResidual<REAL, TFAD>
		       (fN, Np1_FAD, delGamma_FAD, epsRes_FAD, normEpsPErr);   
		  
	tangent.Redim(nVars-6,nVars);
	
	for(i=0; i<nVars-6; i++)
		for(j=0;j<nVars;j++)
		tangent(i,j) = epsRes_FAD[i].dx(j);
      }
      break;
  }
//#ifdef PZ_LOG
//  std::stringstream sout;
//  sout << "Matriz tangent " << tangent;
//  LOGPZ_DEBUG(loggerLadeKim,sout.str().c_str());
//#endif
}

inline void TPZLadeKim::Residual(TPZFMatrix<REAL> &res,int icase)
{
  const int nyield = TPZYCLadeKim::NYield;
  const int nVars = 7+nyield+6; // also includes epsTotal
  int i, j;

  TPZPlasticState<REAL> Np1;
  TPZManVector<REAL, nVars-6> epsRes(nVars-6);
  TPZManVector<REAL, nyield > delGamma(nyield);
  REAL normEpsPErr;

  switch(icase)
  {
    case 0:
      {
	for(i=0;i<6;i++)     Np1.m_eps_p[i] = gRefResidual[i];
	Np1.m_hardening = gRefResidual[i++];
	for(j=0;j<nyield;j++)delGamma[j] = gRefResidual[i++];
	for(j=0;j<6;j++)     Np1.m_eps_t[j] = gRefResidual[i++];
	
	
	//Compute Residual
    PlasticResidual<REAL>(fN, Np1, delGamma, epsRes, normEpsPErr);   

	res.Redim(nVars-6,1);

	for(i=0; i<nVars-6; i++)
		res(i,0) = epsRes[i];
      }
      break;
  }
//#ifdef PZ_LOG
//  std::stringstream sout;
//  sout << "Residual vector " << res;
//  LOGPZ_DEBUG(loggerLadeKim,sout.str().c_str());
//#endif
}




//////////////////CheckConv related methods/////////////////////

#endif //TPZLADEKIM_H
