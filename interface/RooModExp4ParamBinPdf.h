//---------------------------------------------------------------------------
#ifndef HiggsAnalysis_CombinedLimit_RooModExp4ParamBinPdf_h
#define HiggsAnalysis_CombinedLimit_RooModExp4ParamBinPdf_h
//---------------------------------------------------------------------------
#include "RooAbsPdf.h"
#include "RooConstVar.h"
#include "RooRealProxy.h"
//---------------------------------------------------------------------------
class RooRealVar;
class RooAbsReal;

#include "Riostream.h"
#include "TMath.h"
#include <TH1.h>
#include "Math/SpecFuncMathCore.h"
#include "Math/SpecFuncMathMore.h"
#include "Math/Functor.h"
#include "Math/WrappedFunction.h"
#include "Math/IFunction.h"
#include "Math/Integrator.h"

//---------------------------------------------------------------------------
class RooModExp4ParamBinPdf : public RooAbsPdf
{
public:
   RooModExp4ParamBinPdf() {} ;
   RooModExp4ParamBinPdf(const char *name, const char *title,
		    RooAbsReal& _th1x, RooAbsReal& _p1,
		  RooAbsReal& _p2, RooAbsReal& _p3,
		  RooAbsReal& _sqrts, RooAbsReal& _meff, RooAbsReal& _seff);
   RooModExp4ParamBinPdf(const char *name, const char *title,
		    RooAbsReal& _th1x, RooAbsReal& _p1,
		  RooAbsReal& _p2, RooAbsReal& _p3,
		  RooAbsReal& _sqrts);
   RooModExp4ParamBinPdf(const RooModExp4ParamBinPdf& other,
      const char* name = 0);
   void setTH1Binning(TH1* _Hnominal);
   void setAbsTol(double _absTol);
   void setRelTol(double _relTol);
   virtual TObject* clone(const char* newname) const { return new RooModExp4ParamBinPdf(*this,newname); }
   inline virtual ~RooModExp4ParamBinPdf() { }

   Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const;
   Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

protected:   

   RooRealProxy th1x;        // dependent variable
   RooRealProxy p1;       // p1
   RooRealProxy p2;        // p2
   RooRealProxy p3;        // p3
   RooRealProxy sqrts;        // sqrts
   RooRealProxy meff;        // meff
   RooRealProxy seff;        // seff
   Int_t xBins;        // X bins
   Double_t xArray[2000]; // xArray[xBins+1]
   Double_t xMax;        // X max
   Double_t xMin;        // X min
   Double_t relTol;      //relative tolerance for numerical integration
   Double_t absTol;      //absolute tolerance for numerical integration

   Double_t evaluate() const;
private:
   ClassDef(RooModExp4ParamBinPdf,1) // RazorModExp4ParamBinPdf function
    
};
//---------------------------------------------------------------------------
#endif

#include "Math/IFunction.h"
#include "Math/IParamFunction.h"
 
class ModExp4ParamFunction: public ROOT::Math::IParametricFunctionOneDim
{
private:
   const double *pars;
 
public:
   double DoEvalPar(double x,const double* p) const
   {
     double pdf = exp(p[1]*pow(x/p[0],p[2]) + p[1]*pow(1.0-x/p[0],p[3]));
     double eff = 1.;
     if (p[4]>0 && p[5]>0) eff = 0.5 * (1.0 + TMath::Erf((x - p[4])/p[5])) ;
     return pdf*eff;
   }
   
   double DoEval(double x) const
   {
     double pdf = exp(pars[1]*pow(x/pars[0],pars[2]) + pars[1]*pow(1.0-x/pars[0],pars[3]));
     double eff = 1.;
     if (pars[4]>0 && pars[5]>0) eff = 0.5 * (1.0 + TMath::Erf((x - pars[4])/pars[5]));
     return pdf*eff;
   }
 
   ROOT::Math::IBaseFunctionOneDim* Clone() const
   {
      return new ModExp4ParamFunction();
   }
 
   const double* Parameters() const
   {
      return pars;
   }
 
   void SetParameters(const double* p)
   {
      pars = p;
   }
 
   unsigned int NPar() const
   {
      return 6;
   }
};
