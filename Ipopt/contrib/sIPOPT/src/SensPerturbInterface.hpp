// Copyright 2012 Hans Pirnay, Arne Graf
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2012-04-05

#ifndef __SENSPERTURBINTERFACE_HPP__
#define __SENSPERTURBINTERFACE_HPP__

#include "IpReferenced.hpp"
#include "parametricTNLP.hpp"

#include "IpIpoptApplication.hpp"
#include "SensApplication.hpp"
#include "IpPDSearchDirCalc.hpp"
#include "IpIpoptAlg.hpp"
#include "SensRegOp.hpp"
#include "IpOrigIpoptNLP.hpp"
#include "IpDenseVector.hpp"

namespace Ipopt
{

  class SensPerturbInterface : public ReferencedObject
  {
    /** This is the interface for the  */

  public:
    typedef TNLP::LinearityType LinearityType;
    typedef TNLP::IndexStyleEnum IndexStyleEnum;
    /**@name Constructors/Destructors */
    //@{
    SensPerturbInterface()
    {}

    SensPerturbInterface(SmartPtr<IpoptApplication> app_ipopt)
    :
    app_ipopt_(app_ipopt)
    {}

    const SmartPtr<IteratesVector> getParaSensMatrix (const Number* dp_ptrIn)
    {
      SmartPtr<const IteratesVector> curr = app_ipopt_->IpoptDataObject()->curr();
      SmartPtr<const Vector> x = curr->x();
      SmartPtr<const Vector> y_c = curr->y_c();
      SmartPtr<const Vector> y_d = curr->y_d();

      // Get the parameter sensitivity matrices
      SmartPtr<IpoptNLP> ipopt_nlp = app_ipopt_->IpoptNLPObject();
      SmartPtr<const Matrix> opt_jac_c_p = (dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp)))->jac_c_p(*x);
      SmartPtr<const Matrix> opt_jac_d_p = (dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp)))->jac_d_p(*x);
      SmartPtr<const Matrix> opt_h_p     = (dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ipopt_nlp)))->h_p(*x, 1.0, *y_c, *y_d);
      opt_jac_c_p->Print(*app_ipopt_->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "opt_jac_c_p");
      opt_jac_d_p->Print(*app_ipopt_->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "opt_jac_d_p");
          opt_h_p->Print(*app_ipopt_->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "opt_h_p");

      // Set up the perturbed parameters
      SmartPtr<const DenseVector> p0 = dynamic_cast<const DenseVector*>(GetRawPtr(ipopt_nlp->p()));
      SmartPtr<DenseVector>       dp = dynamic_cast<DenseVector*>(p0->MakeNewCopy());
      dp->SetValues(dp_ptrIn); //changed, plz check

      // Set up RHS and LHS for solve
      SmartPtr<IteratesVector> rhs = app_ipopt_->IpoptDataObject()->curr()->MakeNewIteratesVector();
      rhs->Set(0.0);
      SmartPtr<Vector> rhs_x = x->MakeNew();
      opt_h_p->MultVector(1.0, *dp, 0.0, *rhs_x);
      rhs->Set_x_NonConst(*rhs_x);
      SmartPtr<Vector> rhs_c = y_c->MakeNew();
      opt_jac_c_p->MultVector(1.0, *dp, 0.0, *rhs_c);
      rhs->Set_y_c_NonConst(*rhs_c);
      SmartPtr<Vector> rhs_d = y_d->MakeNew();
      opt_jac_d_p->MultVector(1.0, *dp, 0.0, *rhs_d);
      rhs->Set_y_d_NonConst(*rhs_d);

      return rhs;
    }

    const SmartPtr<PDSystemSolver> getKKTmatrix ()
    {
      // Get the (factorized) KKT matrix
      SmartPtr<IpoptAlgorithm> alg = app_ipopt_->AlgorithmObject();
      SmartPtr<PDSearchDirCalculator> pd_search;
      pd_search = dynamic_cast<PDSearchDirCalculator*>(GetRawPtr(alg->SearchDirCalc()));
      return pd_search->PDSolver();
    }

    void solveSens (const SmartPtr<IteratesVector> rhs,
                    const SmartPtr<PDSystemSolver> pd_solver_,
                    const SmartPtr<const IteratesVector> curr)
    {
      SmartPtr<IteratesVector> lhs = rhs->MakeNewIteratesVector();

      pd_solver_->Solve(-1.0, 0.0, *rhs, *lhs, false, false);
      lhs->Axpy(1.0, *curr);
      lhs->Print(*app_ipopt_->Jnlst(), J_INSUPPRESSIBLE, J_DBG, "perturbed_x");
    }




    /** Default destructor */
    virtual ~SensPerturbInterface()
    {}
    //@}

  private:
    SmartPtr<IpoptApplication> app_ipopt_;

    /**@name Default Compiler Generated Methods
    * (Hidden to avoid implicit creation/calling).
    * These methods are not implemented and
    * we do not want the compiler to implement
    * them for us, so we declare them private
    * and do not define them. This ensures that
    * they will not be implicitly created/called. */
    //@{
    /** Copy Constructor */
    SensePerturbInterface(const SensPerturbInterface&);

    /** Overloaded Equals Operator */
    void operator=(const SensPerturbInterface&);
    //@}

  };
}

#endif
