// Copyright 2010, 2011, 2012 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2010-01-05

#include "parametricTNLP.hpp"

#include "IpIpoptApplication.hpp"
#include "SensApplication.hpp"
#include "IpPDSearchDirCalc.hpp"
#include "IpIpoptAlg.hpp"
#include "SensRegOp.hpp"
#include "IpOrigIpoptNLP.hpp"
#include "IpDenseVector.hpp"
#include "SensPerturbInterface.hpp"


int main(int argv, char**argc)
{
  using namespace Ipopt;

  SmartPtr<IpoptApplication> app_ipopt = new IpoptApplication();

  //SmartPtr<SensApplication> app_sens = new SensApplication(app_ipopt->Jnlst(),
  //app_ipopt->Options(),
  //app_ipopt->RegOptions());

  // Register sIPOPT options
  //RegisterOptions_sIPOPT(app_ipopt->RegOptions());
  //app_ipopt->Options()->SetRegisteredOptions(app_ipopt->RegOptions());

  // Call Initialize the first time to create a journalist, but ignore
  // any options file
  ApplicationReturnStatus retval;
  retval = app_ipopt->Initialize("");
  if (retval != Solve_Succeeded) {
    //printf("ampl_ipopt.cpp: Error in first Initialize!!!!\n");
    exit(-100);
  }
  app_ipopt->Initialize();

  // create AmplSensTNLP from argc. This is an nlp because we are using our own TNLP Adapter
  SmartPtr<ParaTNLP> sens_tnlp = new ParametricTNLP();

  retval = app_ipopt->OptimizeTNLP(sens_tnlp);

  SensPerturbInterface sensClass (app_ipopt);

  // Set up the perturbed parameters

  Number* dp_ptr = new Number[2];
  dp_ptr[0] = -0.5;
  dp_ptr[1] = 0.0;
  SmartPtr<IteratesVector> rhs = sensClass.getParaSensMatrix(dp_ptr);

  // Get the (factorized) KKT matrix
  SmartPtr<PDSystemSolver> pd_solver_ = sensClass.getKKTmatrix();
  SmartPtr<IteratesVector> lhs = rhs->MakeNewIteratesVector();
  SmartPtr<const IteratesVector> curr = app_ipopt->IpoptDataObject()->curr();

  sensClass.solveSens(rhs, pd_solver_, curr);
}
