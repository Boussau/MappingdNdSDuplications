/*
 *  mapDnDsDuplicationTree.h
 *  
 *
 *  Created by boussau on 12/05/11.
 *  Copyright 2011 UC Berkeley. All rights reserved.
 *
 */


/*
 Copyright or © or Copr. CNRS
 
 This software is a computer program whose purpose is to describe
 the patterns of substitutions along a phylogeny using substitution mapping.
 
 This software is governed by the CeCILL  license under French law and
 abiding by the rules of distribution of free software.  You can  use, 
 modify and/ or redistribute the software under the terms of the CeCILL
 license as circulated by CEA, CNRS and INRIA at the following URL
 "http://www.cecill.info". 
 
 As a counterpart to the access to the source code and  rights to copy,
 modify and redistribute granted by the license, users are provided only
 with a limited warranty  and the software's author,  the holder of the
 economic rights,  and the successive licensors  have only  limited
 liability. 
 
 In this respect, the user's attention is drawn to the risks associated
 with loading,  using,  modifying and/or developing or reproducing the
 software by the user in light of its specific status of free software,
 that may mean  that it is complicated to manipulate,  and  that  also
 therefore means  that it is reserved for developers  and  experienced
 professionals having in-depth computer knowledge. Users are therefore
 encouraged to load and test the software's suitability as regards their
 requirements in conditions enabling the security of their systems and/or 
 data to be ensured and,  more generally, to use and operate it in the 
 same conditions as regards security. 
 
 The fact that you are presently reading this means that you have had
 knowledge of the CeCILL license and that you accept its terms.
 */

//#include "ChiClustering.h"
//#include "MultinomialClustering.h"

// From the STL:
#include <iostream>
#include <iomanip>
#include <map>
#include <queue>
#include <cmath>


using namespace std;

// From bpp-seq:
#include <Bpp/Seq/Alphabet.all>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From bpp-phyl:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Likelihood.all>
#include <Bpp/Phyl/Mapping.all>
#include <Bpp/Phyl/Simulation.all>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
//#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Io/Nhx.h>
#include <Bpp/Phyl/Model/Nucleotide/JCnuc.h>
#include <Bpp/Phyl/Model/Nucleotide/HKY85.h>
#include <Bpp/Phyl/Model/Protein/JCprot.h>
//#include <Bpp/Phyl/Model/CodonNeutralReversibleSubstitutionModel.h>
#include <Bpp/Phyl/Model/Codon/CodonSubstitutionModel.h>

#include <Bpp/Phyl/Model/Codon/YN98.h>
#include <Bpp/Phyl/Model/SubstitutionModelSet.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>

// From bpp-core:
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Stat/ContingencyTableTest.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>

using namespace bpp;
