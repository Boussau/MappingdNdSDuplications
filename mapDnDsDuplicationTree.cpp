/*
 *  mapDnDsDuplicationTree.cpp
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


#include "mapDnDsDuplicationTree.h"
#include <Bpp/Phyl/Io/BppOFrequenciesSetFormat.h>
#include <Bpp/Phyl/Model/FrequenciesSet/CodonFrequenciesSet.h>

/*************************** COMPILATION

g++  -pipe -o mapDnDsDuplicationTree mapDnDsDuplicationTree.cpp mapDnDsDuplicationTree.h -I/usr/local/include  -L. -L/usr/local/lib  -g -Wall -fopenmp -std=c++0x -lbpp-core -lbpp-seq -lbpp-phyl 
**************************/

/******************************************************************************/

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "mapDnDsDuplicationTree parameter1_name=parameter1_value parameter2_name=parameter2_value").endLine();
  (*ApplicationTools::message << "      ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the package manual for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

/******************************************************************************/

std::map <int, std::vector <int> > breadthFirstreNumber (TreeTemplate<Node> & tree) {
  int index = 0;
  std::map<Node *, int> color ;
  std::map <int, std::vector <int> > DepthToIds; //A std::map where we store the correspondence between the depth of a node (number of branches between the root and the node) and the node id.
  std::map <int, int > IdsToDepths;
  std::vector <Node * > nodes = tree.getNodes();
  //All nodes white
  for (size_t i = 0; i< nodes.size() ; i++) {
    color.insert(std::pair <Node *,int>(nodes[i],0));
  }
  std::queue <Node *> toDo;
  toDo.push(tree.getRootNode());
  color[tree.getRootNode()] = 1;
  tree.getRootNode()->setId(index);
  std::vector <int> v;
  DepthToIds.insert(std::pair <int, std::vector<int> > (0,v));
  DepthToIds[0].push_back(index);
  IdsToDepths[index] = 0;
  index++;
  Node * u;
  while(!toDo.empty()) {
    u = toDo.front();
    toDo.pop();
    int fatherDepth = IdsToDepths[u->getId()];
    std::vector <Node *> sons;
    for (size_t j = 0 ; j< u->getNumberOfSons() ; j++) {
      sons.push_back(u->getSon(j));
    }
    for (size_t j = 0; j< sons.size() ; j++) {
      if (color[sons[j]]==0) {
        color[sons[j]]=1;
        sons[j]->setId(index);
        if (DepthToIds.count(fatherDepth+1)==0) {
          DepthToIds.insert(std::pair <int, std::vector<int> > (fatherDepth+1,v));
        }
        DepthToIds[fatherDepth+1].push_back(index); 
        IdsToDepths[index] = fatherDepth+1;
        index++;
        toDo.push(sons[j]);
      }
    }
    color[u]=2;
  }
  return DepthToIds;
}


/******************************************************************************/


vector< vector<unsigned int> > getCountsPerBranch(
                                                  DRTreeLikelihood& drtl,
                                                  const vector<int>& ids,
                                                  SubstitutionModel* model,
                                                  const SubstitutionRegister& reg,
                                                  bool stationarity = true,
                                                  double threshold = -1)
{
  auto_ptr<SubstitutionCount> count(new UniformizationSubstitutionCount(model, reg.clone()));
  // auto_ptr<SubstitutionCount> count(new DecompositionSubstitutionCount(model, reg.clone()));
  
  //SubstitutionCount* count = new SimpleSubstitutionCount(reg);
  const CategorySubstitutionRegister* creg = 0;
  if (!stationarity) {
    try {
      creg = &dynamic_cast<const CategorySubstitutionRegister&>(reg);
    } catch (Exception& ex) {
      throw Exception("The stationarity option can only be used with a category substitution register.");
    }
  }
 // SubstitutionMappingTools::computeSubstitutionVectors(drtl, ids, *count, false);
  auto_ptr<ProbabilisticSubstitutionMapping> mapping(SubstitutionMappingTools::computeSubstitutionVectors(drtl, ids, *count, false));
  
  vector< vector<unsigned int> > counts(ids.size());
  
  unsigned int nbSites = mapping->getNumberOfSites();
  unsigned int nbTypes = mapping->getNumberOfSubstitutionTypes();
  
  for (size_t k = 0; k < ids.size(); ++k) {
    //vector<double> countsf = SubstitutionMappingTools::computeSumForBranch(*mapping, mapping->getNodeIndex(ids[i]));
    vector<double> countsf(nbTypes, 0);
    vector<double> tmp(nbTypes, 0);
    unsigned int nbIgnored = 0;
    bool error = false;
    for (unsigned int i = 0; !error && i < nbSites; ++i) {
      double s = 0;
      for (unsigned int t = 0; t < nbTypes; ++t) {
        tmp[t] = (*mapping)(k, i, t);
        error = std::isnan((double)tmp[t]);
        if (error)
          goto ERROR;
        s += tmp[t];
      }
      if (threshold >= 0) {
        if (s <= threshold)
          countsf += tmp;
        else {
          nbIgnored++;
        }
      } else {
        countsf += tmp;
      }
    }
    
  ERROR:
    if (error) {
      //We do nothing. This happens for small branches.
      ApplicationTools::displayWarning("On branch " + TextTools::toString(ids[k]) + ", counts could not be computed.");
      for (unsigned int t = 0; t < nbTypes; ++t)
        countsf[t] = 0;
    } else {
      if (nbIgnored > 0) {
        ApplicationTools::displayWarning("On branch " + TextTools::toString(ids[k]) + ", " + TextTools::toString(nbIgnored) + " sites (" + TextTools::toString(ceil(static_cast<double>(nbIgnored * 100) / static_cast<double>(nbSites))) + "%) have been ignored because they are presumably saturated.");
      }
      
      if (!stationarity) {
        vector<double> freqs = DRTreeLikelihoodTools::getPosteriorStateFrequencies(drtl, ids[k]);
        //Compute frequencies for types:
        vector<double> freqsTypes(creg->getNumberOfCategories());
        for (size_t i = 0; i < freqs.size(); ++i) {
          unsigned int c = creg->getCategory(static_cast<int>(i));
          freqsTypes[c - 1] += freqs[i];
        }
        //We divide the counts by the frequencies and rescale:
        double s = VectorTools::sum(countsf);
        for (unsigned int t = 0; t < nbTypes; ++t) {
          countsf[t] /= freqsTypes[creg->getCategoryFrom(t + 1) - 1];
        }
        double s2 = VectorTools::sum(countsf);
        //Scale:
        (countsf / s2) * s;
      }
    }
    
    counts[k].resize(countsf.size());
    for (size_t j = 0; j < countsf.size(); ++j) {
      counts[k][j] = static_cast<unsigned int>(floor(countsf[j] + 0.5)); //Round counts
    }
  }
  return counts;
}



/******************************************************************************/


vector < map< int, vector<unsigned int> > > getCountsPerBranchPerSite(
                                                  DRTreeLikelihood& drtl,
                                                  const vector<int>& ids,
                                                  SubstitutionModel* model,
                                                  const SubstitutionRegister& reg,
                                                  bool stationarity = true,
                                                  double threshold = -1, 
						    const string familyName = "none", 
						    size_t nbTypes = 0
 								    )
{
  unique_ptr<SubstitutionCount> count(new UniformizationSubstitutionCount(model, reg.clone()));
  // auto_ptr<SubstitutionCount> count(new DecompositionSubstitutionCount(model, reg.clone()));
  
  //SubstitutionCount* count = new SimpleSubstitutionCount(reg);
 /* const CategorySubstitutionRegister* creg = 0;
  if (!stationarity) {
    try {
      creg = &dynamic_cast<const CategorySubstitutionRegister&>(reg);
    } catch (Exception& ex) {
      throw Exception("The stationarity option can only be used with a category substitution register.");
    }
  }*/
  //SubstitutionMappingTools::computeSubstitutionVectors(drtl, ids, *count, false);
  unique_ptr<ProbabilisticSubstitutionMapping> mapping(SubstitutionMappingTools::computeSubstitutionVectors(drtl, ids, *count, false));
  
  //vector on branches, map on sites, last vector on substitution types
  vector< map< int, vector<unsigned int> > > counts(ids.size()); 
  
  unsigned int nbSites = mapping->getNumberOfSites();
  nbTypes = mapping->getNumberOfSubstitutionTypes();
  
  for (size_t k = 0; k < ids.size(); ++k) {
    //vector<double> countsf = SubstitutionMappingTools::computeSumForBranch(*mapping, mapping->getNodeIndex(ids[i]));
    vector<double> tmp(nbTypes, 0);
    unsigned int nbIgnored = 0;
   // bool error = false;
    vector<double> countsf(nbTypes, 0);
//    for (unsigned int i = 0; !error && i < nbSites; ++i) {
    for (unsigned int i = 0; i < nbSites; ++i) {
      double s = 0;
      for (unsigned int t = 0; t < nbTypes; ++t) {
        countsf[t] = 0;
        tmp[t] = (*mapping)(k, i, t);
     /*   error = std::isnan((double)tmp[t]);
        if (error)
          goto ERROR;*/
        s += tmp[t];
      }
      if (threshold >= 0) {
        if (s <= threshold)
          countsf += tmp;
        else {
          nbIgnored++;
        }
      } else {
        countsf += tmp;
      }
      //Now we have countsf, which contains substitution counts for the current site.
      //We round the counts to the closest integer.
      for (unsigned int t = 0; t < nbTypes; ++t) {
          countsf[t] = static_cast<unsigned int>(floor(countsf[t] + 0.5)); //Round counts
        //Now we save these counts if they are positive.
        if (countsf[t] > 0) {
            counts[k][i].push_back(t) ;   
        }
      }
    }
    /*
  ERROR:
    if (error) {
      //We do nothing. This happens for small branches.
      ApplicationTools::displayWarning("On branch " + TextTools::toString(ids[k]) + ", some counts could not be computed.");
      for (unsigned int t = 0; t < nbTypes; ++t)
        countsf[t] = 0;
    } else {
      if (nbIgnored > 0) {
        ApplicationTools::displayWarning("On branch " + TextTools::toString(ids[k]) + ", " + TextTools::toString(nbIgnored) + " sites (" + TextTools::toString(ceil(static_cast<double>(nbIgnored * 100) / static_cast<double>(nbSites))) + "%) have been ignored because they are presumably saturated.");
      }
      
    }
    */
  }
  return counts;
}


/******************************************************************************/

vector < map< int, vector<unsigned int> > > getCountsPerBranchPerSiteASR(
                                                  DRTreeLikelihood& tl,
                                                  const vector<int>& ids,
                                                  VectorSiteContainer* sites,
                                                  const SubstitutionRegister& reg,
                                                  bool stationarity = true,
                                                  double threshold = -1, 
                                                  const string familyName = "none", 
                                                  size_t nbTypes = 0
                                                  )

{
  AncestralStateReconstruction *asr = new MarginalAncestralStateReconstruction(&tl);
  //bool  probMethod = true;
 
  TreeTemplate<Node> ttree ( tl.getTree() );
  vector<Node *> nodes = ttree.getNodes();
  size_t nbNodes = nodes.size();
  // Get the rate class with maximum posterior probability:
  vector<size_t> classes = tl.getRateClassWithMaxPostProbOfEachSite();
  // Get the posterior rate, i.e. rate averaged over all posterior probabilities:
  Vdouble rates = tl.getPosteriorRateOfEachSite();
  // Get the ancestral sequences:
  vector<Sequence*> sequences(nbNodes);
  for (size_t i = 0; i < nbNodes; i++) {
    Node *node = nodes[i];
      if (node->isLeaf()) {
        sequences[i] = sites->getSequence(node->getName() ).clone();
      } else {
        sequences[i] = asr->getAncestralSequenceForNode(node->getId() );
      }
  }

      //Now we have all the sequences at all nodes, including leaves
      //We need to compute counts
        //vector on branches, map on sites, last vector on substitution types
  vector< map< int, vector<unsigned int> > > counts(ids.size()); 
  
  unsigned int nbSites = sites->getNumberOfSites();

  nbTypes = reg.getNumberOfSubstitutionTypes();
  const Alphabet *alpha = tl.getAlphabet()->clone();
  /*  
    const SubstitutionModel *model = tl.getSubstitutionModel(0, 0)->clone();*/
  for (size_t k = 0; k < nbNodes; ++k) {
    if (nodes[k]->hasFather() ) {
        //std::cout << "HAHA " << nodes[k]->getId() /*- 1*/ <<std::endl;
      Sequence* seqChild = sequences.at( nodes[k]->getId() /*- 1*/ );

      Sequence* seqFather = sequences.at( nodes[k]->getFather()->getId() /*- 1*/ );

      vector<double> tmp(nbTypes, 0);
     // unsigned int nbIgnored = 0;
      //bool error = false;
      vector<double> countsf(nbTypes, 0.0);

      for (unsigned int i = 0; i < nbSites; ++i) {
        for (unsigned int t = 0; t < nbTypes; ++t) {
          countsf[t] = 0.0;
        }

        try {
        //Now we compare the parent and child sequences
          if ( seqFather->getValue(i) != seqChild->getValue(i) ) { //There is a substitution
          /*  std::cout << "Thre is a sub "<< seqFather->getValue(i) <<" "<< seqChild->getValue(i) << std::endl;
            std::cout << "Thre is a sub "<< seqFather->getChar(i) <<" "<< seqChild->getChar(i) << std::endl;
            std::cout << "Thre is a sub "<< alpha->charToInt( seqFather->getChar(i) ) <<" "<< alpha->charToInt( seqChild->getChar(i) ) << std::endl;
                      std::cout << "Thre is a sub "<< model->getAlphabetStateAsInt( seqFather->getValue(i) ) <<" "<< model->getAlphabetStateAsInt ( seqChild->getValue(i) ) << std::endl;*/

          // size_t type = reg.getType ( alpha->charToInt( seqFather->getChar(i) ), alpha->charToInt( seqChild->getChar(i) ) ) - 1 ;
          
          //std::cout << "seqFather->getValue(i): " << seqFather->getValue(i) << " "<< seqFather->getChar(i) <<std::endl;
          //std::cout << "seqChild->getValue(i): " << seqChild->getValue(i) << " " << seqChild->getChar(i) <<std::endl;
          if (alpha->isUnresolved(seqFather->getValue(i) ) ||  alpha->isUnresolved( seqChild->getValue(i) ) ) {
            //std::cout << "ARE UNRESOLVED" << std::endl;
          }
          else {
            size_t type = reg.getType ( seqFather->getValue(i) , seqChild->getValue(i) ) - 1 ; 
            if (type > nbTypes || type < 0) {
              std::cout << "Substitution out of the admitted alphabet. Not counting."<<std::endl; 
            }
            else {
              countsf[type] = 1;
            }
            //std::cout << "type: "<< type <<std::endl;
            }
          }
        }
        catch (Exception& ex) {
          std::cout <<"Character out of the admitted alphabet. Possibly a substitution involving an unresolved state (e.g. N)... Not counting." << std::endl;
        }
        for (unsigned int t = 0; t < nbTypes; ++t) {
          //std::cout << "k "<< k << " i "<< i << " t "<< t << std::endl;
          if (countsf[t] > 0) {
            counts[k][i].push_back(t) ;               
          }
        }
      }
    }
  }   
  return counts;
}


/******************************************************************************/


void buildCountTree(
                    const vector< vector<unsigned int> >& counts,
                    const vector<int>& ids,
                    Tree* cTree, 
                    unsigned int type)
{
  for (size_t i = 0; i < ids.size(); ++i) {
    if (cTree->hasFather(ids[i])) {
      cTree->setDistanceToFather(ids[i], counts[i][type]);
    }
  }
}

/******************************************************************************/


vector < std::string > buildAnnotatedCountMatrix(
                    SubstitutionRegister* reg,
                    const vector< vector<unsigned int> >& counts,
                    const vector<int>& ids,
                    TreeTemplate<Node> * cTree, 
                    TreeTemplate<Node> * spTree
                                                )
{
  vector < std::string > outputMatrix;
  string head  = string("NodeID") + string("\t") + string("FatherNodeID") + string("\t") + string("SpID") + string("\t") + string("FatherSpID") + string("\t") + string("FatherSon") ;
  for (size_t type = 0; type < counts[0].size(); ++type) {
    head = head + "\t" + reg->getTypeName(type+1) ;
  }
  head = head + "\n";
  outputMatrix.push_back(head);
  string line = "";
  for (size_t i = 0; i < ids.size(); ++i) {
    if (cTree->hasFather(ids[i])) {
      if (cTree->getNode(ids[i])->getFather()->hasNodeProperty("S") && cTree->getNode(ids[i])->hasNodeProperty("S")) {
        line = TextTools::toString(ids[i]) +  "\t" + TextTools::toString(cTree->getNode(ids[i])->getFather()->getId()) + "\t" + TextTools::toString(*(dynamic_cast<const BppString*>(cTree->getNode(ids[i])->getNodeProperty("S")))) + "\t" + TextTools::toString(*(dynamic_cast<const BppString*>(cTree->getNode(ids[i])->getFather()->getNodeProperty("S"))));
        if (spTree->getNode(TextTools::toInt(TextTools::toString(*(dynamic_cast<const BppString*>(cTree->getNode(ids[i])->getNodeProperty("S"))))))->hasFather() && (spTree->getNode(TextTools::toInt(TextTools::toString(*(dynamic_cast<const BppString*>(cTree->getNode(ids[i])->getNodeProperty("S"))))))->getFather()->getId() == TextTools::toInt(TextTools::toString(*(dynamic_cast<const BppString*>(cTree->getNode(ids[i])->getFather()->getNodeProperty("S")))))) ) {
          line = line +"\t" + "Y";
        }
        else {
          line = line +"\t" + "N";
        }
        for (size_t type = 0; type < counts[0].size(); ++type) {
          line = line +"\t" + TextTools::toString(counts[i][type]);
        }
        outputMatrix.push_back(line);
      }
      else {
        std::cout << "No S Node property"<<std::endl;
      }
    }
  }
  return outputMatrix;
}


/******************************************************************************/

      //vector on branches, map on sites, last vector on substitution types
 
vector < std::string > buildAnnotatedSitewiseCountOutput(
                    SubstitutionRegister* reg,
                    const vector< map< int, vector<unsigned int> > > & counts,
                    const vector<int>& ids,
                    TreeTemplate<Node> * cTree, 
                    TreeTemplate<Node> * spTree, 
                    const string &familyName, 
                    const size_t nbSubstitutionTypes, 
                    std::vector<unsigned int>& sumSubst )
{
  vector < std::string > outputMatrix;
  string line = "";
  size_t numNodes = cTree->getNumberOfNodes();
  for (size_t i = 0; i < numNodes; ++i) {
    Node * node = cTree->getNode(ids[i]);
    if (cTree->hasFather(ids[i])) {
      if (node->getFather()->hasNodeProperty("S") && node->hasNodeProperty("S")) {

        if (counts[i].size() > 0) { //There are substitution events
              for ( std::map< int, vector<unsigned int> >::const_iterator it = counts[i].begin(); it != counts[i].end(); it++ ) {
                for (size_t j = 0; j < it->second.size(); ++j) {
                  size_t type = it->second[j];
                  line = "event(" +  TextTools::toString(*(dynamic_cast<const BppString*>(node->getNodeProperty("S")))) + ",\"" + familyName + "\"," + reg->getTypeName(type+1) + "("+ TextTools::toString<int>(it->first + 1) + "))" ; //Adding plus 1, otherwise numbering starts at 0
                  sumSubst[it->first] += 1;
                  outputMatrix.push_back(line);
                }
              }
        }
      }
    }
    
     if ( node->hasNodeProperty("S") ) {
        if (node->getNumberOfSons()>0 ) { // node has children, could be a duplication, and could contain losses
          int fatherSpID = TextTools::toInt ( ( dynamic_cast<const BppString*> ( node->getNodeProperty ( "S" ) ) )->toSTL() );
          std::vector <Node *> sons = node->getSons();
          int a = TextTools::toInt ( ( dynamic_cast<const BppString*> ( sons[0]->getNodeProperty ( "S" ) ) )->toSTL() );
          int b = TextTools::toInt ( ( dynamic_cast<const BppString*> ( sons[1]->getNodeProperty ( "S" ) ) )->toSTL() );

          
          int aold = a;
          int bold = b;
          while ( a!=b ) {
              if ( a>b ) {
                int olda = a;
                Node* nodea = spTree->getNode ( a );
                a = nodea->getFather()->getId();
                if (a == fatherSpID) {  }
                else {
                  std::vector <Node *> sonsA = nodea->getFather()->getSons();
                  int lostBranch;
                  if (sonsA[0]->getId() == olda) {
                    lostBranch = sonsA[1]->getId(); 
                  }
                  else {
                    lostBranch = sonsA[0]->getId();
                  }
                  line = "event(" +  TextTools::toString(lostBranch) + ",\"" + familyName + "\"," + "loss" + ")" ;
                outputMatrix.push_back(line);
                      // lossA = lossA +1;
                }
              }
              else {
                int oldb = b;
                Node* nodeb = spTree->getNode ( b );
                b = nodeb->getFather()->getId();
                if (b == fatherSpID) {   }
                else {
                  std::vector <Node *> sonsb = nodeb->getFather()->getSons();
                  int lostBranch;
                  if (sonsb[0]->getId() == oldb) {
                    lostBranch = sonsb[1]->getId(); 
                  }
                  else {
                    lostBranch = sonsb[0]->getId();
                  }
                  line = "event(" +  TextTools::toString(lostBranch) + ",\"" + familyName + "\"," + "loss" + ")" ;
              outputMatrix.push_back(line);
                }
          //                lossB = lossB + 1;
              }
            }
        /*  sons[0]->setBranchProperty ( "L", BppString ( TextTools::toString ( lossA ) ) );
          sons[1]->setBranchProperty ( "L", BppString ( TextTools::toString ( lossB ) ) );
          node->setNodeProperty ( "S", BppString ( TextTools::toString ( a ) ) );*/
          if ( ( a == aold ) || ( a == bold ) ) { //There was a duplication
            node->setBranchProperty ( "Ev", BppString ( "D" ) );
            int dupBranch = a;
            line = "event(" +  TextTools::toString(dupBranch) + ",\"" + familyName + "\"," + "duplication" + ")" ;
          outputMatrix.push_back(line);
            //We also need to check whether tere have been losses
            if (aold > bold) { //loss in the lineage  leading to a
              a = aold;
                int olda = a;
                Node* nodea = spTree->getNode ( a );
                a = nodea->getFather()->getId();
                  std::vector <Node *> sonsA = nodea->getFather()->getSons();
                  int lostBranch;
                  if (sonsA[0]->getId() == olda) {
                    lostBranch = sonsA[1]->getId(); 
                  }
                  else {
                    lostBranch = sonsA[0]->getId();
                  }
                  line = "event(" +  TextTools::toString(lostBranch) + ",\"" + familyName + "\"," + "loss" + ")" ;
                outputMatrix.push_back(line);
                      // lossA = lossA +1;
                

            }
            else { //loss in the lineage leading to b
              b = bold;
                              int oldb = b;
                Node* nodeb = spTree->getNode ( b );
                b = nodeb->getFather()->getId();
                  std::vector <Node *> sonsb = nodeb->getFather()->getSons();
                  int lostBranch;
                  if (sonsb[0]->getId() == oldb) {
                    lostBranch = sonsb[1]->getId(); 
                  }
                  else {
                    lostBranch = sonsb[0]->getId();
                  }
                  line = "event(" +  TextTools::toString(lostBranch) + ",\"" + familyName + "\"," + "loss" + ")" ;
              outputMatrix.push_back(line);
                
            }
          }
        }
      }
      else {
        std::cout << "No S Node property"<<std::endl;
      }
  }
  return outputMatrix;
}



/*Compilation:
 g++ -lbpp-core -lbpp-seq -lbpp-phyl -o mapDnDsDuplicationTree mapDnDsDuplicationTree.cpp -std=c++0x -I/usr/local/include  -L. -L/usr/local/lib  
 */


int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*                     mapDnDsDuplicationTree, version 0.1.0                      *" << endl;
  cout << "* Authors: B. Boussau                       Created on  12/05/11 *" << endl;
  cout << "*                                           Last Modif. 27/01/13 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;
  
  if (args == 1)
    {
    help();
    exit(0);
    }
  
  try {
    
    BppApplication mapnh(args, argv, "mapDnDsDuplicationTree");
    mapnh.startTimer();

    string countsOut = ApplicationTools::getAFilePath("output.counts.file", mapnh.getParams(), false, false);
    if (FileTools::fileExists(countsOut) ) {
      cout <<"File "<<countsOut<< " already exists, exiting"<<endl;
      exit(-1);
    }

    
    Alphabet* alphabet = SequenceApplicationTools::getAlphabet(mapnh.getParams(), "", false);
    VectorSiteContainer* allSites = SequenceApplicationTools::getSiteContainer(alphabet, mapnh.getParams());
    VectorSiteContainer* sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, mapnh.getParams());
    delete allSites;
    
    ApplicationTools::displayResult("Number of sequences", TextTools::toString(sites->getNumberOfSequences()));
    ApplicationTools::displayResult("Number of sites", TextTools::toString(sites->getNumberOfSites()));
    
    //Get the species tree
    std::string spTreeFile =ApplicationTools::getStringParameter("species.tree.file",mapnh.getParams(),"none");
    if (spTreeFile=="none" )
      {
      std::cout << "\n\nNo Species tree was provided. The option init.species.tree is set to user (by default), which means that the option species.tree.file must be filled with the path of a valid tree file.\n" << std::endl;
      exit(-1);
      }
    ApplicationTools::displayResult("Species Tree file", spTreeFile);
    Newick newick(true);
    TreeTemplate<Node> * spTree = dynamic_cast < TreeTemplate < Node > * > (newick.read(spTreeFile));
    breadthFirstreNumber(*spTree);
    
    string spTreeIdOut = ApplicationTools::getStringParameter("output.species.tree_with_id.file", mapnh.getParams(), spTreeFile + "_withId");
      Nhx nhx(true);
      nhx.write(*spTree, spTreeIdOut);

    //Get the initial tree
//    Tree* tree = PhylogeneticsApplicationTools::getTree(mapnh.getParams());
    string treeIn = ApplicationTools::getAFilePath("input.tree.file", mapnh.getParams(), false, false);

    Tree* tree = nhx.read(treeIn);
    ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree->getNumberOfLeaves()));

    std::cout << nhx.treeToParenthesis (*tree) <<std::endl;
    //Convert to NHX if input tree is newick or nexus?
    string treeIdOut = ApplicationTools::getAFilePath("output.tree_with_id.file", mapnh.getParams(), false, false);
    if (treeIdOut != "none") {
      nhx.write(*tree, treeIdOut);
    }
    
    //Remove sequences not present in the tree
    vector <string> leafNames = tree->getLeavesNames();
    vector <string> sequenceNames = sites->getSequencesNames();
    
    for (unsigned int i = 0 ; i < sequenceNames.size() ; i++) {
      if (! VectorTools::contains(leafNames, sequenceNames[i])){
        sites->removeSequence(sequenceNames[i]);
      }
    }
    
    string familyName = ApplicationTools::getStringParameter("family.name", mapnh.getParams(), "none", "", true, false);

    
    //Perform the mapping:
    SubstitutionRegister* reg = 0;
    string regTypeDesc = ApplicationTools::getStringParameter("map.type", mapnh.getParams(), "All", "", true, false);
    string regType = "";
    map<string, string> regArgs;
    KeyvalTools::parseProcedure(regTypeDesc, regType, regArgs);
    //    auto_ptr<GeneticCode> geneticCode;
    unique_ptr<GeneticCode> geneticCode;
    bool stationarity = true;


      if (AlphabetTools::isCodonAlphabet(alphabet)) {
        string code = regArgs["code"];
        if (TextTools::isEmpty(code)) {
          code = "Standard";
          ApplicationTools::displayWarning("No genetic code provided, standard code used.");
        }
        geneticCode.reset(SequenceApplicationTools::getGeneticCode(dynamic_cast<CodonAlphabet*>(alphabet)->getNucleicAlphabet(), code));
	}


    //Now perform mapping using a JC model for ncleotides, or better model for protens and codons:
    unique_ptr<SubstitutionModel> model ;
    if (AlphabetTools::isNucleicAlphabet(alphabet)) {
//      model.reset( new JCnuc(dynamic_cast<NucleicAlphabet*>(alphabet) ) );
      model.reset( new GTR(dynamic_cast<NucleicAlphabet*>(alphabet) ) );
    } else if (AlphabetTools::isProteicAlphabet(alphabet)) {
//      model.reset( new JCprot(dynamic_cast<ProteicAlphabet*>(alphabet) ) );
//      model.reset( new LLG08_UL3(dynamic_cast<ProteicAlphabet*>(alphabet) ) );
        model.reset( new LG08(dynamic_cast<ProteicAlphabet*>(alphabet) ) );
    } else if (AlphabetTools::isCodonAlphabet(alphabet)) {
      if (ApplicationTools::parameterExists("model", mapnh.getParams())) {
        model.reset(  PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, geneticCode.get(), sites, mapnh.getParams() ) );
      } else {
      //  model = new CodonNeutralReversibleSubstitutionModel(
        //                                                    dynamic_cast<const CodonAlphabet*>(geneticCode->getSourceAlphabet()),
          //                                                  new JCnuc(dynamic_cast<CodonAlphabet*>(alphabet)->getNucleicAlphabet()));
       model.reset( PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, geneticCode.get(), sites, mapnh.getParams() ) );
      }
    }
    else
      throw Exception("Unsupported model!");


    vector <string> lines;
    if (regType == "All"|| regType == "all") {
      stationarity = ApplicationTools::getBooleanParameter("stationarity", regArgs, true);
      reg = new ComprehensiveSubstitutionRegister(model.get(), false);
    }
    else if (regType == "GC") {
      if (AlphabetTools::isNucleicAlphabet(alphabet)) {
        stationarity = ApplicationTools::getBooleanParameter("stationarity", regArgs, true);
        reg = new GCSubstitutionRegister(dynamic_cast<NucleotideSubstitutionModel*>(model.get()), false);
      } else
        throw Exception("GC categorization is only available for nucleotide alphabet!");
    } else if (regType == "TsTv") {
      if (AlphabetTools::isNucleicAlphabet(alphabet))
        reg = new TsTvSubstitutionRegister(dynamic_cast<NucleotideSubstitutionModel*>(model.get()));
      else
        throw Exception("TsTv categorization is only available for nucleotide alphabet!");
    } else if (regType == "DnDs") {
      if (AlphabetTools::isCodonAlphabet(alphabet)) {
     /*   string code = regArgs["code"];
        if (TextTools::isEmpty(code)) {
          code = "Standard";
          ApplicationTools::displayWarning("No genetic code provided, standard code used.");
        }
        geneticCode.reset(SequenceApplicationTools::getGeneticCode(dynamic_cast<CodonAlphabet*>(alphabet)->getNucleicAlphabet(), code));*/

   //By default, we choose F0... string freqOpt = ApplicationTools::getStringParameter("frequencies", mapnh.getParams(), "F0", "", true, 0);
  //  BppOFrequenciesSetFormat freqReader(BppOFrequenciesSetFormat::ALL, 0, 0);
   // freqReader.setGeneticCode(geneticCode.get()); //This uses the same instance as the one that will be used by the model.

	unique_ptr<FrequenciesSet> codonFreqs ;
codonFreqs.reset ( CodonFrequenciesSet::getFrequenciesSetForCodons 	( CodonFrequenciesSet::F0, geneticCode->clone() ) );

//    auto_ptr<FrequenciesSet> codonFreqs(freqReader.read(pCA, freqOpt, data, false));

      model.reset(new YN98( geneticCode.get(), codonFreqs.release()));


        reg = new DnDsSubstitutionRegister( dynamic_cast<CodonSubstitutionModel* >(model.get()), true); //true ensures that we count multiple substitutions as well
      } else
        throw Exception("DnDs categorization is only available for codon alphabet!");
    } else
      throw Exception("Unsupported substitution categorization: " + regType);
    
    //Printing out the alignment in which sites have been selected
    std::string alnNames = ApplicationTools::getAFilePath ("input.sequence.file",
        mapnh.getParams (), true,
        true);
    Fasta* fasta = new Fasta(); 
    fasta->writeSequences( alnNames+".refined", *sites );
    
    DiscreteDistribution* rDist = new ConstantDistribution( 1. );
    
    DRHomogeneousTreeLikelihood drtl(*tree, *sites, model.get(), rDist, false, false);
    drtl.initialize();
    
    //Optimization of parameters:
    PhylogeneticsApplicationTools::optimizeParameters(&drtl, drtl.getParameters(), mapnh.getParams(), "", true, true);
    
    vector<int> ids = drtl.getTree().getNodesId();
    ids.pop_back(); //remove root id.
    double thresholdSat = ApplicationTools::getDoubleParameter("count.max", mapnh.getParams(), -1);
    if (thresholdSat > 0)
      ApplicationTools::displayResult("Saturation threshold used", thresholdSat);

    bool siteBranchWise = ApplicationTools::getBooleanParameter("site.and.branch", mapnh.getParams(), true);

    //We want to get dN/dS ratio, so we need to set a default omega:
    auto_ptr<SubstitutionModel> nullModel(model->clone());
        
        ParameterList pl;
        const ParameterList pl0 = nullModel->getParameters();
        
	vector<string> pn = pl0.getMatchingParameterNames( "YN98.omega*" );
          for (size_t j = 0; j < pn.size(); ++j)
          {
            pl.addParameter(Parameter(pn[j], 1.0) );
          }
	nullModel->matchParametersValues(pl);

    if (siteBranchWise) {
       //vector on branches, map on sites, last vector on substitution types
      vector< map< int, vector<unsigned int> > > counts(ids.size()); 
      size_t nbTypes = 0;
      std::vector<unsigned int> sumSubst (sites->getNumberOfSites(), 0);
      bool asr = ApplicationTools::getBooleanParameter("ancestral.state.reconstruction", mapnh.getParams(), true);
      if (asr) {
        counts = getCountsPerBranchPerSiteASR(drtl, ids, sites, *reg, stationarity, thresholdSat, familyName, nbTypes );
      }
      else {
        counts = getCountsPerBranchPerSite(drtl, ids, model.get(), *reg, stationarity, thresholdSat, familyName, nbTypes);
      }
      lines = buildAnnotatedSitewiseCountOutput(reg, counts, ids, dynamic_cast<TreeTemplate<Node>*>(tree), spTree, familyName, nbTypes, sumSubst );
      //Writing out the site-wise sums of substitutions
      std::string sumSubstFile = alnNames+".sumSubst";
      std::ofstream out (sumSubstFile.c_str(), std::ios::out);
      for (unsigned int i = 0; i < sumSubst.size(); ++i) {
        out << sumSubst[i] <<std::endl;
      }
      out.close(); 
    }
else {
      vector< vector<unsigned int> > counts;
    counts = getCountsPerBranch(drtl, ids, model.get(), *reg, stationarity, thresholdSat);

    //Write count trees:
    string treePathPrefix = ApplicationTools::getStringParameter("output.counts.tree.prefix", mapnh.getParams(), "none");
    if (treePathPrefix != "none") {
      Newick newick;
      for (unsigned int i = 0; i < reg->getNumberOfSubstitutionTypes(); ++i) {
        string path = treePathPrefix + reg->getTypeName(i + 1) + string(".dnd");
        //ApplicationTools::displayResult(string("Output counts of type ") + TextTools::toString(i + 1) + string(" to file"), path);
        ApplicationTools::displayResult(string("Output counts of type ") + reg->getTypeName(i+1) + string(" to file"), path);
        Tree* cTree = tree->clone();
        buildCountTree(counts, ids, cTree, i);
        newick.write(*cTree, path);
        delete cTree;
      }
    }
    
    //Outputs the table containing counts
    lines = buildAnnotatedCountMatrix(reg, counts, ids, dynamic_cast<TreeTemplate<Node>*>(tree), spTree);
}    
    std::ofstream out (countsOut.c_str(), std::ios::out);
    for (unsigned int i = 0; i < lines.size(); ++i) {
      out << lines[i] <<std::endl;
    }
    out.close(); 
    

    //Cleaning up:
    delete alphabet;
    delete sites;
    delete tree;
   // delete model;
    delete rDist;
    delete reg;
    mapnh.done();
    
  }
  catch (exception& e)
  {
  cout << e.what() << endl;
  exit(-1);
  }
  
  return (0);
}
