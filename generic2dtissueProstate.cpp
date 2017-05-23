/**
 * 	
 *	\brief Tissue prostatique
 * Tissue contenant plusieurs cellules prostatique (prostateCell)
 * 
 * 
 * 
 *	generic2dtissueProstate.cpp
 *
 *  Copyright (c) 2001 INSERM. All rights reserved.
 */



#include "generic2dtissueProstate.h"


/**
 * Contructeur.
 * Création d'un vecteur de prostateCell
 */
generic2dtissueProstate::generic2dtissueProstate() : Model(DESS,0,0,0,0,TISSUESIZE_ROW*TISSUESIZE_COL)
{
  // creation of the cells composing the tissu model
  for (int i=0;i<getNumComponents();i++)
  {      
    components->at(i) = new prostateCell();				/**< Création du vecteur de prostateCell */
    numOutputs += (components->at(i))->getNumOutputs();
  }


  ifstream tissueLoad("/home/ciferri/Bureau/stage/branch_github/Ciferri/data/inVes16.dat");
  if(tissueLoad){
	  int i=0;
	  string valeur;
	  while (std::getline(tissueLoad, valeur)){
			istringstream ( valeur ) >> vectTissue[i];
			i++;
		}
  }


  int count = 0;
  for (int i=0;i<TISSUESIZE_ROW;i++)
  {
    for (int j=0;j<TISSUESIZE_COL;j++)
    {
      tissue[i][j] = components->at(count);				/**< tableau du tissue */
      gradTissue[i][j] = gradLin[count];
      count += 1;
    }
  }
  //coupling parameters

  coeff_diff = 1.835*pow(10,-9);
  size = 20;
  time_step_s = 0.07;
  coupParam = (coeff_diff*pow(10,12)*time_step_s)/(size*size);
  
  Vmax = 15.2;
  K_conso = 3.035;
  
  filterIn[0][0]=1; filterIn[0][1]=1; 	filterIn[0][2]=1;
  filterIn[1][0]=1; filterIn[1][1]=-8;	filterIn[1][2]=1;
  filterIn[2][0]=1; filterIn[2][1]=1; 	filterIn[2][2]=1;
  				
  filterOut[0][0]=-1; filterOut[0][1]=-1; filterOut[0][2]=-1;
  filterOut[1][0]=-1; filterOut[1][1]=8;  filterOut[1][2]=-1;
  filterOut[2][0]=-1; filterOut[2][1]=-1; filterOut[2][2]=-1;
  			
}

/**
 * 
 * Destructeur
 * 
 */
generic2dtissueProstate::~generic2dtissueProstate()
{
}

int generic2dtissueProstate::ModelStart()
{

  for (int i=0;i<getNumComponents();i++)
    (components->at(i))->ModelStart();

  return 0;
}

/**
 * 
 * Initialisation du tissue.
 * 
 * */
int generic2dtissueProstate::ModelInitSim()
{

  for (int i=0;i<getNumComponents();i++)
  {
	(components->at(i))->ModelRAZ();	/**< Fonction qui initialise toute les pO2 à zéro */
  }
  return 1;
}

/**
 * 
 * Surcharge de la classe ModelInitSim
 * 
 * */
int generic2dtissueProstate::ModelInitSim(int x, int y,double pO2)
{
  int indice = Coord_XY_to_K(x,y);		/**< Fonction qui prend en paramètre x et y et renvoie un indice */
  (components->at(indice))->ModelRAZ(0.0, 0.0, 0.0, 1.0, pO2); 	/**< Pour un vaisseau, on met la pO2 à 42mmHg */
  return 1;
}

//Fonction qui initialise la matrice d'entrer en fonction de l'image

//~ int generic2dtissueProstate::ModelInitSim(double *tab,int size_x, int size_y)
//~ {
	//~ int count= 0;
  //~ for (int i=0;i<size_x-1;i++)
  //~ {
	//~ for (int j;j<size_y-1;j++)
	//~ {
		//~ //if(valeur de l'image == rouge) { (components->at(indice))->ModelRAZ(0.0, 0.0, 0.0, 1.0, pO2); }
		//~ //else if(valeur de limage == bleu) { (components->at(indice))->ModelRAZ(0.0, 0.0, 1.0, 0.0, pO2); }
		//~ //etc 
	//~ count += 1;
	//~ }
  //~ }
  //~ return 1;
//~ }

int generic2dtissueProstate::NextSampleHit()
{
  for (int i=0;i<getNumComponents();i++)
    (components->at(i))->NextSampleHit();

  return 0;
}

int generic2dtissueProstate::ModelOutputs()
{

  for (int i=0;i<getNumComponents();i++)
    (components->at(i))->ModelOutputs();

  return 0;
}

/**
 * 
 * Mise à jour du tissue 
 * 
 * */
int generic2dtissueProstate::ModelUpdate(double time)
{

  double tempGradIn;
  double tempGradOut;
	//cout << "time : " <<time <<endl;


  //updating of all cells
  for (int i=0;i<getNumComponents();i++)
    (components->at(i))->ModelUpdate(time);
    

  int cpt = 0;
  for (int i=0;i<TISSUESIZE_ROW;i++)
		for (int j=0;j<TISSUESIZE_COL;j++){
			if(vectTissue[cpt] != 0)
				ModelInitSim(j,i,INIT_VASCULAR_PO2);
			cpt += 1;
		}

 
  for (int i=0;i<TISSUESIZE_ROW;i++)
		for (int j=0;j<TISSUESIZE_COL;j++){
			initImageTissue[i][j] = ((prostateCell*)tissue[i][j])->getST_X();
			Conso_PO2[i][j] = (initImageTissue[i][j]*Vmax*time_step_s)/(K_conso + initImageTissue[i][j]);
		}
  /// --Cas particulier--



  /**
   * First row / First Cell
   * */   
  tempGradIn = filterIn[1][2]*((prostateCell*)tissue[1][0])->getST_X() +
    filterIn[2][1]*((prostateCell*)tissue[0][1])->getST_X() +
    filterIn[2][2]*((prostateCell*)tissue[1][1])->getST_X() -
    (filterIn[1][2]+filterIn[2][1]+filterIn[2][2])*((prostateCell*)tissue[0][0])->getST_X();  
  resDiffusionIn[0][0] = tempGradIn*coupParam/3;
  
  tempGradOut = filterOut[1][2]*((prostateCell*)tissue[1][0])->getST_X() +
    filterOut[2][1]*((prostateCell*)tissue[0][1])->getST_X() +
    filterOut[2][2]*((prostateCell*)tissue[1][1])->getST_X() -
    (filterOut[1][2]+filterOut[2][1]+filterOut[2][2])*((prostateCell*)tissue[0][0])->getST_X();  
  resDiffusionOut[0][0] = tempGradOut*coupParam/3;

  /**
   * First row / Last Cell
   * */   
  tempGradIn = filterIn[1][0]*((prostateCell*)tissue[TISSUESIZE_ROW-2][0])->getST_X() +
    filterIn[2][0]*((prostateCell*)tissue[TISSUESIZE_ROW-2][1])->getST_X() +
    filterIn[2][1]*((prostateCell*)tissue[TISSUESIZE_ROW-1][1])->getST_X() -
    (filterIn[1][0]+filterIn[2][0]+filterIn[2][1])*((prostateCell*)tissue[TISSUESIZE_ROW-1][0])->getST_X();  
  resDiffusionIn[TISSUESIZE_ROW-1][0] = tempGradIn*coupParam/3;
  
  tempGradOut = filterOut[1][0]*((prostateCell*)tissue[TISSUESIZE_ROW-2][0])->getST_X() +
    filterOut[2][0]*((prostateCell*)tissue[TISSUESIZE_ROW-2][1])->getST_X() +
    filterOut[2][1]*((prostateCell*)tissue[TISSUESIZE_ROW-1][1])->getST_X() -
    (filterOut[1][0]+filterOut[2][0]+filterOut[2][1])*((prostateCell*)tissue[TISSUESIZE_ROW-1][0])->getST_X();  
  resDiffusionOut[TISSUESIZE_ROW-1][0] = tempGradOut*coupParam/3;

  /**
   * Last row / First Cell
   * */   
  tempGradIn = filterIn[0][1]*((prostateCell*)tissue[0][TISSUESIZE_COL-2])->getST_X() +
    filterIn[0][2]*((prostateCell*)tissue[1][TISSUESIZE_COL-2])->getST_X() +
    filterIn[1][2]*((prostateCell*)tissue[1][TISSUESIZE_COL-1])->getST_X() -
    (filterIn[0][1]+filterIn[0][2]+filterIn[1][2])*((prostateCell*)tissue[0][TISSUESIZE_COL-1])->getST_X();  
  resDiffusionIn[0][TISSUESIZE_COL-1] = tempGradIn*coupParam/3;
  
  tempGradOut = filterOut[0][1]*((prostateCell*)tissue[0][TISSUESIZE_COL-2])->getST_X() +
    filterOut[0][2]*((prostateCell*)tissue[1][TISSUESIZE_COL-2])->getST_X() +
    filterOut[1][2]*((prostateCell*)tissue[1][TISSUESIZE_COL-1])->getST_X() -
    (filterOut[0][1]+filterOut[0][2]+filterOut[1][2])*((prostateCell*)tissue[0][TISSUESIZE_COL-1])->getST_X();  
  resDiffusionOut[0][TISSUESIZE_COL-1] = tempGradOut*coupParam/3;


  /**
   * Last row / Last Cell
   * */   
  tempGradIn = filterIn[0][0]*((prostateCell*)tissue[TISSUESIZE_ROW-2][TISSUESIZE_COL-2])->getST_X() +
    filterIn[0][1]*((prostateCell*)tissue[TISSUESIZE_ROW-1][TISSUESIZE_COL-2])->getST_X() +
    filterIn[1][0]*((prostateCell*)tissue[TISSUESIZE_ROW-2][TISSUESIZE_COL-1])->getST_X() -
    (filterIn[0][0]+filterIn[0][1]+filterIn[1][0])*((prostateCell*)tissue[TISSUESIZE_ROW-1][TISSUESIZE_COL-1])->getST_X();
  resDiffusionIn[TISSUESIZE_ROW-1][TISSUESIZE_COL-1] = tempGradIn*coupParam/3;
  
  tempGradOut = filterOut[0][0]*((prostateCell*)tissue[TISSUESIZE_ROW-2][TISSUESIZE_COL-2])->getST_X() +
    filterOut[0][1]*((prostateCell*)tissue[TISSUESIZE_ROW-1][TISSUESIZE_COL-2])->getST_X() +
    filterOut[1][0]*((prostateCell*)tissue[TISSUESIZE_ROW-2][TISSUESIZE_COL-1])->getST_X() -
    (filterOut[0][0]+filterOut[0][1]+filterOut[1][0])*((prostateCell*)tissue[TISSUESIZE_ROW-1][TISSUESIZE_COL-1])->getST_X();
  resDiffusionOut[TISSUESIZE_ROW-1][TISSUESIZE_COL-1] = tempGradOut*coupParam/3;

  /**
   * First row / Middle cells
   * */
  for (int j=1;j<TISSUESIZE_ROW-1;j++)
  {	
    tempGradIn = filterIn[1][0]*((prostateCell*)tissue[j-1][0])->getST_X() +
      filterIn[1][2]*((prostateCell*)tissue[j+1][0])->getST_X() +
      filterIn[2][0]*((prostateCell*)tissue[j-1][1])->getST_X() +
      filterIn[2][1]*((prostateCell*)tissue[j][1])->getST_X() +
      filterIn[2][2]*((prostateCell*)tissue[j+1][1])->getST_X() -
      (filterIn[0][0]+filterIn[0][1]+filterIn[0][2]+filterIn[1][0]+filterIn[1][2])*((prostateCell*)tissue[j][0])->getST_X();    
    resDiffusionIn[j][0] = tempGradIn*coupParam/5; 

    tempGradOut = filterOut[1][0]*((prostateCell*)tissue[j-1][0])->getST_X() +
      filterOut[1][2]*((prostateCell*)tissue[j+1][0])->getST_X() +
      filterOut[2][0]*((prostateCell*)tissue[j+1][1])->getST_X() +
      filterOut[2][1]*((prostateCell*)tissue[j][1])->getST_X() +
      filterOut[2][2]*((prostateCell*)tissue[j+1][1])->getST_X() -
      (filterOut[0][0]+filterOut[0][1]+filterOut[0][2]+filterOut[1][0]+filterOut[1][2])*((prostateCell*)tissue[j][0])->getST_X();    
    resDiffusionOut[j][0] = tempGradOut*coupParam/5; 
    
  }

  /**
   * Last row / Middle Cells
   * */
  for (int j=1;j<TISSUESIZE_ROW-1;j++)
  {	
    tempGradIn = filterIn[0][0]*((prostateCell*)tissue[j-1][TISSUESIZE_COL-2])->getST_X() +
      filterIn[0][1]*((prostateCell*)tissue[j][TISSUESIZE_COL-2])->getST_X() +
      filterIn[0][2]*((prostateCell*)tissue[j+1][TISSUESIZE_COL-2])->getST_X() +
      filterIn[1][0]*((prostateCell*)tissue[j-1][TISSUESIZE_COL-1])->getST_X() +
      filterIn[1][2]*((prostateCell*)tissue[j+1][TISSUESIZE_COL-1])->getST_X() -
      (filterIn[0][0]+filterIn[0][1]+filterIn[0][2]+filterIn[1][0]+filterIn[1][2])*((prostateCell*)tissue[j][TISSUESIZE_COL-1])->getST_X();    
    resDiffusionIn[j][TISSUESIZE_COL-1] = tempGradIn*coupParam/5;
    
    tempGradOut = filterOut[0][0]*((prostateCell*)tissue[j-1][TISSUESIZE_COL-2])->getST_X() +
      filterOut[0][1]*((prostateCell*)tissue[j][TISSUESIZE_COL-2])->getST_X() +
      filterOut[0][2]*((prostateCell*)tissue[j+1][TISSUESIZE_COL-2])->getST_X() +
      filterOut[1][0]*((prostateCell*)tissue[j-1][TISSUESIZE_COL-1])->getST_X() +
      filterOut[1][2]*((prostateCell*)tissue[j+1][TISSUESIZE_COL-1])->getST_X() -
      (filterOut[0][0]+filterOut[0][1]+filterOut[0][2]+filterOut[1][0]+filterOut[1][2])*((prostateCell*)tissue[j][TISSUESIZE_COL-1])->getST_X();    
    resDiffusionOut[j][TISSUESIZE_COL-1] = tempGradOut*coupParam/5;
    

  }

  /**
   * Last column / Middle Cells  
   * */
  for (int i=1;i<TISSUESIZE_COL-1;i++)
  {	  
    tempGradIn = filterIn[0][0]*((prostateCell*)tissue[TISSUESIZE_ROW-2][i-1])->getST_X() +
      filterIn[0][1]*((prostateCell*)tissue[TISSUESIZE_ROW-1][i-1])->getST_X() +
      filterIn[1][0]*((prostateCell*)tissue[TISSUESIZE_ROW-2][i])->getST_X() +
      filterIn[2][0]*((prostateCell*)tissue[TISSUESIZE_ROW-2][i+1])->getST_X() +
      filterIn[2][1]*((prostateCell*)tissue[TISSUESIZE_ROW-1][i+1])->getST_X() -
      (filterIn[0][0]+filterIn[0][1]+filterIn[1][0]+filterIn[2][0]+filterIn[2][1])*((prostateCell*)tissue[TISSUESIZE_ROW-1][i])->getST_X();
    resDiffusionIn[TISSUESIZE_ROW-1][i] = tempGradIn*coupParam/5;
    
    tempGradOut = filterOut[0][0]*((prostateCell*)tissue[TISSUESIZE_ROW-2][i-1])->getST_X() +
      filterOut[0][1]*((prostateCell*)tissue[TISSUESIZE_ROW-1][i-1])->getST_X() +
      filterOut[1][0]*((prostateCell*)tissue[TISSUESIZE_ROW-2][i])->getST_X() +
      filterOut[2][0]*((prostateCell*)tissue[TISSUESIZE_ROW-2][i+1])->getST_X() +
      filterOut[2][1]*((prostateCell*)tissue[TISSUESIZE_ROW-1][i+1])->getST_X() -
      (filterOut[0][0]+filterOut[0][1]+filterOut[1][0]+filterOut[2][0]+filterOut[2][1])*((prostateCell*)tissue[TISSUESIZE_ROW-1][i])->getST_X();
    resDiffusionOut[TISSUESIZE_ROW-1][i] = tempGradOut*coupParam/5;

  }

  /**
   * First column / Middle Cells
   * */
  for (int i=1;i<TISSUESIZE_COL-1;i++) 
  {	  
		tempGradIn = filterIn[0][1]*((prostateCell*)tissue[0][i-1])->getST_X() +
			filterIn[0][2]*((prostateCell*)tissue[1][i-1])->getST_X() +
			filterIn[1][2]*((prostateCell*)tissue[1][i])->getST_X() +
			filterIn[2][1]*((prostateCell*)tissue[0][i+1])->getST_X() +
			filterIn[2][2]*((prostateCell*)tissue[1][i+1])->getST_X() -
			(filterIn[0][1]+filterIn[0][2]+filterIn[1][2]+filterIn[2][1]+filterIn[2][2])*((prostateCell*)tissue[0][i])->getST_X();
		resDiffusionIn[0][i] = tempGradIn*coupParam/5;
			
		tempGradOut = filterOut[0][1]*((prostateCell*)tissue[0][i-1])->getST_X() +
			filterOut[0][2]*((prostateCell*)tissue[1][i-1])->getST_X() +
			filterOut[1][2]*((prostateCell*)tissue[1][i])->getST_X() +
			filterOut[2][1]*((prostateCell*)tissue[0][i+1])->getST_X() +
			filterOut[2][2]*((prostateCell*)tissue[1][i+1])->getST_X() -
			(filterOut[0][1]+filterOut[0][2]+filterOut[1][2]+filterOut[2][1]+filterOut[2][2])*((prostateCell*)tissue[0][i])->getST_X();
		resDiffusionOut[0][i] = tempGradOut*coupParam/5;

  }

  /**
   * coupling between middle cells
   * */
  for (int i=1;i<TISSUESIZE_ROW-1;i++)
  {
    for (int j=1;j<TISSUESIZE_COL-1;j++)
    {
			tempGradIn = filterIn[0][0]*((prostateCell*)tissue[i-1][j-1])->getST_X() +
				filterIn[0][1]*((prostateCell*)tissue[i-1][j])->getST_X() +
				filterIn[0][2]*((prostateCell*)tissue[i-1][j+1])->getST_X() +
				filterIn[1][0]*((prostateCell*)tissue[i][j-1])->getST_X() +
				filterIn[1][1]*((prostateCell*)tissue[i][j])->getST_X() +
				filterIn[1][2]*((prostateCell*)tissue[i][j+1])->getST_X() +
				filterIn[2][0]*((prostateCell*)tissue[i+1][j-1])->getST_X() +
				filterIn[2][1]*((prostateCell*)tissue[i+1][j])->getST_X() +
				filterIn[2][2]*((prostateCell*)tissue[i+1][j+1])->getST_X();
			resDiffusionIn[i][j] = tempGradIn*coupParam/8;
				
			tempGradOut = filterOut[0][0]*((prostateCell*)tissue[i-1][j-1])->getST_X() +
				filterOut[0][1]*((prostateCell*)tissue[i-1][j])->getST_X() +
				filterOut[0][2]*((prostateCell*)tissue[i-1][j+1])->getST_X() +
				filterOut[1][0]*((prostateCell*)tissue[i][j-1])->getST_X() +
				filterOut[1][1]*((prostateCell*)tissue[i][j])->getST_X() +
				filterOut[1][2]*((prostateCell*)tissue[i][j+1])->getST_X() +
				filterOut[2][0]*((prostateCell*)tissue[i+1][j-1])->getST_X() +
				filterOut[2][1]*((prostateCell*)tissue[i+1][j])->getST_X() +
				filterOut[2][2]*((prostateCell*)tissue[i+1][j+1])->getST_X();
			resDiffusionOut[i][j] = tempGradOut*coupParam/8;
    }
  }


  
  for (int j=0;j<TISSUESIZE_COL;j++)
    for (int i=0;i<TISSUESIZE_ROW;i++){
			((prostateCell*)tissue[i][j])->setIN_Z(initImageTissue[i][j] + resDiffusionIn[i][j] - resDiffusionOut[i][j]);//- Conso_PO2[i][j]);
			((prostateCell*)tissue[i][j])->setIN_X(Conso_PO2[i][j]);
    }
	

		
  
  return 0;
}

int generic2dtissueProstate::ModelDerivatives(double time, ParamVect *variables, ParamVect *derivadas)
{
  for (int i=0;i<getNumComponents();i++)
    (components->at(i))->ModelDerivatives(time,variables,derivadas);

  return 0;
}

int generic2dtissueProstate::ModelTerminate()
{
  for (int i=0;i<getNumComponents();i++)
    (components->at(i))->ModelTerminate();

  return 0;
}

/**
 * Fonction qui prend en paramètre des coordonné x et y et rend une position de vecteur.
 * */
int generic2dtissueProstate::Coord_XY_to_K(int x, int y)
{
		return y*TISSUESIZE_COL + x;

}

/**
 * Fonction qui prend en paramètre une position de vecteur et qui rend un tableau de coordonné x et y.
 * */
//~ int *generic2dtissueProstate::Coord_K_to_XY(int k)
//~ {
	//~ int *tab= new int[2];
	//~ tab[0] = k/TISSUESIZE_ROW;
    //~ tab[1] = k%TISSUESIZE_ROW;
    
    //~ return tab;
    //~ // dans le main : int *array = CoordVectXY(k);	
//~ }
