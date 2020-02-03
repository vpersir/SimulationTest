/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.genetsimulator.entity;

import cern.colt.list.BooleanArrayList;
import cern.colt.list.DoubleArrayList;
import cern.colt.list.IntArrayList;
import cern.colt.map.OpenIntIntHashMap;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.apache.commons.math.stat.correlation.PearsonsCorrelation;
import org.cobi.genetsimulator.controller.LDPair;
import org.cobi.genetsimulator.controller.LDSparseMatrix;
import org.cobi.genetsimulator.controller.SparseLDMatrixElementIndex1Comparator;

/**
 *
 * @author MX Li
 */
public class GenotypeBasedLDSparseMatrix extends LDSparseMatrix {
  // a special matrix to store LD matrix of SNPs
  //it does not store the whole LD data but the genotype data. calcualte the genotype correlation when it is needed  

  List<Individual> indList;
  OpenIntIntHashMap indexGenotypePosMap;
  List<LDPair> ldList = new ArrayList<LDPair>();
  SparseLDMatrixElementIndex1Comparator elementIndex1Comparator = new SparseLDMatrixElementIndex1Comparator();
  LDPair element = new LDPair(0, 0, 0);
  Set<Integer> allUniqueIndexes = new HashSet<Integer>();
  DoubleArrayList genotypeCodeList1 = new DoubleArrayList();
  DoubleArrayList genotypeCodeList2 = new DoubleArrayList();
  PearsonsCorrelation ps = new PearsonsCorrelation();

  public boolean isEmpty() {
    return (indList == null || indList.isEmpty()) || (indexGenotypePosMap == null || indexGenotypePosMap.isEmpty());
  }

  public GenotypeBasedLDSparseMatrix(List<Individual> indList, OpenIntIntHashMap indexGenotypePosMap) {
    this.indList = indList;
    this.indexGenotypePosMap = indexGenotypePosMap;
    IntArrayList temList = indexGenotypePosMap.keys();
    int size = temList.size();
    allUniqueIndexes.clear();
    for (int i = 0; i < size; i++) {
      allUniqueIndexes.add(temList.getQuick(i));
    }
    temList.clear();
  }

  public Set<Integer> getAllUniqueIndexes() {
    return allUniqueIndexes;
  }

  public DoubleMatrix2D subDenseLDMatrix(IntArrayList indexes) throws Exception {
    int dim = indexes.size();
    indexes.quickSort();
    DoubleMatrix2D corrMat = new DenseDoubleMatrix2D(dim, dim);
    double x = 0;
    for (int i = 0; i < dim; i++) {
      corrMat.setQuick(i, i, 1);
      for (int j = i + 1; j < dim; j++) {
        x = getLDAt(indexes.getQuick(i), indexes.getQuick(j));
        corrMat.setQuick(i, j, x);
        corrMat.setQuick(j, i, x);
      }
    }
    // System.out.println(corrMat.toString());
    return corrMat;
  }
 

  
  public DoubleMatrix2D subDenseLDMatrixR2(IntArrayList indexes) throws Exception {
    int dim = indexes.size();
    indexes.quickSort();
    DoubleMatrix2D corrMat = new DenseDoubleMatrix2D(dim, dim);
    double x = 0;
    for (int i = 0; i < dim; i++) {
      corrMat.setQuick(i, i, 1);
      for (int j = i + 1; j < dim; j++) {
        x = getLDAt(indexes.getQuick(i), indexes.getQuick(j));
        //x = x * x;
        corrMat.setQuick(i, j, x);
        corrMat.setQuick(j, i, x);
      }
    }
    // System.out.println(corrMat.toString());
    return corrMat;
  }

  public void mafAndAlleles(IntArrayList indexes, double[] mafs, boolean[] malleles) throws Exception {
    int dim = indexes.size();
    indexes.quickSort();
    Arrays.fill(mafs, 0);
    Arrays.fill(malleles, false);

    int indivSize = indList.size();
    int effectIndivNum = 0;
    for (int i = 0; i < dim; i++) {
      int pos1 = indexGenotypePosMap.get(indexes.getQuick(i));
      effectIndivNum = 0;
      for (int k = 0; k < indivSize; k++) {
        StatusGtySet gty = indList.get(k).markerGtySet;
        if (gty.existence.getQuick(pos1)) {
          if (!gty.paternalChrom.getQuick(pos1) && !gty.maternalChrom.getQuick(pos1)) {
            mafs[i] += 2;
          } else if (gty.paternalChrom.getQuick(pos1) && gty.maternalChrom.getQuick(pos1)) {
            mafs[i] += 0;
          } else {
            mafs[i] += 1;
          }
          effectIndivNum += 2;
        }
      }
      mafs[i] /= effectIndivNum;
      if (mafs[i] > 0.5) {
        mafs[i] = 1 - mafs[i];
        malleles[i] = true;
      }
      System.out.println(i + "\t" + mafs[i] + "\t" + malleles[i]);
    }

  }

  public float getLDAt(int index1, int index2) throws Exception {

    // if (!hasSorted) {
    // throw new Exception("Not a sorted LDSparseMatrix");
    // }
    //ensure curIndex1 is always less than i2, always
    int tmpIndex;
    if (index1 > index2) {
      tmpIndex = index1;
      index1 = index2;
      index2 = tmpIndex;
    } else if (index1 == index2) {
      return 1;
    }
    element.i1 = index1;
    tmpIndex = Collections.binarySearch(ldList, element, elementIndex1Comparator);
    if (tmpIndex < 0) {
      float r = calculateGenotypeCorrelation(index1, index2);
      LDPair ele = new LDPair(index1, index2, r);
      /*
            double x = r * r;
            //when r2 //convert to p-value corrlation coefficent
            //y = 0.7723x6 - 1.5659x5 + 1.201x4 - 0.2355x3 + 0.2184x2 + 0.6086x
            x = (((((0.7723 * x - 1.5659) * x + 1.201) * x - 0.2355) * x + 0.2184) * x + 0.6086) * x;
            
            LDPair ele = new LDPair(index1, index2, (float) x);
             *
       */
      ldList.add(-tmpIndex - 1, ele);
      return (float) r;
    } else {
      int totalSize = ldList.size();
      int i = tmpIndex;
      while (i < totalSize && ldList.get(i).i1 == index1) {
        if (ldList.get(i).i2 == index2) {
          return ldList.get(i).ld;
        }
        i++;
      }
      i = tmpIndex - 1;
      while (i >= 0 && ldList.get(i).i1 == index1) {
        if (ldList.get(i).i2 == index2) {
          return ldList.get(i).ld;
        }
        i--;
      }
      //otherwise has to calcualte the LD
      float r = calculateGenotypeCorrelation(index1, index2);
      LDPair ele = new LDPair(index1, index2, r);
      /*
            double x = r * r;
            //when r2 //convert to p-value corrlation coefficent
            //y = 0.7723x6 - 1.5659x5 + 1.201x4 - 0.2355x3 + 0.2184x2 + 0.6086x
            x = (((((0.7723 * x - 1.5659) * x + 1.201) * x - 0.2355) * x + 0.2184) * x + 0.6086) * x;
            
            LDPair ele = new LDPair(index1, index2, (float) x);
             *
       */
      ldList.add(tmpIndex, ele);
      return (float) r;
    }
  }

  float binarySearch(int tmpIndex2, int left, int right) {
    if (left > right) {
      return 0;
    }
    int middle = (left + right) / 2;

    if (ldList.get(middle).i2 == tmpIndex2) {
      return ldList.get(middle).ld;
    } else if (ldList.get(middle).i2 > tmpIndex2) {
      return binarySearch(tmpIndex2, left, middle - 1);
    } else {
      return binarySearch(tmpIndex2, middle + 1, right);
    }
  }

  public void releaseLDData() {
    ldList.clear();
  }

  public float calculateGenotypeCorrelation(int index1, int index2) throws Exception {
    if (index1 == index2) {
      return 1;
    }
    int indivSize = indList.size();
    genotypeCodeList1.clear();
    genotypeCodeList2.clear();
    int pos1 = indexGenotypePosMap.get(index1);
    int pos2 = indexGenotypePosMap.get(index2);
    if (index1 == 65330464 && index2 == 65428893) {
      int sxss = 0;
    }
    //note the chromsome number is alwalys eaqual to  indivSize * 2 at every locus due to the missing genotypes
    for (int k = 0; k < indivSize; k++) {
      StatusGtySet gty = indList.get(k).markerGtySet;
      if (gty.existence.getQuick(pos1) && gty.existence.getQuick(pos2)) {
        if (!gty.paternalChrom.getQuick(pos1) && !gty.maternalChrom.getQuick(pos1)) {
          genotypeCodeList1.add(2);
        } else if (gty.paternalChrom.getQuick(pos1) && gty.maternalChrom.getQuick(pos1)) {
          genotypeCodeList1.add(0);
        } else {
          genotypeCodeList1.add(1);
        }

        if (!gty.paternalChrom.getQuick(pos2) && !gty.maternalChrom.getQuick(pos2)) {
          genotypeCodeList2.add(2);
        } else if (gty.paternalChrom.getQuick(pos2) && gty.maternalChrom.getQuick(pos2)) {
          genotypeCodeList2.add(0);
        } else {
          genotypeCodeList2.add(1);
        }
      }
    }

    double[] list1 = new double[genotypeCodeList1.size()];
    double[] list2 = new double[genotypeCodeList2.size()];

    for (int i = 0; i < list1.length; i++) {
      list1[i] = genotypeCodeList1.getQuick(i);
    }
    for (int i = 0; i < list2.length; i++) {
      list2[i] = genotypeCodeList2.getQuick(i);
    }
    double r = ps.correlation(list1, list2);

    if (Double.isNaN(r)) {
      return 0;
    }
    return (float) r;
  }
}
