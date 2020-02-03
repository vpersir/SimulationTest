/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.genetsimulator.controller;

import cern.colt.bitvector.BitVector;
import cern.colt.list.IntArrayList;
import cern.jet.stat.Probability;
import java.io.BufferedReader;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.List;
import org.apache.commons.math.stat.inference.ChiSquareTest;
import org.apache.commons.math.stat.inference.ChiSquareTestImpl;
import org.cobi.genetsimulator.entity.Individual;
import org.cobi.genetsimulator.entity.StatusGtySet;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.StringTokenizer;
import org.cobi.genetsimulator.entity.AnnotSNP;
import org.cobi.genetsimulator.entity.AssociationSNP;
import org.cobi.util.stat.FisherExactTest;
import org.cobi.util.stat.LogisticRegression;
import org.cobi.util.stat.SimpleLinearRegression;

/**
 *
 * @author MX Li
 */
public class GenetAssociationAnalyzerP {

  boolean toDebug = false;
  BufferedWriter debugOut = null;

  public boolean isToDebug() {
    return toDebug;
  }

  public void setToDebug(boolean toDebug) {
    this.toDebug = toDebug;
  }

  public GenetAssociationAnalyzerP() {
  }

  public BitVector allelicAssociationTest(List<Individual> indivList, double alpha) throws Exception {

    if (indivList == null || indivList.isEmpty()) {
      return null;
    }

    int snpNum = indivList.get(0).markerGtySet.paternalChrom.size();
    int indivSize = indivList.size();

    long[][] counts = new long[2][2];
    int rowNum = counts.length;
    ChiSquareTest csTest = new ChiSquareTestImpl();
    BitVector hits = new BitVector(snpNum);
    hits.replaceFromToWith(0, snpNum - 1, false);
    double p;
    // table for contigency test
    //           0    1
    // control  counts[0][0]   counts[0][1]
    // case     counts[1][0]   counts[1][1]

    for (int i = 0; i < snpNum; i++) {
      for (int j = 0; j < rowNum; j++) {
        Arrays.fill(counts[j], 0);
      }
      for (int k = 0; k < indivSize; k++) {
        StatusGtySet gty = indivList.get(k).markerGtySet;
        if (indivList.get(k).getAffectedStatus() == 2) {
          if (gty.existence.getQuick(i)) {
            if (gty.paternalChrom.getQuick(i)) {
              counts[1][1] += 1;
            } else {
              counts[1][0] += 1;
            }
            if (gty.maternalChrom.getQuick(i)) {
              counts[1][1] += 1;
            } else {
              counts[1][0] += 1;
            }
          }
        } else if (indivList.get(k).getAffectedStatus() == 1) {
          if (gty.existence.getQuick(i)) {
            if (gty.paternalChrom.getQuick(i)) {
              counts[0][1] += 1;
            } else {
              counts[0][0] += 1;
            }
            if (gty.maternalChrom.getQuick(i)) {
              counts[0][1] += 1;
            } else {
              counts[0][0] += 1;
            }
          }
        }
      }
      p = csTest.chiSquareTest(counts);
      if (toDebug) {
        debugOut.write("SNP " + i + " :" + p + " :" + (counts[0][1] * 1.0) / (counts[0][0] + counts[0][1]) + " :" + (counts[1][1] * 1.0) / (counts[1][0] + counts[1][1]));
      }
      if (p <= alpha) {
        hits.putQuick(i, true);
      }

    }
    return hits;
  }

  public double[] allelicAssociationTest(List<Individual> indivList) throws Exception {
    if (indivList == null || indivList.isEmpty()) {
      return null;
    }
    int snpNum = indivList.get(0).markerGtySet.paternalChrom.size();
    int indivSize = indivList.size();

    long[][] counts = new long[2][2];
    int rowNum = counts.length;
    ChiSquareTest csTest = new ChiSquareTestImpl();
    double[] pValues = new double[snpNum];

    // table for contigency test
    //           0    1
    // control  counts[0][0]   counts[0][1]
    // case     counts[1][0]   counts[1][1]
    for (int i = 0; i < snpNum; i++) {
      for (int j = 0; j < rowNum; j++) {
        Arrays.fill(counts[j], 0);
      }
      for (int k = 0; k < indivSize; k++) {
        StatusGtySet gty = indivList.get(k).markerGtySet;
        if (indivList.get(k).getAffectedStatus() == 2) {
          if (gty.existence.getQuick(i)) {
            if (gty.paternalChrom.getQuick(i)) {
              counts[1][1] += 1;
            } else {
              counts[1][0] += 1;
            }
            if (gty.maternalChrom.getQuick(i)) {
              counts[1][1] += 1;
            } else {
              counts[1][0] += 1;
            }
          }
        } else if (indivList.get(k).getAffectedStatus() == 1) {
          if (gty.existence.getQuick(i)) {
            if (gty.paternalChrom.getQuick(i)) {
              counts[0][1] += 1;
            } else {
              counts[0][0] += 1;
            }
            if (gty.maternalChrom.getQuick(i)) {
              counts[0][1] += 1;
            } else {
              counts[0][0] += 1;
            }
          }
        }
      }
      pValues[i] = Probability.chiSquareComplemented(1, FisherExactTest.pearsonChiSquared22(counts));
      // pValues[i] = csTest.chiSquareTest(counts);
      /*
             if (pValues[snpIndex] <= 0.000001) {
             System.out.println("SNP " + snpIndex + " :" + pValues[snpIndex] + " :" + (counts[0][1] * 1.0) / (counts[0][0] + counts[0][1]) + " :" + (counts[1][1] * 1.0) / (counts[1][0] + counts[1][1]));
             }
             * 
       */

    }
    return pValues;
  }

  public double[] allelicAssociationZScoreTest(List<Individual> indivList) throws Exception {
    if (indivList == null || indivList.isEmpty()) {
      return null;
    }
    int snpNum = indivList.get(0).markerGtySet.paternalChrom.size();
    int indivSize = indivList.size();

    long[][] counts = new long[2][2];
    int rowNum = counts.length;

    double[] zScores = new double[snpNum];

    // table for contigency test
    //           0    1
    // control  counts[0][0]   counts[0][1]
    // case     counts[1][0]   counts[1][1]
    double ncase = 0, ncontrol = 0;
    for (int i = 0; i < snpNum; i++) {
      for (int j = 0; j < rowNum; j++) {
        Arrays.fill(counts[j], 0);
      }
      ncase = 0;
      ncontrol = 0;
      for (int k = 0; k < indivSize; k++) {
        StatusGtySet gty = indivList.get(k).markerGtySet;
        if (indivList.get(k).getAffectedStatus() == 2) {
          if (gty.existence.getQuick(i)) {
            if (gty.paternalChrom.getQuick(i)) {
              counts[1][1] += 1;
            } else {
              counts[1][0] += 1;
            }
            if (gty.maternalChrom.getQuick(i)) {
              counts[1][1] += 1;
            } else {
              counts[1][0] += 1;
            }
            ncase += 1;
          }
        } else if (indivList.get(k).getAffectedStatus() == 1) {
          if (gty.existence.getQuick(i)) {
            if (gty.paternalChrom.getQuick(i)) {
              counts[0][1] += 1;
            } else {
              counts[0][0] += 1;
            }
            if (gty.maternalChrom.getQuick(i)) {
              counts[0][1] += 1;
            } else {
              counts[0][0] += 1;
            }
            ncontrol += 1;
          }
        }
      }
      double pcase = (double) counts[1][1] / (counts[1][0] + counts[1][1]);
      double pcontrol = (double) counts[0][1] / (counts[0][0] + counts[0][1]);
      double pavg = (pcase + pcontrol) / 2;

      zScores[i] = Math.sqrt(ncase * ncontrol / (ncase + ncontrol)) * (pcase - pcontrol) / Math.sqrt(pavg * (1 - pavg));

      // zScores[i] = Probability.chiSquareComplemented(1, zScores[i] * zScores[i]);

      /*
             if (zScores[i] > 0) {
             zScores[i] = 2 * Probability.normal(-zScores[i]);
             } else {
             zScores[i] = 2 * Probability.normal(zScores[i]);
             }
       */
    }
    return zScores;
  }

  public double[] allelicAssociationTestMultTrait(List<Individual> indivList) throws Exception {
    if (indivList == null || indivList.isEmpty()) {
      return null;
    }
    int snpNum = indivList.get(0).markerGtySet.paternalChrom.size();
    int traitNum = indivList.get(0).getMainTrait().length;
    double[] pValues = new double[snpNum * traitNum];

    int indivSize = indivList.size();

    long[][] counts = new long[2][2];
    int rowNum = counts.length;
    ChiSquareTest csTest = new ChiSquareTestImpl();

    // table for contigency test
    //           0    1
    // control  counts[0][0]   counts[0][1]
    // case     counts[1][0]   counts[1][1]
    for (int t = 0; t < traitNum; t++) {
      for (int i = 0; i < snpNum; i++) {
        for (int j = 0; j < rowNum; j++) {
          Arrays.fill(counts[j], 0);
        }
        for (int k = 0; k < indivSize; k++) {
          StatusGtySet gty = indivList.get(k).markerGtySet;
          if (indivList.get(k).getMainTrait()[t] == 2) {
            if (gty.existence.getQuick(i)) {
              if (gty.paternalChrom.getQuick(i)) {
                counts[1][1] += 1;
              } else {
                counts[1][0] += 1;
              }
              if (gty.maternalChrom.getQuick(i)) {
                counts[1][1] += 1;
              } else {
                counts[1][0] += 1;
              }
            }
          } else if (indivList.get(k).getMainTrait()[t] == 1) {
            if (gty.existence.getQuick(i)) {
              if (gty.paternalChrom.getQuick(i)) {
                counts[0][1] += 1;
              } else {
                counts[0][0] += 1;
              }
              if (gty.maternalChrom.getQuick(i)) {
                counts[0][1] += 1;
              } else {
                counts[0][0] += 1;
              }
            }
          } else {
            System.err.println("Unknown trait " + indivList.get(k).getMainTrait()[t]);
          }
        }
        pValues[t * traitNum + i] = csTest.chiSquareTest(counts);
        /*
                 if (pValues[snpIndex] <= 0.000001) {
                 System.out.println("SNP " + snpIndex + " :" + pValues[snpIndex] + " :" + (counts[0][1] * 1.0) / (counts[0][0] + counts[0][1]) + " :" + (counts[1][1] * 1.0) / (counts[1][0] + counts[1][1]));
                 }
                 * 
         */

      }
    }
    return pValues;
  }

  public void correlationOfMultTrait(List<Individual> indivList) throws Exception {
    if (indivList == null || indivList.isEmpty()) {
      return;
    }

    int traitNum = indivList.get(0).getMainTrait().length;
    double[][] correl = new double[traitNum][traitNum];
    for (int t = 0; t < traitNum; t++) {
      Arrays.fill(correl[t], 0);
    }
    int indivSize = indivList.size();

    long[][] counts = new long[2][2];
    int rowNum = counts.length;
    ChiSquareTest csTest = new ChiSquareTestImpl();

    // table for contigency test
    //           0    1
    // control  counts[0][0]   counts[0][1]
    // case     counts[1][0]   counts[1][1]
    for (int t = 0; t < traitNum; t++) {
      for (int i = t; i < traitNum; i++) {
        for (int k = 0; k < indivSize; k++) {
          if (indivList.get(k).getMainTrait()[t] == 2 && indivList.get(k).getMainTrait()[i] == 2) {
            correl[t][i] += 1;
          }
        }
        correl[t][i] /= indivSize;
        System.out.print(correl[t][i] + "\t");
      }
      System.out.println();
    }
    double r = correl[0][1] - correl[0][0] * correl[1][1];
    r = r / (Math.sqrt(correl[0][0] * correl[1][1] * (1 - correl[0][0]) * (1 - correl[1][1])));
    System.out.println(r);
  }

  public double[] allelicAssociationTest1(List<Individual> indivList, List<AnnotSNP> snpList) throws Exception {
    if (indivList == null || indivList.isEmpty()) {
      return null;
    }
    int snpNum = snpList.size();
    int indivSize = indivList.size();

    long[][] counts = new long[2][2];
    int rowNum = counts.length;
    ChiSquareTest csTest = new ChiSquareTestImpl();
    double[] pValues = new double[snpNum];

    // table for contigency test
    //           0    1
    // control  counts[0][0]   counts[0][1]
    // case     counts[1][0]   counts[1][1]
    for (int i = 0; i < snpNum; i++) {
      int snpIndex = snpList.get(i).order;
      for (int j = 0; j < rowNum; j++) {
        Arrays.fill(counts[j], 0);
      }
      for (int k = 0; k < indivSize; k++) {
        StatusGtySet gty = indivList.get(k).markerGtySet;
        if (indivList.get(k).getAffectedStatus() == 2) {
          if (gty.existence.getQuick(snpIndex)) {
            if (gty.paternalChrom.getQuick(snpIndex)) {
              counts[1][1] += 1;
            } else {
              counts[1][0] += 1;
            }
            if (gty.maternalChrom.getQuick(snpIndex)) {
              counts[1][1] += 1;
            } else {
              counts[1][0] += 1;
            }
          }
        } else if (indivList.get(k).getAffectedStatus() == 1) {
          if (gty.existence.getQuick(snpIndex)) {
            if (gty.paternalChrom.getQuick(snpIndex)) {
              counts[0][1] += 1;
            } else {
              counts[0][0] += 1;
            }
            if (gty.maternalChrom.getQuick(snpIndex)) {
              counts[0][1] += 1;
            } else {
              counts[0][0] += 1;
            }
          }
        }
      }
      pValues[i] = csTest.chiSquareTest(counts);
      /*
             if (pValues[snpIndex] <= 0.000001) {
             System.out.println("SNP " + snpIndex + " :" + pValues[snpIndex] + " :" + (counts[0][1] * 1.0) / (counts[0][0] + counts[0][1]) + " :" + (counts[1][1] * 1.0) / (counts[1][0] + counts[1][1]));
             }
             * 
       */

    }
    return pValues;
  }

  public double[][] allelicAssociationTestQuantitativeTrait(List<Individual> indivList) throws Exception {
    if (indivList == null || indivList.isEmpty()) {
      return null;
    }
    int snpNum = indivList.get(0).markerGtySet.paternalChrom.size();
    int traitNum = indivList.get(0).getMainTrait().length;
    double[][] pValues = new double[snpNum * traitNum][2];
    int indivSize = indivList.size();
    double[] Y = new double[indivSize];
    double[] X = new double[indivSize];
    SimpleLinearRegression lr = new SimpleLinearRegression();

    for (int t = 0; t < traitNum; t++) {
      for (int k = 0; k < indivSize; k++) {
        Y[k] = indivList.get(k).getMainTrait()[t];
      }
      lr.setY(Y);

      for (int i = 0; i < snpNum; i++) {
        for (int k = 0; k < indivSize; k++) {
          StatusGtySet gty = indivList.get(k).markerGtySet;
          X[k] = 0;
          if (gty.existence.getQuick(i)) {
            if (gty.paternalChrom.getQuick(i)) {
              X[k] += 1;
            }
            if (gty.maternalChrom.getQuick(i)) {
              X[k] += 1;
            }
          }
        }
        lr.setX(X);
        lr.compute();
        // System.out.println(lr.getRoundedModel());
        pValues[t * snpNum + i][0] = lr.waldTestSlopeChi();
        // pValues[t * snpNum + i][0] = lr.getRSquared();

        pValues[t * snpNum + i][1] = lr.waldTestSlopeP();
        //convert to a  one sided p-value
        /*
                 if (pValues[t * snpNum + i][1] < 0.5) {
                 pValues[t * snpNum + i][1] = -Probability.normalInverse(pValues[t * snpNum + i][1]);
                 } else {
                 pValues[t * snpNum + i][1] = Probability.normalInverse(1 - pValues[t * snpNum + i][1]);
                 }*/

 /*
                 if (pValues[snpIndex] <= 0.000001) {
                 System.out.println("SNP " + snpIndex + " :" + pValues[snpIndex] + " :" + (counts[0][1] * 1.0) / (counts[0][0] + counts[0][1]) + " :" + (counts[1][1] * 1.0) / (counts[1][0] + counts[1][1]));
                 }
                 *
         */
      }
    }
    return pValues;
  }

  public double[][] allelicAssociationTestQuantitativeTraitForGCTA(List<Individual> indivList) throws Exception {
    if (indivList == null || indivList.isEmpty()) {
      return null;
    }
    int snpNum = indivList.get(0).markerGtySet.paternalChrom.size();
    int traitNum = indivList.get(0).getMainTrait().length;
    double[][] pValues = new double[snpNum * traitNum][5];
    int indivSize = indivList.size();
    double[] Y = new double[indivSize];
    double[] X = new double[indivSize];
    SimpleLinearRegression lr = new SimpleLinearRegression();
    double freq = 0;
    for (int t = 0; t < traitNum; t++) {
      for (int k = 0; k < indivSize; k++) {
        Y[k] = indivList.get(k).getMainTrait()[t];
      }
      lr.setY(Y);

      for (int i = 0; i < snpNum; i++) {
        freq = 0;
        for (int k = 0; k < indivSize; k++) {
          StatusGtySet gty = indivList.get(k).markerGtySet;
          X[k] = 0;
          if (gty.existence.getQuick(i)) {
            if (gty.paternalChrom.getQuick(i)) {
              X[k] += 1;
              freq += 1;
            }
            if (gty.maternalChrom.getQuick(i)) {
              X[k] += 1;
              freq += 1;
            }
          }
        }
        freq = freq / (2 * indivSize);
        lr.setX(X);
        lr.compute();
        pValues[t * snpNum + i][0] = freq;
        // System.out.println(lr.getRoundedModel());
        pValues[t * snpNum + i][1] = lr.getSlope();
        pValues[t * snpNum + i][2] = lr.getSlopeSE();
        // pValues[t * snpNum + i][0] = lr.getRSquared();

        pValues[t * snpNum + i][3] = lr.waldTestSlopeP();
        pValues[t * snpNum + i][4] = indivSize;
        //convert to a  one sided p-value
        /*
                 if (pValues[t * snpNum + i][1] < 0.5) {
                 pValues[t * snpNum + i][1] = -Probability.normalInverse(pValues[t * snpNum + i][1]);
                 } else {
                 pValues[t * snpNum + i][1] = Probability.normalInverse(1 - pValues[t * snpNum + i][1]);
                 }*/

 /*
                 if (pValues[snpIndex] <= 0.000001) {
                 System.out.println("SNP " + snpIndex + " :" + pValues[snpIndex] + " :" + (counts[0][1] * 1.0) / (counts[0][0] + counts[0][1]) + " :" + (counts[1][1] * 1.0) / (counts[1][0] + counts[1][1]));
                 }
                 *
         */
      }
    }
    return pValues;
  }

  public double[][] allelicAssociationTestQuantitativeTraitForGCTA(List<Individual> indivList, int traitIndex, IntArrayList snpIndex) throws Exception {
    if (indivList == null || indivList.isEmpty()) {
      return null;
    }
    int snpNum = snpIndex.size();
    int traitNum = 1;
    double[][] pValues = new double[snpNum * traitNum][5];
    int indivSize = indivList.size();
    double[] Y = new double[indivSize];
    double[] X = new double[indivSize];
    SimpleLinearRegression lr = new SimpleLinearRegression();
    double freq = 0;

    for (int k = 0; k < indivSize; k++) {
      Y[k] = indivList.get(k).getMainTrait()[traitIndex];
    }
    lr.setY(Y);

    for (int i = 0; i < snpNum; i++) {
      freq = 0;
      for (int k = 0; k < indivSize; k++) {
        StatusGtySet gty = indivList.get(k).markerGtySet;
        X[k] = 0;
        if (gty.existence.getQuick(snpIndex.getQuick(i))) {
          if (gty.paternalChrom.getQuick(snpIndex.getQuick(i))) {
            X[k] += 1;
            freq += 1;
          }
          if (gty.maternalChrom.getQuick(snpIndex.getQuick(i))) {
            X[k] += 1;
            freq += 1;
          }
        }
      }
      freq = freq / (2 * indivSize);
      lr.setX(X);
      lr.compute();
      pValues[i][0] = freq;
      // System.out.println(lr.getRoundedModel());
      pValues[i][1] = lr.getSlope();
      pValues[i][2] = lr.getSlopeSE();
      // pValues[t * snpNum + i][0] = lr.getRSquared();

      pValues[i][3] = lr.waldTestSlopeP();
      pValues[i][4] = indivSize;
      //convert to a  one sided p-value
      /*
                 if (pValues[t * snpNum + i][1] < 0.5) {
                 pValues[t * snpNum + i][1] = -Probability.normalInverse(pValues[t * snpNum + i][1]);
                 } else {
                 pValues[t * snpNum + i][1] = Probability.normalInverse(1 - pValues[t * snpNum + i][1]);
                 }*/

 /*
                 if (pValues[snpIndex] <= 0.000001) {
                 System.out.println("SNP " + snpIndex + " :" + pValues[snpIndex] + " :" + (counts[0][1] * 1.0) / (counts[0][0] + counts[0][1]) + " :" + (counts[1][1] * 1.0) / (counts[1][0] + counts[1][1]));
                 }
                 *
       */
    }

    return pValues;
  }

  public double[][] allelicAssociationTestQuantitativeTrait(List<Individual> indivList, IntArrayList traitIndex, IntArrayList snpIndex) throws Exception {
    if (indivList == null || indivList.isEmpty()) {
      return null;
    }
    int snpNum = indivList.get(0).markerGtySet.paternalChrom.size();
    int traitNum = indivList.get(0).getMainTrait().length;
    double[][] pValues = new double[snpNum * traitNum][2];
    int indivSize = indivList.size();
    double[] Y = new double[indivSize];
    double[] X = new double[indivSize];
    SimpleLinearRegression lr = new SimpleLinearRegression();

    for (int t = 0; t < traitNum; t++) {
      for (int k = 0; k < indivSize; k++) {
        Y[k] = indivList.get(k).getMainTrait()[t];
      }
      lr.setY(Y);

      for (int i = 0; i < snpNum; i++) {
        for (int k = 0; k < indivSize; k++) {
          StatusGtySet gty = indivList.get(k).markerGtySet;
          X[k] = 0;
          if (gty.existence.getQuick(i)) {
            if (gty.paternalChrom.getQuick(i)) {
              X[k] += 1;
            }
            if (gty.maternalChrom.getQuick(i)) {
              X[k] += 1;
            }
          }
        }
        lr.setX(X);
        lr.compute();
        // System.out.println(lr.getRoundedModel());
        pValues[t * traitNum + i][0] = lr.waldTestSlopeChi();
        pValues[t * traitNum + i][1] = lr.waldTestSlopeP();
        traitIndex.add(t);
        snpIndex.add(i);
        /*
                 if (pValues[snpIndex] <= 0.000001) {
                 System.out.println("SNP " + snpIndex + " :" + pValues[snpIndex] + " :" + (counts[0][1] * 1.0) / (counts[0][0] + counts[0][1]) + " :" + (counts[1][1] * 1.0) / (counts[1][0] + counts[1][1]));
                 }
                 *
         */

      }
    }
    return pValues;
  }

  public ArrayList<double[]> allelicAssociationTestQuantitativeTraitLR(List<Individual> indivList, int traitIndex, IntArrayList snpIndex) throws Exception {
    if (indivList == null || indivList.isEmpty()) {
      return null;
    }
    int snpNum = snpIndex.size();

    double[] pValues = new double[snpNum];
    double[] beta = new double[snpNum];
    double[] se = new double[snpNum];
    int indivSize = indivList.size();
    double[] Y = new double[indivSize];
    double[] X = new double[indivSize];
    SimpleLinearRegression lr = new SimpleLinearRegression();
    for (int k = 0; k < indivSize; k++) {
      Y[k] = indivList.get(k).getMainTrait()[traitIndex];
    }
    lr.setY(Y);

    for (int i = 0; i < snpNum; i++) {
      for (int k = 0; k < indivSize; k++) {
        StatusGtySet gty = indivList.get(k).markerGtySet;
        X[k] = 0;
        if (gty.existence.getQuick(snpIndex.getQuick(i))) {
          if (gty.paternalChrom.getQuick(snpIndex.getQuick(i))) {
            X[k] += 1;
          }
          if (gty.maternalChrom.getQuick(snpIndex.getQuick(i))) {
            X[k] += 1;
          }
        }
      }
      lr.setX(X);
      lr.compute();
      // System.out.println(lr.getRoundedModel());

      // pValues[i] = lr.waldTestSlopeZ();
      pValues[i] = lr.waldTestSlopeP();
      beta[i] = lr.getSlope();
      se[i] = lr.getSlopeSE();

    }

    ArrayList<double[]> lrPara = new ArrayList<>();
    lrPara.add(pValues);
    lrPara.add(beta);
    lrPara.add(se);

    return lrPara;
  }

  public double[] allelicAssociationTestQuantitativeTraitP(List<Individual> indivList, int traitIndex, IntArrayList snpIndex) throws Exception {
    if (indivList == null || indivList.isEmpty()) {
      return null;
    }
    int snpNum = snpIndex.size();

    double[] pValues = new double[snpNum];
    int indivSize = indivList.size();
    double[] Y = new double[indivSize];
    double[] X = new double[indivSize];
    SimpleLinearRegression lr = new SimpleLinearRegression();
    for (int k = 0; k < indivSize; k++) {
      Y[k] = indivList.get(k).getMainTrait()[traitIndex];
    }
    lr.setY(Y);

    for (int i = 0; i < snpNum; i++) {
      for (int k = 0; k < indivSize; k++) {
        StatusGtySet gty = indivList.get(k).markerGtySet;
        X[k] = 0;
        if (gty.existence.getQuick(snpIndex.getQuick(i))) {
          if (gty.paternalChrom.getQuick(snpIndex.getQuick(i))) {
            X[k] += 1;
          }
          if (gty.maternalChrom.getQuick(snpIndex.getQuick(i))) {
            X[k] += 1;
          }
        }
      }
      lr.setX(X);
      lr.compute();
      // System.out.println(lr.getRoundedModel());

     // pValues[i] = lr.waldTestSlopeZ();    
      pValues[i] = lr.waldTestSlopeP();    

    }

    return pValues;
  }

  public double[] allelicAssociationTestQuantitativeTraitZ(List<Individual> indivList, int traitIndex, IntArrayList snpIndex) throws Exception {
    if (indivList == null || indivList.isEmpty()) {
      return null;
    }
    int snpNum = snpIndex.size();

    double[] pValues = new double[snpNum];
    int indivSize = indivList.size();
    double[] Y = new double[indivSize];
    double[] X = new double[indivSize];
    SimpleLinearRegression lr = new SimpleLinearRegression();
    for (int k = 0; k < indivSize; k++) {
      Y[k] = indivList.get(k).getMainTrait()[traitIndex];
    }
    lr.setY(Y);

    for (int i = 0; i < snpNum; i++) {
      for (int k = 0; k < indivSize; k++) {
        StatusGtySet gty = indivList.get(k).markerGtySet;
        X[k] = 0;
        if (gty.existence.getQuick(snpIndex.getQuick(i))) {
          if (gty.paternalChrom.getQuick(snpIndex.getQuick(i))) {
            X[k] += 1;
          }
          if (gty.maternalChrom.getQuick(snpIndex.getQuick(i))) {
            X[k] += 1;
          }
        }
      }
      lr.setX(X);
      lr.compute();
      // System.out.println(lr.getRoundedModel());

      pValues[i] = lr.waldTestSlopeZ();    
     // pValues[i] = lr.waldTestSlopeP();    

    }

    return pValues;
  }

  public double[] allelicAssociationTestQuantitativeTraitChiSquare(List<Individual> indivList, int traitIndex, IntArrayList snpIndex) throws Exception {
    if (indivList == null || indivList.isEmpty()) {
      return null;
    }
    int snpNum = snpIndex.size();

    double[] pValues = new double[snpNum];
    int indivSize = indivList.size();
    double[] Y = new double[indivSize];
    double[] X = new double[indivSize];
    SimpleLinearRegression lr = new SimpleLinearRegression();
    for (int k = 0; k < indivSize; k++) {
      Y[k] = indivList.get(k).getMainTrait()[traitIndex];
    }
    lr.setY(Y);

    for (int i = 0; i < snpNum; i++) {
      for (int k = 0; k < indivSize; k++) {
        StatusGtySet gty = indivList.get(k).markerGtySet;
        X[k] = 0;
        if (gty.existence.getQuick(snpIndex.getQuick(i))) {
          if (gty.paternalChrom.getQuick(snpIndex.getQuick(i))) {
            X[k] += 1;
          }
          if (gty.maternalChrom.getQuick(snpIndex.getQuick(i))) {
            X[k] += 1;
          }
        }
      }
      lr.setX(X);
      lr.compute();
      // System.out.println(lr.getRoundedModel());

      pValues[i] = lr.waldTestSlopeZ();

    }

    return pValues;
  }

  public double[] allelicAssociationTestQuantitativeTrait(List<Individual> indivList, int traitIndex) throws Exception {
    if (indivList == null || indivList.isEmpty()) {
      return null;
    }
    int snpNum = indivList.get(0).markerGtySet.existence.size();

    double[] pValues = new double[snpNum];
    int indivSize = indivList.size();
    double[] Y = new double[indivSize];
    double[] X = new double[indivSize];
    SimpleLinearRegression lr = new SimpleLinearRegression();
    for (int k = 0; k < indivSize; k++) {
      Y[k] = indivList.get(k).getMainTrait()[traitIndex];
    }
    lr.setY(Y);

    for (int i = 0; i < snpNum; i++) {
      for (int k = 0; k < indivSize; k++) {
        StatusGtySet gty = indivList.get(k).markerGtySet;
        X[k] = 0;
        if (gty.existence.getQuick(i)) {
          if (gty.paternalChrom.getQuick(i)) {
            X[k] += 1;
          }
          if (gty.maternalChrom.getQuick(i)) {
            X[k] += 1;
          }
        }
      }
      lr.setX(X);
      lr.compute();
      // System.out.println(lr.getRoundedModel());

     // pValues[i] = lr.waldTestSlopeP();
      pValues[i] = lr.waldTestSlopeZ();

    }

    return pValues;
  }

  public double[][] allelicAssociationTest2QuantitativeTrait(List<Individual> indivList) throws Exception {
    if (indivList == null || indivList.isEmpty()) {
      return null;
    }
    int snpNum = indivList.get(0).markerGtySet.paternalChrom.size();
    int traitNum = indivList.get(0).getMainTrait().length;
    double[][] pValues = new double[snpNum * traitNum][2];
    int indivSize = indivList.size();
    double[] Y1 = new double[indivSize];
    double[] Y2 = new double[indivSize];
    double[] X = new double[indivSize];
    SimpleLinearRegression lr = new SimpleLinearRegression();

    for (int k = 0; k < indivSize; k++) {
      Y1[k] = indivList.get(k).getMainTrait()[0];
    }

    for (int k = 0; k < indivSize; k++) {
      Y2[k] = indivList.get(k).getMainTrait()[1];
    }

    for (int i = 0; i < snpNum; i++) {
      for (int k = 0; k < indivSize; k++) {
        StatusGtySet gty = indivList.get(k).markerGtySet;
        X[k] = 0;
        if (gty.existence.getQuick(i)) {
          if (gty.paternalChrom.getQuick(i)) {
            X[k] += 1;
          }
          if (gty.maternalChrom.getQuick(i)) {
            X[k] += 1;
          }
        }
      }
      lr.setY(Y1);
      lr.setX(X);
      lr.compute();
      pValues[i][0] = lr.waldTestSlopeZ();

      lr.setY(Y2);
      lr.compute();
      pValues[i][1] = lr.waldTestSlopeZ();
    }
    return pValues;
  }

  public double[] allelicAssociationTestSpecifySNP(List<Individual> indivList, List<AnnotSNP> snpList) throws Exception {
    if (indivList == null || indivList.isEmpty()) {
      return null;
    }
    int snpNum = snpList.size();
    int indivSize = indivList.size();

    long[][] counts = new long[2][2];
    int rowNum = counts.length;
    ChiSquareTest csTest = new ChiSquareTestImpl();
    double[] pValues = new double[snpNum];

    // table for contigency test
    //           0    1
    // control  counts[0][0]   counts[0][1]
    // case     counts[1][0]   counts[1][1]
    for (int i = 0; i < snpNum; i++) {
      for (int j = 0; j < rowNum; j++) {
        Arrays.fill(counts[j], 0);
      }
      AnnotSNP snp = snpList.get(i);
      int pos = snp.order;
      for (int k = 0; k < indivSize; k++) {
        StatusGtySet gty = indivList.get(k).markerGtySet;
        if (indivList.get(k).getAffectedStatus() == 2) {
          if (gty.existence.getQuick(pos)) {
            if (gty.paternalChrom.getQuick(pos)) {
              counts[1][1] += 1;
            } else {
              counts[1][0] += 1;
            }
            if (gty.maternalChrom.getQuick(pos)) {
              counts[1][1] += 1;
            } else {
              counts[1][0] += 1;
            }
          }
        } else if (indivList.get(k).getAffectedStatus() == 1) {
          if (gty.existence.getQuick(pos)) {
            if (gty.paternalChrom.getQuick(pos)) {
              counts[0][1] += 1;
            } else {
              counts[0][0] += 1;
            }
            if (gty.maternalChrom.getQuick(pos)) {
              counts[0][1] += 1;
            } else {
              counts[0][0] += 1;
            }
          }
        }
      }
      pValues[i] = csTest.chiSquareTest(counts);
      /*
             if (pValues[snpIndex] <= 0.000001) {
             System.out.println("SNP " + snpIndex + " :" + pValues[snpIndex] + " :" + (counts[0][1] * 1.0) / (counts[0][0] + counts[0][1]) + " :" + (counts[1][1] * 1.0) / (counts[1][0] + counts[1][1]));
             }
             *
       */

    }
    return pValues;
  }

  public double[] allelicAssociationLogTest(List<Individual> indivList) throws Exception {
    if (indivList == null || indivList.isEmpty()) {
      return null;
    }
    int snpNum = indivList.get(0).markerGtySet.paternalChrom.size();
    int indivSize = indivList.size();

    long[][] counts = new long[2][2];
    int rowNum = counts.length;
    ChiSquareTest csTest = new ChiSquareTestImpl();
    double[] pValues = new double[snpNum];

    // table for contigency test
    //           0    1
    // control  counts[0][0]   counts[0][1]
    // case     counts[1][0]   counts[1][1]
    for (int i = 0; i < snpNum; i++) {
      for (int j = 0; j < rowNum; j++) {
        Arrays.fill(counts[j], 0);
      }

      for (int k = 0; k < indivSize; k++) {
        StatusGtySet gty = indivList.get(k).markerGtySet;
        if (indivList.get(k).getAffectedStatus() == 2) {
          if (gty.existence.getQuick(i)) {
            if (gty.paternalChrom.getQuick(i)) {
              counts[1][1] += 1;
            } else {
              counts[1][0] += 1;
            }
            if (gty.maternalChrom.getQuick(i)) {
              counts[1][1] += 1;
            } else {
              counts[1][0] += 1;
            }
          }
        } else if (indivList.get(k).getAffectedStatus() == 1) {
          if (gty.existence.getQuick(i)) {
            if (gty.paternalChrom.getQuick(i)) {
              counts[0][1] += 1;
            } else {
              counts[0][0] += 1;
            }
            if (gty.maternalChrom.getQuick(i)) {
              counts[0][1] += 1;
            } else {
              counts[0][0] += 1;
            }
          }
        }
      }

      pValues[i] = csTest.chiSquareTest(counts);
      pValues[i] = -2 * Math.log(pValues[i]);
    }
    return pValues;
  }

  public double allelicAssociationTestInTotal(List<Individual> indivList) throws Exception {

    if (indivList == null || indivList.isEmpty()) {
      return 0;
    }
    int snpNum = indivList.get(0).markerGtySet.paternalChrom.size();
    int indivSize = indivList.size();
    double totalChi = 0;
    long[][] counts = new long[2][2];
    int rowNum = counts.length;
    ChiSquareTest csTest = new ChiSquareTestImpl();
    double[] pValues = new double[snpNum];

    // table for contigency test
    //           0    1
    // control  counts[0][0]   counts[0][1]
    // case     counts[1][0]   counts[1][1]
    for (int i = 0; i < snpNum; i++) {
      for (int j = 0; j < rowNum; j++) {
        Arrays.fill(counts[j], 0);
      }
      for (int k = 0; k < indivSize; k++) {
        StatusGtySet gty = indivList.get(k).markerGtySet;
        if (indivList.get(k).getAffectedStatus() == 2) {
          if (gty.existence.getQuick(i)) {
            if (gty.paternalChrom.getQuick(i)) {
              counts[1][1] += 1;
            } else {
              counts[1][0] += 1;
            }
            if (gty.maternalChrom.getQuick(i)) {
              counts[1][1] += 1;
            } else {
              counts[1][0] += 1;
            }
          }
        } else if (indivList.get(k).getAffectedStatus() == 1) {
          if (gty.existence.getQuick(i)) {
            if (gty.paternalChrom.getQuick(i)) {
              counts[0][1] += 1;
            } else {
              counts[0][0] += 1;
            }
            if (gty.maternalChrom.getQuick(i)) {
              counts[0][1] += 1;
            } else {
              counts[0][0] += 1;
            }
          }
        }
      }
      totalChi += csTest.chiSquare(counts);
    }
    return totalChi;
  }

  public double allelicAssociationTestInTotal1(List<Individual> indivList, List<AnnotSNP> snpList) throws Exception {
    if (indivList == null || indivList.isEmpty()) {
      return 0;
    }
    int snpNum = snpList.size();
    int indivSize = indivList.size();
    double totalChi = 0;
    long[][] counts = new long[2][2];
    int rowNum = counts.length;
    ChiSquareTest csTest = new ChiSquareTestImpl();

    // table for contigency test
    //           0    1
    // control  counts[0][0]   counts[0][1]
    // case     counts[1][0]   counts[1][1]
    for (int i = 0; i < snpNum; i++) {
      int snpIndex = snpList.get(i).order;
      for (int j = 0; j < rowNum; j++) {
        Arrays.fill(counts[j], 0);
      }
      for (int k = 0; k < indivSize; k++) {
        StatusGtySet gty = indivList.get(k).markerGtySet;
        if (indivList.get(k).getAffectedStatus() == 2) {
          if (gty.existence.getQuick(snpIndex)) {
            if (gty.paternalChrom.getQuick(snpIndex)) {
              counts[1][1] += 1;
            } else {
              counts[1][0] += 1;
            }
            if (gty.maternalChrom.getQuick(snpIndex)) {
              counts[1][1] += 1;
            } else {
              counts[1][0] += 1;
            }
          }
        } else if (indivList.get(k).getAffectedStatus() == 1) {
          if (gty.existence.getQuick(snpIndex)) {
            if (gty.paternalChrom.getQuick(snpIndex)) {
              counts[0][1] += 1;
            } else {
              counts[0][0] += 1;
            }
            if (gty.maternalChrom.getQuick(snpIndex)) {
              counts[0][1] += 1;
            } else {
              counts[0][0] += 1;
            }
          }
        }
      }
      totalChi += csTest.chiSquare(counts);
    }
    return totalChi;
  }

  public double[] allelicAssociationTestChisquare(List<Individual> indivList) throws Exception {
    if (indivList == null || indivList.isEmpty()) {
      return null;
    }
    int snpNum = indivList.get(0).markerGtySet.paternalChrom.size();
    int indivSize = indivList.size();

    long[][] counts = new long[2][2];
    int rowNum = counts.length;
    ChiSquareTest csTest = new ChiSquareTestImpl();
    double[] chisquares = new double[snpNum];

    // table for contigency test
    //           0    1
    // control  counts[0][0]   counts[0][1]
    // case     counts[1][0]   counts[1][1]
    for (int i = 0; i < snpNum; i++) {
      for (int j = 0; j < rowNum; j++) {
        Arrays.fill(counts[j], 0);
      }
      for (int k = 0; k < indivSize; k++) {
        StatusGtySet gty = indivList.get(k).markerGtySet;
        if (indivList.get(k).getAffectedStatus() == 2) {
          if (gty.existence.getQuick(i)) {
            if (gty.paternalChrom.getQuick(i)) {
              counts[1][1] += 1;
            } else {
              counts[1][0] += 1;
            }
            if (gty.maternalChrom.getQuick(i)) {
              counts[1][1] += 1;
            } else {
              counts[1][0] += 1;
            }
          }
        } else if (indivList.get(k).getAffectedStatus() == 1) {
          if (gty.existence.getQuick(i)) {
            if (gty.paternalChrom.getQuick(i)) {
              counts[0][1] += 1;
            } else {
              counts[0][0] += 1;
            }
            if (gty.maternalChrom.getQuick(i)) {
              counts[0][1] += 1;
            } else {
              counts[0][0] += 1;
            }
          }
        }
      }
      chisquares[i] = csTest.chiSquare(counts);
      //chisquares[i] = Math.sqrt(csTest.chiSquare(counts));
    }
    return chisquares;
  }

  public double[] allelicAssociationTestWithAnnote(List<Individual> indivList, List<AnnotSNP> snpList) throws Exception {

    if (indivList == null || indivList.size() == 0) {
      return null;
    }
    if (toDebug) {
      debugOut = new BufferedWriter(new FileWriter("debug.txt"));
    }
    int snpNum = snpList.size();
    int indivSize = indivList.size();

    long[][] counts = new long[2][2];
    int rowNum = counts.length;
    ChiSquareTest csTest = new ChiSquareTestImpl();
    double[] pValues = new double[snpNum];

    // table for contigency test
    //           0    1
    // control  counts[0][0]   counts[0][1]
    // case     counts[1][0]   counts[1][1]
    for (int i = 0; i < snpNum; i++) {
      for (int j = 0; j < rowNum; j++) {
        Arrays.fill(counts[j], 0);
      }
      for (int k = 0; k < indivSize; k++) {
        StatusGtySet gty = indivList.get(k).markerGtySet;
        if (indivList.get(k).getAffectedStatus() == 2) {
          if (gty.existence.getQuick(i)) {
            if (gty.paternalChrom.getQuick(i)) {
              counts[1][1] += 1;
            } else {
              counts[1][0] += 1;
            }
            if (gty.maternalChrom.getQuick(i)) {
              counts[1][1] += 1;
            } else {
              counts[1][0] += 1;
            }
          }
        } else if (indivList.get(k).getAffectedStatus() == 1) {
          if (gty.existence.getQuick(i)) {
            if (gty.paternalChrom.getQuick(i)) {
              counts[0][1] += 1;
            } else {
              counts[0][0] += 1;
            }
            if (gty.maternalChrom.getQuick(i)) {
              counts[0][1] += 1;
            } else {
              counts[0][0] += 1;
            }
          }
        }
      }
      pValues[i] = csTest.chiSquareTest(counts);

      if (toDebug) {
        if (pValues[i] < 0.000001) {
          debugOut.write(snpList.get(i).getRSID() + " :" + pValues[i] + " :" + (counts[0][1] * 1.0) / (counts[0][0] + counts[0][1]) + " :" + (counts[1][1] * 1.0) / (counts[1][0] + counts[1][1]));
          debugOut.write(snpList.get(i).getRSID() + " :" + pValues[i] + " :" + (counts[0][0]) + " " + (counts[1][0]) + " :" + (counts[0][1]) + " " + (counts[1][1]));
          debugOut.newLine();
        }
      }
    }
    debugOut.close();
    return pValues;
  }

  public double[] allelicAssociationTest(List<Individual> genotypeIndivList, List<Individual> matchedPhenotypeIndiviList) throws Exception {

    if (genotypeIndivList == null || genotypeIndivList.size() == 0) {
      return null;
    }
    int snpNum = genotypeIndivList.get(0).markerGtySet.paternalChrom.size();
    int indivSize = genotypeIndivList.size();

    long[][] counts = new long[2][2];
    int rowNum = counts.length;
    ChiSquareTest csTest = new ChiSquareTestImpl();
    double[] pValues = new double[snpNum];

    // table for contigency test
    //           0    1
    // control  counts[0][0]   counts[0][1]
    // case     counts[1][0]   counts[1][1]
    for (int i = 0; i < snpNum; i++) {
      for (int j = 0; j < rowNum; j++) {
        Arrays.fill(counts[j], 0);
      }
      for (int k = 0; k < indivSize; k++) {
        StatusGtySet gty = genotypeIndivList.get(k).markerGtySet;
        if (matchedPhenotypeIndiviList.get(k).getAffectedStatus() == 2) {
          if (gty.existence.getQuick(i)) {
            if (gty.paternalChrom.getQuick(i)) {
              counts[1][1] += 1;
            } else {
              counts[1][0] += 1;
            }
            if (gty.maternalChrom.getQuick(i)) {
              counts[1][1] += 1;
            } else {
              counts[1][0] += 1;
            }
          }
        } else if (matchedPhenotypeIndiviList.get(k).getAffectedStatus() == 1) {
          if (gty.existence.getQuick(i)) {
            if (gty.paternalChrom.getQuick(i)) {
              counts[0][1] += 1;
            } else {
              counts[0][0] += 1;
            }
            if (gty.maternalChrom.getQuick(i)) {
              counts[0][1] += 1;
            } else {
              counts[0][0] += 1;
            }
          }
        }
      }
      pValues[i] = csTest.chiSquareTest(counts);
    }
    return pValues;
  }

  public BitVector genotypicAssociationTest(List<Individual> allIndiv, double alpha) throws Exception {
    if (allIndiv == null || allIndiv.size() == 0) {
      return null;
    }
    int snpNum = allIndiv.get(0).markerGtySet.paternalChrom.size();
    int indivSize = allIndiv.size();

    long[][] counts = new long[2][3];
    int rowNum = counts.length;
    ChiSquareTest csTest = new ChiSquareTestImpl();
    BitVector hits = new BitVector(snpNum);
    hits.replaceFromToWith(0, snpNum - 1, false);
    double p;
    // table for contigency test
    //                 00           01             11
    // control  counts[0][0]   counts[0][1]     counts[0][2]
    // case     counts[1][0]   counts[1][1]     counts[1][2]

    for (int i = 0; i < snpNum; i++) {
      for (int j = 0; j < rowNum; j++) {
        Arrays.fill(counts[j], 0);
      }
      for (int k = 0; k < indivSize; k++) {
        StatusGtySet gty = allIndiv.get(k).markerGtySet;
        if (allIndiv.get(k).getAffectedStatus() == 2) {
          if (gty.existence.getQuick(i)) {
            if (gty.paternalChrom.getQuick(i) || gty.maternalChrom.getQuick(i)) {
              counts[1][1] += 1;
            } else {
              counts[1][0] += 1;
            }
          }
        } else if (allIndiv.get(k).getAffectedStatus() == 1) {
          if (gty.existence.getQuick(i)) {
            if (gty.paternalChrom.getQuick(i) || gty.maternalChrom.getQuick(i)) {
              counts[0][1] += 1;
            } else {
              counts[0][0] += 1;
            }
          }
        }
      }
      p = csTest.chiSquareTest(counts);
      //System.debugOut.println("SNP " + snpIndex + " :" + p + " :" + (counts[0][1] * 1.0) / (counts[0][0] + counts[0][1]) + " :" + (counts[1][1] * 1.0) / (counts[1][0] + counts[1][1]));
      if (p <= alpha) {
        hits.putQuick(i, true);
      }
    }
    return hits;
  }

  public void allelicAssociationLogisticTest(List<Individual> indivList,
      List<AssociationSNP> snpProfile) throws Exception {
    if (snpProfile == null) {
      snpProfile = new ArrayList<AssociationSNP>();
    } else {
      snpProfile.clear();
    }

    if (indivList == null || indivList.size() == 0) {
      System.err.println("No data!");
      return;
    }
    if (toDebug) {
      debugOut = new BufferedWriter(new FileWriter("debug.txt"));
    }
    int indivSize = indivList.size();

    //because all variables are the same exception for the genotypes
    // we could keep them while changing the genotypes
    double[] Y = new double[indivSize];
    List<double[]> tmpFullX = new ArrayList<double[]>();
    int traitN = indivList.get(0).getTraits().size();
    ////  intercept is in the first column
    // for 1 SNPs
    int vacantColNum = 1;
    int nP = 1 + traitN + vacantColNum;
    Individual indiv;
    // assigne scores for individuals
    for (int k = 0; k < indivSize; k++) {
      indiv = indivList.get(k);
      if (indiv.getAffectedStatus() == 2) {
        Y[k] = 1;
      } else if (indiv.getAffectedStatus() == 1) {
        Y[k] = 0;
      } else {
        System.err.println("Invalid disease status for inidividual " + indiv.getIndividualID() + " at family " + indiv.getFamilyID());
      }

      double[] fullXRow = new double[nP];
      fullXRow[0] = 1.0;
      for (int j = 0; j < traitN; j++) {
        fullXRow[j + 1] = Double.parseDouble(indiv.getTraits().get(j));
      }
      tmpFullX.add(fullXRow);
      if (toDebug) {
        debugOut.write(String.valueOf(indiv.getAffectedStatus()));
        debugOut.write("\n");
      }
    }

    LogisticRegression logisticRegresion = new LogisticRegression();
    int gtyColNumIndex = tmpFullX.get(0).length - 1;
    //logistic regression        
    int nonmissingNum = 0;
    int validIndex = 0;
    int snpNum = indivList.get(0).markerGtySet.paternalChrom.size();

    double p = 1;
    double gtyScore = 0;
    double[] localY = null;
    double[][] localX = null;
    for (int snpIndex = 0; snpIndex < snpNum; snpIndex++) {
      //count non missing number
      nonmissingNum = 0;
      for (int k = 0; k < indivSize; k++) {
        StatusGtySet gty = indivList.get(k).markerGtySet;
        if (gty.existence.getQuick(snpIndex)) {
          nonmissingNum++;
        }
      }

      // System.out.println(nonmissingNum);
      localY = new double[nonmissingNum];
      localX = new double[nonmissingNum][0];
      validIndex = 0;

      for (int k = 0; k < indivSize; k++) {
        StatusGtySet gty = indivList.get(k).markerGtySet;
        gtyScore = 0;
        if (gty.existence.getQuick(snpIndex)) {
          localY[validIndex] = Y[k];
          //localX[validIndex] = Arrays.copyOf(fullX[k], gtyColNumIndex+1);
          localX[validIndex] = tmpFullX.get(k);

          if (gty.paternalChrom.getQuick(snpIndex)) {
            gtyScore += 1.0;
          }
          if (gty.maternalChrom.getQuick(snpIndex)) {
            gtyScore += 1.0;
          }
          localX[validIndex][gtyColNumIndex] = gtyScore;
          validIndex++;
        }
      }
      //System.out.println(validIndex);

      logisticRegresion.setY(localY);
      logisticRegresion.setX(localX);
      if (logisticRegresion.fitLM()) {
        p = logisticRegresion.getCoefPValue(gtyColNumIndex);
        //System.out.println(logisticRegresion.toString());
        //System.out.println(p);

        AssociationSNP snp = new AssociationSNP("", -1, snpIndex);
        double od = logisticRegresion.getCoef(gtyColNumIndex);
        if (od < 0) {
          od = -od;
          snp.setHighRiskAllele(false);
        }

        snp.orderInList = (snpIndex);
        //snp.setLogOddsRatio(Math.log(od));
        snp.setLogOddsRatio(od);
        snp.setP(p);
        snpProfile.add(snp);
      }
    }
  }

  public double[] allelicAssociationLogisticTestSimple(List<Individual> indivList) throws Exception {
    if (indivList == null || indivList.isEmpty()) {
      System.err.println("No data!");
      return null;
    }
    if (toDebug) {
      debugOut = new BufferedWriter(new FileWriter("debug.txt"));
    }
    int indivSize = indivList.size();

    //because all variables are the same exception for the genotypes
    // we could keep them while changing the genotypes
    double[] Y = new double[indivSize];
    List<double[]> tmpFullX = new ArrayList<double[]>();
    int traitN = indivList.get(0).getTraits().size();
    ////  intercept is in the first column
    // for 1 SNPs
    int vacantColNum = 1;
    int nP = 1 + traitN + vacantColNum;
    Individual indiv;
    // assigne scores for individuals
    for (int k = 0; k < indivSize; k++) {
      indiv = indivList.get(k);
      if (indiv.getAffectedStatus() == 2) {
        Y[k] = 1;
      } else if (indiv.getAffectedStatus() == 1) {
        Y[k] = 0;
      } else {
        System.err.println("Invalid disease status for inidividual " + indiv.getIndividualID() + " at family " + indiv.getFamilyID());
      }

      double[] fullXRow = new double[nP];
      fullXRow[0] = 1.0;
      for (int j = 0; j < traitN; j++) {
        fullXRow[j + 1] = Double.parseDouble(indiv.getTraits().get(j));
      }
      tmpFullX.add(fullXRow);
      if (toDebug) {
        debugOut.write(String.valueOf(indiv.getAffectedStatus()));
        debugOut.write("\n");
      }
    }

    LogisticRegression logisticRegresion = new LogisticRegression();
    int gtyColNumIndex = tmpFullX.get(0).length - 1;
    //logistic regression
    int nonmissingNum = 0;
    int validIndex = 0;
    int snpNum = indivList.get(0).markerGtySet.paternalChrom.size();
    double[] pValues = new double[snpNum];
    Arrays.fill(pValues, 1);

    double p = 1;
    double gtyScore = 0;
    double[] localY = null;
    double[][] localX = null;
    for (int snpIndex = 0; snpIndex < snpNum; snpIndex++) {
      //count non missing number
      nonmissingNum = 0;
      for (int k = 0; k < indivSize; k++) {
        StatusGtySet gty = indivList.get(k).markerGtySet;
        if (gty.existence.getQuick(snpIndex)) {
          nonmissingNum++;
        }
      }

      // System.out.println(nonmissingNum);
      localY = new double[nonmissingNum];
      localX = new double[nonmissingNum][0];
      validIndex = 0;

      for (int k = 0; k < indivSize; k++) {
        StatusGtySet gty = indivList.get(k).markerGtySet;
        gtyScore = 0;
        if (gty.existence.getQuick(snpIndex)) {
          localY[validIndex] = Y[k];
          //localX[validIndex] = Arrays.copyOf(fullX[k], gtyColNumIndex+1);
          localX[validIndex] = tmpFullX.get(k);

          if (gty.paternalChrom.getQuick(snpIndex)) {
            gtyScore += 1.0;
          }
          if (gty.maternalChrom.getQuick(snpIndex)) {
            gtyScore += 1.0;
          }
          localX[validIndex][gtyColNumIndex] = gtyScore;
          validIndex++;
        }
      }
      //System.out.println(validIndex);

      logisticRegresion.setY(localY);
      logisticRegresion.setX(localX);
      if (logisticRegresion.fitLM()) {
        p = logisticRegresion.getCoefPValue(gtyColNumIndex);
        // System.out.println(logisticRegresion.toString());
        //System.out.println(p);
        pValues[snpIndex] = p;
      }
    }
    return pValues;
  }

  public double[][] allelicAssociationLogisticTestSimpleMainTraits(List<Individual> indivList) throws Exception {
    if (indivList == null || indivList.isEmpty()) {
      System.err.println("No data!");
      return null;
    }
    if (toDebug) {
      debugOut = new BufferedWriter(new FileWriter("debug.txt"));
    }
    int indivSize = indivList.size();
    int traitN = indivList.get(0).getMainTrait().length;
    //because all variables are the same exception for the genotypes
    // we could keep them while changing the genotypes
    double[][] Y = new double[traitN][indivSize];

    ////  intercept is in the first column
    // for 1 SNPs
    int vacantColNum = 1;
    int nP = 1 + traitN + vacantColNum;
    Individual indiv;
    // assigne scores for individuals
    for (int k = 0; k < indivSize; k++) {
      indiv = indivList.get(k);
      for (int i = 0; i < traitN; i++) {
        Y[i][k] = indiv.getMainTrait()[i] - 1;
      }
    }

    LogisticRegression logisticRegresion = new LogisticRegression();

    int snpNum = indivList.get(0).markerGtySet.paternalChrom.size();
    double[][] zScores = new double[snpNum][traitN];

    double p = 1;
    double gtyScore = 0;

    double[][] localX = new double[indivSize][2];
    for (int snpIndex = 0; snpIndex < snpNum; snpIndex++) {

      for (int k = 0; k < indivSize; k++) {
        StatusGtySet gty = indivList.get(k).markerGtySet;
        gtyScore = 0;

        if (gty.paternalChrom.getQuick(snpIndex)) {
          gtyScore += 1.0;
        }
        if (gty.maternalChrom.getQuick(snpIndex)) {
          gtyScore += 1.0;
        }
        localX[k][0] = 1;
        localX[k][1] = gtyScore;
      }
      logisticRegresion.setX(localX);
      //System.out.println(validIndex);
      for (int i = 0; i < traitN; i++) {
        logisticRegresion.setY(Y[i]);
        if (logisticRegresion.fitLM()) {
          p = logisticRegresion.getCoefZ(1);
          //System.out.println(logisticRegresion.toString());
          //System.out.println(p);
          zScores[snpIndex][i] = p;
        }
      }

    }
    return zScores;
  }

  public double[] allelicAssociationLogisticTestSimple(List<Individual> indivList, List<AnnotSNP> snpList) throws Exception {
    if (indivList == null || indivList.size() == 0) {
      System.err.println("No data!");
      return null;
    }
    if (toDebug) {
      debugOut = new BufferedWriter(new FileWriter("debug.txt"));
    }
    int indivSize = indivList.size();

    //because all variables are the same exception for the genotypes
    // we could keep them while changing the genotypes
    double[] Y = new double[indivSize];
    List<double[]> tmpFullX = new ArrayList<double[]>();
    int traitN = indivList.get(0).getTraits().size();
    ////  intercept is in the first column
    // for 1 SNPs
    int vacantColNum = 1;
    int nP = 1 + traitN + vacantColNum;
    Individual indiv;
    // assigne scores for individuals
    for (int k = 0; k < indivSize; k++) {
      indiv = indivList.get(k);
      if (indiv.getAffectedStatus() == 2) {
        Y[k] = 1;
      } else if (indiv.getAffectedStatus() == 1) {
        Y[k] = 0;
      } else {
        System.err.println("Invalid disease status for inidividual " + indiv.getIndividualID() + " at family " + indiv.getFamilyID());
      }

      double[] fullXRow = new double[nP];
      fullXRow[0] = 1.0;
      for (int j = 0; j < traitN; j++) {
        fullXRow[j + 1] = Double.parseDouble(indiv.getTraits().get(j));
      }
      tmpFullX.add(fullXRow);
      if (toDebug) {
        debugOut.write(String.valueOf(indiv.getAffectedStatus()));
        debugOut.write("\n");
      }
    }

    LogisticRegression logisticRegresion = new LogisticRegression();
    int gtyColNumIndex = tmpFullX.get(0).length - 1;
    //logistic regression
    int nonmissingNum = 0;
    int validIndex = 0;
    int snpNum = snpList.size();
    double[] pValues = new double[snpNum];
    Arrays.fill(pValues, 1);

    double p = 1;
    double gtyScore = 0;
    double[] localY = null;
    double[][] localX = null;
    int snpIndex = 0;
    for (int i = 0; i < snpNum; i++) {
      snpIndex = snpList.get(i).order;
      //count non missing number
      nonmissingNum = 0;
      for (int k = 0; k < indivSize; k++) {
        StatusGtySet gty = indivList.get(k).markerGtySet;
        if (gty.existence.getQuick(snpIndex)) {
          nonmissingNum++;
        }
      }
      // System.out.println(nonmissingNum);
      localY = new double[nonmissingNum];
      localX = new double[nonmissingNum][0];
      validIndex = 0;

      for (int k = 0; k < indivSize; k++) {
        StatusGtySet gty = indivList.get(k).markerGtySet;
        gtyScore = 0;
        if (gty.existence.getQuick(snpIndex)) {
          localY[validIndex] = Y[k];
          //localX[validIndex] = Arrays.copyOf(fullX[k], gtyColNumIndex+1);
          localX[validIndex] = tmpFullX.get(k);

          if (gty.paternalChrom.getQuick(snpIndex)) {
            gtyScore += 1.0;
          }
          if (gty.maternalChrom.getQuick(snpIndex)) {
            gtyScore += 1.0;
          }
          localX[validIndex][gtyColNumIndex] = gtyScore;
          validIndex++;
        }
      }
      //System.out.println(validIndex);

      logisticRegresion.setY(localY);
      logisticRegresion.setX(localX);
      if (logisticRegresion.fitLM()) {
        p = logisticRegresion.getCoefPValue(gtyColNumIndex);
        //System.out.println(logisticRegresion.toString());
        //System.out.println(p);
        pValues[i] = p;
      }
    }
    return pValues;
  }

  public void readCovariables(String pathName, List<Individual> indivList) throws Exception {
    File path = new File(pathName);
    BufferedReader br = new BufferedReader(new FileReader(path));
    String line = null;

    String delmilit = ", \t";
    Map<String, String> indivMap = new HashMap<String, String>();
    int lineCounter = 0;
    br.readLine(); //skip title line
    while ((line = br.readLine()) != null) {
      if (line.trim().length() == 0) {
        continue;
      }
      lineCounter++;
      StringTokenizer tokenizer = new StringTokenizer(line, delmilit);

      String famID = tokenizer.nextToken();
      String indivID = tokenizer.nextToken();
      indivMap.put(famID + indivID, line);

    }
    br.close();

    int size = indivList.size();
    //map the index of the individuals in the list
    for (int i = 0; i < size; i++) {
      Individual indiv = indivList.get(i);

      String fileLine = indivMap.get(indiv.getFamilyID() + indiv.getIndividualID());
      if (fileLine == null) {
        System.err.println("Failed to find individual " + indiv.getIndividualID() + " at family " + indiv.getFamilyID());
        continue;
      }

      StringTokenizer tokenizer = new StringTokenizer(fileLine, delmilit);
      String famID = tokenizer.nextToken();
      String indivID = tokenizer.nextToken();

      //System.out.println(indiv.getFamilyID() + " -> " + indiv.getIndividualID());
      while (tokenizer.hasMoreTokens()) {
        indiv.addTrait(tokenizer.nextToken());
      }

      //System.out.println(indiv.getFamilyID() +" -> "+ indiv.getIndividualID());
    }
  }
}
