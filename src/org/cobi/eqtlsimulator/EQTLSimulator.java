/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.eqtlsimulator;

import cern.colt.list.BooleanArrayList;
import cern.colt.list.DoubleArrayList;
import cern.colt.list.IntArrayList;
import cern.colt.map.OpenIntIntHashMap;
import cern.colt.matrix.DoubleMatrix2D;
import cern.jet.stat.Descriptive;
import cern.jet.stat.Gamma;

import java.io.*;
import java.lang.reflect.Array;
import java.util.*;

import jsc.independentsamples.MannWhitneyTest;
import jsc.tests.H1;
import org.apache.commons.math.distribution.ChiSquaredDistribution;
import org.apache.commons.math.random.MersenneTwister;
import org.apache.commons.math.random.RandomData;
import org.apache.commons.math.random.RandomDataImpl;

import org.cobi.genetsimulator.controller.GeneInforProcessor;
import org.cobi.genetsimulator.controller.GenetAssociationAnalyzerP;
import org.cobi.genetsimulator.controller.GenotypeQC;
import org.cobi.genetsimulator.controller.PValuePainter;
import org.cobi.genetsimulator.entity.AnnotSNP;
import org.cobi.genetsimulator.entity.GenotypeBasedLDSparseMatrix;
import org.cobi.genetsimulator.entity.Individual;
import org.cobi.genetsimulator.entity.PlinkDataset;
import org.cobi.genetsimulator.entity.Population;
import org.cobi.genetsimulator.entity.SNPPosiComparator;
import org.cobi.genetsimulator.entity.StatusGtySet;
import org.cobi.util.stat.LogisticRegression;
import org.cobi.util.stat.MultipleTestingMethod;
import org.cobi.util.stat.SimpleLinearRegression;
import org.ejml.data.DenseMatrix64F;
import org.rosuda.REngine.Rserve.RConnection;
import org.cobi.eqtlsimulator.DoubleArrayListComparatorR;
import org.cobi.eqtlsimulator.DoubleArrayListComparator;

/**
 *
 * @author limx54
 */
public class EQTLSimulator {

  /**
   * @param args the command line arguments
   */
  public static void main(String[] args) {
    RandomParam param = new RandomParam("1",1,2, 1000, 0.005);
    // TODO code application logic here
    Options option = new Options();
    EQTLSimulator simulator = new EQTLSimulator();
    List<List<AnnotSNP>> geneFullSnpList = new ArrayList<List<AnnotSNP>>();
    List<DoubleMatrix2D> ldMatrix = new ArrayList<DoubleMatrix2D>();

    List<List<Individual>> geneIndList = new ArrayList<List<Individual>>();
    List<List<AnnotSNP>> geneIndependentSnpList = new ArrayList<List<AnnotSNP>>();

    try {
      option.readOptions("F:\\java\\SimulationTestCode\\param.txt");
      option.parseOptions();
      int simulationTime = option.permutationTime;
      simulator.loadGenotypes(option, geneFullSnpList, geneIndependentSnpList, geneIndList, ldMatrix);

      int geneNum = geneFullSnpList.size();

      DoubleArrayList pvalues = new DoubleArrayList();
      BufferedWriter debugOut = new BufferedWriter(new FileWriter("debug.txt"));
      //BufferedWriter results = new BufferedWriter(new FileWriter("results.txt"));
      String[] type = {"positive", "negative", "random"};
      PrintStream ps = new PrintStream(new FileOutputStream("results.txt"));
      for (int k = 0; k < type.length; k++) {
        double errors = 0.0;
//        List<int[]> nmlPowerCounts = new ArrayList<int[]>();
//        List<int[]> orgPowerCounts = new ArrayList<int[]>();
//        //List<int[]> fmCounts = new ArrayList<int[]>();
//        for (int i = 0; i < geneNum; i++) {
//          nmlPowerCounts.add(new int[4]);
//          orgPowerCounts.add(new int[2]);
//          //fmCounts.add(new int[4]);
//        }
        System.setOut(new PrintStream(new FileOutputStream(FileDescriptor.out)));
        simulationTime = 3;
        for (int i = 0; i < simulationTime; i++) {
          //simulator.simulateExpression(option, geneFullSnpList, geneIndependentSnpList, geneIndList, ldMatrix);
          errors += simulator.expressionQTLFineMapping(option, geneFullSnpList, geneIndependentSnpList, geneIndList, type[k]);
          //simulator.simpleExpressionQTLRelationCheck(option, geneFullSnpList, geneIndependentSnpList, geneIndList, ldMatrix, nmlPowerCounts, orgPowerCounts, pvalues, debugOut);
          //to do some analysis with the simulated phenotypes  mIndivi.mainTrait[1]
          //  ....
        }
        simulator.plotPValues(new DoubleArrayList[]{pvalues});
        System.setOut(ps);
//        System.out.println("weights  " + type[k] + ":");
//        for (int i = 0; i < geneNum; i++) {
//          int[] counts = orgPowerCounts.get(i);
//          System.out.print(i);
//          for (int j = 0; j < counts.length; j++) {
//            System.out.print("\t" + counts[j]);
//          }
//          System.out.println();
//        }
//
//        for (int i = 0; i < geneNum; i++) {
//          int[] counts0 = orgPowerCounts.get(i);
//          int[] counts = nmlPowerCounts.get(i);
//          System.out.print(i);
//          for (int j = 0; j < counts.length; j++) {
//            if (j == 0) {
//              System.out.print("\t" + counts[j] + "(" + String.format("%.2f", ((double) counts[j]) / counts0[0]) + ")");
//            } else {
//              System.out.print("\t" + counts[j]);
//            }
//          }
//          System.out.println();
//        }
//
//        System.out.println("Average MSE:" + errors / simulationTime);
      }
      debugOut.close();

//      for(int i = 0; i < geneNum; i++){
//        int[] counts0 = orgPowerCounts.get(i);
//        int[] counts = fmCounts.get(i);
//        System.out.print(i);
//        for (int j = 0; j < counts.length; j++) {
//          if (j == 0) {
//            System.out.print("\t" + counts[j] + "(" + String.format("%.2f", ((double) counts[j]) / counts0[0]) + ")");
//          } else {
//            System.out.print("\t" + counts[j]);
//          }
//        }
//        System.out.println();
//      }

    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }

  public void plotPValues(DoubleArrayList[] pValueArray) throws Exception {
    PValuePainter pvPainter = new PValuePainter(450, 450);
    List<DoubleArrayList> pvalueList = new ArrayList<DoubleArrayList>();
    List<String> nameList = new ArrayList<String>();
    for (int i = 0; i < pValueArray.length; i++) {
      pvalueList.add(pValueArray[i]);
      nameList.add("GeneP" + i);
    }
    File plotFile2 = new File("ecsvar.qq.png");

    pvPainter.drawMultipleQQPlot(pvalueList, nameList, "Test", plotFile2.getCanonicalPath(), 1E-10);
    String info = "The QQ plot saved in " + plotFile2.getCanonicalPath();

  }

  public void loadGenotypes(Options option, List<List<AnnotSNP>> geneFullSnpList, List<List<AnnotSNP>> geneIndependentSnpList,
                            List<List<Individual>> geneIndList, List<DoubleMatrix2D> geneLDMatrix) throws Exception {
    String[] genes = option.geneList;
    String chroName = option.chromosome;
    String geneS;
    GenotypeQC qc = new GenotypeQC();
    int snpSize;
    String plinkFileName = option.plinFileName;
    PlinkDataset plinkDataset = new PlinkDataset(plinkFileName + ".fam", plinkFileName + ".bim", plinkFileName + ".bed");

    GeneInforProcessor gip = new GeneInforProcessor();
    Map<String, int[]> geneRegions = gip.readGenePos("F:\\java\\SimulationTestCode\\SeqGeneB36.txt", 100000, chroName);
    System.out.println("---------------Start to read genotypes from files----------------!");
    OpenIntIntHashMap indexGenotypePosMap = new OpenIntIntHashMap();
    //read genotypes from files
    for (int g = 0; g < genes.length; g++) {
      geneS = genes[g];
      int[] poss = geneRegions.get(geneS + "," + chroName);

      indexGenotypePosMap.clear();

      List<Individual> indListTmp = new ArrayList<Individual>();
      List<AnnotSNP> snpList = new ArrayList<AnnotSNP>();
      if (plinkDataset.readPlinkBinaryFormatMapByPysic(snpList, chroName, poss) <= 0) {
        continue;
        //throw new Exception("No SNP data!");
      }

      plinkDataset.readPlinkBinaryFormatPedigreeGenotype(indListTmp, snpList);
      qc.removeByMAF(snpList, indListTmp, 0.05);
      Collections.sort(snpList, new SNPPosiComparator());

      GenotypeBasedLDSparseMatrix gtyLDMatrix = new GenotypeBasedLDSparseMatrix(indListTmp, indexGenotypePosMap);
      if (gtyLDMatrix == null) {
        throw new Exception("LD data!");
      }
      IntArrayList indexPoses = new IntArrayList();
      snpSize = snpList.size();

      for (int j = 0; j < snpSize; j++) {
        indexGenotypePosMap.put(snpList.get(j).getPhysicalPosition(), snpList.get(j).order);
        indexPoses.add(snpList.get(j).getPhysicalPosition());
      }
      DoubleMatrix2D ldCorr = gtyLDMatrix.subDenseLDMatrix(indexPoses);

      int testSNPNum = ldCorr.columns();
      if (testSNPNum < 2) {
        throw new Exception("Only one SNP left!");
      }

      String info = qc.ldPruning(snpList, gtyLDMatrix, 0.975, true);
      System.out.println(info);

      if (snpList.size() < 3 || snpList.size() > 1500) {
        continue;
      }

      //update the SNP list after LD pruning
      snpSize = snpList.size();
      indexGenotypePosMap.clear();
      indexPoses.clear();
      gtyLDMatrix.releaseLDData();
      for (int j = 0; j < snpSize; j++) {
        indexGenotypePosMap.put(snpList.get(j).getPhysicalPosition(), snpList.get(j).order);
        indexPoses.add(snpList.get(j).getPhysicalPosition());
      }
      ldCorr = gtyLDMatrix.subDenseLDMatrixR2(indexPoses);
      System.out.println(ldCorr.toString());
      geneLDMatrix.add(ldCorr);

      qc.calculateMAF(snpList, indListTmp);

      geneFullSnpList.add(snpList);
      geneIndList.add(indListTmp);

      indexGenotypePosMap.clear();

      List<AnnotSNP> independentSnpList = new ArrayList<AnnotSNP>();
      independentSnpList.addAll(snpList);
      qc.ldPruning(independentSnpList, gtyLDMatrix, 0.01, true);
      geneIndependentSnpList.add(independentSnpList);

    }
  }

  public void simulateExpression(Options option, List<List<AnnotSNP>> geneFullSnpList, List<List<AnnotSNP>> geneIndependentSnpList,
                                 List<List<Individual>> geneIndList, List<DoubleMatrix2D> geneLDMatrix) throws Exception {
    double varPercentage = option.genetVarPercentage;
    double[][] determineCoeffs = option.determineCoeffs;
    int maxQTLNum = option.qtlNum;
    RandomData randGenerator = new RandomDataImpl(new MersenneTwister());

    List<AnnotSNP> qtlSnpList = new ArrayList<AnnotSNP>();

    int indivSize;
    double gtyScore = 0;

    DoubleArrayList trats = new DoubleArrayList();
    DoubleArrayList tratEnviroment = new DoubleArrayList();
    GenetAssociationAnalyzerP analyer = new GenetAssociationAnalyzerP();
    //   analyer.readCovariables("PCA_as_covariates_TW_SC.txt", indListSubSNP);

    int selectedNum = 0;
    int independentSNPNum = 0;
    int allSNPNum = 0;
    int selectedIndex;
    Set<Integer> avaialbeIndex = new HashSet<Integer>();
    int geneSize = geneFullSnpList.size();

    Population popu = new Population();
    indivSize = geneIndList.get(0).size();
    double[][] phenos = new double[indivSize][geneSize];
    IntArrayList indexOrders = new IntArrayList();
    IntArrayList QTLindexOrders = new IntArrayList();
    IntArrayList QTLindexOrders1 = new IntArrayList();
    RConnection rcon = new RConnection();
    rcon.eval("library(NNLM)");

    System.out.println("---------------Start to simulate gene expression by genotypes----------------!");
    //simulate gene expression for each gene
    for (int g = 0; g < geneSize; g++) {
      List<Individual> indListTmp = geneIndList.get(g);
      gtyScore = 0;
      popu.setAllIndiv(indListTmp);

      List<AnnotSNP> independentSnpList = geneIndependentSnpList.get(g);
      independentSNPNum = independentSnpList.size();
      while (true) {
        trats.clear();
        tratEnviroment.clear();
        qtlSnpList.clear();

        if (independentSNPNum < maxQTLNum) {
          qtlSnpList.addAll(independentSnpList);
        } else {
          selectedNum = 0;
          avaialbeIndex.clear();
          while (selectedNum < maxQTLNum) {
            selectedIndex = randGenerator.nextInt(0, independentSNPNum - 1);
            if (avaialbeIndex.contains(selectedIndex)) {
              continue;
            }
            avaialbeIndex.add(selectedIndex);
            AnnotSNP selVar = independentSnpList.get(selectedIndex);
            qtlSnpList.add(selVar);
            selectedNum++;
          }
        }
        Collections.sort(qtlSnpList, new SNPPosiComparator());

        QTLindexOrders.clear();
        int qtlNum = qtlSnpList.size();
        for (int i = 0; i < qtlNum; i++) {
          QTLindexOrders.add(qtlSnpList.get(i).order);
        }
        System.out.println();
        for (int i = 0; i < indivSize; i++) {
          Individual mIndivi = indListTmp.get(i);
          gtyScore = 0;

          StatusGtySet gtySet = mIndivi.markerGtySet;
          for (int snpPos = 0; snpPos < qtlNum; snpPos++) {
            int pos = qtlSnpList.get(snpPos).order;
            if (gtySet.existence.getQuick(pos)) {
              if (gtySet.paternalChrom.getQuick(pos)) {
                gtyScore += 1;
              }
              if (gtySet.maternalChrom.getQuick(pos)) {
                gtyScore += 1;
              }
            }
          }
          // mIndivi.setMainTrait(mIndivi.getMainTrait() + gtyScore);
          indListTmp.get(i).getMainTrait()[0] = gtyScore;
        }

        double acculated2PQ = 0;
        for (int i = 0; i < qtlNum; i++) {
          acculated2PQ += (2 * qtlSnpList.get(i).getAAlleleFreq() * (1 - qtlSnpList.get(i).getAAlleleFreq()));
        }

        trats.clear();
        tratEnviroment.clear();
        double vg = varPercentage;
        double alaph = Math.sqrt(vg / acculated2PQ);
        double ve = Math.sqrt(1 - vg);
        double ve1;

        for (int i = 0; i < indivSize; i++) {
          Individual mIndivi = indListTmp.get(i);
          ve1 = randGenerator.nextGaussian(0, ve);
          mIndivi.getMainTrait()[1] = mIndivi.getMainTrait()[0] * alaph + ve1;
          trats.add(mIndivi.getMainTrait()[1]);
          tratEnviroment.add(ve1);
          phenos[i][g] = mIndivi.getMainTrait()[1];
        }

        double meanT = Descriptive.mean(trats);
        double varT = Descriptive.sampleVariance(trats, meanT);

        if (Math.abs(varT - 1) > 0.01) {
          continue;
        }
        System.out.println("Overall mean and SD " + meanT + " " + varT);

        meanT = Descriptive.mean(tratEnviroment);
        varT = Descriptive.sampleVariance(tratEnviroment, meanT);
        System.out.println("Environmental  mean and SD " + meanT + " " + varT);

        indexOrders.clear();
        allSNPNum = geneFullSnpList.get(g).size();
        for (int i = 0; i < allSNPNum; i++) {
          indexOrders.add(geneFullSnpList.get(g).get(i).order);
        }
        double[] pValues = analyer.allelicAssociationTestQuantitativeTraitP(popu.getAllIndiv(), 1, indexOrders);
        //double[] chiseqr = analyer.allelicAssociationTestQuantitativeTraitP(popu.getAllIndiv(), 1);

        QTLindexOrders1.clear();
        for (int t = 0; t < indexOrders.size(); t++) {
          System.out.println(t + "\t" + geneFullSnpList.get(g).get(t).order + "\t" + geneFullSnpList.get(g).get(t).physicalPosition
                  + "\t" + geneFullSnpList.get(g).get(t).getaAlleleFreq() + "\t" + pValues[t] + "\t" + (QTLindexOrders.contains(geneFullSnpList.get(g).get(t).order) ? "QTL" : ""));
          // System.out.println(qtlSnpList.get(t).physicalPosition + "\t" + chiseqr[t]);
          QTLindexOrders1.add(t);
        }

        break;

        //  MultipleTestingMethod
      }
    }

    rcon.close();

    System.out.println("---------------Start to operate gene expression by genotypes----------------!");

    DoubleArrayList aList = new DoubleArrayList();
    DoubleArrayList bList = new DoubleArrayList();
    double mean1, sd1, mean2, sd2;
    double coeff;
    for (int g = 0; g < geneSize; g++) {
      aList.clear();
      for (int i = 0; i < indivSize; i++) {
        aList.add(phenos[i][g]);
      }
      mean1 = Descriptive.mean(aList);
      sd1 = Descriptive.sampleVariance(aList, mean1);

      for (int k = g + 1; k < geneSize; k++) {
        if (determineCoeffs[g][k] != 0) {
          bList.clear();
          for (int i = 0; i < indivSize; i++) {
            bList.add(phenos[i][k]);
          }
          mean2 = Descriptive.mean(bList);
          sd2 = Descriptive.sampleVariance(bList, mean2);
          coeff = sd2 * determineCoeffs[g][k] / sd1 / (1 - determineCoeffs[g][k]);
          coeff = Math.sqrt(coeff);
          for (int i = 0; i < indivSize; i++) {
            phenos[i][k] += (phenos[i][g] * coeff);
          }
        }
      }
    }

//update the phenotype scores at each individual
    for (int g = 0; g < geneSize; g++) {
      List<Individual> indListTmp = geneIndList.get(g);
      for (int i = 0; i < indivSize; i++) {
        Individual mIndivi = indListTmp.get(i);
        mIndivi.mainTrait[1] = phenos[i][g];

      }
    }

    calculateCorrelation(phenos);

  }

  public double expressionQTLFineMapping(Options option, List<List<AnnotSNP>> geneFullSnpList, List<List<AnnotSNP>> geneIndependentSnpList,
                                         List<List<Individual>> geneIndList, List<DoubleMatrix2D> geneLDMatrix, List<int[]> nnlPowerCounts, List<int[]> powerCounts, DoubleArrayList nhpvalues, List<int[]> fmCounts, BufferedWriter debugOut, double weight) throws Exception {
    double varPercentage = option.genetVarPercentage;
    double[][] determineCoeffs = option.determineCoeffs;
    int maxQTLNum = option.qtlNum;
    RandomData randGenerator = new RandomDataImpl(new MersenneTwister());

    List<AnnotSNP> qtlSnpList = new ArrayList<AnnotSNP>();

    double errors = 0.0;
    int indivSize;
    double gtyScore = 0;

    DoubleArrayList trats = new DoubleArrayList();
    DoubleArrayList tratEnviroment = new DoubleArrayList();
    GenetAssociationAnalyzerP analyer = new GenetAssociationAnalyzerP();
    //   analyer.readCovariables("PCA_as_covariates_TW_SC.txt", indListSubSNP);

    int selectedNum = 0;
    int independentSNPNum = 0;
    int allSNPNum = 0;
    int selectedIndex;
    Set<Integer> avaialbeIndex = new HashSet<Integer>();
    int geneSize = geneFullSnpList.size();

    Population popu = new Population();
    indivSize = geneIndList.get(0).size();
    double[][] phenos = new double[indivSize][geneSize];
    IntArrayList indexOrders = new IntArrayList();
    IntArrayList qtlPosOrders = new IntArrayList();
    IntArrayList qtlIndexOrders = new IntArrayList();
    IntArrayList testedIndexOrders = new IntArrayList();
    RConnection rcon = new RConnection();
    rcon.eval("library(NNLM)");
    //rcon.eval("library(limSolve)");
    //rcon.eval("library(bvls)");

    System.out.println("---------------Start to simulate gene expression by genotypes----------------!");
    //simulate gene expression for each gene
    for (int g = 0; g < geneSize; g++) {
      if (g != 1) {
        // continue;
      }
      List<Individual> indListTmp = geneIndList.get(g);
      gtyScore = 0;
      popu.setAllIndiv(indListTmp);

      List<AnnotSNP> independentSnpList = geneIndependentSnpList.get(g);
      independentSNPNum = independentSnpList.size();

      while (true) {
        trats.clear();
        tratEnviroment.clear();
        qtlSnpList.clear();

        if (independentSNPNum <= maxQTLNum) {
          qtlSnpList.addAll(independentSnpList);
        } else {
          selectedNum = 0;
          avaialbeIndex.clear();
          while (selectedNum < maxQTLNum) {
            selectedIndex = randGenerator.nextInt(0, independentSNPNum - 1);
            if (avaialbeIndex.contains(selectedIndex)) {
              continue;
            }
            avaialbeIndex.add(selectedIndex);
            AnnotSNP selVar = independentSnpList.get(selectedIndex);
            qtlSnpList.add(selVar);
            selectedNum++;
          }
        }
        Collections.sort(qtlSnpList, new SNPPosiComparator());

        qtlPosOrders.clear();
        int qtlNum = qtlSnpList.size();
        for (int i = 0; i < qtlNum; i++) {
          qtlPosOrders.add(qtlSnpList.get(i).order);
        }

        for (int i = 0; i < indivSize; i++) {
          Individual mIndivi = indListTmp.get(i);
          gtyScore = 0;

          StatusGtySet gtySet = mIndivi.markerGtySet;
          for (int snpPos = 0; snpPos < qtlNum; snpPos++) {
            int pos = qtlSnpList.get(snpPos).order;
            if (gtySet.existence.getQuick(pos)) {
              if (gtySet.paternalChrom.getQuick(pos)) {
                gtyScore += 1;
              }
              if (gtySet.maternalChrom.getQuick(pos)) {
                gtyScore += 1;
              }
            }
          }
          // mIndivi.setMainTrait(mIndivi.getMainTrait() + gtyScore);
          indListTmp.get(i).getMainTrait()[0] = gtyScore;
        }

        double acculated2PQ = 0;
        for (int i = 0; i < qtlNum; i++) {
          acculated2PQ += (2 * qtlSnpList.get(i).getAAlleleFreq() * (1 - qtlSnpList.get(i).getAAlleleFreq()));
        }

        trats.clear();
        tratEnviroment.clear();
        double vg = varPercentage;
        double alaph = Math.sqrt(vg / acculated2PQ);
        double ve = Math.sqrt(1 - vg);
        double ve1;

        for (int i = 0; i < indivSize; i++) {
          Individual mIndivi = indListTmp.get(i);
          ve1 = randGenerator.nextGaussian(0, ve);
          mIndivi.getMainTrait()[1] = mIndivi.getMainTrait()[0] * alaph + ve1;
          trats.add(mIndivi.getMainTrait()[1]);
          tratEnviroment.add(ve1);
          phenos[i][g] = mIndivi.getMainTrait()[1];
        }

        double meanT = Descriptive.mean(trats);
        double varT = Descriptive.sampleVariance(trats, meanT);

        if (Math.abs(varT - 1) > 0.01) {
          continue;
        }
        System.out.println("\nGene " + g);
        System.out.println("Overall mean and SD " + meanT + " " + varT);

        meanT = Descriptive.mean(tratEnviroment);
        varT = Descriptive.sampleVariance(tratEnviroment, meanT);
        System.out.println("Environmental  mean and SD " + meanT + " " + varT);

        indexOrders.clear();
        allSNPNum = geneFullSnpList.get(g).size();
        for (int i = 0; i < allSNPNum; i++) {
          indexOrders.add(geneFullSnpList.get(g).get(i).order);
        }
        ArrayList<double[]> lrPara = analyer.allelicAssociationTestQuantitativeTraitLR(popu.getAllIndiv(), 1, indexOrders);
        double[] pValues = lrPara.get(0);
        double[] beta = lrPara.get(1);
        double[] se = lrPara.get(2);

        //FINEMAP_FILE finemap = new FINEMAP_FILE();
        //double[] pValues = analyer.allelicAssociationTestQuantitativeTraitP(popu.getAllIndiv(), 1, indexOrders);
        //double[] zscores = analyer.allelicAssociationTestQuantitativeTraitZ(popu.getAllIndiv(), 1, indexOrders);
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < pValues.length; i++) {
          sb.append(pValues[i] + "\t");
        }
        sb.setCharAt(sb.length() - 1, '\n');
        //debugOut.write(sb.toString());
        //double[] chiseqr = analyer.allelicAssociationTestQuantitativeTraitP(popu.getAllIndiv(), 1);
        qtlIndexOrders.clear();
        testedIndexOrders.clear();
        int hitQTL1 = 0, hitQTL2 = 0, hitNonQTL = 0;
        double pt = 0.05 / allSNPNum;
        pt = 0.001;
        IntArrayList indexPoses = new IntArrayList();
        OpenIntIntHashMap indexGenotypePosMap = new OpenIntIntHashMap();
        //make the LD matrix is consistent to the positve z score
        for (int j = 0; j < pValues.length; j++) {
          int pos = geneFullSnpList.get(g).get(j).order;
          indexPoses.add(geneFullSnpList.get(g).get(j).getPhysicalPosition());
          indexGenotypePosMap.put(geneFullSnpList.get(g).get(j).getPhysicalPosition(), pos);

          /*
          if (pValues[j] < 0) {
            for (int i = 0; i < indivSize; i++) {
              Individual mIndivi = indListTmp.get(i);

              StatusGtySet gtySet = mIndivi.markerGtySet;
              if (gtySet.existence.getQuick(pos)) {
                gtySet.paternalChrom.putQuick(pos, !gtySet.paternalChrom.getQuick(pos));
                gtySet.maternalChrom.putQuick(pos, !gtySet.maternalChrom.getQuick(pos));
              }
            }
          }
           */
        }
        GenotypeBasedLDSparseMatrix gtyLDMatrix = new GenotypeBasedLDSparseMatrix(indListTmp, indexGenotypePosMap);
        if (gtyLDMatrix == null) {
          throw new Exception("LD data!");
        }
        DoubleMatrix2D ldCorr = gtyLDMatrix.subDenseLDMatrix(indexPoses);

        for (int t = 0; t < indexOrders.size(); t++) {
          //System.out.println(t + "\t" + geneFullSnpList.get(g).get(t).order + "\t" + geneFullSnpList.get(g).get(t).physicalPosition + "\t" + geneFullSnpList.get(g).get(t).getaAlleleFreq() + "\t" + chiseqr[t] + "\t" + (QTLindexOrders.contains(geneFullSnpList.get(g).get(t).order) ? "QTL" : ""));
          // System.out.println(qtlSnpList.get(t).physicalPosition + "\t" + chiseqr[t]);
          if (pValues[t] < pt) {
            testedIndexOrders.add(t);
          }

          if (qtlPosOrders.contains(geneFullSnpList.get(g).get(t).order)) {
            qtlIndexOrders.add(t);
            if (pValues[t] <= pt) {
              hitQTL1++;
            }
          }
        }
        if (testedIndexOrders.isEmpty()) {
          continue;
        }

        IntArrayList testedIndexOdersTmp = new IntArrayList();
        double[] pdcValues;
        int snpNum = pValues.length;

        // System.out.println(geneLDMatrix.get(g).toString());

        do {
          snpNum = testedIndexOrders.size();
          List<double[]> nnlsList = calculateEffectiveChiW(rcon, pValues, geneLDMatrix.get(g), testedIndexOrders, qtlIndexOrders, weight);
          pdcValues = nnlsList.get(0);
          errors = nnlsList.get(1)[0];
          // pdcValues = calculateEffectiveZ(rcon, zscores,pValues, ldCorr, testedIndexOrders, qtlIndexOrders, debugOut);
          //finemap.master_file(geneFullSnpList.get(g), testedIndexOrders, geneLDMatrix.get(g), beta, se);
          //finemap.fineMap();

          testedIndexOdersTmp.clear();
          for (int i = 0; i < snpNum; i++) {
            if (pdcValues[i] < pt) {
              testedIndexOdersTmp.add(testedIndexOrders.getQuick(i));
            }
          }
          if (testedIndexOdersTmp.isEmpty()) {
            break;
          }
          if (testedIndexOdersTmp.size() == testedIndexOrders.size()) {
            break;
          } else {
            testedIndexOrders.clear();
            testedIndexOrders.addAllOf(testedIndexOdersTmp);
          }

        } while (true);

        int snpSize = testedIndexOrders.size();
        int[] counts = new int[4];
        Arrays.fill(counts, 0);
        double r;
        /*
        for (int i = 0; i < snpSize; i++) {
          System.out.print(i);
          for (int t = 0; t <= i; t++) {
            System.out.print("\t-");
          }
          for (int j = i + 1; j < snpSize; j++) {
            r = geneLDMatrix.get(g).getQuick(testedIndexOrders.getQuick(i), testedIndexOrders.getQuick(j));
            System.out.print("\t" + r);
          }
          System.out.println();
        }
         */
        double chi;
        for (int i = 0; i < snpSize; i++) {
          chi = pValues[testedIndexOrders.getQuick(i)];
          // + "\t" + chi * chi
          System.out.println(testedIndexOrders.getQuick(i) + "\t" + pValues[testedIndexOrders.getQuick(i)] + "\t" + pdcValues[i] + "\t" + pdcValues[i + snpSize] + "\t" + pdcValues[i + snpSize + snpSize]
                  + "\t" + (qtlIndexOrders.contains(testedIndexOrders.getQuick(i)) ? "QTL" : "") + "\t" + (pdcValues[i] <= pt ? "Sig" : ""));
          if (pdcValues[i] <= pt) {
            if (qtlIndexOrders.contains(testedIndexOrders.getQuick(i))) {
              hitQTL2++;
            } else {
              hitNonQTL++;
            }
          }
          if (!qtlIndexOrders.contains(testedIndexOrders.getQuick(i)) && pdcValues[i + snpSize + snpSize] > 0) {
            nhpvalues.add(pdcValues[i]);
          }
        }
        if (hitQTL2 != 0) {
          if (hitQTL2 == qtlIndexOrders.size()) {
            counts[0]++;
          }
          if (hitQTL2 >= (qtlIndexOrders.size() - 1)) {
            counts[1]++;
          }
        }
        if (hitNonQTL != 0) {
          if (hitNonQTL > 0) {
            counts[2]++;
          }
          if (hitNonQTL > 1) {
            counts[3]++;
          }
        }

        int[] counts0 = nnlPowerCounts.get(g);
        for (int i = 0; i < counts.length; i++) {
          counts0[i] += counts[i];
        }

        //no fine mapping
        counts0 = powerCounts.get(g);
        if (hitQTL1 == qtlIndexOrders.size()) {
          counts0[0]++;
        }
        if (hitQTL1 >= (qtlIndexOrders.size() - 1)) {
          counts0[1]++;
        }


        //Finemap processing
//        int hitQTL3 = 0, hitNonQTL3 = 0;
//        int[] finemapIndex = finemap.snp_file();
//        for(int i = 0; i < finemapIndex.length; i++){
//          if(qtlIndexOrders.contains(finemapIndex[i] - 1)){
//            hitQTL3++;
//          }else {
//            hitNonQTL3++;
//          }
//        }
//        counts = fmCounts.get(g);
//        if (hitQTL3 != 0) {
//          if (hitQTL3 == qtlIndexOrders.size()) {
//            counts[0]++;
//          }
//          if (hitQTL3 >= (qtlIndexOrders.size() - 1)) {
//            counts[1]++;
//          }
//        }
//        if (hitNonQTL3 != 0) {
//          if (hitNonQTL3 > 0) {
//            counts[2]++;
//          }
//          if (hitNonQTL3 > 1) {
//            counts[3]++;
//          }
//        }

        break;
        //  MultipleTestingMethod
      }
    }

    rcon.close();
    return errors;
  }

  public void simpleExpressionQTLRelationCheck(Options option, List<List<AnnotSNP>> geneFullSnpList, List<List<AnnotSNP>> geneIndependentSnpList, List<List<Individual>> geneIndList,
                                               List<DoubleMatrix2D> geneLDMatrix, List<int[]> nnlPowerCounts, List<int[]> powerCounts, DoubleArrayList nhpvalues, BufferedWriter debugOut, double weight) throws Exception {
    double varPercentage = option.genetVarPercentage;
    double[][] determineCoeffs = option.determineCoeffs;
    int maxQTLNum = option.qtlNum;
    RandomData randGenerator = new RandomDataImpl(new MersenneTwister());

    List<AnnotSNP> qtlSnpList = new ArrayList<AnnotSNP>();

    int indivSize;
    double gtyScore = 0;

    DoubleArrayList trats = new DoubleArrayList();
    DoubleArrayList tratEnviroment = new DoubleArrayList();
    GenetAssociationAnalyzerP analyer = new GenetAssociationAnalyzerP();
    //   analyer.readCovariables("PCA_as_covariates_TW_SC.txt", indListSubSNP);

    int selectedNum = 0;
    int independentSNPNum = 0;
    int allSNPNum = 0;
    int selectedIndex;
    Set<Integer> avaialbeIndex = new HashSet<Integer>();
    int geneSize = geneFullSnpList.size();

    Population popu = new Population();
    indivSize = geneIndList.get(0).size();
    double[][] phenos = new double[indivSize][geneSize];
    IntArrayList indexOrders = new IntArrayList();
    IntArrayList qtlPosOrders = new IntArrayList();
    IntArrayList qtlIndexOrders = new IntArrayList();
    IntArrayList testedIndexOrders = new IntArrayList();
    RConnection rcon = new RConnection();
    rcon.eval("library(NNLM)");

    System.out.println("---------------Start to simulate gene expression by genotypes----------------!");
    int g = 0;
    List<Individual> indListTmp = geneIndList.get(g);
    gtyScore = 0;
    popu.setAllIndiv(indListTmp);

    List<AnnotSNP> independentSnpList = geneIndependentSnpList.get(g);
    independentSNPNum = independentSnpList.size();
    int sampleNum = 0;
    while (true) {
      trats.clear();
      tratEnviroment.clear();
      qtlSnpList.clear();

      if (independentSNPNum <= maxQTLNum) {
        qtlSnpList.addAll(independentSnpList);
      } else {
        selectedNum = 0;
        avaialbeIndex.clear();
        while (selectedNum < maxQTLNum) {
          selectedIndex = randGenerator.nextInt(0, independentSNPNum - 1);
          selectedIndex = 0;
          if (avaialbeIndex.contains(selectedIndex)) {
            continue;
          }
          avaialbeIndex.add(selectedIndex);
          AnnotSNP selVar = independentSnpList.get(selectedIndex);
          qtlSnpList.add(selVar);
          selectedNum++;
        }
      }
      Collections.sort(qtlSnpList, new SNPPosiComparator());

      qtlPosOrders.clear();
      int qtlNum = qtlSnpList.size();
      for (int i = 0; i < qtlNum; i++) {
        qtlPosOrders.add(qtlSnpList.get(i).order);
      }

      for (int i = 0; i < indivSize; i++) {
        Individual mIndivi = indListTmp.get(i);
        gtyScore = 0;

        StatusGtySet gtySet = mIndivi.markerGtySet;
        for (int snpPos = 0; snpPos < qtlNum; snpPos++) {
          int pos = qtlSnpList.get(snpPos).order;
          if (gtySet.existence.getQuick(pos)) {
            if (gtySet.paternalChrom.getQuick(pos)) {
              gtyScore += 1;
            }
            if (gtySet.maternalChrom.getQuick(pos)) {
              gtyScore += 1;
            }
          }
        }
        // mIndivi.setMainTrait(mIndivi.getMainTrait() + gtyScore);
        indListTmp.get(i).getMainTrait()[0] = gtyScore;
      }

      double acculated2PQ = 0;
      for (int i = 0; i < qtlNum; i++) {
        acculated2PQ += (2 * qtlSnpList.get(i).getAAlleleFreq() * (1 - qtlSnpList.get(i).getAAlleleFreq()));
      }

      trats.clear();
      tratEnviroment.clear();
      double vg = varPercentage;
      double alaph = Math.sqrt(vg / acculated2PQ);
      double ve = Math.sqrt(1 - vg);
      double ve1;

      for (int i = 0; i < indivSize; i++) {
        Individual mIndivi = indListTmp.get(i);
        ve1 = randGenerator.nextGaussian(0, ve);
        mIndivi.getMainTrait()[1] = mIndivi.getMainTrait()[0] * alaph + ve1;
        trats.add(mIndivi.getMainTrait()[1]);
        tratEnviroment.add(ve1);
        phenos[i][g] = mIndivi.getMainTrait()[1];
      }

      double meanT = Descriptive.mean(trats);
      double varT = Descriptive.sampleVariance(trats, meanT);

      if (Math.abs(varT - 1) > 0.01) {
        continue;
      }
      System.out.println("\nGene " + g);
      System.out.println("Overall mean and SD " + meanT + " " + varT);

      meanT = Descriptive.mean(tratEnviroment);
      varT = Descriptive.sampleVariance(tratEnviroment, meanT);
      System.out.println("Environmental  mean and SD " + meanT + " " + varT);

      indexOrders.clear();
      allSNPNum = geneFullSnpList.get(g).size();
      for (int i = 0; i < allSNPNum; i++) {
        indexOrders.add(geneFullSnpList.get(g).get(i).order);
      }
      // double[] chiseqr = analyer.allelicAssociationTestQuantitativeTraitChiSquare(popu.getAllIndiv(), 1, indexOrders);
      double[] chiseqr = analyer.allelicAssociationTestQuantitativeTraitChiSquare(popu.getAllIndiv(), 1, indexOrders);

      qtlIndexOrders.clear();
      testedIndexOrders.clear();
      int hitQTL1 = 0, hitQTL2 = 0, hitNonQTL = 0;
      double pt = 2;
      sampleNum++;
      debugOut.write(String.valueOf(chiseqr[0]));
      for (int t = 1; t < indexOrders.size(); t++) {
        //System.out.println(t + "\t" + geneFullSnpList.get(g).get(t).order + "\t" + geneFullSnpList.get(g).get(t).physicalPosition + "\t" + geneFullSnpList.get(g).get(t).getaAlleleFreq() + "\t" + chiseqr[t] + "\t" + (QTLindexOrders.contains(geneFullSnpList.get(g).get(t).order) ? "QTL" : ""));
        // System.out.println(qtlSnpList.get(t).physicalPosition + "\t" + chiseqr[t]);
        if (chiseqr[t] >= 0) {
          testedIndexOrders.add(t);
        }
        debugOut.write("\t" + chiseqr[t]);

        if (qtlPosOrders.contains(geneFullSnpList.get(g).get(t).order)) {
          qtlIndexOrders.add(t);
          if (chiseqr[t] <= pt) {
            hitQTL1++;
          }
        }
      }
      debugOut.write("\n");
      if (sampleNum > 10) {
        continue;
      }
      if (testedIndexOrders.isEmpty()) {
        continue;
      }

      IntArrayList testedIndexOdersTmp = new IntArrayList();
      double[] pdcValues;
      int snpNum = 0;
      // System.out.println(geneLDMatrix.get(g).toString());

      do {
        snpNum = testedIndexOrders.size();
        List<double[]> nnlsList = calculateEffectiveChiW(rcon, chiseqr, geneLDMatrix.get(g), testedIndexOrders, qtlIndexOrders, weight);
        pdcValues = nnlsList.get(0);

        testedIndexOdersTmp.clear();
        for (int i = 0; i < snpNum; i++) {
          if (pdcValues[i] < 2) {
            testedIndexOdersTmp.add(testedIndexOrders.getQuick(i));
          }
        }
        if (testedIndexOdersTmp.isEmpty()) {
          break;
        }
        if (testedIndexOdersTmp.size() == testedIndexOrders.size()) {
          break;
        } else {
          testedIndexOrders.clear();
          testedIndexOrders.addAllOf(testedIndexOdersTmp);
        }
      } while (true);

      int snpSize = testedIndexOrders.size();
      int[] counts = new int[4];
      Arrays.fill(counts, 0);
      double r;
      /*
        for (int i = 0; i < snpSize; i++) {
          System.out.print(i);
          for (int t = 0; t <= i; t++) {
            System.out.print("\t-");
          }
          for (int j = i + 1; j < snpSize; j++) {
            r = geneLDMatrix.get(g).getQuick(testedIndexOrders.getQuick(i), testedIndexOrders.getQuick(j));
            System.out.print("\t" + r);
          }
          System.out.println();
        }
       */
 /*
      double chi;
      for (int i = 0; i < snpSize; i++) {
        chi = MultipleTestingMethod.zScore(chiseqr[testedIndexOrders.getQuick(i)] / 2);
        System.out.println(i + "\t" + chiseqr[testedIndexOrders.getQuick(i)] + "\t" + chi * chi + "\t" + pdcValues[i] + "\t" + pdcValues[i + snpSize] + "\t" + pdcValues[i + snpSize + snpSize]
            + "\t" + (qtlIndexOrders.contains(testedIndexOrders.getQuick(i)) ? "QTL" : "") + "\t" + (pdcValues[i] <= pt ? "Sig" : ""));
        if (pdcValues[i] <= pt) {
          if (qtlIndexOrders.contains(testedIndexOrders.getQuick(i))) {
            hitQTL2++;
          } else {
            hitNonQTL++;
          }
        }
        if (!qtlIndexOrders.contains(testedIndexOrders.getQuick(i))) {
          nhpvalues.add(pdcValues[i]);
        }
      }
      if (hitQTL2 != 0) {
        if (hitQTL2 == qtlIndexOrders.size()) {
          counts[0]++;
        }
        if (hitQTL2 >= (qtlIndexOrders.size() - 1)) {
          counts[1]++;
        }
      }
      if (hitNonQTL != 0) {
        if (hitNonQTL > 0) {
          counts[2]++;
        }
        if (hitNonQTL > 1) {
          counts[3]++;
        }
      }

      int[] counts0 = nnlPowerCounts.get(g);
      for (int i = 0; i < counts.length; i++) {
        counts0[i] += counts[i];
      }
      counts0 = powerCounts.get(g);
      if (hitQTL1 == qtlIndexOrders.size()) {
        counts0[0]++;
      }
      if (hitQTL1 >= (qtlIndexOrders.size() - 1)) {
        counts0[1]++;
      }
       */
      break;
      //  MultipleTestingMethod
    }

    rcon.close();

  }

  public void nnlsSolver(RConnection rcon, double[] A, double[] b, DenseMatrix64F finalX1, DenseMatrix64F finalX2) throws Exception {
    int size = b.length / 2;
    rcon.assign("A", A);
    rcon.voidEval("A<-matrix(A, nrow=" + size + ", ncol=" + size + ", byrow = TRUE)");
    rcon.assign("b", b);
    rcon.voidEval("b<-matrix(b, nrow=" + size + ", ncol=" + 2 + ", byrow = TRUE)");
    //rcon.voidEval("beta <- nnlm(A, b, loss = 'mkl')");
    //rcon.voidEval("beta <- nnlm(A, b, loss = 'mse')");

    rcon.voidEval("beta <- nnlm(A, b)");
    double[][] coeffs = rcon.eval("beta$coefficients").asDoubleMatrix();
    for (int i = 0; i < size; i++) {
      finalX1.set(i, 0, coeffs[i][0]);
      finalX2.set(i, 0, coeffs[i][1]);
    }

  }

  public double[] nnlsSolver(RConnection rcon, double[] A, double[] b, double[] w) throws Exception {
    int size = b.length;
    rcon.assign("A", A);
    rcon.voidEval("A<-matrix(A, nrow=" + size + ", ncol=" + size + ", byrow = TRUE)");
    rcon.assign("b", b);
    rcon.voidEval("b<-matrix(b, nrow=" + size + ", ncol=" + 1 + ", byrow = TRUE)");
    rcon.assign("w", w);
    rcon.voidEval("w<-diag(w)");
    rcon.voidEval("w<-sqrt(w)");

    //rcon.voidEval("beta <- nnlm(A, b, loss = 'mkl')");
    rcon.voidEval("beta <- nnlm(A%*%w, w%*%b)");
    double[] coeffs = rcon.eval("beta$coefficients").asDoubles();

    return coeffs;
  }

  public double[] nnlsSolver0(RConnection rcon, double[] A, double[] b, double[] w, DenseMatrix64F finalX1, DenseMatrix64F finalX2) throws Exception {
    int size = b.length / 2;
    double[] A1 = new double[A.length];
    for (int i = 0; i < w.length; i++) {
      A1[i * w.length + i] = 1;
      for (int j = i + 1; j < w.length; j++) {
        double r = Math.random() * 10.0;
        A1[i * w.length + j] = r;
        A1[j * w.length + i] = r;
      }
    }
    rcon.assign("A", A);
    rcon.voidEval("A<-matrix(A, nrow=" + size + ", ncol=" + size + ", byrow = TRUE)");
    rcon.assign("A1", A1);
    rcon.voidEval("A1<-matrix(A1, nrow=" + size + ", ncol=" + size + ", byrow = TRUE)");
    rcon.assign("b", b);
    rcon.voidEval("b<-matrix(b, nrow=" + size + ", ncol=" + 2 + ", byrow = TRUE)");
    rcon.assign("w", w);
    rcon.voidEval("w<-diag(w)");
    // rcon.voidEval("w<-sqrt(w)");
    // rcon.voidEval("beta <- nnlm(A%*%w, w%*%b)");


    //rcon.voidEval("beta <- nnlm(A, b, loss = 'mkl')");
    //rcon.voidEval("beta <- nnlm(A, b, loss = 'mse')");
    rcon.voidEval("beta <- nnlm(A1%*%A, A1%*%b, max.iter = 1000000)");

    //rcon.voidEval("beta <- nnlm(A%*%w, w%*%b, loss = 'mse')");
    double[][] coeffs = rcon.eval("beta$coefficients").asDoubleMatrix();
    double[] errors = rcon.eval("beta$error").asDoubles();
    for (int i = 0; i < size; i++) {
      finalX1.set(i, 0, coeffs[i][0]);
      finalX2.set(i, 0, coeffs[i][1]);
    }
    /*
    double[][] out = rcon.eval("A%*%beta$coefficients").asDoubleMatrix();
    for (int i = 0; i < size; i++) {
      System.out.println((out[i][1] - b[i * 2 + 1]) + "\t" + (out[i][0] - b[i * 2]));

    }
    //method = \"lee\",
    rcon.voidEval("beta <- nnlm(A%*%w, w%*%b[,1],loss =\"mkl\",max.iter = 100000L)");
    //rcon.voidEval("beta <- nnlm(A, b, loss = 'mkl')");
    double[] coeffs1 = rcon.eval("beta$coefficients").asDoubles();

    double[] out1 = rcon.eval("A%*%beta$coefficients").asDoubles();
    for (int i = 0; i < size; i++) {
      System.out.println(out1[i] - b[i * 2]);

    }
     */
    int ssss = 0;
    return errors;
  }

  public void nnlsSolver(RConnection rcon, double[] A, double[] b, double[] w, DenseMatrix64F finalX1, DenseMatrix64F finalX2) throws Exception {
    int size = b.length / 2;
    rcon.assign("A", A);
    rcon.voidEval("A<-matrix(A, nrow=" + size + ", ncol=" + size + ", byrow = TRUE)");
    rcon.assign("b", b);
    rcon.voidEval("b<-matrix(b, nrow=" + size + ", ncol=" + 2 + ", byrow = TRUE)");
    rcon.voidEval("beta0 <- nnlm(A, b)");
    double[][] coeffs0 = rcon.eval("beta0$coefficients").asDoubleMatrix();
    rcon.eval("limit<-sum(beta0$coefficients[,1])");
    rcon.eval("F <- matrix(limit)");
    rcon.eval("n<-" + size);
    rcon.eval("E <- matrix(rep(1, n), ncol = n)");
    rcon.eval("F <- matrix(limit)");
    rcon.eval("G <- diag(rep(1, n), n)");
    rcon.eval("H <- matrix(rep(0, n), nrow = n)");
    rcon.assign("w", w);
    double[] bb = new double[size];
    for (int i = 0; i < size; i++) {
      bb[i] = b[i * 2];
    }
    rcon.assign("b", bb);
    rcon.eval("beta<-lsei(A, b, E, F, G, H,Wa = w)");

    double[] coeffs = rcon.eval("beta$X").asDoubles();
    for (int i = 0; i < size; i++) {
      finalX1.set(i, 0, coeffs[i]);
      finalX2.set(i, 0, coeffs0[i][1]);
    }

  }

  public double[] lmSolver(RConnection rcon, double[] A, double[] b, double[] w) throws Exception {
    int size = b.length;
    rcon.assign("A", A);
    rcon.voidEval("A<-matrix(A, nrow=" + size + ", ncol=" + size + ", byrow = TRUE)");
    rcon.assign("b", b);
    rcon.voidEval("b<-matrix(b, nrow=" + size + ", ncol=" + 1 + ", byrow = TRUE)");
    rcon.assign("w", w);
    rcon.voidEval("w<-diag(w)");
    rcon.voidEval("w<-sqrt(w)");

    //rcon.voidEval("beta <- nnlm(A, b, loss = 'mkl')");
    //rcon.voidEval("beta <- nnlm(A, b, loss = 'mse')");
    rcon.voidEval("beta <- Solve(A%*%w, w%*%b)");
    double[] coeffs = rcon.eval("beta").asDoubles();

    return coeffs;
  }

  public double[] bvlsSolver(RConnection rcon, double[] A, double[] b, double[] w) throws Exception {
    int size = b.length;
    rcon.assign("A", A);
    rcon.voidEval("A<-matrix(A, nrow=" + size + ", ncol=" + size + ", byrow = TRUE)");
    rcon.assign("b", b);
    rcon.voidEval("b<-matrix(b, nrow=" + size + ", ncol=" + 1 + ", byrow = TRUE)");
    rcon.assign("w", w);
    rcon.voidEval("w<-diag(w)");
    rcon.voidEval("w<-sqrt(w)");

    rcon.voidEval("bl <- b");
    rcon.voidEval("bl[bl>0] <- 0");
    rcon.voidEval("bu <- b");
    rcon.voidEval("bu[bu<0] <- 0");

    //rcon.voidEval("bl <- rep(min(b), " + size + ")");
    //rcon.voidEval("bl <- rep(0, " + size + ")");
    // rcon.voidEval("bu <- rep(max(b), " + size + ")");
    //rcon.voidEval("beta <- nnlm(A, b, loss = 'mkl')");
    //rcon.voidEval("beta <- nnlm(A, b, loss = 'mse')");
    //bl <- rep(min(A), 10)
//bu <- rep(max(A), 10)
//x<-bvls(A=A, b=b, bl=bl, bu=bu)
    rcon.voidEval(" tfit<-bvls(A=A%*%w, b=w%*%b, bl=bl, bu=bu)");

    double[] coeffs = rcon.eval("(tfit$x)").asDoubles();

    return coeffs;
  }

  public int[] calculateEffectiveChi(RConnection rcon, double[] pValueArray, DoubleMatrix2D ldCorr, IntArrayList QTLindexOrders, IntArrayList testedIndexOrders) throws Exception {
    int snpSize = testedIndexOrders.size();
    if (snpSize == 0) {
      return null;
    } else if (snpSize == 1) {
      return null;
    }

    double[] chisquares1 = new double[snpSize];
    for (int i = 0; i < snpSize; i++) {
      chisquares1[i] = pValueArray[testedIndexOrders.getQuick(i)] / 2;
    }
    double[] chisquares = MultipleTestingMethod.zScores(chisquares1);

    double[] arrayA = new double[snpSize * snpSize];
    double[] arrayB = new double[snpSize * 2];
    DenseMatrix64F xx1 = new DenseMatrix64F(snpSize, 1);
    DenseMatrix64F xx2 = new DenseMatrix64F(snpSize, 1);
    double r;
    int adjIndivSize = 2507;
    for (int i = 0; i < snpSize; i++) {
      arrayB[i * 2] = chisquares[i] * chisquares[i];
      arrayB[i * 2 + 1] = 1;
      arrayA[i * snpSize + i] = 1;
      for (int j = i + 1; j < snpSize; j++) {
        r = ldCorr.getQuick(testedIndexOrders.getQuick(i), testedIndexOrders.getQuick(j));
        r = r * r;
        //r = 1 - (adjIndivSize - 3.0) / (adjIndivSize - 2.0) * (1.0 - r) * (1 + 2.0 * (1 - r) / (adjIndivSize - 3.3));
        // r=Math.pow(r, 0.9);
        arrayA[i * snpSize + j] = r;
        arrayA[j * snpSize + i] = r;
      }
    }
    double df1 = 0;
    double Y1 = 0;
    double total = 0;

    int hitQTL = 0, hitNonQTL = 0;
    double p = 0, pt = 0.05 / pValueArray.length;
    int effectiveQTLNum = QTLindexOrders.size();
    Set<Integer> sigSetList1 = new HashSet<Integer>();
    Set<Integer> sigSetList2 = new HashSet<Integer>();

    nnlsSolver(rcon, arrayA, arrayB, xx1, xx2);
    int[] counts = new int[4];
    for (int i = 0; i < snpSize; i++) {
      if (QTLindexOrders.contains(i) && pValueArray[testedIndexOrders.getQuick(i)] >= pt) {
        //  effectiveQTLNum--;
      }
      Y1 += (xx1.get(i, 0));
      //  Y1 += chisquareArray[t][0];
      //   pValueArray[t].var *
      df1 += (xx2.get(i, 0));
      total += arrayB[i * 2];
      p = Gamma.incompleteGammaComplement(xx2.get(i, 0) / 2, xx1.get(i, 0) / 2);
      System.out.println(i + "\t" + pValueArray[testedIndexOrders.getQuick(i)] + "\t" + arrayB[i * 2] + "\t"
              + xx2.get(i, 0) + "\t" + xx1.get(i, 0) + "\t" + p + "\t" + (QTLindexOrders.contains(i) ? "QTL" : "") + "\t" + (p <= pt ? "Sig" : ""));
      if (p <= pt) {
        if (QTLindexOrders.contains(i)) {
          hitQTL++;
        } else {
          hitNonQTL++;
        }
      }
    }

    if (hitQTL == effectiveQTLNum) {
      counts[0]++;
    }
    if (hitQTL >= (effectiveQTLNum - 1)) {
      counts[1]++;
    }
    if (hitNonQTL > 0) {
      counts[2]++;
    }
    if (hitNonQTL > 1) {
      counts[3]++;
    }

    return counts;
  }

  public List<double[]> calculateEffectiveChiW(RConnection rcon, double[] pValueArray, DoubleMatrix2D ldCorr,
                                               IntArrayList testedIndexOrders, IntArrayList qtlIndex, double weight) throws Exception {
    int snpSize = testedIndexOrders.size();
    if (snpSize == 0) {
      return null;
    }

    List<double[]> nnlsList = new ArrayList<double[]>();
    double[] chisquares1 = new double[snpSize];
    for (int i = 0; i < snpSize; i++) {
      chisquares1[i] = pValueArray[testedIndexOrders.getQuick(i)] / 2;
    }
    double[] chisquares = MultipleTestingMethod.zScores(chisquares1);

    if (snpSize == 1) {
      nnlsList.add(new double[]{pValueArray[testedIndexOrders.getQuick(0)], 1, chisquares[0] * chisquares[0]});
      nnlsList.add(new double[]{0.0});
      return nnlsList;
    }

    double[] arrayA = new double[snpSize * snpSize];
    double[] arrayB = new double[snpSize * 2];
    double[] arrayW = new double[snpSize];

    DenseMatrix64F xx1 = new DenseMatrix64F(snpSize, 1);
    DenseMatrix64F xx2 = new DenseMatrix64F(snpSize, 1);
    double r;
    Arrays.fill(arrayW, 1);
    for (int i = 0; i < snpSize; i++) {
      if (qtlIndex.contains(testedIndexOrders.getQuick(i))) {
        arrayW[i] = 1;
      }
      arrayB[i * 2] = chisquares[i] * chisquares[i];
      arrayB[i * 2 + 1] = 1;
      arrayA[i * snpSize + i] = 1;

      for (int j = i + 1; j < snpSize; j++) {
        r = ldCorr.getQuick(testedIndexOrders.getQuick(i), testedIndexOrders.getQuick(j));
        // r = r * r;
        r = Math.abs(r);
        //r = 1 - (adjIndivSize - 3.0) / (adjIndivSize - 2.0) * (1.0 - r) * (1 + 2.0 * (1 - r) / (adjIndivSize - 3.3));
        // r=Math.pow(r, 0.9);
        arrayA[i * snpSize + j] = r;
        arrayA[j * snpSize + i] = r;

      }

    }

    BufferedWriter br = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("matrixA.txt")));
    for (int i = 0; i < arrayA.length; i++) {
      br.write(Double.toString(arrayA[i]));
      if (i % snpSize == snpSize - 1) {
        br.write("\n");
      } else {
        br.write("\t");
      }
    }
    br.flush();
    br.close();

    br = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("matrixB.txt")));
    for (int i = 0; i < arrayB.length; i++) {
      if (i % 2 == 0) {
        br.write(Double.toString(arrayB[i]));
      } else {
        br.write("\n");
      }
    }
    br.flush();
    br.close();

    //new FINEMAP_FILE().k_file(arrayW);

    double df1 = 0;
    double Y1 = 0;
    double total = 0;

    double p = 0;

    double[] values = new double[snpSize * 3];
    Arrays.fill(values, 1);
    //nnlsSolver(rcon, corrMatrix, arrayB, xx1, xx2);
    double[] errors = nnlsSolver0(rcon, arrayA, arrayB, arrayW, xx1, xx2);
    // nnlsSolver(rcon, arrayA, arrayB, arrayW, xx1, xx2);

    for (int i = 0; i < snpSize; i++) {
      Y1 += (xx1.get(i, 0));
      //  Y1 += chisquareArray[t][0];
      //   pValueArray[t].var *
      df1 += (xx2.get(i, 0));
      total += arrayB[i * 2];
      if (xx2.get(i, 0) < 0.001) {
        int ssss = 0;
      }
      p = Gamma.incompleteGammaComplement(xx2.get(i, 0) / 2, xx1.get(i, 0) / 2);
      if (xx2.get(i, 0) < 0.001 && p < 0.001) {
        int ssss = 0;
      }

      if (p < 0.001 && !qtlIndex.contains(testedIndexOrders.getQuick(i))) {
        int sss = 0;
      }
      values[i] = p;
      values[i + snpSize] = xx2.get(i, 0);
      values[i + snpSize + snpSize] = xx1.get(i, 0);
    }

    nnlsList.clear();
    nnlsList.add(values);
    nnlsList.add(errors);
    return nnlsList;
  }

  public double[] calculateEffectiveZ(RConnection rcon, double[] zscores, double[] pvalues, DoubleMatrix2D ldCorr,
                                      IntArrayList testedIndexOrders, IntArrayList qtlIndex, BufferedWriter debugOut) throws Exception {
    int snpSize = testedIndexOrders.size();
    if (snpSize == 0) {
      return null;
    }

    double[] zscores1 = new double[snpSize];
    for (int i = 0; i < snpSize; i++) {
      zscores1[i] = zscores[testedIndexOrders.getQuick(i)];
    }
    /*
    double[] halfp = new double[snpSize];
    for (int i = 0; i < snpSize; i++) {
      halfp[i] = zscores[testedIndexOrders.getQuick(i)] / 2;
    }
    double[] zscores = MultipleTestingMethod.zScores(halfp);
    double[] zscores1 = new double[snpSize];
    System.arraycopy(zscores, 0, zscores1, 0, snpSize);
     */
    if (snpSize == 1) {
      return new double[]{pvalues[testedIndexOrders.getQuick(0)], 1, zscores1[0]};
    }

    double[] corrMatrix = new double[snpSize * snpSize];

    double[] arrayW = new double[snpSize];

    double r;
    Arrays.fill(arrayW, 1);
    for (int i = 0; i < snpSize; i++) {
      if (qtlIndex.contains(testedIndexOrders.getQuick(i))) {
        // arrayW[i] = 0.2;
      }
      // zscores1[i] = Math.abs(zscores1[i]);

      corrMatrix[i * snpSize + i] = 1;

      for (int j = i + 1; j < snpSize; j++) {
        r = ldCorr.getQuick(testedIndexOrders.getQuick(i), testedIndexOrders.getQuick(j));
        //r = Math.abs(r);
        // r = r * r;
        //r = 1 - (adjIndivSize - 3.0) / (adjIndivSize - 2.0) * (1.0 - r) * (1 + 2.0 * (1 - r) / (adjIndivSize - 3.3));
        // r=Math.pow(r, 0.9);
        corrMatrix[i * snpSize + j] = r;
        corrMatrix[j * snpSize + i] = r;
      }

    }
    double df1 = 0;
    double Y1 = 0;
    double total = 0;

    double p = 0;

    double[] values = new double[snpSize * 3];
    Arrays.fill(values, 1);

    //nnlsSolver(rcon, corrMatrix, arrayB, xx1, xx2);
    // double[] effecZ = lmSolver(rcon, corrMatrix, zscores1, arrayW);
    // p=Gamma.incompleteGammaComplement(0.5, zscores1[testedIndexOrders.getQuick(0)]*zscores1[testedIndexOrders.getQuick(0)] / 2);
    //  double[] effecZ = nnlsSolver(rcon, corrMatrix, zscores1, arrayW);
    double[] effecZ = bvlsSolver(rcon, corrMatrix, zscores1, arrayW);

    double chi;
    StringBuilder sb = new StringBuilder();
    for (int i = 0; i < snpSize; i++) {
      sb.append((zscores1[i]) + "\t");

      chi = effecZ[i] * effecZ[i];
      Y1 += chi;
      //  Y1 += chisquareArray[t][0];
      //   zscores[t].var *
      df1 += 1;

      p = Gamma.incompleteGammaComplement(0.5, chi / 2);

      if (p < 1E-8) {
        int ssss = 0;
      }
      values[i] = p;
      values[i + snpSize] = 1;
      values[i + snpSize + snpSize] = chi;
    }
    sb.setCharAt(sb.length() - 1, '\n');
    debugOut.write(sb.toString());
    return values;
  }

  public void calculateCorrelation(double[][] expressoinValues) {
    DoubleArrayList aList = new DoubleArrayList();
    DoubleArrayList bList = new DoubleArrayList();
    int indivSize = expressoinValues.length, geneSize = expressoinValues[0].length;
    double mean1, sd1, mean2, sd2;
    double corr;
    for (int g = 0; g < geneSize; g++) {
      aList.clear();
      for (int i = 0; i < indivSize; i++) {
        aList.add(expressoinValues[i][g]);
      }
      mean1 = Descriptive.mean(aList);
      sd1 = Descriptive.sampleVariance(aList, mean1);
      for (int k = g + 1; k < geneSize; k++) {
        bList.clear();
        for (int i = 0; i < indivSize; i++) {
          bList.add(expressoinValues[i][k]);
        }
        mean2 = Descriptive.mean(bList);
        sd2 = Descriptive.sampleVariance(bList, mean2);
        corr = Descriptive.correlation(aList, Math.sqrt(sd1), bList, Math.sqrt(sd2));
        System.out.print("\t" + corr);
      }
      System.out.println();
    }
  }


  private int[] partitionPartMatrixLDBlocks(DoubleMatrix2D ldMatrix, int start0, int end0, int windowSize, double maxR2) throws Exception {
    //System.out.println("partition start:");
    if (start0 > ldMatrix.rows()) {
      System.out.println("Error: start index is over the matrix row");
      return null;
    }
    if (end0 > ldMatrix.rows()) {
      end0 = ldMatrix.rows();
    }
    int snpNum = end0;
    if (end0 - start0 <= 0) {
      System.out.println("Error: null matrix");
      return null;
    }

    if (end0 - start0 <= 1) {
      return new int[]{0, snpNum};
    }

    windowSize = Math.min(windowSize, snpNum - 1);
    double tolrateRatio = 0.03;
    int tolerateNum = (int) Math.rint(tolrateRatio * windowSize);
    IntArrayList boundaries = new IntArrayList();
    boundaries.add(start0);

    int overNum = 0, curRow = start0, curCol = start0 + 1;
    int start = curCol, end = curCol + windowSize;
    double maxR = Math.sqrt(maxR2);
    //double maxR = 0.3;
    for (int i = start; i < end; ++i) {
      overNum += ldMatrix.getQuick(curRow, i) > maxR ? 1 : 0;
    }
    while (true) {
//            System.out.println("row: " + curRow);
//            System.out.println("col: " + curCol);
//            System.out.println("overNum: " + overNum);
      if (overNum <= tolerateNum) {
        if (curCol == curRow) {
//                    overNum = overNum == 0 ? 1 : overNum;
//                    curRow += overNum == -1 ? 0 : overNum;
          //System.out.println("boundary: " + curRow);
          if (curRow >= snpNum) {
            boundaries.add(snpNum);
            break;
          }
          int bound = boundaries.getQuick(boundaries.size() - 1);
          overNum = 0;
          for (int i = bound; i < start; ++i) {
            overNum += ldMatrix.getQuick(i, curRow) > maxR ? 1 : 0;
          }
          if (overNum > tolrateRatio * (start - bound)) {
            curRow = start;
          } else {
            boundaries.add(curRow);
          }
          curCol = curRow + 1;
          start = curCol;
          if (curCol + windowSize > snpNum) {
            windowSize = snpNum - curCol;
            tolerateNum = (int) Math.rint(tolrateRatio * windowSize);
          }
          end = curCol + windowSize;
          overNum = 0;
          for (int i = start; i < end; ++i) {
            overNum += ldMatrix.getQuick(curRow, i) > maxR ? 1 : 0;
          }
        } else {
          while (overNum <= tolerateNum && curRow < curCol) {
            curRow++;
            if (curCol == curRow) {
              break;
            }
            overNum = 0;
            for (int i = start; i < end; ++i) {
              overNum += ldMatrix.getQuick(curRow, i) > maxR ? 1 : 0;
            }
          }
        }
      } else {
        if (curCol >= snpNum - 1) {
          boundaries.add(snpNum);
          break;
        }

        if (ldMatrix.getQuick(curRow, curCol) > maxR) {
          overNum--;
        }
        if (curCol + windowSize >= snpNum) {
          windowSize = snpNum - curCol - 1;
          tolerateNum = (int) Math.rint(tolrateRatio * windowSize);
        } else if (ldMatrix.getQuick(curRow, curCol + windowSize) > maxR) {
          overNum++;
        }
        curCol++;
        start = curCol;
        end = curCol + windowSize;
      }
    }

    int[] boundaryArray = new int[boundaries.size()];
    for (int i = 0; i < boundaryArray.length; i++) {
      boundaryArray[i] = boundaries.getQuick(i);
    }
    return boundaryArray;
  }

  private void recursivePartitionMatrixLDBlocks(DoubleMatrix2D ldMatrix, int start0, int end0, int windowSize, double maxR2, IntArrayList boundaries) throws Exception {
    int[] bds = partitionPartMatrixLDBlocks(ldMatrix, start0, end0, windowSize, maxR2);
    int checkingLen = 150;
    double checkingMaxR2 = 0.3;
    double checkingR2Inc = 0.025;
    for (int i = 1; i < bds.length; i++) {
      if (bds[i] - bds[i - 1] > checkingLen) {
        if (maxR2 <= checkingMaxR2) {
          recursivePartitionMatrixLDBlocks(ldMatrix, bds[i - 1], bds[i], windowSize, maxR2 + checkingR2Inc, boundaries);
        } else {
          boundaries.add(bds[i]);
        }
      } else {
        boundaries.add(bds[i]);
      }
    }
  }

  private void slidingWindowPartion(List<AnnotSNP> selectedVariantList, List<List<AnnotSNP>> blockedVariants, DoubleMatrix2D ldmatrix, int windowWidth, double maxR2) throws Exception {
    int totalSize = selectedVariantList.size();
    //System.out.println(totalSize);
    int startIndex = 0, startPos, endIndex = 0, endPos = 0;
    //List<Variant> independentList = new ArrayList<Variant>(selectedVariantList);
    List<AnnotSNP> selectedVariantListTmp = new ArrayList<>();
    int windowWithTmp = windowWidth / 10;
    startPos = selectedVariantList.get(0).physicalPosition;
    endIndex++;
    //List<Variant> independentListTmp = new ArrayList<Variant>();
    IntArrayList boundaries = new IntArrayList();

    while (endIndex < totalSize) {
      //endPos = independentList.get(endIndex).refStartPosition;
      endPos = selectedVariantList.get(endIndex).physicalPosition;
      if (endPos - startPos < windowWithTmp) {
        endIndex++;
        if (endIndex >= totalSize) {
          break;
        }
      } else {
//                selectedVariantListTmp.clear();
//                independentListTmp.clear();
//                boundaries.clear();
//                boundaries.add(0);
//                independentListTmp.addAll(independentList.subList(startIndex, endIndex));
//                ldmatrix = ldC.obtainLD(independentListTmp, totalPedSubjectNum, isPhased, threadNum, selectedVariantListTmp);
        boundaries.clear();
        boundaries.add(startIndex);
        recursivePartitionMatrixLDBlocks(ldmatrix, startIndex, endIndex, windowWidth, maxR2, boundaries);
        //int[] bounderies = partitionMatrixLDBlocks(ldmatrix, 50, maxR2);
        int length = boundaries.size();
        for (int i = 0; i < length - 2; ++i) {
          List<AnnotSNP> blockedVarList = new ArrayList<>();
          for (int j = boundaries.getQuick(i); j < boundaries.getQuick(i + 1); ++j) {
            //blockedVarList.add(independentListTmp.get(j));
            blockedVarList.add(selectedVariantList.get(j));
          }
          //System.out.println("block size:    " + blockedVarList.size());
          blockedVariants.add(blockedVarList);
        }
        //startIndex += boundaries.getQuick(length - 2);
        //startPos = independentList.get(startIndex).refStartPosition;
        startIndex = boundaries.getQuick(length - 2);
        startPos = selectedVariantList.get(startIndex).physicalPosition;
        endIndex++;
      }
    }

//        independentListTmp.clear();
//        selectedVariantListTmp.clear();
//        boundaries.clear();
//        boundaries.add(0);
//        independentListTmp.addAll(independentList.subList(startIndex, endIndex));
//        ldmatrix = ldC.obtainLD(independentListTmp, totalPedSubjectNum, isPhased, threadNum, selectedVariantListTmp);
    boundaries.clear();
    boundaries.add(startIndex);
    recursivePartitionMatrixLDBlocks(ldmatrix, startIndex, endIndex, windowWidth, maxR2, boundaries);
    //int[] bounderies = partitionMatrixLDBlocks(ldmatrix, 50, maxR2);
    int length = boundaries.size();
    for (int i = 0; i < length - 1; ++i) {
      List<AnnotSNP> blockedVarList = new ArrayList<>();
      for (int j = boundaries.getQuick(i); j < boundaries.getQuick(i + 1); ++j) {
        //blockedVarList.add(independentListTmp.get(j));
        blockedVarList.add(selectedVariantList.get(j));
      }
      //System.out.println("block size:    " + blockedVarList.size());
      blockedVariants.add(blockedVarList);
    }
  }

  public void nnlmSolver(RConnection rcon, double[] A, double[] b, DenseMatrix64F finalX1, DenseMatrix64F finalX2) throws Exception {
    int size = b.length / 2;
    double[] inis = new double[size];
    Arrays.fill(inis, 0);
    rcon.assign("inis", inis);
    rcon.voidEval("inis<-matrix(inis, nrow=" + size + ", ncol=" + 2 + ", byrow = TRUE)");

    rcon.assign("b", b);
    rcon.voidEval("b<-matrix(b, nrow=" + size + ", ncol=" + 2 + ", byrow = TRUE)");

    rcon.assign("A", A);
    rcon.voidEval("A<-matrix(A, nrow=" + size + ", ncol=" + size + ", byrow = TRUE)");

    double minSquareDiff = Double.MAX_VALUE;
    double squareDiff = 0;

    //n.threads
    //should set a high rel.tol or low max.iter=50,   loss = 'mkl' , n.threads=0 init=inits,
    rcon.voidEval("beta <- nnlm(A, b, max.iter=1000000, n.threads=0 )");
    double[][] coeffs = rcon.eval("beta$coefficients").asDoubleMatrix();
    for (int i = 0; i < size; i++) {
      finalX1.set(i, 0, coeffs[i][0]);
      finalX2.set(i, 0, coeffs[i][1]);
    }
    squareDiff = 0;
    double[][] out = rcon.eval("A%*%beta$coefficients").asDoubleMatrix();
    for (int i = 0; i < size; i++) {
      squareDiff += ((out[i][0] - b[i * 2]) * (out[i][0] - b[i * 2]));
      //System.out.println((out[t][1] - b[t * 2 + 1]) + "\t" + (out[t][0] - b[t * 2]));
    }
    if (squareDiff < minSquareDiff) {
      minSquareDiff = squareDiff;
    }

    // System.out.println("MinSquareDiff\t" + minSquareDiff);

        /*

    //method = \"lee\",
    rcon.voidEval("beta <- nnlm(A%*%w, w%*%b[,1],loss =\"mkl\",max.iter = 100000L)");
    //rcon.voidEval("beta <- nnlm(A, b, loss = 'mkl')");
    double[] coeffs1 = rcon.eval("beta$coefficients").asDoubles();

    double[] out1 = rcon.eval("A%*%beta$coefficients").asDoubles();
    for (int t = 0; t < size; t++) {
      System.out.println(out1[t] - b[t * 2]);

    }
         */
  }

  public double[] calculateEffectiveChiW(RConnection rcon, double[] arrayA, double[] arrayB, int snpSize) throws Exception {
    if (snpSize == 0) {
      return null;
    }

    if (snpSize == 1) {
      return new double[]{1, arrayB[0]};
    }

//        double[] arrayA = new double[snpSize * snpSize];
//        double[] arrayB = new double[snpSize * 2];
    DenseMatrix64F xx1 = new DenseMatrix64F(snpSize, 1);
    DenseMatrix64F xx2 = new DenseMatrix64F(snpSize, 1);
//        double r;

//        // build an exploration way of searching optimal weighting scale
//        //
//        for (int i = 0; i < snpSize; i++) {
//
//            arrayB[i * 2] = chisquares1.getQuick(i);
//            arrayB[i * 2 + 1] = 1;
//            arrayA[i * snpSize + i] = 1;
//
//            for (int j = i + 1; j < snpSize; j++) {
//                r = ldCorr.getQuick(i, j);
//                // r = r * r;
//                r = Math.abs(r);
//                //r = 1 - (adjIndivSize - 3.0) / (adjIndivSize - 2.0) * (1.0 - r) * (1 + 2.0 * (1 - r) / (adjIndivSize - 3.3));
//                // r=Math.pow(r, 0.9);
//                arrayA[i * snpSize + j] = r;
//                arrayA[j * snpSize + i] = r;
//
//            }
//
//        }
    double df1 = 0;
    double Y1 = 0;
    double total = 0;

    double p = 0;

    double[] values = new double[3];
    Arrays.fill(values, 1);
    //nnlsSolver(rcon, corrMatrix, arrayB, xx1, xx2);
    nnlmSolver(rcon, arrayA, arrayB, xx1, xx2);
    // nnlsSolver(rcon, arrayA, arrayB, xx1, xx2);

    for (int i = 0; i < snpSize; i++) {
      Y1 += (xx1.get(i, 0));
      //  Y1 += chisquareArray[t][0];
      //   pValueArrayECS[t].var *
      df1 += (xx2.get(i, 0));
      total += arrayB[i * 2];
      if (xx2.get(i, 0) < 0.001) {
        int ssss = 0;
      }

    }

    p = Gamma.incompleteGammaComplement(df1 / 2, Y1 / 2);

    values[0] = df1;
    values[1] = Y1;
    values[2] = p;

    return values;
  }

  //vary with type
  /*
  private void cellLineEnrichmentAnalysis(String[] labels, IntArrayList cellLines, DoubleArrayList idWeights, List<AnnotSNP> selectedVariantList, List<AnnotSNP> contrastVariantList, double pCut) {
    if (selectedVariantList == null || selectedVariantList.size() == 0) {
      System.out.println("Error: selectedVariantList is null!");
    }

    int cellLineNums = labels.length;
    if (cellLineNums == 0) {
      System.out.println("Error: No cellLines!");
    }

    String[] valueStr;
    DoubleArrayList values = new DoubleArrayList();
    int ii;
    for (int i = 0; i < cellLineNums; i++) {
      ii = 5 + i * 2;
      int selectedVarNums = selectedVariantList.size();
      for (int j = 0; j < selectedVarNums; ++j) {
        valueStr = selectedVariantList.get(j).getFeatureValues();
        if (valueStr[ii] == null) {
          continue;
        }
        values.add(Double.parseDouble(valueStr[ii]));
      }
      double[] score1Array = new double[values.size()];
      for (int j = 0; j < score1Array.length; j++) {
        score1Array[j] = values.getQuick(j);
      }
      values.clear();

      int contrastVarNums = contrastVariantList.size();
      for (int j = 0; j < contrastVarNums; ++j) {
        valueStr = contrastVariantList.get(j).getFeatureValues();
        if (valueStr[ii] == null) {
          continue;
        }
        values.add(Double.parseDouble(valueStr[ii]));
      }
      double[] score2Array = new double[values.size()];
      for (int j = 0; j < score2Array.length; j++) {
        score2Array[j] = values.getQuick(j);
      }
      values.clear();
      double pValue = new MannWhitneyTest(score1Array, score2Array, H1.GREATER_THAN).getSP();
      if (pValue < pCut) {
        cellLines.add(ii);
        idWeights.add(-Math.log10(pValue));

      }
      System.out.println(labels[(ii - 4) / 2] + "\t" + pValue);
    }

  }

  private List<double[]> sortOrgVariantList(List<AnnotSNP> orgVariantList, IntArrayList cellLines, DoubleArrayList idWeights) {
    List<double[]> order1 = new ArrayList<>();
    List<double[]> order2 = new ArrayList<>();
    List<double[]> order3 = new ArrayList<>();
    int orgVarNums = orgVariantList.size();
    int cellLineNums = cellLines.size();
    for (int i = 0; i < orgVarNums; ++i) {
      AnnotSNP var = orgVariantList.get(i);
      double score = 0.0;
      for (int j = 0; j < cellLineNums; ++j) {
        score += Double.parseDouble(var.getFeatureValues()[cellLines.getQuick(j)]) * idWeights.getQuick(j);
      }
      order1.add(new double[]{i, score});
      order2.add(new double[]{i, -Math.log10(var.getpValue())});
      order3.add(new double[]{i, -1});
    }

    Collections.sort(order1, new DoubleArrayListComparatorR(1));
    Collections.sort(order2, new DoubleArrayListComparatorR(1));
    double[] items;
    for (int i = 0; i < orgVarNums; ++i) {
      items = order3.get((int) (order1.get(i)[0]));
      items[1] = (i + 1);
    }
    for (int i = 0; i < orgVarNums; ++i) {
      items = order3.get((int) (order2.get(i)[0]));
      items[1] += (i + 1);
    }
    Collections.sort(order3, new DoubleArrayListComparator(1));
    return order3;
  }
  */

  private void conditionalPcut(RConnection rcon, double[] pValues, double[] chisquares, DoubleMatrix2D orgMatrix,
                               List<AnnotSNP> orgVariantList, List<AnnotSNP> selectedVariantList, List<AnnotSNP> contrastVariantList,
                               List<double[]> order, double pCut) throws Exception {
    int snpNum = order.size();
    double[] results0 = new double[2];
    double[] results1 = new double[2];
//        DoubleArrayList chiSquares = new DoubleArrayList();
    double df, chiq, p, r;
    AnnotSNP var;
    int index;
    double[] condiPvalues = new double[pValues.length];

//        DoubleMatrix2D matrixForCalc = null;
    double[] arrayA;
    double[] arrayB;
    List<Integer> indexForCalc = new ArrayList<>();



/*
        if (true) {
            //for testing
            MatrixBoxPainter test = new MatrixBoxPainter(800, 800);
            String outputPath = "./MatrixImageTest.png";
            test.setMaxValue(1);
            test.setMinValue(0);
            List<Variant> selectedVariantList11 = new ArrayList<Variant>();

            //System.out.println(orgMatrix.toString());
            String[] myIndex = new String[selectedVariantList11.size()];
            double[][] matrix = orgMatrix.toArray();
            for (int end = 0; end < myIndex.length; end++) {
                myIndex[end] = selectedVariantList11.get(end).label;
                for (int j = 0; j < myIndex.length; j++) {
                    matrix[end][j] = matrix[end][j] * matrix[end][j];

                }
            }
            test.drawMatrixBox(matrix, myIndex, null, outputPath);
        }
*/
//        int num1 = 0, num2 = 0;
    double ignoredLD = 0.98;
    for (int i = 0; i < snpNum; i++) {
      index = ((int) order.get(i)[0]);

      if (indexForCalc.size() == 0) {
        var = orgVariantList.get(index);
        //var.setCondiPValue(var.getpValue());
        condiPvalues[var.index] = pValues[var.index];
        results0[0] = 1;
        results0[1] = chisquares[var.index];
//        results0[1] = var.getChiSquare();
        indexForCalc.add(index);
//                chiSquares.clear();
//                chiSquares.add(results0[1]);
        results1[0] = results0[0];
        results1[1] = results0[1];
        continue;
      }

      indexForCalc.add(index);
      int snpNumTmp = indexForCalc.size();

      arrayA = new double[snpNumTmp * snpNumTmp];
      arrayB = new double[snpNumTmp * 2];
      boolean breakflag = false;
      for (int row = 0; row < snpNumTmp; ++row) {
        if (breakflag) {
          break;
        }

        //arrayB[row * 2] = orgVariantList.get(indexForCalc.get(row)).getChiSquare();
        arrayB[row * 2] = chisquares[orgVariantList.get(indexForCalc.get(row)).index];
        arrayB[row * 2 + 1] = 1;
        arrayA[row * snpNumTmp + row] = 1;

        for (int col = row + 1; col < snpNumTmp; ++col) {
          r = orgMatrix.getQuick(indexForCalc.get(row), indexForCalc.get(col));
          r = Math.abs(r);
          arrayA[row * snpNumTmp + col] = r;
          arrayA[col * snpNumTmp + row] = r;
          if (col == snpNumTmp - 1 && r > ignoredLD) {
            breakflag = true;
            break;
          }
        }
      }

      var = orgVariantList.get(index);

      if (breakflag) {
        //num1++;
        p = 1.0;
      } else {
        results1 = calculateEffectiveChiW(rcon, arrayA, arrayB, snpNumTmp);
        df = results1[0] - results0[0];
        chiq = results1[1] - results0[1];
        p = Gamma.incompleteGammaComplement(df / 2, chiq / 2);
      }

      if (p >=pValues[var.index]) {
        //var.setCondiPValue(p);
        condiPvalues[var.index] = pValues[var.index];
        results0[0] = results1[0];
        results0[1] = results1[1];
      } else {
        //num2++;
        //var.setCondiPValue(var.getpValue());
        condiPvalues[var.index] = pValues[var.index];
        //System.out.println(var.getCondiPValue());
        //varListForCalc.remove(snpNumTmp - 1);
        //chiSquares.remove(snpNumTmp - 1);
        indexForCalc.remove(snpNumTmp - 1);
      }
    }
    for (int i = 0; i < snpNum; ++i) {
      var = orgVariantList.get(i);
      if (condiPvalues[var.index] < pCut) {
        //selectedVarNums++;
        selectedVariantList.add(var);
      } else {
        contrastVariantList.add(var);
      }
    }
  }

  public List<AnnotSNP> iterativeConditionalECSVar(RConnection rcon, double[] pValues, double[] chisquares, DoubleMatrix2D orgMatrix,
                                                  List<List<AnnotSNP>> blockedVariantLists, List<AnnotSNP> qtlSnpList, double pCut, String type) throws Exception {

    //double pCut0 = 1.0E-5;
    //double cellLinepCut = 1.0E-3;
    // double pCut2 = 1.0E-5;

    List<AnnotSNP> selectedVariantList = new ArrayList<>();
    List<AnnotSNP> contrastVariantList = new ArrayList<>();


    int varListNums = blockedVariantLists.size();
    for (int i = 0; i < varListNums; ++i) {
      List<AnnotSNP> blockedVariantList = blockedVariantLists.get(i);
      int varListLength = blockedVariantList.size();
      List<double[]> order = new ArrayList<>();
      int[] views = new int[varListLength];
      //System.out.println("length:" + varListLength);
      for (int j = 0; j < varListLength; ++j) {
        AnnotSNP var = blockedVariantList.get(j);
        views[j] = var.index;
        //order.add(new double[]{j, -Math.log10(pValues[var.index])});
        if (qtlSnpList.contains(var)) {
          order.add(new double[]{j, 10});
        } else {
          order.add(new double[]{j, 1});
        }
      }
      Collections.sort(order, new DoubleArrayListComparatorR(1));

      DoubleMatrix2D blockedMatrix = orgMatrix.viewSelection(views, views);
      conditionalPcut(rcon, pValues, chisquares, blockedMatrix, blockedVariantList, selectedVariantList, contrastVariantList, order, pCut);
      //System.out.println(blockedVariantList.size() + "\t" + selectedVariantList.size());
      order.clear();
    }

    int size = selectedVariantList.size();
    List<AnnotSNP> results = new ArrayList<>();
    for (int i = 0; i < size; ++i) {
      results.add(selectedVariantList.get(i));
    }
//    double[] scores= new double[size];
//    Random rand = new Random();
//    for (int i = 0; i < size; i++) {
//      switch (type) {
//        case "positive":
//          scores[i] = -Math.log10(pValues[selectedVariantList.get(i).index]);
//          break;
//        case "negative":
//          scores[i] = Math.log10(pValues[selectedVariantList.get(i).index]);
//          break;
//        case "random":
//          scores[i] = rand.nextDouble();
//          break;
//      }
//    }
//
//    List<double[]> order1 = new ArrayList<>();
//    List<double[]> order2 = new ArrayList<>();
//    List<double[]> order3 = new ArrayList<>();
//    for (int i = 0; i < size; ++i) {
//      order1.add(new double[]{i, scores[i]});
//      order2.add(new double[]{i, -Math.log10(pValues[selectedVariantList.get(i).index])});
//      order3.add(new double[]{i, -1});
//    }
//
//    Collections.sort(order1, new DoubleArrayListComparatorR(1));
//    Collections.sort(order2, new DoubleArrayListComparatorR(1));
//    double[] items;
//    for (int i = 0; i < size; ++i) {
//      items = order3.get((int) (order1.get(i)[0]));
//      items[1] = (i + 1);
//    }
//    for (int i = 0; i < size; ++i) {
//      items = order3.get((int) (order2.get(i)[0]));
//      items[1] += (i + 1);
//    }
//    Collections.sort(order3, new DoubleArrayListComparator(1));
//
//    List<AnnotSNP> results = new ArrayList<>();
//    for (int i = 0; i < size; ++i) {
//      results.add(selectedVariantList.get((int) (order3.get(i)[0])));
//    }

    return results;
  }

  public double expressionQTLFineMapping(Options option, List<List<AnnotSNP>> geneFullSnpList, List<List<AnnotSNP>> geneIndependentSnpList,
                                      List<List<Individual>> geneIndList, String type) throws Exception {
    double varPercentage = option.genetVarPercentage;
    double[][] determineCoeffs = option.determineCoeffs;
    int maxQTLNum = option.qtlNum;
    RandomData randGenerator = new RandomDataImpl(new MersenneTwister());

    List<AnnotSNP> qtlSnpList = new ArrayList<AnnotSNP>();

    double errors = 0.0;
    int indivSize;
    double gtyScore = 0;

    RConnection rcon = new RConnection();
    rcon.eval("library(NNLM)");

    DoubleArrayList trats = new DoubleArrayList();
    DoubleArrayList tratEnviroment = new DoubleArrayList();
    GenetAssociationAnalyzerP analyer = new GenetAssociationAnalyzerP();

    int selectedNum = 0;
    int independentSNPNum = 0;
    int allSNPNum = 0;
    int selectedIndex;
    Set<Integer> avaialbeIndex = new HashSet<Integer>();
    int geneSize = geneFullSnpList.size();

    Population popu = new Population();
    indivSize = geneIndList.get(0).size();
    double[][] phenos = new double[indivSize][geneSize];
    IntArrayList indexOrders = new IntArrayList();
    IntArrayList qtlPosOrders = new IntArrayList();
    IntArrayList qtlIndexOrders = new IntArrayList();
    IntArrayList testedIndexOrders = new IntArrayList();


    System.out.println("---------------Start to simulate gene expression by genotypes----------------!");
    //simulate gene expression for each gene
    for (int g = 0; g < geneSize; g++) {
      if (g != 1) {
        // continue;
      }
      List<Individual> indListTmp = geneIndList.get(g);
      popu.setAllIndiv(indListTmp);

      List<AnnotSNP> independentSnpList = geneIndependentSnpList.get(g);
      independentSNPNum = independentSnpList.size();

      while (true) {
        trats.clear();
        tratEnviroment.clear();
        qtlSnpList.clear();

        if (independentSNPNum <= maxQTLNum) {
          qtlSnpList.addAll(independentSnpList);
        } else {
          selectedNum = 0;
          avaialbeIndex.clear();
          while (selectedNum < maxQTLNum) {
            selectedIndex = randGenerator.nextInt(0, independentSNPNum - 1);
            if (avaialbeIndex.contains(selectedIndex)) {
              continue;
            }
            avaialbeIndex.add(selectedIndex);
            AnnotSNP selVar = independentSnpList.get(selectedIndex);
            qtlSnpList.add(selVar);
            selectedNum++;
          }
        }
        Collections.sort(qtlSnpList, new SNPPosiComparator());

        qtlPosOrders.clear();
        int qtlNum = qtlSnpList.size();
        for (int i = 0; i < qtlNum; i++) {
          qtlPosOrders.add(qtlSnpList.get(i).order);
        }

        for (int i = 0; i < indivSize; i++) {
          Individual mIndivi = indListTmp.get(i);
          gtyScore = 0;

          StatusGtySet gtySet = mIndivi.markerGtySet;
          for (int snpPos = 0; snpPos < qtlNum; snpPos++) {
            int pos = qtlSnpList.get(snpPos).order;
            if (gtySet.existence.getQuick(pos)) {
              if (gtySet.paternalChrom.getQuick(pos)) {
                gtyScore += 1;
              }
              if (gtySet.maternalChrom.getQuick(pos)) {
                gtyScore += 1;
              }
            }
          }
          // mIndivi.setMainTrait(mIndivi.getMainTrait() + gtyScore);
          indListTmp.get(i).getMainTrait()[0] = gtyScore;
        }

        double acculated2PQ = 0;
        for (int i = 0; i < qtlNum; i++) {
          acculated2PQ += (2 * qtlSnpList.get(i).getAAlleleFreq() * (1 - qtlSnpList.get(i).getAAlleleFreq()));
        }

        trats.clear();
        tratEnviroment.clear();
        double vg = varPercentage;
        double alaph = Math.sqrt(vg / acculated2PQ);
        double ve = Math.sqrt(1 - vg);
        double ve1;

        for (int i = 0; i < indivSize; i++) {
          Individual mIndivi = indListTmp.get(i);
          ve1 = randGenerator.nextGaussian(0, ve);
          mIndivi.getMainTrait()[1] = mIndivi.getMainTrait()[0] * alaph + ve1;
          trats.add(mIndivi.getMainTrait()[1]);
          tratEnviroment.add(ve1);
          phenos[i][g] = mIndivi.getMainTrait()[1];
        }

        double meanT = Descriptive.mean(trats);
        double varT = Descriptive.sampleVariance(trats, meanT);

        if (Math.abs(varT - 1) > 0.01) {
          continue;
        }
        System.out.println("\nGene " + g);
        System.out.println("Overall mean and SD " + meanT + " " + varT);

        meanT = Descriptive.mean(tratEnviroment);
        varT = Descriptive.sampleVariance(tratEnviroment, meanT);
        System.out.println("Environmental  mean and SD " + meanT + " " + varT);

        indexOrders.clear();
        allSNPNum = geneFullSnpList.get(g).size();
        for (int i = 0; i < allSNPNum; i++) {
          indexOrders.add(geneFullSnpList.get(g).get(i).order);
        }
        ArrayList<double[]> lrPara = analyer.allelicAssociationTestQuantitativeTraitLR(popu.getAllIndiv(), 1, indexOrders);
        double[] pValues = lrPara.get(0);

//        StringBuilder sb = new StringBuilder();
//        for (int i = 0; i < pValues.length; i++) {
//          sb.append(pValues[i] + "\t");
//        }
//        sb.setCharAt(sb.length() - 1, '\n');
        //debugOut.write(sb.toString());
        //double[] chiseqr = analyer.allelicAssociationTestQuantitativeTraitP(popu.getAllIndiv(), 1);
        qtlIndexOrders.clear();
        testedIndexOrders.clear();
        int hitQTL1 = 0, hitQTL2 = 0, hitNonQTL = 0;
        double pt = 0.05 / allSNPNum;
        pt = 0.001;
        IntArrayList indexPoses = new IntArrayList();
        OpenIntIntHashMap indexGenotypePosMap = new OpenIntIntHashMap();
        //make the LD matrix is consistent to the positve z score
        for (int j = 0; j < pValues.length; j++) {
          int pos = geneFullSnpList.get(g).get(j).order;
          indexPoses.add(geneFullSnpList.get(g).get(j).getPhysicalPosition());
          indexGenotypePosMap.put(geneFullSnpList.get(g).get(j).getPhysicalPosition(), pos);

          /*
          if (pValues[j] < 0) {
            for (int i = 0; i < indivSize; i++) {
              Individual mIndivi = indListTmp.get(i);

              StatusGtySet gtySet = mIndivi.markerGtySet;
              if (gtySet.existence.getQuick(pos)) {
                gtySet.paternalChrom.putQuick(pos, !gtySet.paternalChrom.getQuick(pos));
                gtySet.maternalChrom.putQuick(pos, !gtySet.maternalChrom.getQuick(pos));
              }
            }
          }
           */
        }
        GenotypeBasedLDSparseMatrix gtyLDMatrix = new GenotypeBasedLDSparseMatrix(indListTmp, indexGenotypePosMap);
        if (gtyLDMatrix == null) {
          throw new Exception("LD data!");
        }
        DoubleMatrix2D ldCorr = gtyLDMatrix.subDenseLDMatrix(indexPoses);

        for (int t = 0; t < indexOrders.size(); t++) {
          //System.out.println(t + "\t" + geneFullSnpList.get(g).get(t).order + "\t" + geneFullSnpList.get(g).get(t).physicalPosition + "\t" + geneFullSnpList.get(g).get(t).getaAlleleFreq() + "\t" + (qtlIndexOrders.contains(geneFullSnpList.get(g).get(t).order) ? "QTL" : ""));
          // System.out.println(qtlSnpList.get(t).physicalPosition + "\t" + chiseqr[t]);
          //System.out.println(pValues[t]);
          if (pValues[t] < pt) {
            testedIndexOrders.add(t);
          }

          if (qtlPosOrders.contains(geneFullSnpList.get(g).get(t).order)) {
            qtlIndexOrders.add(t);
            if (pValues[t] <= pt) {
              hitQTL1++;
            }
          }
        }
        if (testedIndexOrders.isEmpty()) {
          continue;
        }

        IntArrayList testedIndexOdersTmp = new IntArrayList();
        double[] pdcValues;
        int snpNum = pValues.length;

        // System.out.println(geneLDMatrix.get(g).toString());

        System.out.println("qtl rsid:");
        for (int i = 0; i < qtlSnpList.size(); ++i) {
          System.out.println(qtlSnpList.get(i).getRSID());
        }

        int snpSize = pValues.length;
        //int[] counts = new int[4];
        //Arrays.fill(counts, 0);
        //double r;
//        System.out.println();
//        System.out.println("size:" + geneFullSnpList.get(g).size());
//        for (int i = 0; i < geneFullSnpList.get(g).size(); ++i) {
//          System.out.println(geneFullSnpList.get(g).get(i).probeSetID);
//        }
//        System.out.println();
        double[] chisquares1 = new double[snpSize];
        for (int i = 0; i < snpSize; i++) {
          chisquares1[i] = pValues[i] / 2;
        }
        double[] chisquares = MultipleTestingMethod.zScores(chisquares1);

        for (int i = 0; i < snpSize; ++i) {
          geneFullSnpList.get(g).get(i).index = i;
        }

        List<List<AnnotSNP>> blockedVariants = new ArrayList<>();
        // r or r-square?
        slidingWindowPartion(geneFullSnpList.get(g), blockedVariants, ldCorr, 50, 0.15);
        //System.out.println(blockedVariants.size());
        List<AnnotSNP> results = iterativeConditionalECSVar(rcon, pValues, chisquares, ldCorr, blockedVariants, qtlSnpList, 5 * Math.pow(10, -3), type);

        System.out.println("selected rsid(" + type + "): ");
        for (int i = 0; i < results.size(); ++i) {
          System.out.println("order " + i + ":" + results.get(i).getRSID());
        }
        break;
        /*
        for (int i = 0; i < snpSize; i++) {
          System.out.print(i);
          for (int t = 0; t <= i; t++) {
            System.out.print("\t-");
          }
          for (int j = i + 1; j < snpSize; j++) {
            r = geneLDMatrix.get(g).getQuick(testedIndexOrders.getQuick(i), testedIndexOrders.getQuick(j));
            System.out.print("\t" + r);
          }
          System.out.println();
        }
         */
//        double chi;
//        for (int i = 0; i < snpSize; i++) {
//          chi = pValues[testedIndexOrders.getQuick(i)];
//          // + "\t" + chi * chi
//          System.out.println(testedIndexOrders.getQuick(i) + "\t" + pValues[testedIndexOrders.getQuick(i)] + "\t" + pdcValues[i] + "\t" + pdcValues[i + snpSize] + "\t" + pdcValues[i + snpSize + snpSize]
//                  + "\t" + (qtlIndexOrders.contains(testedIndexOrders.getQuick(i)) ? "QTL" : "") + "\t" + (pdcValues[i] <= pt ? "Sig" : ""));
//          if (pdcValues[i] <= pt) {
//            if (qtlIndexOrders.contains(testedIndexOrders.getQuick(i))) {
//              hitQTL2++;
//            } else {
//              hitNonQTL++;
//            }
//          }
//          if (!qtlIndexOrders.contains(testedIndexOrders.getQuick(i)) && pdcValues[i + snpSize + snpSize] > 0) {
//            nhpvalues.add(pdcValues[i]);
//          }
//        }
//        if (hitQTL2 != 0) {
//          if (hitQTL2 == qtlIndexOrders.size()) {
//            counts[0]++;
//          }
//          if (hitQTL2 >= (qtlIndexOrders.size() - 1)) {
//            counts[1]++;
//          }
//        }
//        if (hitNonQTL != 0) {
//          if (hitNonQTL > 0) {
//            counts[2]++;
//          }
//          if (hitNonQTL > 1) {
//            counts[3]++;
//          }
//        }
//
//        int[] counts0 = nnlPowerCounts.get(g);
//        for (int i = 0; i < counts.length; i++) {
//          counts0[i] += counts[i];
//        }
//
//        //no fine mapping
//        counts0 = powerCounts.get(g);
//        if (hitQTL1 == qtlIndexOrders.size()) {
//          counts0[0]++;
//        }
//        if (hitQTL1 >= (qtlIndexOrders.size() - 1)) {
//          counts0[1]++;
//        }
//
//        break;
//        //  MultipleTestingMethod
      }
    }

    rcon.close();
    return errors;
  }



}
