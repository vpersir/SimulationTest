/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.genetsimulator.controller;

import cern.colt.list.DoubleArrayList;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.CholeskyDecomposition;
import cern.jet.stat.Descriptive;
import cern.jet.stat.Gamma;
import cern.jet.stat.Probability;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.StringTokenizer;
import org.apache.commons.math.distribution.ChiSquaredDistributionImpl;
import org.apache.commons.math.distribution.GammaDistribution;
import org.apache.commons.math.distribution.GammaDistributionImpl;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.random.CorrelatedRandomVectorGenerator;
import org.apache.commons.math.random.GaussianRandomGenerator;
import org.apache.commons.math.random.MersenneTwister;
 
 
import org.cobi.eqtlsimulator.Constants;
import static org.cobi.eqtlsimulator.Constants.FILE_READER_BUFFER_SIZE;
import org.cobi.genetsimulator.entity.Individual;
import org.cobi.genetsimulator.entity.PValueSet;
import org.cobi.genetsimulator.entity.PValueWeight;
import org.cobi.genetsimulator.entity.PValueWeightPComparator;
 
import org.cobi.util.stat.MultipleTestingMethod;
import umontreal.ssj.probdistmulti.BiNormalDist;
import umontreal.ssj.randvar.NormalGen;
import umontreal.ssj.randvarmulti.MultinormalCholeskyGen;
import umontreal.ssj.rng.MT19937;
import umontreal.ssj.rng.WELL607;
 

/**
 *
 * @author mxli
 */
public class PValuesAnalyzer implements Constants {

    private static List<Individual> indList = Collections.synchronizedList(new ArrayList<Individual>());
    private static List<PValueSet[]> permuHistogramSets = Collections.synchronizedList(new ArrayList<PValueSet[]>());

    public PValueSet[] histogram(double[] pValues, double[][] range) {
        int rowNum = range.length;
        int pNum = pValues.length;
        double effectiveNum = 0;
        int count = 0;
        for (int j = 0; j < pNum; j++) {
            if (pValues[j] >= 0 && pValues[j] <= 1.0) {
                effectiveNum += 1.0;
            } else {
                System.err.println(pValues[j]);
            }

        }
        PValueSet[] pvalueCount = new PValueSet[rowNum];

        for (int i = 0; i < rowNum; i++) {
            count = 0;
            for (int j = 0; j < pNum; j++) {
                if (pValues[j] >= range[i][0] && pValues[j] < range[i][1]) {
                    count++;
                }
            }

            pvalueCount[i] = new PValueSet(count, count / effectiveNum);
        }
        return pvalueCount;
    }

    
    public double[] readPValues(String pValueFileName) throws Exception {
        File pValueFile = new File(pValueFileName);
        BufferedReader br = new BufferedReader(new FileReader(pValueFile), FILE_READER_BUFFER_SIZE);
        String line = null;
        String delmilit = ", \t";
        int pValueIndex = 8;

        DoubleArrayList pValueList = new DoubleArrayList();

        //skip the head line
        br.readLine();
        String tmpStr;
        int index;
        int lineCounter = -1;
        StringBuffer tmpBuffer = new StringBuffer();
        try {
            while ((line = br.readLine()) != null) {
                if (line.trim().length() == 0) {
                    continue;
                }
                lineCounter++;
                StringTokenizer tokenizer = new StringTokenizer(line, delmilit);
                index = 0;

                while (tokenizer.hasMoreTokens()) {
                    tmpBuffer.delete(0, tmpBuffer.length());
                    tmpBuffer.append(tokenizer.nextToken().trim());
                    tmpStr = tmpBuffer.toString();
                    if (index == pValueIndex) {
                        double p = Double.parseDouble(tmpStr);
                        pValueList.add(p);
                    }
                    if (index == pValueIndex) {
                        break;
                    }
                    index++;
                }
            }

            StringBuffer runningInfo = new StringBuffer();
            runningInfo.append("The number of SNPs  in map file ");
            runningInfo.append(pValueFile.getName());
            runningInfo.append(" is ");
            runningInfo.append(pValueList.size());
            runningInfo.append(".");

          //  GlobalVariables.addInforLog(runningInfo.toString());
            int size = pValueList.size();
            double[] pValues = new double[size];
            for (int i = 0; i < size; i++) {
                pValues[i] = pValueList.get(i);
            }
            return (pValues);
        } finally {
            br.close();
        }

    }

    public void pValueDistribution(double[] pValues) throws Exception {
        double[][] thresholds = {{0, 0.01}, {0, 0.05}, {0, 0.1}, {0, 0.2}, {0, 0.3}, {0, 0.4}, {0, 0.5}, {0.01, 0.2}, {0.05, 0.2}, {0.05, 0.5}, {0.2, 0.5}};
        //double[][] thresholds = {{0, 0.01}, {0.01, 0.02}, {0.02, 0.03}, {0.03, 0.04}, {0.04, 0.05}, {0.05, 0.06}, {0.06, 0.07}, {0.07, 0.08}, {0.08, 0.09}, {0.09, 0.1}, {0.1, 0.11}, {0.11, 0.12}, {0.12, 0.13}, {0.13, 0.14}, {0.14, 0.15}, {0.15, 0.16}, {0.16, 0.17}, {0.17, 0.18}, {0.18, 0.19}, {0.19, 0.2}, {0.2, 0.21}, {0.21, 0.22}, {0.22, 0.23}, {0.23, 0.24}, {0.24, 0.25}, {0.25, 0.26}, {0.26, 0.27}, {0.27, 0.28}, {0.28, 0.29}, {0.29, 0.3}, {0.3, 0.31}, {0.31, 0.32}, {0.32, 0.33}, {0.33, 0.34}, {0.34, 0.35}, {0.35, 0.36}, {0.36, 0.37}, {0.37, 0.38}, {0.38, 0.39}, {0.39, 0.4}, {0.4, 0.41}, {0.41, 0.42}, {0.42, 0.43}, {0.43, 0.44}, {0.44, 0.45}, {0.45, 0.46}, {0.46, 0.47}, {0.47, 0.48}, {0.48, 0.49}, {0.49, 0.5}, {0.5, 0.51}, {0.51, 0.52}, {0.52, 0.53}, {0.53, 0.54}, {0.54, 0.55}, {0.55, 0.56}, {0.56, 0.57}, {0.57, 0.58}, {0.58, 0.59}, {0.59, 0.6}, {0.6, 0.61}, {0.61, 0.62}, {0.62, 0.63}, {0.63, 0.64}, {0.64, 0.65}, {0.65, 0.66}, {0.66, 0.67}, {0.67, 0.68}, {0.68, 0.69}, {0.69, 0.7}, {0.7, 0.71}, {0.71, 0.72}, {0.72, 0.73}, {0.73, 0.74}, {0.74, 0.75}, {0.75, 0.76}, {0.76, 0.77}, {0.77, 0.78}, {0.78, 0.79}, {0.79, 0.8}, {0.8, 0.81}, {0.81, 0.82}, {0.82, 0.83}, {0.83, 0.84}, {0.84, 0.85}, {0.85, 0.86}, {0.86, 0.87}, {0.87, 0.88}, {0.88, 0.89}, {0.89, 0.9}, {0.9, 0.91}, {0.91, 0.92}, {0.92, 0.93}, {0.93, 0.94}, {0.94, 0.95}, {0.95, 0.96}, {0.96, 0.97}, {0.97, 0.98}, {0.98, 0.99}, {0.99, 1}};

        PValueSet[] observedSet = histogram(pValues, thresholds);
        for (int j = 0; j < thresholds.length; j++) {
            System.out.println(observedSet[j].getCount() + " " + observedSet[j].getProportion());
        }
    }

    // a method proposed by Zaykin et al (2007) Combining minP-values in large-scale genomics experiments.
    public double combinePValuebyGammaMethod(double[] pValueArray) throws Exception {
        //a =0.0137 gives ��soft truncation threshold�� (STT) =0.05; a =0.0383 gives STT=0.1; and a =1 gives STT =1/e
        double a = 0.0137;
        double b = 1;
        int snpSize = pValueArray.length;
        if (snpSize == 0) {
            return -1;
        } else if (snpSize == 1) {
            return pValueArray[0];
        }
        GammaDistribution gameDis = new GammaDistributionImpl(a, b);

        int testNum = 0;
        double Y = 0;
        for (int i = 0; i < snpSize; i++) {
            double p = pValueArray[i];
            if (p < 1e-16) {
                p = 1e-16;
                Y += gameDis.inverseCumulativeProbability(1 - p);
                testNum++;
            } else if (p < 1.0) {
                Y += gameDis.inverseCumulativeProbability(1 - p);
                testNum++;
            }
        }
        if (testNum > 0) {
            double p = 1 - gameDis.cumulativeProbability(Y);
            if (p < 0.0) {
                p = 0.0;
            }
            return p;
        } else {
            return 1.0;
        }
    }

    public double combinePValuebyFisherCombinationTest(double[] pValueArray) throws Exception {
        int snpSize = pValueArray.length;
        if (snpSize == 0) {
            return -1;
        } else if (snpSize == 1) {
            return pValueArray[0];
        }

        double Y = 0;
        for (int i = 0; i < snpSize; i++) {
            double p = pValueArray[i];
            Y += (-2 * Math.log(p));
        }
        double p = Probability.chiSquareComplemented(snpSize * 2, Y);
        if (p < 0.0) {
            p = 0.0;
        }
        return p;

    }

    public double combinePValuebyFisherCombinationTest(double[] pValueArray, int snpSize) throws Exception {
        if (snpSize == 0) {
            return -1;
        } else if (snpSize == 1) {
            return pValueArray[0];
        }

        double Y = 0;
        for (int i = 0; i < snpSize; i++) {
            double p = pValueArray[i];
            Y += (-2 * Math.log(p));
        }
        double p = Probability.chiSquareComplemented(snpSize * 2, Y);
        if (p < 0.0) {
            p = 0.0;
        }
        return p;
    }

    public double combinePValuebyScaleedFisherCombinationTest(double[] pValueArray, DoubleMatrix2D ldCorr) throws Exception {
        int snpSize = pValueArray.length;
        if (snpSize == 0) {
            return -1;
        } else if (snpSize == 1) {
            return pValueArray[0];
        }

        double Y = 0;
        for (int i = 0; i < snpSize; i++) {
            double p = pValueArray[i];
            Y += (-2 * Math.log(p));
        }
        //calcualte the scalled chi
        double varTD = 4 * snpSize;
        double a1 = 3.263119, a2 = 0.709866, a3 = 0.026589, a4 = 0.709866;
        //double a1 = 3.263, a2 = 0.710, a3 = 0.027, a4 = 0.709866;
        for (int i = 0; i < snpSize; i++) {
            for (int j = i + 1; j < snpSize; j++) {
                varTD += 2 * (a1 * ldCorr.get(i, j) + a2 * ldCorr.get(i, j) * ldCorr.get(i, j) + a3 * ldCorr.get(i, j) * ldCorr.get(i, j) * ldCorr.get(i, j));
            }
        }
        // varTD = 8 + 2 * 0.98;
        double c = varTD / (4 * snpSize);
        double f = 8 * snpSize * snpSize / varTD;
        double p = Probability.chiSquareComplemented(f, Y / c);
        if (p < 0.0) {
            p = 0.0;
        }
        return p;

    }

    public double combinePValuebyStouffer(PValueWeight[] pValueArray, DoubleMatrix2D ldCorr) throws Exception {
        int snpSize = pValueArray.length;
        if (snpSize == 0) {
            return -1;
        } else if (snpSize == 1) {
            return pValueArray[0].pValue;
        }

        Set<Integer> selectedIndex = new HashSet<Integer>();

        DoubleMatrix2D covP = new DenseDoubleMatrix2D(snpSize, snpSize);
        for (int i = 0; i < snpSize; i++) {
            covP.setQuick(i, i, 1);
            selectedIndex.add(pValueArray[i].index);
            for (int j = i + 1; j < snpSize; j++) {
                if (i == j) {
                    continue;
                }
                //poweredcorrMat.setQuick(i, j, Math.pow(corrMat.getQuick(i, j), power));
                double x = ldCorr.getQuick(pValueArray[i].index, pValueArray[j].index);
                x = x * x;
                //when r2                 
                //I do not know why it seems if I use x1=x1*x1  it woks better in terms of type 1 error
                x = (((((0.7723 * x - 1.5659) * x + 1.201) * x - 0.2355) * x + 0.2184) * x + 0.6086) * x;
                //x = x1 * x1;
                covP.setQuick(pValueArray[i].index, pValueArray[j].index, x);
                covP.setQuick(pValueArray[j].index, pValueArray[i].index, x);
            }
        }

        DoubleMatrix1D pvs = new DenseDoubleMatrix1D(snpSize);
        CholeskyDecomposition cd = new CholeskyDecomposition(covP);

        Algebra algb = new Algebra();

        DoubleMatrix2D lowInversTriangle = algb.inverse(cd.getL());

        // System.out.println(cd.getL().zMult(lowInversTriangle,null).toString());
        double q = 0;
        for (int i = 0; i < snpSize; i++) {
            if (pValueArray[i].pValue > 0.5) {
                if (pValueArray[i].pValue >= 1) {
                    //the Probability.normalInverse could handle 1E-323 but cannot handle 1-(1E-323)
                    q = Probability.normalInverse(1E-323);
                } else {
                    q = Probability.normalInverse(1 - pValueArray[i].pValue);
                }
            } else {
                if (pValueArray[i].pValue < 1E-320) {
                    pValueArray[i].pValue = 1E-320;
                }
                q = -Probability.normalInverse(pValueArray[i].pValue);
            }
            pvs.setQuick(i, q);
        }

        DoubleMatrix1D transvertP = lowInversTriangle.zMult(pvs, null);
        double Y = 0;
        for (int i = 0; i < snpSize; i++) {
            Y += transvertP.getQuick(i);
        }

        Y /= Math.sqrt(snpSize);
        double p = 1 - Probability.normal(Y);
        if (p < 0.0) {
            p = 0.0;
        }
        return p;
    }

    public double combinePValuebyStoufferTwo(PValueWeight[] pValueArray, DoubleMatrix2D ldCorr) throws Exception {
        int snpSize = pValueArray.length;
        if (snpSize == 0) {
            return -1;
        } else if (snpSize == 1) {
            return pValueArray[0].pValue;
        }

        double x = ldCorr.getQuick(pValueArray[0].index, pValueArray[1].index);
        x = x * x;
        //when r2                 
        //I do not know why it seems if I use x1=x1*x1  it woks better in terms of type 1 error
        x = (((((0.7723 * x - 1.5659) * x + 1.201) * x - 0.2355) * x + 0.2184) * x + 0.6086) * x;

        double p1 = pValueArray[0].pValue;
        double p2 = pValueArray[1].pValue;

        double q1 = 0, q2 = 0;

        if (p1 > 0.5) {
            if (p1 >= 1) {
                //the Probability.normalInverse could handle 1E-323 but cannot handle 1-(1E-323)
                q1 = Probability.normalInverse(1E-323);
            } else {
                q1 = Probability.normalInverse(1 - p1);
            }
        } else {
            if (p1 < 1E-320) {
                q1 = 1E-320;
            }
            q1 = -Probability.normalInverse(p1);
        }
        if (p2 > 0.5) {
            if (p2 >= 1) {
                //the Probability.normalInverse could handle 1E-323 but cannot handle 1-(1E-323)
                q2 = Probability.normalInverse(1E-323);
            } else {
                q2 = Probability.normalInverse(1 - p2);
            }
        } else {
            if (p2 < 1E-320) {
                q2 = 1E-320;
            }
            q2 = -Probability.normalInverse(p2);
        }

        p2 = q2 - x * q1;
        p2 = p2 / Math.sqrt(1 - x * x);
        //p2 = 1 - Probability.normal(p2);
        double Y = 0;
        Y = pValueArray[0].var * q1 + pValueArray[1].var * p2;
        Y /= Math.sqrt(pValueArray[0].var * pValueArray[0].var + pValueArray[1].var * pValueArray[1].var);
        double p = 1 - Probability.normal(Y);

        if (p < 0.0) {
            p = 0.0;
        }
        return p;
    }

    public double combinePValuebyScaleedFisherCombinationTestCovLogP(double[] pValueArray, DoubleMatrix2D ldCorr) throws Exception {
        int snpSize = pValueArray.length;
        if (snpSize == 0) {
            return -1;
        } else if (snpSize == 1) {
            return pValueArray[0];
        }

        DoubleMatrix2D covLogP = new DenseDoubleMatrix2D(snpSize, snpSize);
        for (int i = 0; i < snpSize; i++) {
            for (int j = i + 1; j < snpSize; j++) {
                if (i == j) {
                    continue;
                }
                double x1 = ldCorr.getQuick(i, j);

                //for r
                //y = 0.0079x3 + 3.9459x2 - 0.0024x ; y = 0.0331x3 + 3.9551x2 - 0.0156x
                //x = 0.0331 * (Math.pow(x1, 3)) + 3.9551 * (Math.pow(x1, 2)) - 0.0156 * x1;
                if (x1 < 0) {
                    x1 = -x1;
                }
                x1 = (0.75 * x1 + 3.25) * x1;

                /*
                 //the formular I used in KGG
                 if (x > 0.5) {
                 x = 0.75 * x + 3.25 * Math.sqrt(x);
                 } else {
                 x = 8.6 * x;
                 //  x1 = 0.75 * x1 + 3.25 * Math.sqrt(x1);
                 }
                 */
                //for r2
                // x1 = x1 * x1;
                // x1 = 3.9471 * x1;
                // x1 = 4 * x1;
                covLogP.setQuick(i, j, x1);
            }
        }

        double Y = 0;
        for (int i = 0; i < snpSize; i++) {
            double p = pValueArray[i];
            Y += (-2 * Math.log(p));
        }
        //calcualte the scalled chi
        double varTD = 4 * snpSize;
        for (int i = 0; i < snpSize; i++) {
            for (int j = i + 1; j < snpSize; j++) {
                varTD += 2 * covLogP.get(i, j);
            }
        }
        double c = varTD / (4 * snpSize);
        double f = 8 * snpSize * snpSize / varTD;
        double p = Probability.chiSquareComplemented(f, Y / c);
        if (p < 0.0) {
            p = 0.0;
        }
        return p;
    }

    public double combinePValuebyCorrectedChiFisherCombinationTestMXLi(double[] pValueArray, DoubleMatrix2D ldCorr) throws Exception {
        int snpSize = pValueArray.length;
        if (snpSize == 0) {
            return -1;
        } else if (snpSize == 1) {
            return pValueArray[0];
        }

        double[] chisquares1 = new double[snpSize];
        for (int i = 0; i < snpSize; i++) {
            chisquares1[i] = pValueArray[i] / 2;
        }
        double[] chisquares = MultipleTestingMethod.zScores(chisquares1);
        DoubleMatrix2D A = new DenseDoubleMatrix2D(snpSize, snpSize);
        for (int i = 0; i < snpSize; i++) {
            A.set(i, i, 1);
            for (int j = i + 1; j < snpSize; j++) {
                double r = ldCorr.getQuick(i, j);
                r = r * r;
                r = Math.pow(r, 0.75);
                A.set(i, j, r);
                A.set(j, i, r);
            }
        }

        double threshold = 0.05;
        int columnNum = 997;
        double scale = 10;
        double pairThreshold = 0.3;
        int maxBlock = 5;

        //  System.out.println(A.toString());
        /// ArrayList<Integer[]> blockCluster = new MarixDensity().getCluster(A, threshold, maxBlock, scale, pairThreshold);
        return mySudoSVDSolver(A, chisquares);

    }

    
   
    
     
    public double mySudoSVDSolver(DoubleMatrix2D A, double[] chisquares) throws Exception {
        cern.colt.matrix.linalg.SingularValueDecomposition svd = new cern.colt.matrix.linalg.SingularValueDecomposition(A);
        int size = A.columns();
        DoubleMatrix2D x1 = null;
        DoubleMatrix2D x2 = null;
        DoubleMatrix2D b1 = new DenseDoubleMatrix2D(size, 1);
        DoubleMatrix2D b2 = new DenseDoubleMatrix2D(size, 1);

        for (int i = 0; i < size; i++) {
            b1.set(i, 0, chisquares[i] * chisquares[i]);
            b2.set(i, 0, 1);
        }

        Algebra alg = new Algebra();
        DoubleMatrix2D U = alg.transpose(svd.getU());
        DoubleMatrix2D W = svd.getS();
        // System.out.println(W.toString());
        DoubleMatrix2D V = svd.getV();

        double threshlod = 1e-6;
        double maxX1 = 0, maxB1 = 0;
        double maxX2 = 0, maxB2 = 0;
        int trancatedID = size;

        DoubleMatrix2D W1 = null;
        DoubleMatrix2D tb1 = null;
        DoubleMatrix2D tb2 = null;
        for (int t = 0; t < size; t++) {
            if (maxB1 < Math.abs(b1.get(t, 0))) {
                maxB1 = Math.abs(b1.get(t, 0));
            }
            if (maxB2 < Math.abs(b2.get(t, 0))) {
                maxB2 = Math.abs(b2.get(t, 0));
            }
            if (W.get(t, t) < threshlod) {
                W.set(t, t, 0);
            } else {
                W.set(t, t, 1 / W.get(t, t));
            }
        }
        double minDiff = 0.001 * size;
        maxB1 = maxB1 * 5;
        maxB2 = maxB2 * 5;
        trancatedID--;
        boolean needContinue = true;
        double diff = 0;
        do {
            W1 = alg.mult(alg.mult(V, W), U);
            x1 = alg.mult(W1, b1);
            maxX1 = 0;
            maxX2 = 0;
            for (int t = 0; t < size; t++) {
                if (maxX1 < Math.abs(x1.get(t, 0))) {
                    maxX1 = Math.abs(x1.get(t, 0));
                }
            }
            if (maxX1 <= maxB1) {
                x2 = alg.mult(W1, b2);
                for (int t = 0; t < size; t++) {
                    if (maxX2 < Math.abs(x2.get(t, 0))) {
                        maxX2 = Math.abs(x2.get(t, 0));
                    }
                }
                if (maxX2 <= maxB2) {
                    tb1 = alg.mult(A, x1);
                    diff = 0;
                    for (int t = 0; t < size; t++) {
                        diff += Math.abs(tb1.getQuick(t, 0) - b1.getQuick(t, 0));
                    }
                    if (diff < minDiff) {
                        tb2 = alg.mult(A, x2);
                        diff = 0;
                        for (int t = 0; t < size; t++) {
                            diff += Math.abs(tb2.getQuick(t, 0) - b2.getQuick(t, 0));
                        }
                        if (diff < minDiff) {
                            needContinue = false;
                            break;
                        }
                    }
                }
            }
            W.set(trancatedID, trancatedID, W.get(trancatedID, trancatedID) / 2);
            if (W.get(trancatedID, trancatedID) <= threshlod) {
                W.set(trancatedID, trancatedID, 0);
                trancatedID--;
                if (trancatedID < 0) {
                    break;
                }
            }
        } while (needContinue);

        double df = 0;
        double Y = 0;
        for (int i = 0; i < size; i++) {
            Y += (x1.get(i, 0));
            df += (x2.get(i, 0));
        }
        if (Y <= 0 || df <= 0) {
            return 1;
        }
        if (df < 1) {
            df = 1;
        }
        double p1 = Gamma.incompleteGammaComplement(df / 2, Y / 2);
        if (p1 < 1E-8) {
            int sss = 0;
        }
        return p1;
    }

  
  
    public static int[] partitionEvenBlock(int startIndex, int endIndex, int intervalLen) {
        int totalSnpSize = endIndex - startIndex;
        int blockNum = totalSnpSize / intervalLen;
        if (blockNum <= 1) {
            int[] bigBlockIndexes = new int[2];
            bigBlockIndexes[0] = startIndex;
            bigBlockIndexes[1] = endIndex;
            return bigBlockIndexes;
        }

        int[] bigBlockIndexes = new int[blockNum + blockNum + 1];
        bigBlockIndexes[0] = startIndex;
        blockNum++;
        for (int i = 1; i < blockNum; i++) {
            bigBlockIndexes[i * 2] = startIndex + i * intervalLen;
            bigBlockIndexes[i * 2 - 1] = (bigBlockIndexes[i * 2] + bigBlockIndexes[i * 2 - 2]) / 2;
        }
        blockNum--;
        if (bigBlockIndexes[blockNum * 2] < endIndex) {
            bigBlockIndexes[blockNum * 2] = endIndex;
            bigBlockIndexes[blockNum * 2 - 1] = (bigBlockIndexes[blockNum * 2] + bigBlockIndexes[blockNum * 2 - 2]) / 2;
        }
        return bigBlockIndexes;
    }
 
    public double combinePValuebyCorrectedChiFisherCombinationTestMXLiS(double[] pValueArray, DoubleMatrix2D ldCorr) throws Exception {
        int snpSize = pValueArray.length;
        if (snpSize == 0) {
            return -1;
        } else if (snpSize == 1) {
            return pValueArray[0];
        }

        double[] chisquares1 = new double[snpSize];
        for (int i = 0; i < snpSize; i++) {
            chisquares1[i] = pValueArray[i] / 2;
        }
        double[] chisquares = MultipleTestingMethod.zScores(chisquares1);
        double[] dfs = new double[snpSize];
        for (int i = 0; i < snpSize; i++) {
            chisquares[i] = chisquares[i] * chisquares[i];
            chisquares1[i] = chisquares[i];
            dfs[i] = 1;
        }
        double r;

        double df1 = 0;
        double Y1 = 0;
        for (int t = 0; t < snpSize; t++) {
            for (int k = t + 1; k < snpSize; k++) {
                r = ldCorr.getQuick(t, k);
                r = r * r;
                r = Math.pow(r, 0.75);
                chisquares[k] -= (chisquares[t] * r);
                dfs[k] -= (dfs[t] * r);
            }
        }

        for (int i = 0; i < snpSize; i++) {
            Y1 += chisquares[i];
            df1 += dfs[i];
        }

        if (Y1 <= 0 || df1 <= 0) {
            return 1;
        }

        //calcualte the scalled chi 
        // p1 = Probability.chiSquareComplemented((df1 + df2) / 2, (Y1 + Y2) / 2);
        double p1 = Probability.chiSquareComplemented(df1, Y1);
        if (p1 < 0.0) {
            p1 = 0.0;
        }

        return p1;
    }

    public double combinePValuebyCorrectedChiFisherCombinationTestJohnnyDf(double[] pValueArray, DoubleMatrix2D ldCorr) throws Exception {
        int snpSize = pValueArray.length;
        if (snpSize == 0) {
            return -1;
        } else if (snpSize == 1) {
            return pValueArray[0];
        }

        double[] chisquares1 = new double[snpSize];
        for (int i = 0; i < snpSize; i++) {
            chisquares1[i] = pValueArray[i] / 2;
        }
        double[] chisquares = MultipleTestingMethod.zScores(chisquares1);
        double[][] chisquareArray = new double[snpSize][1];
        for (int i = 0; i < snpSize; i++) {
            chisquareArray[i][0] = chisquares[i] * chisquares[i];
        }

        DoubleMatrix2D chisequreMat = new DenseDoubleMatrix2D(chisquareArray);
        Algebra alg = new Algebra();
        DoubleMatrix2D inv = alg.inverse(ldCorr);
        DoubleMatrix2D indChisequreMat = alg.mult(inv, chisequreMat);
        for (int i = 0; i < snpSize; i++) {
            chisequreMat.setQuick(i, 0, 1);
        }
        double df = 0;
        double Y = 0;
        for (int i = 0; i < snpSize; i++) {
            Y += indChisequreMat.getQuick(0, i);
        }

        df = 0;
        for (int i = 0; i < snpSize; i++) {
            for (int j = i + 1; j < snpSize; j++) {
                df += ldCorr.getQuick(i, j);
            }
        }
        df += df;
        df += snpSize;
        df = snpSize * snpSize / df;

        //calcualte the scalled chi 
        double p = Probability.chiSquareComplemented(df, Y);
        if (p < 0.0) {
            p = 0.0;
        }
        return p;
    }

    public double combinePValuebyCorrectedChiFisherCombinationTestJohnnyAll(double[] pValueArray, DoubleMatrix2D ldCorr) throws Exception {
        int snpSize = pValueArray.length;
        if (snpSize == 0) {
            return -1;
        } else if (snpSize == 1) {
            return pValueArray[0];
        }

        double[] chisquares1 = new double[snpSize];
        for (int i = 0; i < snpSize; i++) {
            chisquares1[i] = pValueArray[i] / 2;
        }
        double[] chisquares = MultipleTestingMethod.zScores(chisquares1);
        double df = 0;
        double Y = 0;
        for (int i = 0; i < snpSize; i++) {
            Y += chisquares[i] * chisquares[i];
        }

        df = 0;
        for (int i = 0; i < snpSize; i++) {
            for (int j = i + 1; j < snpSize; j++) {
                df += ldCorr.getQuick(i, j);
            }
        }
        df += df;
        df += snpSize;
        Y /= df;
        Y *= snpSize;
        df = snpSize * snpSize / df;

        //calcualte the scalled chi 
        double p = Probability.chiSquareComplemented(df, Y);
        if (p < 0.0) {
            p = 0.0;
        }
        return p;
    }

    public double combinePValuebyCorrectedChiFisherCombinationTestWeightJohnnyAll(PValueWeight[] pValueArray, DoubleMatrix2D ldCorr) throws Exception {
        int snpSize = pValueArray.length;
        if (snpSize == 0) {
            return -1;
        } else if (snpSize == 1) {
            return pValueArray[0].pValue;
        }

        double[] chisquares1 = new double[snpSize];
        for (int i = 0; i < snpSize; i++) {
            chisquares1[i] = pValueArray[i].pValue / 2;
        }
        double[] chisquares = MultipleTestingMethod.zScores(chisquares1);
        double df = 0;
        double Y = 0;
        double allWeight = 0;
        for (int i = 0; i < snpSize; i++) {
            Y += (pValueArray[i].var * chisquares[i] * chisquares[i]);
            allWeight += pValueArray[i].var;
        }

        df = 0;
        for (int i = 0; i < snpSize; i++) {
            df += (pValueArray[i].var * pValueArray[i].var);
            for (int j = i + 1; j < snpSize; j++) {
                df += (2 * pValueArray[i].var * pValueArray[j].var * ldCorr.getQuick(i, j));
            }
        }
        Y /= df;
        Y *= allWeight;
        df = allWeight * allWeight / df;

        //calcualte the scalled chi 
        double p = Probability.chiSquareComplemented(df, Y);
        if (p < 0.0) {
            p = 0.0;
        }
        return p;
    }

    public double combinePValuebyScaleedFisherCombinationTestCovLogP(PValueWeight[] pValueArray, DoubleMatrix2D ldCorr) throws Exception {
        int snpSize = pValueArray.length;
        if (snpSize == 0) {
            return -1;
        } else if (snpSize == 1) {
            return pValueArray[0].pValue;
        }
        DoubleMatrix2D covLogP = new DenseDoubleMatrix2D(snpSize, snpSize);
        for (int i = 0; i < snpSize; i++) {
            for (int j = i + 1; j < snpSize; j++) {
                if (i == j) {
                    continue;
                }
                double x = ldCorr.getQuick(pValueArray[i].index, pValueArray[j].index);

                //for r
                //y = 0.0079x3 + 3.9459x2 - 0.0024x ; y = 0.0331x3 + 3.9551x2 - 0.0156x
                //x = 0.0331 * (Math.pow(x1, 3)) + 3.9551 * (Math.pow(x1, 2)) - 0.0156 * x1;
                if (x < 0) {
                    x = -x;
                }
                // x1 = (0.75 * x1 + 3.25) * x1;

                //for r2                          
                //the new approximation
                // x1 = 8.595 * x1 * x1;
                //the formular I used in KGG
                if (x > 0.5) {
                    x = 0.75 * x + 3.25 * Math.sqrt(x);
                } else {
                    x = 8.6 * x;
                    //  x1 = 0.75 * x1 + 3.25 * Math.sqrt(x1);
                }

                covLogP.setQuick(i, j, x);
            }
        }

        double Y = 0;
        for (int i = 0; i < snpSize; i++) {
            Y += (-2 * Math.log(pValueArray[i].pValue) * pValueArray[i].var);
        }
        /*
        
         //calcualte the scalled chi
         double varTD = 4 * snpSize;
         for (int i = 0; i < snpSize; i++) {
         for (int j = i + 1; j < snpSize; j++) {
         varTD += 2 * covLogP.get(i, j);
         }
         }
         double c = varTD / (4 * snpSize);
         double f = 8 * snpSize * snpSize / varTD;
         double p1 = Probability.chiSquareComplemented(f, Y1 / c);
         */
        double sumWeight = 0;
        double sumWeight2 = 0;
        double sumCov = 0;
        for (int i = 0; i < snpSize; i++) {
            sumWeight += pValueArray[i].var;
        }
        //http://www.sciencedirect.com/science/article/pii/S016771520500009X 
        //A simple approximation for the distribution of the weighted combination of non-independent or independent probabilities Chia-Ding Hou

        for (int i = 0; i < snpSize; i++) {
            sumWeight2 += pValueArray[i].var * pValueArray[i].var;
        }

        for (int i = 0; i < snpSize; i++) {
            for (int j = i + 1; j < snpSize; j++) {
                sumCov += covLogP.get(i, j) * pValueArray[i].var * pValueArray[j].var;
            }
        }

        double c = (sumWeight2 + sumCov / 2) / sumWeight;
        double f = 4 * sumWeight * sumWeight / (2 * sumWeight2 + sumCov);
        double p = Probability.chiSquareComplemented(f, Y / c);
        if (p < 0.0) {
            p = 0.0;
        }
        return p;
    }

    public double combinePValuebyScaleedFisherCombinationTestCovLogPChisquare(PValueWeight[] pValueArray, DoubleMatrix2D ldCorr, double scale) throws Exception {
        int snpSize = pValueArray.length;
        if (snpSize == 0) {
            return -1;
        } else if (snpSize == 1) {
            return pValueArray[0].pValue;
        }

        DoubleMatrix2D covLogP = new DenseDoubleMatrix2D(snpSize, snpSize);
        for (int i = 0; i < snpSize; i++) {
            for (int j = i + 1; j < snpSize; j++) {
                if (i == j) {
                    continue;
                }
                double x = ldCorr.getQuick(pValueArray[i].index, pValueArray[j].index);

                //for r
                //y = 0.0079x3 + 3.9459x2 - 0.0024x ; y = 0.0331x3 + 3.9551x2 - 0.0156x
                //x = 0.0331 * (Math.pow(x1, 3)) + 3.9551 * (Math.pow(x1, 2)) - 0.0156 * x1;
                /*
                 if (x1 < 0) {
                 x1 = -x1;
                 }
                 x1 = (0.75 * x1 + 3.25) * x1;
                 * 
                 */
                //for r2
                x = x * x;
                x = 3.9471 * x;
                // x1 = 4 * x1;
                covLogP.setQuick(i, j, x);
            }
        }

        double Y = 0;
        for (int i = 0; i < snpSize; i++) {
            Y += (-2 * Math.log(pValueArray[i].pValue) * pValueArray[i].var);
        }
        /*
        
         //calcualte the scalled chi
         double varTD = 4 * snpSize;
         for (int i = 0; i < snpSize; i++) {
         for (int j = i + 1; j < snpSize; j++) {
         varTD += 2 * covLogP.get(i, j);
         }
         }
         double c = varTD / (4 * snpSize);
         double f = 8 * snpSize * snpSize / varTD;
         double p1 = Probability.chiSquareComplemented(f, Y1 / c);
         */
        double sumWeight = 0;
        double sumWeight2 = 0;
        double sumCov = 0;
        for (int i = 0; i < snpSize; i++) {
            sumWeight += pValueArray[i].var;
        }
        //http://www.sciencedirect.com/science/article/pii/S016771520500009X A simple approximation for the distribution of the weighted combination of non-independent or independent probabilities Chia-Ding Hou

        for (int i = 0; i < snpSize; i++) {
            sumWeight2 += pValueArray[i].var * pValueArray[i].var;
        }

        for (int i = 0; i < snpSize; i++) {
            for (int j = i + 1; j < snpSize; j++) {
                sumCov += covLogP.get(i, j) * pValueArray[i].var * pValueArray[j].var;
            }
        }

        double c = (sumWeight2 + sumCov / 2) / sumWeight;
        double f = 4 * sumWeight * sumWeight / (2 * sumWeight2 + sumCov);
        // Probability.chiSquareComplemented(f, Y1 / c);
        double chi = Math.pow(Y / c / f, 1 / scale) - (1 - 2 / (scale * scale * f));
        chi = scale * chi / Math.sqrt(2 / f);

        return chi;
    }

    public double combinePValuebyFisherCombinationTestChisquare(double[] pValueArray) throws Exception {
        int snpSize = pValueArray.length;
        if (snpSize == 0) {
            return -1;
        } else if (snpSize == 1) {
            return pValueArray[0];
        }

        double Y = 0;
        for (int i = 0; i < snpSize; i++) {
            double p = pValueArray[i];
            Y += (-2 * Math.log(p));
        }
        return Y;
    }

    public double combinePValuebyBonferroniCombinationTest(double[] pValueArray) throws Exception {
        int snpSize = pValueArray.length;
        if (snpSize == 0) {
            return 1;
        } else if (snpSize == 1) {
            return pValueArray[0];
        }

        double p = pValueArray[0];
        for (int i = 1; i < snpSize; i++) {
            if (pValueArray[i] < p) {
                p = pValueArray[i];
            }
        }
        return p * snpSize;
    }

    public double combinePValuebySidakCombinationTest(double[] pValueArray, double effectiveSampleSize) throws Exception {
        int snpSize = pValueArray.length;
        if (snpSize == 0) {
            return 1;
        } else if (snpSize == 1) {
            return pValueArray[0];
        }

        double p = pValueArray[0];
        for (int i = 1; i < snpSize; i++) {
            if (pValueArray[i] < p) {
                p = pValueArray[i];
            }
        }
        p = 1 - (Math.pow(1.0 - p, effectiveSampleSize));
        return p;
    }

    public double combinePValuebySidakCombinationTest(double[] pValueArray) throws Exception {
        int snpSize = pValueArray.length;
        if (snpSize == 0) {
            return 1;
        } else if (snpSize == 1) {
            return pValueArray[0];
        }

        double p = pValueArray[0];
        for (int i = 1; i < snpSize; i++) {
            if (pValueArray[i] < p) {
                p = pValueArray[i];
            }
        }
        p = 1 - (Math.pow(1.0 - p, snpSize));
        return p;
    }

    public double combinePValuebyWeightedSimeCombinationTest(PValueWeight[] pValueArray) throws Exception {
        int snpSize = pValueArray.length;
        if (snpSize == 0) {
            return 1;
        } else if (snpSize == 1) {
            return pValueArray[0].pValue;
        }
        Arrays.sort(pValueArray, new PValueWeightPComparator());
        double accumulatedWeight = pValueArray[0].var;
        double minP = snpSize * pValueArray[0].pValue / accumulatedWeight;
        double p;

        for (int i = 1; i < snpSize; i++) {
            accumulatedWeight += pValueArray[i].var;
            p = (snpSize * pValueArray[i].pValue / accumulatedWeight);
            if (p < minP) {
                minP = p;
            }
        }
        return minP;
    }

    public double combinePValuebySimeCombinationTest(double[] pValueArray) throws Exception {
        int snpSize = pValueArray.length;
        if (snpSize == 0) {
            return 1;
        } else if (snpSize == 1) {
            return pValueArray[0];
        }
        Arrays.sort(pValueArray);

        double minP = snpSize * pValueArray[0];
        double p;

        for (int i = 1; i < snpSize; i++) {
            p = snpSize * pValueArray[i] / (i + 1);
            if (p < minP) {
                minP = p;
            }
        }
        return minP;
    }

    public double combinePValuebyBiNormalTest(double[] pValueArray, DoubleMatrix2D ldCorr) throws Exception {
        int snpSize = pValueArray.length;
        if (snpSize == 0) {
            return 1;
        } else if (snpSize == 1) {
            return pValueArray[0];
        }
        Arrays.sort(pValueArray);

        //note only valid for 2 SNPs
        BiNormalDist bionrmal = new BiNormalDist(ldCorr.getQuick(0, 1));
        //two-tails
        double q1 = Probability.normalInverse(pValueArray[0] / 2);
        double q2 = Probability.normalInverse(pValueArray[1] / 2);
        q2 = q1;
        double p = bionrmal.cdf(q1, q2, -q1, -q2);

        return (1 - p);
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        PValuesAnalyzer pa = new PValuesAnalyzer();
        try {
            int lociNum = 3;

            double[][] ldr2 = new double[][]{{1, 0.110236, 0.055118, 0.200787, 0.728346, 0.507874, 0.311024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {
                0.110236, 1, 0, 0, 0.110236, 0, 0.098425, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {
                0.055118, 0, 1, 0, 0.055118, 0, 0.051181, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {
                0.200787, 0, 0, 1, 0.200787, 0.094488, 0.732283, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {
                0.728346, 0.110236, 0.055118, 0.200787, 1, 0.507874, 0.80315, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {
                0.507874, 0, 0, 0.094488, 0.507874, 1, 0.437008, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {
                0.311024, 0.098425, 0.051181, 0.732283, 0.80315, 0.437008, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {
                0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0.011811, 0, 0, 0, 0, 0, 0, 0.011811, 0, 0, 0, 0, 0, 0, 0.007874}, {
                0, 0, 0, 0, 0, 0, 0, 0, 1, 0.059055, 0.059055, 0, 0, 0, 0.051181, 0, 0.047244, 0.047244, 0, 0, 0, 0, 0, 0.047244, 0.047244, 0.047244, 0.047244, 0.035433, 0.03937, 0}, {
                0, 0, 0, 0, 0, 0, 0, 0, 0.059055, 1, 0.055118, 0, 0, 0, 0.047244, 0, 0.043307, 0.043307, 0, 0, 0, 0, 0, 0.043307, 0.043307, 0.043307, 0.043307, 0.035433, 0.035433, 0}, {
                0, 0, 0, 0, 0, 0, 0, 0, 0.059055, 0.055118, 1, 0, 0, 0, 0.051181, 0, 0.043307, 0.047244, 0, 0, 0, 0, 0, 0.043307, 0.047244, 0.043307, 0.043307, 0.035433, 0.03937, 0}, {
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.007874, 0.011811, 0, 0.015748, 0, 0, 0.011811, 0.011811, 0.011811, 0.011811, 0.015748, 0, 0, 0, 0, 0.007874, 0, 0.011811}, {
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.007874, 1, 0, 0, 0.007874, 0, 0, 0, 0, 0, 0, 0.007874, 0, 0, 0, 0, 0, 0, 0}, {
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.011811, 0, 1, 0, 0.007874, 0, 0, 0, 0, 0, 0, 0.007874, 0, 0, 0, 0, 0, 0, 0.007874}, {
                0, 0, 0, 0, 0, 0, 0, 0, 0.051181, 0.047244, 0.051181, 0, 0, 0, 1, 0, 0.055118, 0.059055, 0, 0, 0, 0, 0, 0.051181, 0.055118, 0.051181, 0.051181, 0.043307, 0.043307, 0}, {
                0, 0, 0, 0, 0, 0, 0, 0.011811, 0, 0, 0, 0.015748, 0.007874, 0.007874, 0, 1, 0, 0, 0.007874, 0.007874, 0.007874, 0.007874, 0.011811, 0, 0, 0, 0, 0, 0, 0.011811}, {
                0, 0, 0, 0, 0, 0, 0, 0, 0.047244, 0.043307, 0.043307, 0, 0, 0, 0.055118, 0, 1, 0.055118, 0, 0, 0, 0, 0, 0.051181, 0.051181, 0.051181, 0.051181, 0.03937, 0.03937, 0}, {
                0, 0, 0, 0, 0, 0, 0, 0, 0.047244, 0.043307, 0.047244, 0, 0, 0, 0.059055, 0, 0.055118, 1, 0, 0, 0, 0, 0, 0.051181, 0.055118, 0.051181, 0.051181, 0.043307, 0.043307, 0}, {
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.011811, 0, 0, 0, 0.007874, 0, 0, 1, 0, 0, 0, 0.007874, 0, 0, 0, 0, 0, 0, 0}, {
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.011811, 0, 0, 0, 0.007874, 0, 0, 0, 1, 0, 0, 0.007874, 0, 0, 0, 0, 0, 0, 0}, {
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.011811, 0, 0, 0, 0.007874, 0, 0, 0, 0, 1, 0, 0.007874, 0, 0, 0, 0, 0, 0, 0}, {
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.011811, 0, 0, 0, 0.007874, 0, 0, 0, 0, 0, 1, 0.007874, 0, 0, 0, 0, 0, 0, 0}, {
                0, 0, 0, 0, 0, 0, 0, 0.011811, 0, 0, 0, 0.015748, 0.007874, 0.007874, 0, 0.011811, 0, 0, 0.007874, 0.007874, 0.007874, 0.007874, 1, 0, 0, 0, 0, 0, 0, 0.011811}, {
                0, 0, 0, 0, 0, 0, 0, 0, 0.047244, 0.043307, 0.043307, 0, 0, 0, 0.051181, 0, 0.051181, 0.051181, 0, 0, 0, 0, 0, 1, 0.059055, 0.055118, 0.055118, 0.047244, 0.047244, 0}, {
                0, 0, 0, 0, 0, 0, 0, 0, 0.047244, 0.043307, 0.047244, 0, 0, 0, 0.055118, 0, 0.051181, 0.055118, 0, 0, 0, 0, 0, 0.059055, 1, 0.059055, 0.059055, 0.051181, 0.051181, 0}, {
                0, 0, 0, 0, 0, 0, 0, 0, 0.047244, 0.043307, 0.043307, 0, 0, 0, 0.051181, 0, 0.051181, 0.051181, 0, 0, 0, 0, 0, 0.055118, 0.059055, 1, 0.059055, 0.047244, 0.047244, 0}, {
                0, 0, 0, 0, 0, 0, 0, 0, 0.047244, 0.043307, 0.043307, 0, 0, 0, 0.051181, 0, 0.051181, 0.051181, 0, 0, 0, 0, 0, 0.055118, 0.059055, 0.059055, 1, 0.047244, 0.047244, 0}, {
                0, 0, 0, 0, 0, 0, 0, 0, 0.035433, 0.035433, 0.035433, 0.007874, 0, 0, 0.043307, 0, 0.03937, 0.043307, 0, 0, 0, 0, 0, 0.047244, 0.051181, 0.047244, 0.047244, 1, 0.043307, 0}, {
                0, 0, 0, 0, 0, 0, 0, 0, 0.03937, 0.035433, 0.03937, 0, 0, 0, 0.043307, 0, 0.03937, 0.043307, 0, 0, 0, 0, 0, 0.047244, 0.051181, 0.047244, 0.047244, 0.043307, 1, 0}, {
                0, 0, 0, 0, 0, 0, 0, 0.007874, 0, 0, 0, 0.011811, 0, 0.007874, 0, 0.011811, 0, 0, 0, 0, 0, 0, 0.011811, 0, 0, 0, 0, 0, 0, 1}
            };
            DoubleMatrix2D ldCorr = new DenseDoubleMatrix2D(ldr2);

            double[] pvalues = new double[]{0.124844607, 0.110661309, 0.822461041, 0.127453284, 0.417854711, 0.868397721, 0.067436906, 0.020542319, 0.009341871, 0.015883614, 0.012985321, 0.025744335, 0.012498142, 0.019169445, 0.015311359, 0.052160798, 0.054297026, 0.030886073, 0.021841838, 0.024432454, 0.013435841, 0.014925982, 0.029405817, 0.033456928, 0.047903982, 0.023374816, 0.017867578, 0.039527871, 0.104978108, 0.019891388};
            int startIndex = 7;
            int endIndex = 30;
            PValueWeight[] pvalueWeights = new PValueWeight[endIndex - startIndex];
            for (int i = startIndex; i < endIndex; i++) {
                pvalueWeights[i - startIndex] = new PValueWeight();
                pvalueWeights[i - startIndex].pValue = pvalues[i];
                pvalueWeights[i - startIndex].index = i - startIndex;
                pvalueWeights[i - startIndex].var = 1;
            }

            ldCorr = ldCorr.viewPart(startIndex, startIndex, endIndex - startIndex, endIndex - startIndex);
            // System.out.println(ldCorr.toString());
            System.out.println(pa.combinePValuebyScaleedFisherCombinationTestCovLogP(pvalueWeights, ldCorr));
            //for Pitman, Roger Keith analysis 
           

        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    public void testScaledChi() throws Exception {
        int lociNum = 10;
        final double PRECISION = 1.0e-8;
        double[] mean = new double[lociNum];
        Arrays.fill(mean, 0.0);

        RealMatrix covarianceMatrix = new Array2DRowRealMatrix(lociNum, lociNum);

        double roh = 0.5;

        for (int i = 0; i < lociNum; i++) {
            for (int j = 0; j < lociNum; j++) {
                if (i == j) {
                    covarianceMatrix.setEntry(i, j, 1);
                } else {
                    covarianceMatrix.setEntry(i, j, roh);
                }
            }
        }

        CorrelatedRandomVectorGenerator sg = new CorrelatedRandomVectorGenerator(mean, covarianceMatrix, PRECISION, new GaussianRandomGenerator(new MersenneTwister()));

        DoubleArrayList[] chisList = new DoubleArrayList[lociNum];
        DoubleArrayList[] pvList = new DoubleArrayList[lociNum];
        double[][] covLogP = new double[lociNum][lociNum];
        double[][] correPV = new double[lociNum][lociNum];

        for (int i = 0; i < lociNum; i++) {
            Arrays.fill(covLogP[i], 0);
            Arrays.fill(correPV[i], 0);
            chisList[i] = new DoubleArrayList();
            pvList[i] = new DoubleArrayList();
        }

        int popuSize = 500000;
        double currentP = 0.01;

        double adjP = 0;

        double varTD = 2 * lociNum;
        double meanTD = lociNum;
        double a1 = 3.263119, a2 = 0.709866, a3 = 0.026589, a4 = 0.709866;

        // double a1 = 0.0037, a2 = 3.9457, a3 = 0.0004, a4 = 0.709866;
        double X = 0;
        for (int i = 0; i < lociNum; i++) {
            for (int j = i + 1; j < lociNum; j++) {
                //y = 0.0037x3 + 3.9457x2 + 0.0004x
                //when here is r
                //x = 0.026589 * (Math.pow(x1, 3)) + 0.709866* (Math.pow(x1, 2)) + 3.263119 * x1;
                //x = 0.0037 * (Math.pow(x1, 3)) + 3.9457 * (Math.pow(x1, 2)) + 0.0004 * x1;
                double r = covarianceMatrix.getEntry(i, j);
                //r = Math.abs(r);
                // X += 2 * (a1 * r + a2 * r * r + a3 * r * r * r);
                //X += 2 * 2 * r * r;
                X += 0.180;
                //X += 2 * (4.0 * r * r);
                //y = 0.0565x3 + 0.0542x2 + 3.8891x
                //when here is r2
                //r = r * r;
                //X +=  2*(0.0565 * (Math.pow(r, 3)) + 0.0542 * (Math.pow(r, 2)) + 3.8891 * r);

                // varTD += 2 * (a1 * ldCorr.get(i, j) + a2 * ldCorr.get(i, j) * ldCorr.get(i, j) + a3 * ldCorr.get(i, j) * ldCorr.get(i, j) * ldCorr.get(i, j));
            }
        }
//roh
        varTD = X + varTD;
        double c = varTD / (2 * meanTD);
        double f = 2 * meanTD * meanTD / varTD;
        double sum1 = 0, sum2 = 0;

        double adjP1 = 0;
        DoubleMatrix2D ldCorr = new DenseDoubleMatrix2D(2, 2);
        double[] pvs = new double[lociNum];
        ChiSquaredDistributionImpl chisqare = new ChiSquaredDistributionImpl(1);
        DoubleArrayList list1 = new DoubleArrayList();
        DoubleArrayList list2 = new DoubleArrayList();

        /*
         pp <- scan("D:/home/mxli/MyJava/GenetSimulator/test.txt", quiet= TRUE);
         hist(pp, breaks=100);
         */
        BufferedWriter debugOut = new BufferedWriter(new FileWriter("test.txt"));

        for (int i = 0; i < popuSize; i++) {

            double[] samples = sg.nextVector();
            sum1 = 0;

            for (int j = 0; j < lociNum; j++) {
                // System.out.print(samples[j] + " ");
                double pt = 1 - Probability.normal(samples[j]);
                //System.out.println(samples[j]);

                sum1 += -2 * Math.log(pt);
                chisList[j].add(-2 * Math.log(pt));
            }
            list1.add(sum1);
            //System.out.println();

            double pt = 1;
            sum2 = 0;
            for (int j = 0; j < lociNum; j++) {
                // pvList[j].add(samples[j]);
                // sum2 +=samples[j];
                samples[j] = samples[j] * samples[j];
                pt = Probability.chiSquareComplemented(1, samples[j]);
                sum2 += Probability.normalInverse(1 - pt);

                //sum2 += -2 * Math.log(2 - 2*Probability.normal(Math.abs(samples[j])));
                //pvList[j].add(Probability.normalInverse(1-pt));
                pvList[j].add(Probability.normalInverse(1 - pt));
            }
            debugOut.write(String.valueOf(sum2));
            debugOut.newLine();
            list2.add(sum2);
            double p = 0;

            p = Probability.chiSquareComplemented(f, sum1 / c);
            if (p <= currentP) {
                adjP1 += 1;
            }
            //p = Probability.chiSquareComplemented(f, sum2 / c);
            p = Probability.normal(0, lociNum + 2 * X, sum2);
            if (p <= currentP) {
                adjP += 1;
            }
        }
        debugOut.close();

        System.out.println(adjP1 / popuSize);
        System.out.println(adjP / popuSize);

        System.out.println(Descriptive.mean(list1));
        System.out.println(Descriptive.sampleVariance(list1, Descriptive.mean(list1)));

        System.out.println(Descriptive.mean(list2));
        System.out.println(Descriptive.sampleVariance(list2, Descriptive.mean(list2)));

        for (int j = 0; j < lociNum; j++) {
            for (int s = j + 1; s < lociNum; s++) {
                double mean1 = Descriptive.mean(chisList[j]);
                double mean2 = Descriptive.mean(chisList[s]);
                double sd1 = Descriptive.sampleVariance(chisList[j], mean1);
                double sd2 = Descriptive.sampleVariance(chisList[s], mean2);

                //covLogP[j][s] += Descriptive.correlation(chisList[j], Math.sqrt(sd1), chisList[s], Math.sqrt(sd2));
                covLogP[j][s] += Descriptive.covariance(chisList[j], chisList[s]);
                mean1 = Descriptive.mean(pvList[j]);
                mean2 = Descriptive.mean(pvList[s]);
                sd1 = Descriptive.sampleVariance(pvList[j], mean1);
                sd2 = Descriptive.sampleVariance(pvList[s], mean2);
                //correPV[j][s] += Descriptive.correlation(pvList[j], Math.sqrt(sd1), pvList[s], Math.sqrt(sd2));
                correPV[j][s] += Descriptive.covariance(pvList[j], pvList[s]);
            }
        }
        for (int j = 0; j < lociNum; j++) {
            for (int s = 0; s < lociNum; s++) {
                System.out.print(String.valueOf(covLogP[j][s]));
                System.out.print("\t");
            }
            System.out.println();
        }

        for (int j = 0; j < lociNum; j++) {
            for (int s = 0; s < lociNum; s++) {
                System.out.print(String.valueOf(correPV[j][s]));
                System.out.print("\t");
            }
            System.out.println();
        }

        // System.out.println(crudeEp / simulationTime);
        //System.out.println(adjP / simulationTime);
    }

    public double veagsPNormal(double[] vars, DoubleMatrix2D ldCorr) throws Exception {
        //VEAGA orginal implemenation
        int lociNum = ldCorr.rows();
        CholeskyDecomposition cd = new CholeskyDecomposition(ldCorr);
        DoubleMatrix2D lowTriangle = cd.getL();
        //System.out.println(ltriangle.toString());
        final double PRECISION = 1.0e-8;
        double[] mean = new double[lociNum];
        Arrays.fill(mean, 0.0);

        RealMatrix covarianceMatrix = new Array2DRowRealMatrix(lociNum, lociNum);
        for (int i = 0; i < lociNum; i++) {
            for (int j = 0; j < lociNum; j++) {
                if (i != j) {
                    covarianceMatrix.setEntry(i, j, 0);
                } else {
                    covarianceMatrix.setEntry(i, j, 1);
                }
            }
        }

        double sumAB = 0;
        for (int i = 0; i < lociNum; i++) {
            sumAB += (vars[i] * vars[i]);
        }
        CorrelatedRandomVectorGenerator sg = new CorrelatedRandomVectorGenerator(mean, covarianceMatrix, PRECISION, new GaussianRandomGenerator(new MersenneTwister()));
        double crudeEp = 0;
        int popuSize = 10000;
        DoubleMatrix2D vector = new DenseDoubleMatrix2D(1, lociNum);
        DoubleMatrix2D tmpMatrix = new DenseDoubleMatrix2D(1, lociNum);
        for (int i = 0; i < popuSize; i++) {
            double[] samples = sg.nextVector();
            for (int j = 0; j < lociNum; j++) {
                vector.setQuick(0, j, samples[j]);
            }
            vector.zMult(lowTriangle, tmpMatrix);
            double sum = 0;
            for (int j = 0; j < lociNum; j++) {
                sum += (tmpMatrix.getQuick(0, j) * tmpMatrix.getQuick(0, j));
            }
            if (sumAB <= sum) {
                crudeEp += 1;
            }
        }
        //System.out.println(crudeEp / simulationTime);
        return crudeEp / popuSize;
    }

    /*
     library(mvtnorm)
     m <- 3
     sigma <- diag(3)
     sigma[2, 1] <- 3/5
     sigma[3, 1] <- 1/3
     sigma[3, 2] <- 11/15
     pmvnorm(mean = rep(0, m), sigma, lower = rep(-Inf, m), upper = c(1, 4, 2))
     *
     *
     *
     *
     */
    public double combinePValuebyMultiNormalTest(double[] pValueArray, DoubleMatrix2D ldCorr) throws Exception {
        int snpSize = pValueArray.length;
        if (snpSize == 0) {
            return 1;
        } else if (snpSize == 1) {
            return pValueArray[0];
        }
        Arrays.sort(pValueArray);
        //two-tails
        double[] a = new double[snpSize];
        for (int i = 0; i < snpSize; i++) {
            a[i] = Probability.normalInverse(pValueArray[i] / 2);
        }

        CholeskyDecomposition cd = new CholeskyDecomposition(ldCorr);
        DoubleMatrix2D ltriangle = cd.getL();
        System.out.println(ltriangle.toString());

        final double PRECISION = 1.0e-8;
        double alpha = 2.5;
        int iterMax = 100;
        double intSum = 0, varSum = 0;
        double d0 = Probability.normal(a[0] / ltriangle.getQuick(0, 0));
        double e0 = Probability.normal(-a[0] / ltriangle.getQuick(0, 0));
        double f0 = e0 - d0;
        double error = 1;

        double uniVariable = 0;
        cern.jet.random.engine.RandomEngine twister = new cern.jet.random.engine.MersenneTwister(new java.util.Date());
        double[] y = new double[snpSize];
        int iterNum = 0;
        double tmpD = 0;
        while (error >= PRECISION && iterNum < iterMax) {
            for (int i = 1; i < snpSize; ++i) {
                uniVariable = twister.raw();
                tmpD = d0 + uniVariable * f0;
                System.out.println(tmpD);
                y[i - 1] = Probability.normalInverse(tmpD);
                double sum = 0;
                for (int j = 0; j < i; j++) {
                    sum += ltriangle.getQuick(i, j) * y[j];
                }
                d0 = Probability.normal(a[i] - sum) / ltriangle.getQuick(i, i);
                e0 = Probability.normal(-a[i] - sum) / ltriangle.getQuick(i, i);
                f0 = (e0 - d0) * f0;
            }
            intSum += f0;
            varSum += f0 * f0;
            iterNum++;
            if (iterNum > 2) {
                error = alpha * Math.sqrt(((varSum / iterNum - (intSum * intSum / iterNum / iterNum)) / iterNum));
            }
        }

        double p = intSum / iterNum;
        return (1 - p) / 2;
    }

  
    public double combinePValuebyCorrectedChiFisherCombinationTestMXLiWeight(PValueWeight[] pValueArray, DoubleMatrix2D ldRMatrix) throws Exception {
        int snpSize = pValueArray.length;
        if (snpSize == 0) {
            return 1;
        } else if (snpSize == 1) {
            return pValueArray[0].pValue;
        }

        double[] chisquares1 = new double[snpSize];
        for (int i = 0; i < snpSize; i++) {
            chisquares1[i] = pValueArray[i].pValue / 2;
        }
        double[] chisquares = MultipleTestingMethod.zScores(chisquares1);
        double[][] chisquareArray = new double[snpSize][1];
        for (int i = 0; i < snpSize; i++) {
            chisquareArray[i][0] = chisquares[i] * chisquares[i];
        }

        DoubleMatrix2D chisequreMat = new DenseDoubleMatrix2D(chisquareArray);
        Algebra alg = new Algebra();
        for (int i = 0; i < snpSize; i++) {
            //* pValueArray[i].var
            ldRMatrix.setQuick(i, i, pValueArray[i].var * pValueArray[i].var);
            for (int j = i + 1; j < snpSize; j++) {
                //Math.sqrt(pValueArray[i].var * pValueArray[j].var) 
                ldRMatrix.setQuick(i, j, pValueArray[i].var * pValueArray[j].var * ldRMatrix.getQuick(i, j));
                ldRMatrix.setQuick(j, i, ldRMatrix.getQuick(i, j));
            }
        }

        DoubleMatrix2D inv = alg.inverse(ldRMatrix);
        DoubleMatrix2D indChisequreMat = alg.mult(inv, chisequreMat);

        for (int i = 0; i < snpSize; i++) {
            chisequreMat.setQuick(i, 0, 1);
        }
        DoubleMatrix2D indChiDfMat = alg.mult(inv, chisequreMat);
        double df = 0;
        double Y = 0;
        for (int i = 0; i < snpSize; i++) {
            Y += (pValueArray[i].var * indChisequreMat.getQuick(0, i));
            // Y1 += chisquareArray[i][0];
            // pValueArray[i].var *
            df += (pValueArray[i].var * indChiDfMat.getQuick(0, i));
        }

        //calcualte the scalled chi 
        double p = Probability.chiSquareComplemented(df, Y);
        if (p < 0.0) {
            p = 0.0;
        }
        return p;
    }

   
     
    public double[] combinePValuebyVEGAS(double totalTestStatistics, double bestStatistics, DoubleMatrix2D ldCorr) throws Exception {
        int snpSize = ldCorr.columns();
        if (snpSize == 0) {
            return new double[]{1, 1};
        }

        /*
         for (int i = 0; i < snpSize; i++) {
         for (int j = 0; j < snpSize; j++) {
         ldCorr.setQuick(i, j, ldCorr.getQuick(i, j) * ldCorr.getQuick(i, j));
         }
         }
         *
         */
        CholeskyDecomposition cd = new CholeskyDecomposition(ldCorr);
        DoubleMatrix2D lowTriangle = cd.getL();
        //System.out.println(ltriangle.toString());

        final double PRECISION = 1.0e-8;
        double[] mean = new double[snpSize];
        Arrays.fill(mean, 0.0);
        int simulationTime = 1000000;
        RealMatrix covarianceMatrix = new Array2DRowRealMatrix(snpSize, snpSize);
        DoubleMatrix2D vector = new DenseDoubleMatrix2D(1, snpSize);
        DoubleMatrix2D tmpMatrix = new DenseDoubleMatrix2D(1, snpSize);

        double empiricalP1 = 0.0;
        double empiricalP2 = 0.0;
        for (int i = 0; i < snpSize; i++) {
            covarianceMatrix.setEntry(i, i, 1);
        }
        double tmp = 0;
        double betStat = 0;

        CorrelatedRandomVectorGenerator sg = new CorrelatedRandomVectorGenerator(mean, covarianceMatrix, PRECISION, new GaussianRandomGenerator(new MersenneTwister()));
        for (int i = 0; i < simulationTime; i++) {
            double[] samples = sg.nextVector();
            for (int j = 0; j < snpSize; j++) {
                vector.setQuick(0, j, samples[j]);
            }

            vector.zMult(lowTriangle, tmpMatrix);

            double sum = 0;
            betStat = -1;
            for (int j = 0; j < snpSize; j++) {
                if (Double.isInfinite(tmpMatrix.getQuick(0, j)) || Double.isNaN(tmpMatrix.getQuick(0, j))) {
                    continue;
                }
                tmp = (tmpMatrix.getQuick(0, j) * tmpMatrix.getQuick(0, j));
                if (tmp >= betStat) {
                    betStat = tmp;
                }
                sum += tmp;
            }

            if (sum >= totalTestStatistics) {
                empiricalP1 += 1;
            }

            if (betStat >= bestStatistics) {
                empiricalP2 += 1;
            }
        }
        empiricalP1 += 1;
        empiricalP2 += 1;
        simulationTime += 1;
        return new double[]{empiricalP1 / simulationTime, empiricalP2 / simulationTime};
    }

    public double[] combinePValuebyVEGAS1(double totalTestStatistics, double bestStatistics, DoubleMatrix2D ldCorr) throws Exception {
        int snpSize = ldCorr.columns();
        if (snpSize == 0) {
            return new double[]{1, 1};
        }

        /*
        for (int i = 0; i < snpSize; i++) {
            for (int j = 0; j < snpSize; j++) {
                ldCorr.setQuick(i, j, Math.abs(ldCorr.getQuick(i, j)));
            }
        }
         */
        //CholeskyDecomposition cd = new CholeskyDecomposition(ldCorr);
        //DoubleMatrix2D lowTriangle = cd.getL();
        //System.out.println(ltriangle.toString());
        final double PRECISION = 1.0e-8;
        double[] mean = new double[snpSize];
        Arrays.fill(mean, 0.0);
        int simulationTime = 10000;

        double empiricalP1 = 0.0;
        double empiricalP2 = 0.0;

        double betStat = 0;
        double sum = 0;
        int[] seeds = new int[19];
        for (int i = 0; i < seeds.length; i++) {
            seeds[i] = (int) (Math.random() * 1000000);
            // System.out.println(seeds[i]);
        }
        WELL607 we = new WELL607();
        we.setSeed(seeds);
        NormalGen ng = new NormalGen(new MT19937(we));

        MultinormalCholeskyGen sg = new MultinormalCholeskyGen(ng, mean, ldCorr); //How about the alternative one MultinormalPCAGen?
        double[] samples = new double[snpSize];

        for (int i = 0; i < simulationTime; i++) {
            sg.nextPoint(samples);
            sum = 0;
            betStat = -1;
            for (int j = 0; j < snpSize; j++) {
                samples[j] = samples[j] * samples[j];
                if (samples[j] > betStat) {
                    betStat = samples[j];
                }
                sum += samples[j];
            }

            if (sum >= totalTestStatistics) {
                empiricalP1 += 1;
            }

            if (betStat >= bestStatistics) {
                empiricalP2 += 1;
            }
        }
        empiricalP1 += 1;
        empiricalP2 += 1;
        simulationTime += 1;
        return new double[]{empiricalP1 / simulationTime, empiricalP2 / simulationTime};
    }

    //Genetic Epidemiology 22:170?185 (2002)
    public double combinePValuebyZaykin(double[] pvalues, DoubleMatrix2D ldCorr) throws Exception {
        int snpSize = ldCorr.columns();
        if (snpSize == 0) {
            return 1;
        }
        DoubleMatrix2D zs = new DenseDoubleMatrix2D(pvalues.length, 1);
        for (int i = 0; i < snpSize; i++) {
            zs.setQuick(i, 0, MultipleTestingMethod.zScore(pvalues[i]));
        }

        CholeskyDecomposition cd = new CholeskyDecomposition(ldCorr);
        DoubleMatrix2D lowTriangle = cd.getL();

        DoubleMatrix2D zs1 = new DenseDoubleMatrix2D(pvalues.length, 1);
        lowTriangle.zMult(zs, zs1);
        double[] pp = new double[snpSize];
        for (int i = 0; i < snpSize; i++) {
            pp[i] = 1 - Probability.normal(zs1.getQuick(i, 0));
        }
        return combinePValuebyFisherCombinationTest(pp);
    }

    public double combinePValuebyBestSVEGAS1(double[] pValueArray, DoubleMatrix2D ldCorr) throws Exception {
        int snpSize = ldCorr.columns();
        if (snpSize == 0) {
            return 1;
        }
        double bestStatistics = 0;
        for (int i = 0; i < snpSize; i++) {
            double a = Probability.normalInverse(pValueArray[i] / 2);
            if (a * a > bestStatistics) {
                bestStatistics = (a * a);
            }
        }

        double betStat;
        CholeskyDecomposition cd = new CholeskyDecomposition(ldCorr);
        DoubleMatrix2D lowTriangle = cd.getL();
        //System.out.println(ltriangle.toString());

        final double PRECISION = 1.0e-8;
        double[] mean = new double[snpSize];
        Arrays.fill(mean, 0.0);
        int simulationTime = 10000;
        RealMatrix covarianceMatrix = new Array2DRowRealMatrix(snpSize, snpSize);
        DoubleMatrix2D vector = new DenseDoubleMatrix2D(1, snpSize);
        DoubleMatrix2D tmpMatrix = new DenseDoubleMatrix2D(1, snpSize);

        double empiricalP = 0.0;
        for (int i = 0; i < snpSize; i++) {
            covarianceMatrix.setEntry(i, i, 1);
        }
        CorrelatedRandomVectorGenerator sg = new CorrelatedRandomVectorGenerator(mean, covarianceMatrix, PRECISION, new GaussianRandomGenerator(new MersenneTwister()));
        for (int i = 0; i < simulationTime; i++) {
            double[] samples = sg.nextVector();

            for (int j = 0; j < snpSize; j++) {
                vector.setQuick(0, j, samples[j]);
            }
            vector.zMult(lowTriangle, tmpMatrix);
            betStat = 0;
            for (int j = 0; j < snpSize; j++) {
                if (Double.isInfinite(tmpMatrix.getQuick(0, j)) || Double.isNaN(tmpMatrix.getQuick(0, j))) {
                    continue;
                }
                if (tmpMatrix.getQuick(0, j) * tmpMatrix.getQuick(0, j) > betStat) {
                    betStat = tmpMatrix.getQuick(0, j) * tmpMatrix.getQuick(0, j);
                }
            }

            if (betStat >= bestStatistics) {
                empiricalP += 1;
            }
            //System.out.println(tmpMatrix.toString());
        }
        return empiricalP / simulationTime;
    }

    //note this function is very sensitive to the non positive definite covariance matrix
    public double combinePValuebyBestSVEGAS(double[] pValueArray, DoubleMatrix2D ldCorr) throws Exception {
        int snpSize = ldCorr.columns();
        if (snpSize == 0) {
            return 1;
        }
        double bestStatistics = 0;
        for (int i = 0; i < snpSize; i++) {
            double a = Probability.normalInverse(pValueArray[i] / 2);
            if (a * a > bestStatistics) {
                bestStatistics = (a * a);
            }
        }

        final double PRECISION = 1.0e-8;
        double[] mean = new double[snpSize];
        Arrays.fill(mean, 0.0);
        int simulationTime = 10000;
        RealMatrix covarianceMatrix = new Array2DRowRealMatrix(snpSize, snpSize);

        double empiricalP = 0.0;

        for (int i = 0; i < snpSize; i++) {
            for (int j = 0; j < snpSize; j++) {
                covarianceMatrix.setEntry(i, j, ldCorr.getQuick(i, j));
            }
        }

        double betStat = 0;
        CorrelatedRandomVectorGenerator sg = new CorrelatedRandomVectorGenerator(mean, covarianceMatrix, PRECISION, new GaussianRandomGenerator(new MersenneTwister()));
        for (int i = 0; i < simulationTime; i++) {
            double[] samples = sg.nextVector();
            betStat = samples[0] * samples[0];
            for (int j = 1; j < snpSize; j++) {
                samples[j] = samples[j] * samples[j];
                if (samples[j] > betStat) {
                    betStat = samples[j];
                }
            }
            if (betStat >= bestStatistics) {
                empiricalP += 1;
            }
            //System.out.println(tmpMatrix.toString());
        }
        return empiricalP / simulationTime;
    }

    public double combinePValuebyMultiChi(double[] pValueArray, DoubleMatrix2D ldCorr) throws Exception {
        int snpSize = pValueArray.length;
        if (snpSize == 0) {
            return 1;
        } else if (snpSize == 1) {
            return pValueArray[0];
        }

        //two-tails
        double[] a = new double[snpSize];
        for (int i = 0; i < snpSize; i++) {
            a[i] = Probability.normalInverse(pValueArray[i] / 2);
        }

        for (int i = 0; i < snpSize; i++) {
            for (int j = 0; j < snpSize; j++) {
                //  ldCorr.setQuick(j, j, ldCorr.getQuick(j, j) * ldCorr.getQuick(j, j));
            }
        }

        Algebra algebra = new Algebra();
        DoubleMatrix2D inverseCorr = algebra.inverse(ldCorr);
        //System.out.println(inverseCorr.toString());
        DoubleMatrix2D vector1 = new DenseDoubleMatrix2D(1, snpSize);
        DoubleMatrix2D vector2 = new DenseDoubleMatrix2D(snpSize, 1);
        for (int j = 0; j < snpSize; j++) {
            vector1.setQuick(0, j, a[j]);
            vector2.setQuick(j, 0, a[j]);
        }

        DoubleMatrix2D middle = new DenseDoubleMatrix2D(1, snpSize);
        vector1.zMult(inverseCorr, middle);
        //System.out.println(middle.toString());
        DoubleMatrix2D result = new DenseDoubleMatrix2D(1, 1);
        middle.zMult(vector2, result);
        //System.out.println(result.toString());
        //System.out.println(result.getQuick(0, 0));
        double p = Probability.chiSquareComplemented(snpSize, (result.getQuick(0, 0)));

        /*
         CholeskyDecomposition cd = new CholeskyDecomposition(ldCorr);
         DoubleMatrix2D lowTriangle = cd.getL();
         DoubleMatrix2D tmpMatrix = new DenseDoubleMatrix2D(1, snpSize);
         */
        return p;
    }

    public double combineNormalByMultiChi(double[] normalArray, DoubleMatrix2D ldCorr) throws Exception {
        int snpSize = normalArray.length;
        if (snpSize == 0) {
            return 1;
        }

        for (int i = 0; i < snpSize; i++) {
            for (int j = 0; j < snpSize; j++) {
                // ldCorr.setQuick(j, j, ldCorr.getQuick(j, j) * ldCorr.getQuick(j, j));
            }
        }

        Algebra algebra = new Algebra();
        DoubleMatrix2D inverseCorr = algebra.inverse(ldCorr);
        //System.out.println(inverseCorr.toString());
        DoubleMatrix2D vector1 = new DenseDoubleMatrix2D(1, snpSize);
        DoubleMatrix2D vector2 = new DenseDoubleMatrix2D(snpSize, 1);
        for (int j = 0; j < snpSize; j++) {
            vector1.setQuick(0, j, normalArray[j]);
            vector2.setQuick(j, 0, normalArray[j]);
        }

        DoubleMatrix2D middle = new DenseDoubleMatrix2D(1, snpSize);
        vector1.zMult(inverseCorr, middle);
        //System.out.println(middle.toString());
        DoubleMatrix2D result = new DenseDoubleMatrix2D(1, 1);
        middle.zMult(vector2, result);
        //System.out.println(result.toString());
        // System.out.println(result.getQuick(0, 0));
        double p = Probability.chiSquareComplemented(snpSize, (result.getQuick(0, 0)));

        /*
         CholeskyDecomposition cd = new CholeskyDecomposition(ldCorr);
         DoubleMatrix2D lowTriangle = cd.getL();
         DoubleMatrix2D tmpMatrix = new DenseDoubleMatrix2D(1, snpSize);
         */
        return p;
    }
}
