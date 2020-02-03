/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.genetsimulator.controller;

import cern.colt.list.DoubleArrayList;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import cern.jet.stat.Descriptive;
import cern.jet.stat.Probability;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.random.CorrelatedRandomVectorGenerator;
import org.apache.commons.math.random.GaussianRandomGenerator;
import org.apache.commons.math.random.MersenneTwister;
import umontreal.ssj.probdistmulti.BiNormalDist;
 

/**
 *
 * @author mxli
 */
public class EffectiveNumberStudy {

    public double calculateEffectSampleSize(DoubleMatrix2D covarianceMatrix) throws Exception {
        //I found this function is less error-prone  than the  EigenDecompositionImpl 2.0 and slightly faster
        //System.out.println(poweredCorrMat.toString());
        EigenvalueDecomposition ed = new EigenvalueDecomposition(covarianceMatrix);

        DoubleMatrix1D eVR = ed.getRealEigenvalues();

        double effectSampleSize = covarianceMatrix.columns();
        int realSampleSize = covarianceMatrix.columns();
        for (int i = 0; i < realSampleSize; i++) {
            if (Double.isNaN(eVR.get(i))) {
                System.err.println("NaN error for eigen values!");
            }
            if (eVR.getQuick(i) > 1) {
                effectSampleSize -= (eVR.getQuick(i) - 1);
            }
        }
        return (effectSampleSize);
    }

    /*
     * path<-"D:/home/mxli/MyJava/GenetSimulator/debug.txt";
    dat<-read.table(path,sep = "\t",header=TRUE);

    dat <- scan(path, quiet= TRUE);
    hist(dat,breaks=1000)
     */
    public void lookintoProperties(RealMatrix covarianceMatrix, double alpha) throws Exception {
        int lociNum = covarianceMatrix.getColumnDimension();
        final double PRECISION = 1.0e-8;
        double[] mean = new double[lociNum];
        Arrays.fill(mean, 0.0);
        int popuSize = 5000000;

        double max = -Probability.normalInverse(alpha / 2);
        BufferedWriter debugOut = new BufferedWriter(new FileWriter("debug.txt"));
        CorrelatedRandomVectorGenerator sg = new CorrelatedRandomVectorGenerator(mean, covarianceMatrix, PRECISION, new GaussianRandomGenerator(new MersenneTwister()));
        for (int i = 0; i < popuSize; i++) {
            double[] samples = sg.nextVector();
            if (Math.abs(samples[0]) < max) {
                debugOut.write(String.valueOf(samples[1]));
                debugOut.newLine();
            }
        }
        debugOut.close();

    }

    public EffectiveNumberStudy() {
    }

    public void validateVaryAlpha() throws Exception {
        double rho = -0.8;
        double x = rho;
        //when r2
        //y = y = 2.3894x6 - 6.6051x5 + 6.5648x4 - 2.6244x3 - 0.3279x2 + 1.605x
        // x = (((((2.3894 * x - 6.6051) * x + 6.5648) * x - 2.6244) * x - 0.3279) * x + 1.605) * x;

        BiNormalDist binorD = new BiNormalDist(x);

        double apha = .05 / 2;
        int grid = 100;
        double delta = apha / grid;
        double inc = delta;
        while (inc <= apha) {
            double x1 = Probability.normalInverse(inc);
            double y1 = Probability.normalInverse(inc);
            double cdf = binorD.cdf(-x1, -y1) - binorD.cdf(x1, -y1) - binorD.cdf(-x1, y1) + binorD.cdf(x1, y1);
            System.out.println(2 * inc + "\t" + (1 - cdf) + "\t" + Math.log(cdf) / Math.log(1 - 2 * inc));
            //System.out.println(inc+"\t"+cdf + "\t" +Math.pow(1 - 2 * inc,2));
            //cdf = 1 - cdf;
            //System.out.println(cdf + " " + cdf / (2 * apha1));
            inc += delta;
        }
        x = rho * rho;
        //when r2
        //y = 0.7723x6 - 1.5659x5 + 1.201x4 - 0.2355x3 + 0.2184x2 + 0.6086x
        x = (((((0.7723 * x - 1.5659) * x + 1.201) * x - 0.2355) * x + 0.2184) * x + 0.6086) * x;
        pValueCorrelation(rho);
        System.out.println(x);

        DoubleMatrix2D corrMat = new DenseDoubleMatrix2D(2, 2);

        corrMat.setQuick(0, 0, 1);
        corrMat.setQuick(1, 1, 1);
        corrMat.setQuick(0, 1, x);
        corrMat.setQuick(1, 0, x);
        EigenvalueDecomposition ed = new EigenvalueDecomposition(corrMat);
        // System.out.println(corrMat.toString());
        DoubleMatrix1D eVR = ed.getRealEigenvalues();

        double effectSampleSize = 2;

        for (int i = 0; i < 2; i++) {
            if (Double.isNaN(eVR.get(i))) {
                System.err.println("NaN error for eigen values!");
            }
            if (eVR.getQuick(i) > 1) {
                effectSampleSize -= (eVR.getQuick(i) - 1);//(eVR.getQuick(j));
            }
        }
        System.out.println("Estimate: " + effectSampleSize);
    }

    public void pValueCorrelation(double rho) {
        double[] mean = new double[2];
        Arrays.fill(mean, 0.0);
        int num = 1000000;
        final double PRECISION = 1.0e-8;
        RealMatrix covarianceMatrix = new Array2DRowRealMatrix(2, 2);
        covarianceMatrix.setEntry(0, 0, 1);
        covarianceMatrix.setEntry(1, 1, 1);
        covarianceMatrix.setEntry(0, 1, rho);
        covarianceMatrix.setEntry(1, 0, rho);
        DoubleArrayList var1 = new DoubleArrayList();
        DoubleArrayList var2 = new DoubleArrayList();
        double corre = 0;
        try {
            CorrelatedRandomVectorGenerator sg = new CorrelatedRandomVectorGenerator(mean, covarianceMatrix, PRECISION, new GaussianRandomGenerator(new MersenneTwister()));
            for (int i = 0; i < num; i++) {
                double[] samples = sg.nextVector();
                var1.add(Probability.chiSquareComplemented(1, samples[0] * samples[0]));
                var2.add(Probability.chiSquareComplemented(1, samples[1] * samples[1]));
                //note: samples[0] can be negative
                // var1.add(2 - 2*Probability.normal(samples[0]));
                //var2.add(2 - 2*Probability.normal(samples[1]));
            }
            double mean1 = Descriptive.mean(var1);
            double mean2 = Descriptive.mean(var2);
            double sd1 = Descriptive.sampleVariance(var1, mean1);
            double sd2 = Descriptive.sampleVariance(var2, mean2);

            corre = Descriptive.correlation(var1, Math.sqrt(sd1), var2, Math.sqrt(sd2));
            System.out.println("The correlation is " + corre);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(-1);
        }
    }

    public void validateVaryRho() throws Exception {
        double alpha = 1E-8;
        double apha1 = 0;
        int grid = 10000;
        double delta = 1.0 / grid;
        double rho = 0;
        BiNormalDist binorD = null;
        DoubleMatrix2D corrMat = new DenseDoubleMatrix2D(2, 2);
        alpha = 0.1;
        List<StringBuilder> results = new ArrayList<StringBuilder>();
        for (int j = 0; j <= grid; j++) {
            StringBuilder sb = new StringBuilder();
            sb.append(rho);
            results.add(sb);
            rho += delta;
        }

        StringBuilder head = new StringBuilder();
        head.append('r');
        while (alpha >= 1E-10) {
            alpha *= 0.5;
            head.append('\t');
            head.append(alpha);
            apha1 = alpha / 2;
            rho = 0;
            for (int j = 0; j <= grid; j++) {
                if (rho > 1) {
                    rho = 1;
                }
                double x = rho;
                binorD = new BiNormalDist(x);
                double x1 = Probability.normalInverse(apha1);
                double y1 = Probability.normalInverse(apha1);
                double cdf = binorD.cdf(-x1, -y1) - binorD.cdf(x1, -y1) - binorD.cdf(-x1, y1) + binorD.cdf(x1, y1);
                //two tailed
                // System.out.println(cdf + " " + cdf / (2* apha1));
                x = rho * rho;
                //when r2
                //y = 0.7723x6 - 1.5659x5 + 1.201x4 - 0.2355x3 + 0.2184x2 + 0.6086x
                //x = (((((0.7723 * x - 1.5659) * x + 1.201) * x - 0.2355) * x + 0.2184) * x + 0.6086) * x;
                corrMat.setQuick(0, 0, 1);
                corrMat.setQuick(1, 1, 1);
                corrMat.setQuick(0, 1, x);
                corrMat.setQuick(1, 0, x);
                EigenvalueDecomposition ed = new EigenvalueDecomposition(corrMat);
                //System.out.println(corrMat.toString());
                DoubleMatrix1D eVR = ed.getRealEigenvalues();

                double effectSampleSize = 2;
                for (int i = 0; i < 2; i++) {
                    if (Double.isNaN(eVR.get(i))) {
                        System.err.println("NaN error for eigen values!");
                    }
                    if (eVR.getQuick(i) > 1) {
                        effectSampleSize -= (eVR.getQuick(i) - 1);//(eVR.getQuick(j));
                    }
                }
                results.get(j).append('\t');
                results.get(j).append(Math.log(cdf) / Math.log(1 - 2 * apha1) / effectSampleSize);
                //System.out.println(rho + "\t" + Math.log(cdf) / Math.log(1 - 2 * apha1) / effectSampleSize);
                //System.out.println(rho + "\t" + cdf / (1 - Math.pow(1 - 2 * apha1, effectSampleSize)));
                rho += delta;
            }
        }

        BufferedWriter debugOut = new BufferedWriter(new FileWriter("debug.txt"));
        //System.out.println(head);
        debugOut.write(head.toString());
        debugOut.newLine();
        for (int j = 0; j <= grid; j++) {
            //System.out.println(results.get(j));
            debugOut.write(results.get(j).toString());
            debugOut.newLine();
        }
        debugOut.close();
    }

    public static void main(String[] args) {
        try {
            int originalSampleSize = 7;
          // RealMatrix corMat = ApacheMatrixBasic.readUpperTriangleMatrixFromFile("test.txt", originalSampleSize, originalSampleSize);
            EffectiveNumberStudy enf = new EffectiveNumberStudy();
            double rho = 0.5;
            double fwer = 0.05;
            int locusNum = 2;
            DoubleMatrix2D corrMat = new DenseDoubleMatrix2D(locusNum, locusNum);


            int rowNum = corrMat.rows(), colNum = corrMat.columns();
            for (int i = 0; i < rowNum; i++) {
                corrMat.setQuick(i, i, 1);
                for (int j = i + 1; j < colNum; j++) {
                    corrMat.setQuick(i, j, rho);
                    corrMat.setQuick(j, i, rho);
                }
            }


          //  double realEffectNum1 = enf.iterativeEstimator(corrMat, fwer);
            enf.validateVaryRho();
            RealMatrix corMat=new Array2DRowRealMatrix(2, 2);
            corMat.setEntry(0, 0, 1);
            corMat.setEntry(1, 1, 1);
            corMat.setEntry(0, 1, rho);
            corMat.setEntry(1, 0, rho);
            enf.lookintoProperties(corMat, fwer);
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    /*
     * RealMatrix mat = new Array2DRowRealMatrix(2, 2);
    mat.setEntry(0, 0, 1);
    mat.setEntry(1, 1, 1);
    mat.setEntry(0, 1, rho);
    mat.setEntry(1, 0, rho);
    enf.lookintoProperties(mat, fwer);
     */
    public double iterativeEstimator(DoubleMatrix2D corrMat, double fwer) throws Exception {
        DoubleMatrix2D uniqueCorrMat = removeRedundantItems(corrMat, 1);
        double effectiveNum1 = uniqueCorrMat.columns();
        double effectiveNum2 = effectiveNum1;
        double threshold;
        double tolerate = 0.0001;
        double cdf = 0;
        BiNormalDist binorD = new BiNormalDist(corrMat.getQuick(0, 1));
        do {
            effectiveNum1 = effectiveNum2;
            //System.out.println(effectiveNum1);
            threshold = 1 - Math.pow(1 - fwer, 1 / effectiveNum1);
            double a1 = Probability.normalInverse(threshold / 2);
            cdf = 1 - (binorD.cdf(-a1, -a1) - binorD.cdf(a1, -a1) - binorD.cdf(-a1, a1) + binorD.cdf(a1, a1));
            System.out.println(cdf);
            DoubleMatrix2D adjustedCorrMat = adjustPairwiseCorrelation(uniqueCorrMat, threshold);
            effectiveNum2 = calculateEffectSampleSize(adjustedCorrMat);
        } while (Math.abs(effectiveNum1 - effectiveNum2) > tolerate);

        return effectiveNum2;
    }

    public DoubleMatrix2D adjustPairwiseCorrelation(DoubleMatrix2D corMat, double alpha) throws Exception {
        BiNormalDist binorD = new BiNormalDist(0);
        DoubleMatrix2D adjustedMatrix = corMat.copy();
        int rowNum = corMat.rows(), colNum = corMat.columns();
        double apha1 = alpha / 2;
        double a1 = Probability.normalInverse(apha1);
        double realEffectNum;
        double cdf;
        double rho;
        for (int i = 0; i < rowNum; i++) {
            for (int j = i + 1; j < colNum; j++) {
                binorD = new BiNormalDist(corMat.getQuick(i, j));
                cdf = binorD.cdf(-a1, -a1) - binorD.cdf(a1, -a1) - binorD.cdf(-a1, a1) + binorD.cdf(a1, a1);
                //realEffectNum = Math.log(binorD.cdf(a1, -a1, a1, -a1)) / Math.log(1 - apha1);
                realEffectNum = Math.log(cdf) / Math.log(1 - alpha);
                rho = 2 - realEffectNum;
                //rho = rho * rho;
                adjustedMatrix.setQuick(i, j, rho);
                adjustedMatrix.setQuick(j, i, rho);
            }
        }

        return adjustedMatrix;
    }

    public static DoubleMatrix2D removeRedundantItems(DoubleMatrix2D corrMat, double maxCorr) {
        int originalSampleSize = corrMat.columns();
        int newSampleSize = originalSampleSize;
        Set<Integer> highlyCorrIndexes = new HashSet<Integer>();

        for (int i = 0; i < originalSampleSize; i++) {
            for (int j = i + 1; j < originalSampleSize; j++) {
                if (Math.abs(corrMat.getQuick(i, j)) >= maxCorr) {
                    if (!highlyCorrIndexes.contains(j) && !highlyCorrIndexes.contains(i)) {
                        highlyCorrIndexes.add(j);
                        //  System.out.println(i + " <-> " + j);
                    }
                }
            }
        }

        if (highlyCorrIndexes.size() > 0) {
            // System.out.println("Removed columns and rows: " + highlyCorrIndexes.toString());
            newSampleSize = originalSampleSize - highlyCorrIndexes.size();

            DoubleMatrix2D poweredCorrMat = new DenseDoubleMatrix2D(newSampleSize, newSampleSize);
            int incRow = 0;
            int incCol = 0;
            for (int i = 0; i < originalSampleSize; i++) {
                if (highlyCorrIndexes.contains(i)) {
                    continue;
                }
                incCol = 0;
                for (int j = 0; j < originalSampleSize; j++) {
                    if (highlyCorrIndexes.contains(j)) {
                        continue;
                    }
                    poweredCorrMat.setQuick(incRow, incCol, corrMat.getQuick(i, j));
                    incCol++;
                }
                incRow++;
            }

            // System.out.println(corrMat.toString());
            return poweredCorrMat;
        } else {
            return corrMat.copy();
        }

    }
}
