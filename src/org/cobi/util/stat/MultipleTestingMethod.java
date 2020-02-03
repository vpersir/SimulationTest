/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.stat;

import cern.colt.list.DoubleArrayList;
import cern.jet.stat.Descriptive;
import cern.jet.stat.Probability;
import org.rosuda.REngine.Rserve.RConnection;

/**
 *
 * @author mxli
 */
public class MultipleTestingMethod {

    public void nnlsSolver(RConnection rcon, double[] A, double[] b, double[] finalX1, double[] finalX2) throws Exception {
        int size = b.length / 2;
        rcon.assign("A", A);
        rcon.voidEval("A<-matrix(A, nrow=" + size + ", ncol=" + size + ", byrow = TRUE)");
        rcon.assign("b", b);
        rcon.voidEval("b<-matrix(b, nrow=" + size + ", ncol=" + 2 + ", byrow = TRUE)");
        rcon.voidEval("beta <- nnlm(A, b, loss = 'mkl')");
        double[][] coeffs = rcon.eval("beta$coefficients").asDoubleMatrix();
        for (int i = 0; i < size; i++) {
            finalX1[i] = coeffs[i][0];
            finalX2[i] = coeffs[i][1];
        }

    }

    //would be very slow
    public static double iterativeChisquareInverse(double df, double p) {
        double chil = 0;
        //this is the maximal value for  df 2 under the computer precise 
        double chih = 1432;
        double precise = p / 100000000;
        double p1 = 1;
        double chi = 0;
        do {
            chi = (chil + chih) / 2;
            p1 = Probability.chiSquareComplemented(df, chi);
            if (p1 < p) {
                chih = chi;
            } else if (p1 > p) {
                chil = chi;
            }
        } while (Math.abs(p1 - p) > precise);
        return chi;
    }

    public static double combinationHeterogeneityCochranQTest(final double[] pValues) {
        double p = 1;
        double Q = 0;
        double q = 1;
        double meanQ = 0;
        int size = pValues.length;
        double[] qValues = new double[size];
        // assume they are two-tailed I2-values
        for (int i = 0; i < size; i++) {
            q = 1 - pValues[i];
            //the Probability.normalInverse could handle 1E-323 but cannot handle 1-(1E-323)
            if (q > 0.5) {
                q = 1 - q;
                if (q < 1E-323) {
                    q = 1E-323;
                }
                q = Probability.normalInverse(q);
                qValues[i] = -q;
            } else {
                if (q < 1E-323) {
                    q = 1E-323;
                }
                q = Probability.normalInverse(q);
                qValues[i] = q;
            }

            //Cochran's Q statistic  http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0000841
            // meanQ += qValues[i];
        }

        Q = qValues[0] - qValues[1];
        Q = Q * Q / 2;
        p = Probability.chiSquareComplemented(1, Q);
        //Quantifying heterogeneity in a meta-analysis Julian P. T. Higgins?; ? and Simon G. Thompson  I 2 = 100%��(Q ? df)/Q http://onlinelibrary.wiley.com/doi/10.1002/sim.1186/pdf       
        if (p < 0) {
            p = 0;
        }
        return p;
    }

    public static double[] zScores(final double[] pValues) {
        double pValue = 1;

        int size = pValues.length;
        double[] qValues = new double[size];
        // assume they are two-tailed I2-values
        for (int i = 0; i < size; i++) {
            pValue = pValues[i];
            //the Probability.normalInverse could handle 1E-323 but cannot handle 1-(1E-323)
            if (pValue > 0.5) {
                pValue = 1 - pValue;
                if (pValue < 1E-323) {
                    pValue = 1E-323;
                }
                qValues[i] = Probability.normalInverse(pValue);
            } else {
                if (pValue < 1E-323) {
                    pValue = 1E-323;
                }
                qValues[i] = -Probability.normalInverse(pValue);
            }
            // double x=  Probability.chiSquareComplemented(1, qValues[i]* qValues[i]);
            // int sss=0;
        }

        return qValues;
    }

    public static void correctInflationFactor(DoubleArrayList pValues) {
        double pValue = 1;

        int size = pValues.size();
        DoubleArrayList qValues = new DoubleArrayList();
        // assume they are two-tailed I2-values
        for (int i = 0; i < size; i++) {
            pValue = pValues.getQuick(i) / 2;
            //the Probability.normalInverse could handle 1E-323 but cannot handle 1-(1E-323)
            if (pValue > 0.5) {
                pValue = 1 - pValue;
                if (pValue < 1E-323) {
                    pValue = 1E-323;
                }
                pValue = Probability.normalInverse(pValue);
                qValues.add(pValue * pValue);
            } else {
                if (pValue < 1E-323) {
                    pValue = 1E-323;
                }
                pValue = -Probability.normalInverse(pValue);
                qValues.add(pValue * pValue);
            }
        }
        qValues.quickSort();
        double median = Descriptive.median(qValues);
        double expectedMedian = 0.456;
        median = median / expectedMedian;
        System.out.println("Lambda:" + median);
        median = 1.25;
        pValues.clear();
        for (int i = 0; i < size; i++) {
            pValue = qValues.getQuick(i);
            pValue = pValue / median;
            pValues.add(Probability.chiSquareComplemented(1, pValue));
        }

    }

    public static double zScore(double pValue) {
        //the Probability.normalInverse could handle 1E-323 but cannot handle 1-(1E-323)
        if (pValue < 0.5) {
            if (pValue < 1E-323) {
                pValue = 1E-323;
            }
            pValue = -Probability.normalInverse(pValue);
        } else {
            pValue = 1 - pValue;
            pValue = Probability.normalInverse(pValue);
        }
        return pValue;
    }

    public static double combinationHeterogeneityI2(final double[] pValues) {
        double I2 = 1;
        double Q = 0;
        double q = 1;
        double meanQ = 0;
        int size = pValues.length;
        double[] qValues = new double[size];
        // assume they are two-tailed I2-values
        for (int i = 0; i < size; i++) {
            q = 1 - pValues[i];
            //the Probability.normalInverse could handle 1E-323 but cannot handle 1-(1E-323)
            if (q > 0.5) {
                q = 1 - q;
                if (q < 1E-323) {
                    q = 1E-323;
                }
                q = Probability.normalInverse(q);
                qValues[i] = -q;
            } else {
                if (q < 1E-323) {
                    q = 1E-323;
                }
                q = Probability.normalInverse(q);
                qValues[i] = q;
            }
            //Cochran's Q statistic  http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0000841
            // meanQ += qValues[i];
        }

        Q = qValues[0] - qValues[1];
        Q = Q * Q / 2;

        //Quantifying heterogeneity in a meta-analysis Julian P. T. Higgins?; ? and Simon G. Thompson  I 2 = 100%��(Q ? df)/Q http://onlinelibrary.wiley.com/doi/10.1002/sim.1186/pdf
        I2 = 1 - 1 / Q;
        if (I2 < 0) {
            I2 = 0;
        }
        return I2;
    }
}
