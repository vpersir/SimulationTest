/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.stat;

import cern.colt.bitvector.BitVector;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.CholeskyDecomposition;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.StringTokenizer;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.probdistmulti.BiNormalDist;
import umontreal.ssj.randvar.NormalGen;
import umontreal.ssj.randvarmulti.MultinormalGen;
import umontreal.ssj.rng.MT19937;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.rng.WELL607;
 
/**
 *
 * @author mxli
 */
public class MultiNormalRandGenerator {

    final double PRECISION = 1.0e-8;
    DoubleMatrix2D jointProb = null;
    DoubleMatrix2D covarianceMatrix = null;
    double[] quantileProb = null;

    public DoubleMatrix2D getJointProb() {
        return jointProb;
    }

    public void setJointProb(DoubleMatrix2D jointProb) {
        this.jointProb = jointProb;
    }

    public DoubleMatrix2D getCovarianceMatrix() {
        return covarianceMatrix;
    }

    public void setCovarianceMatrix(DoubleMatrix2D covarianceMatrix) {
        this.covarianceMatrix = covarianceMatrix;
    }

    public double[] getQuantileProb() {
        return quantileProb;
    }

    public void setQuantileProb(double[] quantileProb) {
        this.quantileProb = quantileProb;
    }

    public MultiNormalRandGenerator() {
    }

    void stat(BitVector[] s) {
        DoubleMatrix2D testJointProb = new DenseDoubleMatrix2D(3, 3);//How about the initialization
        int size = 0;
        double conditionalProb12 = 0.0;
        double conditionalProb1 = 0.0;
        double conditionalProb2 = 0.0;
        double count = 0.0;
        for (int r = 0; r < s.length; r++) {
            size = s[r].size();
            for (int c = 0; c < size; c++) {
                for (int d = c; d < size; d++) {
                    if (s[r].getQuick(c) && s[r].getQuick(d)) {
                        testJointProb.setQuick(c, d, testJointProb.getQuick(c, d) + 1);
                    }
                }
            }
            if (s[r].getQuick(0)) {
                if (s[r].getQuick(1) && s[r].getQuick(2)) {
                    conditionalProb12 += 1.0;
                }
                if (s[r].getQuick(1)) {
                    conditionalProb1 += 1.0;
                }
                if (s[r].getQuick(2)) {
                    conditionalProb2 += 1.0;
                }
                count += 1.0;
            }
        }

        for (int c = 0; c < size; c++) {
            for (int d = c; d < size; d++) {
                //testJointProb.multiplyEntry(c, d, 1.0 / s.length);
                testJointProb.setQuick(c, d, testJointProb.getQuick(c, d) / s.length);//Right?
            }
        }
        // LocalString.print2DRealMatrix(testJointProb);//No such function in LocalString Class!!!

        System.out.println(conditionalProb12 / count - conditionalProb1 * conditionalProb2 / count / count);
        System.out.println(testJointProb.getQuick(1, 2) - testJointProb.getQuick(1, 1) * testJointProb.getQuick(2, 2));

    }

    public void generateMultiNormalRandomVector() throws Exception {
        double[] mean = new double[3];
        Arrays.fill(mean, 0.0);
        int num = 1000000;
        BitVector[] genotypes = new BitVector[num];
        RandomStream rs = new MT19937(new WELL607());//What's about the other random stream?Such as  F2NL607, GenF2w32, LFSR113, LFSR258, MRG31k3p, MRG32k3a, MRG32k3aL, MT19937, RandMrg?

        try {
            MultinormalGen sg = new MultinormalGen(new NormalGen(rs), mean.length);
            for (int i = 0; i < num; i++) {
                //double[] samples = sg.nextVector();
                double[] samples = new double[mean.length];
                sg.nextPoint(samples);
                genotypes[i] = new BitVector(3);
                for (int j = 0; j < 3; j++) {
                    if (samples[j] <= quantileProb[j]) {
                        genotypes[i].putQuick(j, true);
                    }
                }
            }
            stat(genotypes);
        } catch (Exception e) {
            System.exit(-1);
        }
    }

    public void readJointProb(String filePath) throws Exception {
        String delmiliter = "\t  ,";
        BufferedReader br = new BufferedReader(new FileReader(filePath));
        String line = null;
        int rowNum = 0;
        int colNum = -1;
        List<Double> tmpVector = new ArrayList<Double>();
        while ((line = br.readLine()) != null) {
            StringTokenizer tokenizer = new StringTokenizer(line, delmiliter);
            if (colNum < 0) {
                colNum = tokenizer.countTokens();
            }
            while (tokenizer.hasMoreTokens()) {
                tmpVector.add(new Double(tokenizer.nextToken().trim()));
            }
            rowNum++;
        }
        br.close();
        int i = 0;
        int len = tmpVector.size();
        jointProb = new DenseDoubleMatrix2D(rowNum, colNum);
        for (i = 0; i < len; i++) {
            jointProb.setQuick(i / colNum, i % colNum, tmpVector.get(i));
        }

        tmpVector.clear();
        tmpVector = null;
        System.out.println(" A matrix with " + rowNum + " rows and " + colNum + " columns has been read from " + filePath);
        // LocalString.print2DRealMatrix(jointProb);

    }

    public boolean exploreCovarMatrix() throws Exception {
        //int rowNum = jointProb.getRowDimension();
        int rowNum = jointProb.rows();
        if (rowNum != jointProb.columns()) {
            throw (new Exception("unequal number of row and column for the joint probability matrix!"));
        }
        covarianceMatrix = new DenseDoubleMatrix2D(rowNum, rowNum);
        NormalDist normalDis = new NormalDist();
        quantileProb = new double[rowNum];
        for (int i = 0; i < rowNum; i++) {
            covarianceMatrix.setQuick(i, i, 1.0);//Why the diagonal elements are set 1? 
            quantileProb[i] = normalDis.inverseF(jointProb.getQuick(i, i));
        }

        for (int i = 0; i < rowNum; i++) {
            for (int j = i + 1; j < rowNum; j++) {
                double element = findRho(quantileProb[i], quantileProb[j], jointProb.getQuick(i, j));
                covarianceMatrix.setQuick(i, j, element);
                covarianceMatrix.setQuick(j, i, element);
            }
        }
        // System.out.println("The covariance matrix with " + rowNum + " rows and " + rowNum + " columns has been generated!");
        CholeskyDecomposition cdMatrix = new CholeskyDecomposition(covarianceMatrix);
        if (!cdMatrix.isSymmetricPositiveDefinite()) {
            // LocalString.print2DRealMatrix(covarianceMatrix);
            System.out.println("covarianceMatrix is not PositiveDefinite! Approximating an PositiveDefinite Matrix!");
            covarianceMatrix = approximatePositiveDefiniteMatrix(covarianceMatrix);
            return false;
        }
        return true;
        // System.out.println(covarianceMatrix.toString());
    }

    private DoubleMatrix2D covarianceMat2correlationMat(final DoubleMatrix2D covMat, DoubleMatrix2D sd) throws Exception { //This method is weird! 
        /*
        "covarianceMat2correlationMat" <-
        function(cov.mat){
        dc <- diag(cov.mat)
        So <- sqrt(dc)
        ds <- diag(1/So)
        Mo <- ds %*% cov.mat %*% ds
        Mo <- Mo - 0.5*(Mo-t(Mo)) # Correct any round-off error
        return(list(mat=Mo,sd=So))
        }
         */
        //DoubleMatrix2D dc = ApacheMatrixBasic.getDiagonal(covMat);

        DoubleMatrix2D dc = SSJAFC.getDiagonal(covMat);
//        for (int i = 0; i < dc.rows(); i++) {
//            sd.setRowVector(i, dc.viewRow(i).mapSqrt());
//        }
        sd = SSJAFC.getSqrt(dc);
        DoubleMatrix2D ds = SSJAFC.getDiagonal(SSJAFC.getInv(sd));
        //DoubleMatrix2D mO = ds.multiply(covMat).multiply(ds);
        //DoubleMatrix2D m0=new DenseDoubleMatrix2D(ds.rows(),ds.columns());
        DoubleMatrix2D m0 = ds.zMult(covMat, null).zMult(ds, null);
        //m0.zMult(ds, m0);
        //m0 = m0.subtract((m0.subtract(m0.viewDice())).scalarMultiply(0.5));
        DoubleMatrix2D result = SSJAFC.subtract(m0, SSJAFC.scalarMultiply(SSJAFC.subtract(m0, m0.viewDice()), 0.5));
        return result;
    }

    private DoubleMatrix2D sumsqscale(final DoubleMatrix2D mat) throws Exception {
        /*
        "sumsqscale" <-
        function(mat){
        s <- sign(mat)
        sq <- mat^2
        ssq <- apply(sq, 2, sum)
        for (i in 1:ncol(mat)) sq[,i] <- sq[,i] / ssq[i]
        return (s*sqrt(sq))
        }
         */

        //DoubleMatrix2D s = ApacheMatrixBasic.mapSignum(mat);
        DoubleMatrix2D s = SSJAFC.signum(mat);
        //DoubleMatrix2D sq = ApacheMatrixBasic.mapPow(mat, 2);
        DoubleMatrix2D sq = SSJAFC.pow(mat, 2);
        //RealVector ssq = ApacheMatrixBasic.sumColumn(sq);
        // DoubleMatrix1D ssq=SSJAFC.sumColumn(sq);
        //int num = sq.getColumnDimension();
//        int num=sq.columns();
//        for (int i = 0; i < num; i++) {
//            sq=SSJAFC.setColumnVector(sq, SSJAFC.divide(sq.viewColumn(i), ssq.getQuick(i)),i);
//        }
        sq = SSJAFC.divideColumnSum(sq);  //sq is not changed!
        //DoubleMatrix2D sumsqscale = ApacheMatrixBasic.ebeMultiply(s, SSJAFC.getSqrt(sq));
        DoubleMatrix2D sumsqscale = SSJAFC.ebeMultiply(s, SSJAFC.getSqrt(sq));
        return sumsqscale;
    }

    private DoubleMatrix2D correlationMat2covarianceMat(final DoubleMatrix2D covMat, final DoubleMatrix2D sd) throws Exception {
        /*
        "correlationMat2covarianceMat" <-
        function(cor.mat, sd){
        dc <- diag(cor.mat)
        ds <- diag(sd)
        Mo <- ds %*% cor.mat %*% ds
        Mo <- Mo - 0.5 * (Mo-t(Mo)) # Correct any round-off error
        Mo <- Mo - diag(diag(Mo)-sd * sd)
        return(Mo)
        }
         */
        // RealMatrix dc = ApacheMatrixBasic.getDiagonal(covMat);
        DoubleMatrix2D ds = SSJAFC.getDiagonal(sd);

        // DoubleMatrix2D m0 = ds.multiply(covMat).multiply(ds);
        DoubleMatrix2D m0 = ds.zMult(covMat, null).zMult(ds, null);//Should make sure the expression. 
        //mO = mO.subtract((mO.subtract(mO.transpose())).scalarMultiply(0.5));
        m0 = SSJAFC.subtract(m0, SSJAFC.scalarMultiply(SSJAFC.subtract(m0, m0.viewDice()), 0.5));
        //mO = mO.subtract(ApacheMatrixBasic.getDiagonal(ApacheMatrixBasic.getDiagonal(mO).subtract(ApacheMatrixBasic.ebeMultiply(sd, sd))));
        m0 = SSJAFC.subtract(m0, SSJAFC.getDiagonal(SSJAFC.subtract(SSJAFC.getDiagonal(m0), SSJAFC.ebeMultiply(sd, sd))));
        return m0;
    }

    public DoubleMatrix2D approximatePositiveDefiniteMatrix(DoubleMatrix2D mat) throws Exception {
        /*
        "makepd" <-
        function(mat, eig.tol=1.0000e-06){
        cor.mat <- covarianceMat2correlationMat(mat)
        mat <- cor.mat$mat
        eig <- eigenValues(mat,sym = TRUE, EISPACK = TRUE)
        D <- sort(eig$values)
        S <- eig$vectors[,order(sort(D, decreasing=TRUE ))]
        D[which(D<eig.tol)] <- eig.tol
        for (i in 1:ncol(S)) S[,i] <- S[,i] * D[i]
        B <- t(sumsqscale(t(S)))
        A <- B %*% t(B)
        return (correlationMat2covarianceMat(A, cor.mat$sd))
        }
         */

        int minDim = Math.min(mat.rows(), mat.columns());
        DoubleMatrix2D sd = new DenseDoubleMatrix2D(minDim, minDim);
        for (int i = 0; i < sd.rows(); i++) {
            sd.setQuick(i, i, mat.getQuick(i, i));
        }
        DoubleMatrix2D corMat = covarianceMat2correlationMat(mat, sd);

//eigenValues
        EigenvalueDecomposition eigenDeco = new EigenvalueDecomposition(corMat);
        double[] eigenValues = eigenDeco.getRealEigenvalues().toArray();
        //sort eigenValues value and vectors
        double tmpMax = 0;
        int tmpMaxIndex = 0;
        int num = eigenValues.length;
        DoubleMatrix2D eigenVectors = new DenseDoubleMatrix2D(num, num);
        for (int i = 0; i < num; i++) {
            //eigenVectors.setColumnVector(i, eigenDeco.getEigenvector(i));
            SSJAFC.setColumnVector(eigenVectors, eigenDeco.getV().viewColumn(i), i);
        }
        eigenDeco = null;//?
        for (int i = num - 1; i >= 0; i--) {
            tmpMax = eigenValues[i];
            tmpMaxIndex = i;
            for (int j = i - 1; j >= 0; j--) {
                if (eigenValues[j] > tmpMax) {
                    tmpMax = eigenValues[j];
                    tmpMaxIndex = j;
                }
            }
            if (tmpMaxIndex != i) {
                eigenValues[tmpMaxIndex] = eigenValues[i];
                eigenValues[i] = tmpMax;
                DoubleMatrix1D tmpMatrix = eigenVectors.viewColumn(i);
                SSJAFC.setColumnVector(eigenVectors, eigenVectors.viewColumn(tmpMaxIndex), i);
                SSJAFC.setColumnVector(eigenVectors, tmpMatrix, tmpMaxIndex);
            }
            if (eigenValues[i] < PRECISION) {
                eigenValues[i] = PRECISION;
            }
            for (int k = 0; k < num; k++) {
                //eigenVectors.multiplyEntry(k, i, eigenValues[i]);
                SSJAFC.multiplyEntry(eigenVectors, k, i, eigenValues[i]);
            }
        }

        DoubleMatrix2D B = sumsqscale(eigenVectors.viewDice()).viewDice();
        DoubleMatrix2D A = B.zMult(B.viewDice(), null);
        //LocalString.print2DRealMatrix(A);
        //System.out.println();

        DoubleMatrix2D result = correlationMat2covarianceMat(A, sd);
        //LocalString.print2DRealMatrix(result);
        //System.out.println();
        return result;
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here

        try {
            double[][] data = {{1.0000000, 0.0569067, -0.0504999, 0.0612860},
            {0.0569067, 1.0000000, 0.9999990, 0.9993868},
            {-0.0504999, 0.9999990, 1.0000000, 0.9999990},
            {0.0612860, 0.9993868, 0.9999990, 1.0000000}};
            DoubleMatrix2D mat = new DenseDoubleMatrix2D(data);
            EigenvalueDecomposition eigenDeco = new EigenvalueDecomposition(mat);
            double[][] d = {{1, 0, 0, 0},
            {0, 1, 0, 0},
            {0, 0, 1, 0},
            {0, 0, 0, 1}};

            DoubleMatrix2D E = new DenseDoubleMatrix2D(d);

            double[] eigenValues = eigenDeco.getRealEigenvalues().toArray();
            //DoubleMatrix2D eigenVector = new DenseDoubleMatrix2D(eigenDeco.getEigenvector(0).getData());
            DoubleMatrix2D eigenVector = eigenDeco.getV();
            // LocalString.print2DRealMatrix(eigenVector);
            //eigenVector = ApacheMatrixBasic.mapMultiply(eigenVector, -1);
            // LocalString.print2DRealMatrix(eigenVector);
            //DoubleMatrix2D re = E.scalarMultiply(eigenValues[0]).subtract(mat).multiply(eigenVector);
            DoubleMatrix2D re = SSJAFC.subtract(SSJAFC.scalarMultiply(E, eigenValues[0]), mat).zMult(eigenVector, null);
            // LocalString.print2DRealMatrix(re);
        } catch (Exception ex) {
            ex.printStackTrace();
        }

    }

    /* solves the equation for rho */
    double findRho(double x, double y, double c) throws Exception {
        double a = -1;
        double b = 1;
        double m = 0;
        BiNormalDist biNormDist = new BiNormalDist(0.0);
        double w = biNormDist.cdf(x, y, 1);
        while ((b - a) > PRECISION) {
            m = (a + b) / 2;
            if (biNormDist.cdf(x, y, m) == c) {
                break;
            }
            if (biNormDist.cdf(x, y, m) < c) {
                a = m;
            } else {
                b = m;
            }
        }
        return m;
    }

}
