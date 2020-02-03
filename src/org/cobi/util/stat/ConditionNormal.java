/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.stat;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.jet.stat.Probability;

/**
 *
 * @author mxli
 */
public class ConditionNormal {

    Algebra algebra = new Algebra();

    public ConditionNormal() {
    }

    /**
    ondNormal <- function(x.given, mu, sigma, given.ind, req.ind){
    # Returns conditional mean and variance of x[req.ind] 
    # Given x[given.ind] = x.given
    # where X is multivariate Normal with
    # mean = mu and covariance = sigma
    # 
    B <- sigma[req.ind, req.ind]
    C <- sigma[req.ind, given.ind, drop=FALSE]
    D <- sigma[given.ind, given.ind]
    CDinv <- C %*% solve(D)
    cMu <- c(mu[req.ind] + CDinv %*% (x.given - mu[given.ind]))
    cVar <- B - CDinv %*% t(C)
    list(condMean=cMu, condVar=cVar)
    }
    
    n <- 10
    E <- matrix(0.8,n^2,n,n)
    diag(E) <- 1
    condNormal(x=rep(2.3,n-1), mu=rep(0,n), sigma=E, req=n, given=seq(1,n-1))
     */
    public static void main(String[] args) {
        double[][] cov = {
            {1, 0.75, 0.1},
            {0.75, 1, 0.6},
            {0.1, 0.6, 1}};

        double[] x = {1E-10, 1E-1, 1E-1};
        int partCol = 1;

        ConditionNormal condNor = new ConditionNormal();
        DoubleMatrix2D condMean = new DenseDoubleMatrix2D(1, partCol);
        DoubleMatrix2D condVar = new DenseDoubleMatrix2D(partCol, partCol);
        condNor.conditionMeanCov(cov, x, partCol, condMean, condVar);
        System.out.println(condMean);
        System.out.println(condVar);
    }

    public void conditionMeanCov(double[][] cov, double[] x, int partCol, DoubleMatrix2D condMean, DoubleMatrix2D condVar) {
        x = MultipleTestingMethod.zScores(x);
        //transform r2 to coeff of z
        int totalCol = cov.length;
        //y = 0.622x6 - 1.3728x5 + 1.3263x4 - 0.5189x3 + 0.3495x2 + 0.5912x
        for (int i = 0; i < totalCol; i++) {
            cov[i][i] = 1;
            for (int j = i + 1; j < totalCol; j++) {
                cov[i][j] = (((((0.622 * cov[i][j] - 1.3728) * cov[i][j] + 1.3263) * cov[i][j] - 0.5189) * cov[i][j] + 0.3495) * cov[i][j] + 0.5912) * cov[i][j];
                cov[j][i] = cov[i][j];
            }
        }

        DoubleMatrix2D covarianceMatrix = new DenseDoubleMatrix2D(cov);
        System.out.println(covarianceMatrix);


        DoubleMatrix2D X = new DenseDoubleMatrix2D(x.length - partCol, 1);
        for (int i = partCol; i < x.length; i++) {
            X.setQuick(i - partCol, 0, x[i]);
        }

        DoubleMatrix2D part1 = covarianceMatrix.viewPart(0, 0, partCol, partCol);
        DoubleMatrix2D part2 = covarianceMatrix.viewPart(partCol, partCol, totalCol - partCol, totalCol - partCol);
        DoubleMatrix2D part12 = covarianceMatrix.viewPart(0, partCol, partCol, totalCol - partCol);


        DoubleMatrix2D ipart2 = algebra.inverse(part2);
        DoubleMatrix2D interM = algebra.mult(part12, ipart2);
        DoubleMatrix2D condVar1 = minus(part1, algebra.mult(interM, part12.viewDice()));
        condVar.assign(condVar1);

        DoubleMatrix2D condMean1 = algebra.mult(interM, X);
        condMean.assign(condMean1);


        double sd = Math.sqrt(condVar.getQuick(0, 0));
        double q5 = Probability.normalInverse(0.025);
        double q95 = Probability.normalInverse(0.975);
        q5 = q5 * sd + condMean.getQuick(0, 0);
        q95 = q95 * sd + condMean.getQuick(0, 0);
        q5 = 1 - Probability.normal(q5);
        q95 = 1 - Probability.normal(q95);

        System.out.println(q5 + ", " + q95 + ", " + (1 - Probability.normal(condMean.getQuick(0, 0), condVar.getQuick(0, 0), x[0])));

        condMean.setQuick(0, 0, 1 - Probability.normal(condMean.getQuick(0, 0)));
    }

    public double conditionP(DoubleMatrix2D cov, double[] x, int partCol, double[] ci, double[] normalParams) {
        x = MultipleTestingMethod.zScores(x);
        DoubleMatrix2D covarianceMatrix = cov.like();

        int totalCol = cov.columns();
        double v = 0;
        //transform r2 to coeff of z
        //y = 0.622x6 - 1.3728x5 + 1.3263x4 - 0.5189x3 + 0.3495x2 + 0.5912x
        for (int i = 0; i < totalCol; i++) {
            covarianceMatrix.setQuick(i, i, 1);
            for (int j = i + 1; j < totalCol; j++) {
                v = cov.getQuick(i, j);
                v = (((((0.5863 * v - 1.1934) * v + 1.036) * v - 0.32) * v + 0.2906) * v + 0.5982) * v;
                //v = (((((0.622 * v - 1.3728) * v + 1.3263) * v - 0.5189) * v + 0.3495) * v + 0.5912) * v;
                covarianceMatrix.setQuick(i, j, v);
                covarianceMatrix.setQuick(j, j, v);
            }
        }


        // System.out.println(covarianceMatrix);


        DoubleMatrix2D X = new DenseDoubleMatrix2D(x.length - 1, 1);
        for (int i = partCol; i < x.length; i++) {
            X.setQuick(i - partCol, 0, x[i]);
        }

        DoubleMatrix2D part1 = covarianceMatrix.viewPart(0, 0, partCol, partCol);
        DoubleMatrix2D part2 = covarianceMatrix.viewPart(partCol, partCol, totalCol - partCol, totalCol - partCol);
        DoubleMatrix2D part12 = covarianceMatrix.viewPart(0, partCol, partCol, totalCol - partCol);


        DoubleMatrix2D ipart2 = algebra.inverse(part2);
        DoubleMatrix2D interM = algebra.mult(part12, ipart2);
        DoubleMatrix2D condVar = minus(part1, algebra.mult(interM, part12.viewDice()));


        DoubleMatrix2D condMean = algebra.mult(interM, X);


        double sd = Math.sqrt(condVar.getQuick(0, 0));
        ci[0] = Probability.normalInverse(0.025);
        ci[1] = Probability.normalInverse(0.975);

        ci[0] = ci[0] * sd + condMean.getQuick(0, 0);
        ci[1] = ci[1] * sd + condMean.getQuick(0, 0);

        ci[0] = 1 - Probability.normal(ci[0]);
        ci[1] = 1 - Probability.normal(ci[1]);

        v = 1 - Probability.normal(condMean.getQuick(0, 0));
        double x1 = condVar.getQuick(0, 0);
        if (x1 < 0.6) {
            //use slope between actual p and imputed p to adjust the imputed p
            //y = 0.4024x3 - 0.8002x2 - 0.5379x + 1.0037
            //x = ((0.4024 * x - 0.8002) * x - 0.5379) * x + 1.0037;

            //y = 0.5068x3 - 0.9796x2 - 0.4623x + 0.9996
            x1 = ((0.5068 * x1 - 0.9796) * x1 - 0.4623) * x1 + 0.9996;
            v = Math.log10(v) / x1;
            v = Math.pow(10, v);
            normalParams[0] = MultipleTestingMethod.zScore(v);
        } else {
            normalParams[0] = condMean.getQuick(0, 0);
        }
        normalParams[1] = sd;
        return (v);

    }

    public double conditionPQC(DoubleMatrix2D cov, double[] x, int partCol, double[] ci) {
        x = MultipleTestingMethod.zScores(x);
        DoubleMatrix2D covarianceMatrix = cov.like();

        int totalCol = cov.columns();
        double v = 0;
        //transform r2 to coeff of z
        //y = 0.622x6 - 1.3728x5 + 1.3263x4 - 0.5189x3 + 0.3495x2 + 0.5912x
        for (int i = 0; i < totalCol; i++) {
            covarianceMatrix.setQuick(i, i, 1);
            for (int j = i + 1; j < totalCol; j++) {
                v = cov.getQuick(i, j);
                v = (((((0.5863 * v - 1.1934) * v + 1.036) * v - 0.32) * v + 0.2906) * v + 0.5982) * v;
                //v = (((((0.622 * v - 1.3728) * v + 1.3263) * v - 0.5189) * v + 0.3495) * v + 0.5912) * v;
                covarianceMatrix.setQuick(i, j, v);
                covarianceMatrix.setQuick(j, j, v);
            }
        }


        // System.out.println(covarianceMatrix);


        DoubleMatrix2D X = new DenseDoubleMatrix2D(x.length, 1);
        for (int i = 0; i < x.length; i++) {
            X.setQuick(i, 0, x[i]);
        }

        DoubleMatrix2D part1 = covarianceMatrix.viewPart(0, 0, partCol, partCol);
        DoubleMatrix2D part2 = covarianceMatrix.viewPart(partCol, partCol, totalCol - partCol, totalCol - partCol);
        DoubleMatrix2D part12 = covarianceMatrix.viewPart(0, partCol, partCol, totalCol - partCol);


        DoubleMatrix2D ipart2 = algebra.inverse(part2);
        DoubleMatrix2D interM = algebra.mult(part12, ipart2);
        DoubleMatrix2D condVar = minus(part1, algebra.mult(interM, part12.viewDice()));


        DoubleMatrix2D condMean = algebra.mult(interM, X);


        double sd = Math.sqrt(condVar.getQuick(0, 0));
        ci[0] = Probability.normalInverse(0.025);
        ci[1] = Probability.normalInverse(0.975);

        ci[0] = ci[0] * sd + condMean.getQuick(0, 0);
        ci[1] = ci[1] * sd + condMean.getQuick(0, 0);

        ci[0] = 1 - Probability.normal(ci[0]);
        ci[1] = 1 - Probability.normal(ci[1]);





        return (1 - Probability.normal(condMean.getQuick(0, 0)));

    }

    public double conditionVar(DoubleMatrix2D cov, int partCol) {

        DoubleMatrix2D covarianceMatrix = cov.like();
        int totalCol = cov.columns();
        double v = 0;
        //transform r2 to coeff of z
        //y = 0.622x6 - 1.3728x5 + 1.3263x4 - 0.5189x3 + 0.3495x2 + 0.5912x
        for (int i = 0; i < totalCol; i++) {
            covarianceMatrix.setQuick(i, i, 1);
            for (int j = i + 1; j < totalCol; j++) {
                v = cov.getQuick(i, j);
                v = (((((0.5863 * v - 1.1934) * v + 1.036) * v - 0.32) * v + 0.2906) * v + 0.5982) * v;
                //v = (((((0.622 * v - 1.3728) * v + 1.3263) * v - 0.5189) * v + 0.3495) * v + 0.5912) * v;
                covarianceMatrix.setQuick(i, j, v);
                covarianceMatrix.setQuick(j, j, v);
            }
        }


        // System.out.println(covarianceMatrix);

        DoubleMatrix2D part1 = covarianceMatrix.viewPart(0, 0, partCol, partCol);
        DoubleMatrix2D part2 = covarianceMatrix.viewPart(partCol, partCol, totalCol - partCol, totalCol - partCol);
        DoubleMatrix2D part12 = covarianceMatrix.viewPart(0, partCol, partCol, totalCol - partCol);


        DoubleMatrix2D ipart2 = algebra.inverse(part2);
        DoubleMatrix2D interM = algebra.mult(part12, ipart2);
        DoubleMatrix2D condVar = minus(part1, algebra.mult(interM, part12.viewDice()));


        return (condVar.getQuick(0, 0));

    }

    /**
     * C = A + B
     *
     * @param B another matrix
     * @return A + B
     */
    public DoubleMatrix2D plus(DoubleMatrix2D A, DoubleMatrix2D B) {
        int rowNum = A.rows();
        int colNum = A.columns();
        DoubleMatrix2D C = A.like(rowNum, colNum);
        for (int i = 0; i < rowNum; i++) {
            for (int j = 0; j < rowNum; j++) {
                C.setQuick(i, j, A.getQuick(i, j) + B.getQuick(i, j));
            }
        }
        return C;
    }

    /**
     * C = A - B
     *
     * @param B another matrix
     * @return A - B
     */
    public DoubleMatrix2D minus(DoubleMatrix2D A, DoubleMatrix2D B) {
        int rowNum = A.rows();
        int colNum = A.columns();
        DoubleMatrix2D C = A.like(rowNum, colNum);
        for (int i = 0; i < rowNum; i++) {
            for (int j = 0; j < rowNum; j++) {
                C.setQuick(i, j, A.getQuick(i, j) - B.getQuick(i, j));
            }
        }
        return C;
    }
}
