/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.util.stat;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import java.util.Arrays;

/**
 *
 * @author mxli
 */
public class SSJAFC {

    public static DoubleMatrix2D getDiagonal(DoubleMatrix2D dm) {
        double[][] dc = new double[dm.rows()][dm.columns()];
        for (int i = 0; i < dm.rows(); i++) {
            Arrays.fill(dc[i], 0);    //Should be tested! Are all the other elements set as 0?
            dc[i][i] = dm.getQuick(i, i);
        }
        DoubleMatrix2D result = new SparseDoubleMatrix2D(dc);
        return result;
    }

    public static DoubleMatrix1D getSqrt(DoubleMatrix1D dm) {
        DoubleMatrix1D result = new DenseDoubleMatrix1D(dm.size());
        for (int i = 0; i < dm.size(); i++) {
            result.setQuick(i, Math.sqrt(dm.getQuick(i)));
        }
        return result;
    }

    public static DoubleMatrix2D getSqrt(DoubleMatrix2D dm) {
        DoubleMatrix2D result = new DenseDoubleMatrix2D(dm.rows(), dm.columns());
        for (int i = 0; i < dm.rows(); i++) {
            for (int j = 0; j < dm.columns(); j++) {
                result.setQuick(i, j, Math.sqrt(dm.getQuick(i, j)));
            }
        }
        return result;
    }

    public static DoubleMatrix2D getInv(DoubleMatrix2D dm) {
        DoubleMatrix2D result = new DenseDoubleMatrix2D(dm.rows(), dm.columns());
        for (int i = 0; i < dm.rows(); i++) {
            for (int j = 0; j < dm.columns(); j++) {
                if (dm.getQuick(i, j) == 0) {
                    result.setQuick(i, j, 0); //This should be discussed!
                } else {
                    result.setQuick(i, j, 1 / Math.sqrt(dm.getQuick(i, j)));
                }
            }
        }
        return result;
    }

    public static DoubleMatrix2D subtract(DoubleMatrix2D dm1, DoubleMatrix2D dm2) {
        DoubleMatrix2D dm0 = new DenseDoubleMatrix2D(dm1.rows(), dm1.columns());
        for (int i = 0; i < dm0.rows(); i++) {
            for (int j = 0; j < dm0.columns(); j++) {
                dm0.setQuick(i, j, dm1.get(i, j) - dm2.getQuick(i, j));
            }
        }
        return dm0;
    }

    public static void subtract(DoubleMatrix2D dm1, DoubleMatrix2D dm2, DoubleMatrix2D dm3) {
        for (int i = 0; i < dm3.rows(); i++) {
            for (int j = 0; j < dm3.columns(); j++) {
                dm3.setQuick(i, j, dm1.get(i, j) - dm2.getQuick(i, j));
            }
        }
    }

    public static DoubleMatrix2D scalarMultiply(DoubleMatrix2D dm, double dblNum) {
        DoubleMatrix2D result = new DenseDoubleMatrix2D(dm.rows(), dm.columns());
        for (int i = 0; i < dm.rows(); i++) {
            for (int j = 0; j < dm.columns(); j++) {
                result.setQuick(i, j, dm.getQuick(i, j) * dblNum);
            }
        }
        return result;
    }

    public static DoubleMatrix2D signum(DoubleMatrix2D dm) {
        DoubleMatrix2D result = new DenseDoubleMatrix2D(dm.rows(), dm.columns());
        for (int i = 0; i < dm.rows(); i++) {
            for (int j = 0; j < dm.columns(); j++) {
                result.setQuick(i, j, Math.signum(dm.getQuick(i, j)));
            }
        }
        return result;
    }

    public static DoubleMatrix2D pow(DoubleMatrix2D dm, double dblPow) {
        DoubleMatrix2D result = new DenseDoubleMatrix2D(dm.rows(), dm.columns());
        for (int i = 0; i < dm.rows(); i++) {
            for (int j = 0; j < dm.columns(); j++) {
                result.setQuick(i, j, Math.pow(dm.getQuick(i, j), dblPow));
            }
        }
        return result;
    }

    public static DoubleMatrix1D sumColumn(DoubleMatrix2D dm) {
        DoubleMatrix1D dm1 = new DenseDoubleMatrix1D(dm.columns());
        for (int j = 0; j < dm1.size(); j++) {
            dm1.setQuick(j, dm.viewColumn(j).zSum());
        }
        return dm1;
    }

    public static double sumAllAbs(DoubleMatrix2D dm) {
        double sum = 0;
        for (int i = 0; i < dm.rows(); i++) {
            for (int j = 0; j < dm.columns(); j++) {
                sum += Math.abs(dm.getQuick(i, j));
            }
        }
        return sum;
    }

    public static DoubleMatrix1D divide(DoubleMatrix1D dm, double dblNum) {
        DoubleMatrix1D result = new DenseDoubleMatrix1D(dm.size());
        for (int i = 0; i < dm.size(); i++) {
            result.setQuick(i, dm.getQuick(i) / dblNum);
        }
        return result;
    }

    public static DoubleMatrix2D setColumnVector(DoubleMatrix2D dm, DoubleMatrix1D dm1, int intIndex) {
        DoubleMatrix2D result = dm;
        for (int i = 0; i < result.rows(); i++) {
            result.setQuick(i, intIndex, dm1.getQuick(i));
        }
        return result;
    }

    public static DoubleMatrix2D divideColumnSum(DoubleMatrix2D dm) {
        DoubleMatrix2D result = new DenseDoubleMatrix2D(dm.rows(), dm.columns());
        for (int j = 0; j < dm.columns(); j++) {
            DoubleMatrix1D dm1 = dm.viewColumn(j);
            dm1 = divide(dm1, dm1.zSum());  //dm1.zSum=1, why?
            result = setColumnVector(result, dm1, j);
        }
        return result;
    }

    public static DoubleMatrix2D ebeMultiply(DoubleMatrix2D dm1, DoubleMatrix2D dm2) {
        DoubleMatrix2D result = new DenseDoubleMatrix2D(dm1.rows(), dm1.columns());
        for (int i = 0; i < dm1.rows(); i++) {
            for (int j = 0; j < dm1.columns(); j++) {
                result.setQuick(i, j, dm1.getQuick(i, j) * dm2.getQuick(i, j));
            }
        }
        return result;
    }

    public static DoubleMatrix2D multiplyEntry(DoubleMatrix2D dm, int x, int y, double dblFactor) {
        DoubleMatrix2D result = new DenseDoubleMatrix2D(dm.rows(), dm.columns());
        result.setQuick(x, y, dm.getQuick(x, y) * dblFactor);
        return result;
    }
}
