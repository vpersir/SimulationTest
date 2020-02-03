/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.genetsimulator.entity;

import cern.colt.bitvector.BitVector;
import cern.colt.list.FloatArrayList;
import cern.colt.list.IntArrayList;
import cern.colt.map.OpenIntIntHashMap;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

/**
 *
 * @author mxli
 */
public class GentypeBasedFloatFastLDSparseMatrix extends LDSparseMatrix {

    final static long m1 = 0x5555555555555555L; //binary: 0101...  
    final static long m2 = 0x3333333333333333L; //binary: 00110011..  
    final static long m4 = 0x0f0f0f0f0f0f0f0fL; //binary:  4 zeros,  4 ones ...  
    final static long m8 = 0x00ff00ff00ff00ffL; //binary:  8 zeros,  8 ones ...  
    final static long m16 = 0x0000ffff0000ffffL; //binary: 16 zeros, 16 ones ...  
    final static long m32 = 0x00000000ffffffffL; //binary: 32 zeros, 32 ones ...  
    final static long hff = 0xffffffffffffffffL; //binary: all ones  
    final static long h01 = 0x0101010101010101L; //the sum of 256 to the power of 0,1,2,3...  
    // a specil matrix to store pair-wsie LD information
    // it is designed for the triganle matrix
    private static final long serialVersionUID = 100L;
    //it need high resoluation LD, otherwise acculaatively it intorduce a lot noise
    List<FloatArrayList> leadingSiteLDList = new ArrayList<FloatArrayList>();
    OpenIntIntHashMap lowerIndexMap = new OpenIntIntHashMap();
    int ldNum = 0;
    boolean hasSorted = false;
    //1/254
    public boolean isPhased = false;
    public long[][] bits1;
    public long[][] bits2;
    public long[][] bits3;
    public long[][] bits4;
    public int indivSize;
    double douIndivSize;
    public int varSize;
    public int unitNum;
    public double[] sum1;
    public double[] sum12;
    public boolean[] hasMissingGty;
    int adjIndivSize;

    public void setIsPhased(boolean isPhased) {
        this.isPhased = isPhased;
    }

    public GentypeBasedFloatFastLDSparseMatrix(List<AnnotSNP> snpList, List<Individual> indListForLD, boolean isPhased) {
        this.isPhased = isPhased;
        varSize = snpList.size();
        indivSize = indListForLD.size();
        if (isPhased) {
            adjIndivSize = indivSize * 2;
            bits1 = new long[varSize][];
            bits2 = new long[varSize][];
            bits3 = new long[varSize][];
        } else {
            adjIndivSize = indivSize;
            bits1 = new long[varSize][];
            bits2 = new long[varSize][];
            bits3 = new long[varSize][];
            bits4 = new long[varSize][];
        }

        hasMissingGty = new boolean[varSize];
        BitVector temp1BV = new BitVector(indivSize);
        BitVector temp2BV = new BitVector(indivSize);
        BitVector temp3BV = new BitVector(indivSize);
        BitVector temp4BV = new BitVector(indivSize);

        sum1 = new double[varSize];
        Arrays.fill(sum1, 0);
        sum12 = new double[varSize];
        Arrays.fill(sum12, 0);
        Arrays.fill(hasMissingGty, false);
        unitNum = 0;

        boolean tmpBool1, tmpBool2;
        int missingGtyNum = 0;
        int nonNumIndiv = 0;
        //note the oder of SNPs cannot be changed any more
        if (isPhased) {
            for (int j = 0; j < varSize; j++) {
                AnnotSNP snp = snpList.get(j);
                if (snp.order < 0) {
                    //System.out.println(snp1.rsID+"  "+snp1.physicalPosition);
                    continue;
                }
                temp1BV.replaceFromToWith(0, indivSize - 1, false);
                temp2BV.replaceFromToWith(0, indivSize - 1, false);
                temp3BV.replaceFromToWith(0, indivSize - 1, false);

                nonNumIndiv = 0;
                for (int k = 0; k < indivSize; k++) {
                    StatusGtySet gty = indListForLD.get(k).markerGtySet;
                    if (gty.existence.getQuick(snp.order)) {
                        tmpBool1 = gty.paternalChrom.get(snp.order);
                        temp1BV.putQuick(k, tmpBool1);
                        if (tmpBool1) {
                            sum1[j] += 1;
                        }
                        tmpBool1 = gty.maternalChrom.get(snp.order);
                        temp2BV.putQuick(k, tmpBool1);
                        if (tmpBool1) {
                            sum1[j] += 1;
                        }
                        nonNumIndiv += 2;
                    } else {
                        temp3BV.putQuick(k, true);
                        hasMissingGty[j] = true;
                    }
                }
                sum1[j] /= nonNumIndiv;
                sum12[j] = (1 - sum1[j]) * sum1[j];
                long[] tempLong = temp1BV.elements();
                bits1[j] = new long[tempLong.length];
                System.arraycopy(tempLong, 0, bits1[j], 0, tempLong.length);
                tempLong = temp2BV.elements();
                bits2[j] = new long[tempLong.length];
                System.arraycopy(tempLong, 0, bits2[j], 0, tempLong.length);
                tempLong = temp3BV.elements();
                bits3[j] = new long[tempLong.length];
                System.arraycopy(tempLong, 0, bits3[j], 0, tempLong.length);
                if (unitNum == 0) {
                    unitNum = tempLong.length;
                }
                //change the ID finally
                snp.order = j;
            }
        } else {
            for (int i = 0; i < varSize; i++) {
                AnnotSNP snp = snpList.get(i);
                if (snp.order < 0) {
                    //System.out.println(snp1.rsID+"  "+snp1.physicalPosition);
                    continue;
                }
                nonNumIndiv = 0;
                temp1BV.replaceFromToWith(0, indivSize - 1, false);
                temp2BV.replaceFromToWith(0, indivSize - 1, false);
                temp3BV.replaceFromToWith(0, indivSize - 1, false);
                temp4BV.replaceFromToWith(0, indivSize - 1, true);
                for (int k = 0; k < indivSize; k++) {
                    StatusGtySet gty = indListForLD.get(k).markerGtySet;
                    if (gty.existence.getQuick(snp.order)) {
                        tmpBool1 = gty.paternalChrom.get(snp.order);
                        tmpBool2 = gty.maternalChrom.get(snp.order);
                        if (tmpBool1 && tmpBool2) {
                            temp1BV.putQuick(k, true);
                            temp2BV.putQuick(k, true);
                            temp3BV.putQuick(k, false);
                            sum1[i] += 2;
                            sum12[i] += 4;
                        } else if (tmpBool1 || tmpBool2) {
                            temp1BV.putQuick(k, false);
                            temp2BV.putQuick(k, true);
                            temp3BV.putQuick(k, true);
                            sum1[i] += 1;
                            sum12[i] += 1;
                        } else {
                            temp1BV.putQuick(k, false);
                            temp2BV.putQuick(k, false);
                            temp3BV.putQuick(k, false);
                        }
                        nonNumIndiv += 1;
                    } else {
                        temp4BV.putQuick(k, false);
                        if (!hasMissingGty[i]) {
                            missingGtyNum++;
                            hasMissingGty[i] = true;
                        }

                        ///break;
                    }
                }
                sum1[i] /= Math.sqrt(nonNumIndiv);
                //sum12[i] = Math.sqrt(sum12[i] - sum1[i] * sum1[i]);
                long[] tempLong = temp1BV.elements();
                bits1[i] = new long[tempLong.length];
                System.arraycopy(tempLong, 0, bits1[i], 0, tempLong.length);
                tempLong = temp2BV.elements();
                bits2[i] = new long[tempLong.length];
                System.arraycopy(tempLong, 0, bits2[i], 0, tempLong.length);
                tempLong = temp3BV.elements();
                bits3[i] = new long[tempLong.length];
                System.arraycopy(tempLong, 0, bits3[i], 0, tempLong.length);
                tempLong = temp4BV.elements();
                bits4[i] = new long[tempLong.length];
                System.arraycopy(tempLong, 0, bits4[i], 0, tempLong.length);
                //change the ID finally
                if (unitNum == 0) {
                    unitNum = tempLong.length;
                }
                snp.order = i;
            }
        }

        /*
         //code for  testing
         int index1 = 11, index2 = 12;
         DoubleArrayList a1 = new DoubleArrayList();
         DoubleArrayList a2 = new DoubleArrayList();
         for (int k = 0; k < indivSize; k++) {
         StatusGtySet gty = indListForLD.get(k).markerGtySet;
         if (gty.existence.getQuick(index1) && gty.existence.getQuick(index2)) {
         tmpBool1 = gty.paternalChrom.get(index1);
         tmpBool2 = gty.maternalChrom.get(index1);
         if (tmpBool1 && tmpBool2) {
         a1.add(2);
         } else if (tmpBool1 || tmpBool2) {
         a1.add(1);
         } else {
         a1.add(0);
         }
        
         tmpBool1 = gty.paternalChrom.get(index2);
         tmpBool2 = gty.maternalChrom.get(index2);
         if (tmpBool1 && tmpBool2) {
         a2.add(2);
         } else if (tmpBool1 || tmpBool2) {
         a2.add(1);
         } else {
         a2.add(0);
         }
         }
         }
         double mean1 = Descriptive.mean(a1);
         double mean2 = Descriptive.mean(a2);
         double sd1 = Descriptive.sampleVariance(a1, mean1);
         double sd2 = Descriptive.sampleVariance(a2, mean2);
         double r = Descriptive.correlation(a1, Math.sqrt(sd1), a2, Math.sqrt(sd2));
         r=r*r;
         double r1 = 0;
         try {
         r1 = this.getLDAt(index1, index2);
         } catch (Exception ex) {
         }
         System.out.println(r + "\t" + r1);
         */
    }

    public GentypeBasedFloatFastLDSparseMatrix(long[][] bits1, long[][] bits2, long[][] bits3, double[] sum1, double[] sum12, int indivSize) {
        this.bits1 = bits1;
        this.bits2 = bits2;
        this.bits3 = bits3;
        this.sum1 = sum1;
        this.sum12 = sum12;
        unitNum = bits1[0].length;
        varSize = sum1.length;
        this.indivSize = indivSize;
        this.douIndivSize = indivSize * 2;
    }

    public GentypeBasedFloatFastLDSparseMatrix(GentypeBasedFloatFastLDSparseMatrix ldSpareMatrix) {
        this.bits1 = ldSpareMatrix.bits1;
        this.bits2 = ldSpareMatrix.bits2;
        this.bits3 = ldSpareMatrix.bits3;
        this.sum1 = ldSpareMatrix.sum1;
        this.sum12 = ldSpareMatrix.sum12;
        this.unitNum = ldSpareMatrix.unitNum;
        this.varSize = ldSpareMatrix.varSize;
        this.indivSize = ldSpareMatrix.indivSize;
        this.douIndivSize = ldSpareMatrix.douIndivSize;
        this.isPhased = ldSpareMatrix.isPhased;

    }

    @Override
    public boolean isEmpty() {
        return leadingSiteLDList.isEmpty();
    }

    public int size() {
        return ldNum;
    }

    @Override
    public Set<Integer> getAllUniqueIndexes() {
        return null;
    }

    public boolean addLDAt(int index1, int index2, double ld) throws Exception {

        int acturalPos = -1;
        //ensure curIndex1 is less than index2 always
        if (index1 < index2) {
            if (lowerIndexMap.containsKey(index1)) {
                acturalPos = lowerIndexMap.get(index1);
            } else {
                acturalPos = leadingSiteLDList.size();
                lowerIndexMap.put(index1, acturalPos);
                leadingSiteLDList.add(new FloatArrayList());
            }
            int offSet = index2 - index1;
            FloatArrayList subLDList = leadingSiteLDList.get(acturalPos);
            if (subLDList == null) {
                subLDList = new FloatArrayList();
            }
            for (int i = subLDList.size(); i < offSet; i++) {
                subLDList.add(Float.NaN);
            }
            offSet--;
            subLDList.setQuick(offSet, (float) ld);
            ldNum++;
        } else if (index1 > index2) {
            if (lowerIndexMap.containsKey(index2)) {
                acturalPos = lowerIndexMap.get(index2);
            } else {
                acturalPos = leadingSiteLDList.size();
                lowerIndexMap.put(index2, acturalPos);
                leadingSiteLDList.add(new FloatArrayList());
            }

            int offSet = index1 - index2;
            FloatArrayList subLDList = leadingSiteLDList.get(acturalPos);
            if (subLDList == null) {
                subLDList = new FloatArrayList();
            }

            for (int i = subLDList.size(); i < offSet; i++) {
                subLDList.add(Float.NaN);
            }
            offSet--;
            subLDList.setQuick(offSet, (float) ld);
            ldNum++;
        } else {
            //  throw new Exception("No need to store this LD!");
            System.err.println("No need to store this LD for " + index1 + " " + index2 + " " + ld);
            return false;
        }

        return true;
    }

    @Override
    public boolean calculateGenotypeCorrelationBlock(IntArrayList indexes) throws Exception {

        return true;
    }

    @Override
    public DoubleMatrix2D subDenseLDMatrix(IntArrayList poss) throws Exception {
        int dim = poss.size();
        poss.quickSort();
        DoubleMatrix2D corrMat = new DenseDoubleMatrix2D(dim, dim);
        double x = 0;
        for (int i = 0; i < dim; i++) {
            corrMat.setQuick(i, i, 1);
            for (int j = i + 1; j < dim; j++) {
                x = getLDAt(poss.getQuick(i), poss.getQuick(j));
                corrMat.setQuick(i, j, x);
                corrMat.setQuick(j, i, x);
            }
        }
        // System.out.println(corrMat.toString());
        return corrMat;
    }

    @Override
    public void releaseLDData() {
        if (lowerIndexMap.size() > 1000) {
            leadingSiteLDList.clear();
            lowerIndexMap.clear();
        }

        System.gc();
    }

    public double calculateGenotypeCorrelation(int index1, int index2) throws Exception {
        if (index1 == index2) {
            return 1;
        }
        double r = 0;
        int freqAB = 0;
        long x;
        int tempInt;

        for (int k = 0; k < unitNum; k++) {
            x = bits1[index1][k] & bits1[index2][k];
            freqAB += (Long.bitCount(x) << 1);  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
            x = bits2[index1][k] & bits2[index2][k];
            freqAB += (Long.bitCount(x) << 1);
            x = bits3[index1][k] & bits3[index2][k];
            freqAB -= (Long.bitCount(x));
        }
        if (index1 == 11 && index2 == 12) {
            int sss = 0;
        }

        if (hasMissingGty[index1] || hasMissingGty[index2]) {
            int sum1i = 0;
            int sum1j = 0;
            int sum12i = 0;
            int sum12j = 0;
            int nij = 0;
            long exsit;
            for (int k = 0; k < unitNum; k++) {
                exsit = bits4[index1][k] & bits4[index2][k];
                x = exsit;
                nij += (Long.bitCount(x));

                x = bits1[index1][k] & bits1[index1][k] & exsit;
                tempInt = Long.bitCount(x);
                sum1i += (tempInt);
                sum12i += (tempInt << 1);

                x = bits2[index1][k] & bits2[index1][k] & exsit;
                tempInt = Long.bitCount(x);
                sum1i += (tempInt);
                sum12i += (tempInt << 1);

                x = bits3[index1][k] & bits3[index1][k] & exsit;
                sum12i -= (Long.bitCount(x));

                x = bits1[index2][k] & bits1[index2][k] & exsit;
                tempInt = Long.bitCount(x);
                sum1j += (tempInt);
                sum12j += (tempInt << 1);
                x = bits2[index2][k] & bits2[index2][k] & exsit;
                tempInt = Long.bitCount(x);
                sum1j += (tempInt);
                sum12j += (tempInt << 1);

                x = bits3[index2][k] & bits3[index2][k] & exsit;
                sum12j -= (Long.bitCount(x));
            }

            //Strange! this formula is even faster
            r = (freqAB - sum1i * sum1j / ((double) nij)) / Math.sqrt((sum12i - sum1i * sum1i / ((double) nij)) * (sum12j - sum1j * sum1j / ((double) nij)));

            //r = (freqAB - sum1[index1] * sum1[index2]) / (sum12[index1] * sum12[index2]);
        } else {
            //Strange! this formula is even faster
            r = (freqAB - sum1[index1] * sum1[index2]) / Math.sqrt((sum12[index1] - sum1[index1] * sum1[index1]) * (sum12[index2] - sum1[index2] * sum1[index2]));
            //r = (freqAB - sum1[index1] * sum1[index2]) / (sum12[index1] * sum12[index2]);
        }

      //  r = r * r;
        if (r > 1) {
            int sss = 0;
        }
        //correction
        // r = (douSize * r - 1) / (douSize - 3);
        //johny's correction                
       // r = 1 - (adjIndivSize - 3.0) / (adjIndivSize - 2.0) * (1.0 - r) * (1.0 + 2.0 * (1.0 - r) / (adjIndivSize - 3.3));
        return (float) (r);
    }
    //assume the genotypes are phased 

    public double calculateLDR(int index1, int index2) throws Exception {
        if (index1 == index2) {
            return 1;
        }
        //note the chromsome number is alwalys eaqual to  indivSize * 2 at every locus due to the missing genotypes      
        int douSize = indivSize * 2;
        long x;
        double r;
        int freqAB = 0;
        int missingNum = 0;
        //refecen http://blog.csdn.net/hitwhylz/article/details/10122617
        for (int k = 0; k < unitNum; k++) {
            x = bits1[index1][k] & bits1[index2][k];
            x -= (x >>> 1) & m1;             //put count of each 2 bits into those 2 bits  
            x = (x & m2) + ((x >>> 2) & m2); //put count of each 4 bits into those 4 bits   
            x = (x + (x >>> 4)) & m4;        //put count of each 8 bits into those 8 bits   
            freqAB += ((x * h01) >>> 56);  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
            x = bits2[index1][k] & bits2[index2][k];
            x -= (x >>> 1) & m1;             //put count of each 2 bits into those 2 bits  
            x = (x & m2) + ((x >>> 2) & m2); //put count of each 4 bits into those 4 bits   
            x = (x + (x >>> 4)) & m4;        //put count of each 8 bits into those 8 bits   
            freqAB += ((x * h01) >>> 56);
        }

        r = ((double) freqAB / douSize - sum1[index1] * sum1[index2]);
        r = r / Math.sqrt(sum12[index1] * sum12[index2]);
        //r = r * r / (sum12[index1] * sum12[index2]);

        if (Double.isNaN(r)) {
            return 0;
        }

        //correction
        // r = (douSize * r - 1) / (douSize - 3);
        //johny's correction                
        //r = 1 - (adjIndivSize - 3.0) / (adjIndivSize - 2.0) * (1.0 - r) * (1.0 + 2 * (1 - r) / (adjIndivSize - 3.3));
        return (double) (r);
    }

    @Override
    public double getLDAt(int index1, int index2) throws Exception {

        int acturalPos = -1;
        double r2 = 0;
        if (index1 < index2) {
            if (lowerIndexMap.containsKey(index1)) {
                acturalPos = lowerIndexMap.get(index1);
            } else {
                if (isPhased) {
                    r2 = calculateLDR(index1, index2);
                } else {
                    r2 = calculateGenotypeCorrelation(index1, index2);
                }
                addLDAt(index1, index2, r2);
                return r2;
            }

            FloatArrayList subLDList = leadingSiteLDList.get(acturalPos);
            if (subLDList == null) {
                if (isPhased) {
                    r2 = calculateLDR(index1, index2);
                } else {
                    r2 = calculateGenotypeCorrelation(index1, index2);
                }
                addLDAt(index1, index2, r2);
                return r2;
            }

            int offSet = index2 - index1;
            if (offSet > subLDList.size()) {
                if (isPhased) {
                    r2 = calculateLDR(index1, index2);
                } else {
                    r2 = calculateGenotypeCorrelation(index1, index2);
                }
                addLDAt(index1, index2, r2);
                return r2;
            }
            offSet--;
            double val = subLDList.get(offSet);
            if (Double.isNaN(val)) {
                if (isPhased) {
                    r2 = calculateLDR(index1, index2);
                } else {
                    r2 = calculateGenotypeCorrelation(index1, index2);
                }
                addLDAt(index1, index2, r2);
                return r2;
            } else if (val == -127) {
                return 0;
            } else if (val == 127) {
                return 1;
            }
            return val;
        } else if (index1 > index2) {
            if (lowerIndexMap.containsKey(index2)) {
                acturalPos = lowerIndexMap.get(index2);
            } else {
                if (isPhased) {
                    r2 = calculateLDR(index1, index2);
                } else {
                    r2 = calculateGenotypeCorrelation(index1, index2);
                }
                addLDAt(index1, index2, r2);
                return r2;
            }

            FloatArrayList subLDList = leadingSiteLDList.get(acturalPos);
            if (subLDList == null) {
                if (isPhased) {
                    r2 = calculateLDR(index1, index2);
                } else {
                    r2 = calculateGenotypeCorrelation(index1, index2);
                }
                addLDAt(index1, index2, r2);
                return r2;
            }

            int offSet = index1 - index2;
            if (offSet > subLDList.size()) {
                if (isPhased) {
                    r2 = calculateLDR(index1, index2);
                } else {
                    r2 = calculateGenotypeCorrelation(index1, index2);
                }
                addLDAt(index1, index2, r2);
                return r2;
            }
            offSet--;
            double val = subLDList.get(offSet);
            if (Double.isNaN(val)) {
                if (isPhased) {
                    r2 = calculateLDR(index1, index2);
                } else {
                    r2 = calculateGenotypeCorrelation(index1, index2);
                }
                addLDAt(index1, index2, r2);
                return r2;
            }
            return val;
        }
        return 1;
    }
}
