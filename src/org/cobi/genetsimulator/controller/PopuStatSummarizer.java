/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.genetsimulator.controller;

import cern.colt.list.DoubleArrayList;
import cern.colt.list.IntArrayList;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.jet.stat.Descriptive;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.RealMatrix;
import org.cobi.genetsimulator.controller.MultiLocusLiabilityThresholdModel.LiabilityDiseaseSNP;
import org.cobi.genetsimulator.controller.MultilocusLogitModel.LogitDiseaseSNP;
import org.cobi.genetsimulator.entity.AnnotSNP;
import org.cobi.genetsimulator.entity.DiseaseSNP;
import org.cobi.genetsimulator.entity.Individual;
import org.cobi.genetsimulator.entity.StatusGtySet;
import org.cobi.util.text.LocalNumber;
 
import org.cobi.util.text.LocalString;
 

/**
 *
 * @author mxli
 */
public class PopuStatSummarizer {

    public RealMatrix calculate2DProbability(List<Individual> allIndiv, int[] indexes) throws Exception {
        int len = indexes.length;
        int indivSize = allIndiv.size();
        int chromSize = indivSize * 2;
        RealMatrix twoDProb = new Array2DRowRealMatrix(len, len);
        for (int i = 0; i < len; i++) {
            for (int j = i; j < len; j++) {
                twoDProb.setEntry(i, j, 0.0);
            }
        }
        for (int k = 0; k < indivSize; k++) {
            StatusGtySet gty = allIndiv.get(k).markerGtySet;
            for (int i = 0; i < len; i++) {
                for (int j = i; j < len; j++) {
                    if (gty.existence.getQuick(indexes[i]) && gty.existence.getQuick(indexes[j])) {
                        if (gty.paternalChrom.getQuick(indexes[i]) && gty.paternalChrom.getQuick(indexes[j])) {
                            twoDProb.setEntry(i, j, twoDProb.getEntry(i, j) + 1.0);
                        }
                        if (gty.maternalChrom.getQuick(indexes[i]) && gty.maternalChrom.getQuick(indexes[j])) {
                            twoDProb.setEntry(i, j, twoDProb.getEntry(i, j) + 1.0);
                        }
                    }
                }
            }
        }

        for (int i = 0; i < len; i++) {
            for (int j = i; j < len; j++) {
                twoDProb.setEntry(i, j, twoDProb.getEntry(i, j) / chromSize);
            }
        }
        LocalString.print2DRealMatrix(twoDProb);
        return twoDProb;
    }

    public DoubleMatrix2D convert2LDIndex(DoubleMatrix2D twoDProb) throws Exception {
        int i, j;
        int rowNum = twoDProb.rows(), colNum = twoDProb.columns();
        DoubleMatrix2D dMatrix = twoDProb.copy();
        double[] alleleFreq = new double[rowNum];
        System.out.println("Allele Frequencies:");
        for (i = 0; i < rowNum; i++) {
            alleleFreq[i] = dMatrix.getQuick(i, i);
            System.out.println(alleleFreq[i]);
        }

        System.out.println("D Matrix:");
        for (i = 0; i < rowNum; i++) {
            for (j = i; j < colNum; j++) {
                dMatrix.setQuick(i, j, dMatrix.getQuick(i, j) - alleleFreq[i] * alleleFreq[j]);
            }
        }
        // LocalString.print2DRealMatrix(dMatrix);

        System.out.println("R Matrix:");
        DoubleMatrix2D rMatrix = dMatrix.copy();
        DoubleMatrix2D covP = new DenseDoubleMatrix2D(rowNum, rowNum);

        for (i = 0; i < rowNum; i++) {
            for (j = i; j < colNum; j++) {
                rMatrix.setQuick(i, j, dMatrix.getQuick(i, j) / Math.sqrt(alleleFreq[i] * alleleFreq[j] * (1 - alleleFreq[i]) * (1 - alleleFreq[j])));
                covP.setQuick(i, j, dMatrix.getQuick(i, j) / Math.sqrt(alleleFreq[i] * alleleFreq[j] * (1 - alleleFreq[i]) * (1 - alleleFreq[j])));
            }
        }
        // LocalString.print2DRealMatrix(rMatrix);

        System.out.println("D' Matrix:");
        for (i = 0; i < rowNum; i++) {
            for (j = i; j < colNum; j++) {
                if (dMatrix.getQuick(i, j) >= 0) {
                    dMatrix.setQuick(i, j, dMatrix.getQuick(i, j) / Math.max(alleleFreq[i] * (1 - alleleFreq[j]), alleleFreq[j] * (1 - alleleFreq[i])));
                } else {
                    dMatrix.setQuick(i, j, dMatrix.getQuick(i, j) / Math.min(-alleleFreq[i] * alleleFreq[j], -(1 - alleleFreq[i]) * (1 - alleleFreq[j])));
                }
            }
        }
        // LocalString.print2DRealMatrix(dMatrix);
        return covP;
    }

    public DoubleMatrix2D convertLDR2JointProbability(final double[][] ldR) throws Exception {
        int rowNum = ldR.length, colNum = ldR[0].length;
        int i, j;
        double tmp;
        DoubleMatrix2D probMatrix = new DenseDoubleMatrix2D(ldR);

       // System.out.println("Joint probabilies:");
        double[] alleleFreq = new double[rowNum];
        for (i = 0; i < rowNum; i++) {
            alleleFreq[i] = ldR[i][i];
        }
        for (i = 0; i < rowNum; i++) {
            for (j = i + 1; j < colNum; j++) {
                tmp = ldR[i][j] * Math.sqrt(alleleFreq[i] * alleleFreq[j] * (1 - alleleFreq[i]) * (1 - alleleFreq[j]));
                tmp = tmp + alleleFreq[i] * alleleFreq[j];
                probMatrix.setQuick(i, j, tmp);
            }
        }
        //System.out.println(probMatrix.toString());
        return probMatrix;
    }

    public DoubleMatrix2D convertLDR2JointProbability(final DoubleMatrix2D ldR, double[] alleleFreq) throws Exception {
        int rowNum = alleleFreq.length, colNum = alleleFreq.length;
        int i, j;
        double tmp;
        DoubleMatrix2D probMatrix = new DenseDoubleMatrix2D(rowNum, colNum);

       // System.out.println("Joint probabilies:");

        for (i = 0; i < rowNum; i++) {
            probMatrix.setQuick(i, i, alleleFreq[i]);
            for (j = i + 1; j < colNum; j++) {
                tmp = ldR.getQuick(i, j) * Math.sqrt(alleleFreq[i] * alleleFreq[j] * (1 - alleleFreq[i]) * (1 - alleleFreq[j]));
                tmp = tmp + alleleFreq[i] * alleleFreq[j];
                if (tmp < 0) {
                    tmp = 0;
                }
                probMatrix.setQuick(i, j, tmp);
            }
        }
        //LocalString.print2DRealMatrix(probMatrix);
        return probMatrix;
    }

    public double[] calculateJointSusceptibGtyProbability(List<Individual> allIndiv, List<LogitDiseaseSNP> susbSnpList) throws Exception {
        int len = susbSnpList.size();
        int[] indexes = new int[len];
        int indivSize = allIndiv.size();
        for (int i = 0; i < len; i++) {
            indexes[i] = susbSnpList.get(i).getLDMarkerPosition();
        }
        int times = (int) Math.pow(3, len);
        double[] freq = new double[times];
        Arrays.fill(freq, 0.0);
        StringBuffer gtyLabel = new StringBuffer();
        gtyLabel.setLength(len);
        boolean noMiss;
        int effectIndivSize = 0;
        boolean riskAllele;
        for (int k = 0; k < indivSize; k++) {
            StatusGtySet gty = allIndiv.get(k).markerGtySet;
            noMiss = true;
            for (int i = 0; i < len; i++) {
                riskAllele = susbSnpList.get(i).riskAlleleLable;
                if (gty.existence.getQuick(indexes[i])) {
                    //ensure 2 denote 2 risk alleles and 1 1 risk allele and 0 no risk allele
                    if (gty.paternalChrom.getQuick(indexes[i]) == riskAllele && gty.maternalChrom.getQuick(indexes[i]) == riskAllele) {
                        gtyLabel.setCharAt(i, '2');
                    } else if (gty.paternalChrom.getQuick(indexes[i]) != riskAllele && gty.maternalChrom.getQuick(indexes[i]) != riskAllele) {
                        gtyLabel.setCharAt(i, '0');
                    } else {
                        gtyLabel.setCharAt(i, '1');
                    }
                } else {
                    noMiss = false;
                    break;
                }
            }
            if (noMiss) {
                freq[LocalNumber.xSystem2decimal(gtyLabel, 3)] += 1;
                effectIndivSize++;
            }
        }
        for (int i = 0; i < times; i++) {
            freq[i] /= effectIndivSize;
            System.out.print(LocalNumber.decimal2XSystem(i, 3, 2));
            System.out.print("\t");
            System.out.println(freq[i]);
        }
        return freq;
    }

    public DoubleMatrix2D calculate2DProbability(List<Individual> allIndiv, int snpNum) throws Exception {
        int indivSize = allIndiv.size();
        int chromSize = indivSize * 2;

        DoubleMatrix2D twoDProb = new DenseDoubleMatrix2D(snpNum, snpNum);
        for (int i = 0; i < snpNum; i++) {
            for (int j = i; j < snpNum; j++) {
                twoDProb.setQuick(i, j, 0.0);
            }
        }

        for (int k = 0; k < indivSize; k++) {
            StatusGtySet gty = allIndiv.get(k).markerGtySet;
            for (int i = 0; i < snpNum; i++) {
                for (int j = i; j < snpNum; j++) {
                    if (gty.existence.getQuick(i) && gty.existence.getQuick(j)) {
                        if (gty.paternalChrom.getQuick(i) && gty.paternalChrom.getQuick(j)) {
                            twoDProb.setQuick(i, j, twoDProb.getQuick(i, j) + 1.0);
                        }
                        if (gty.maternalChrom.getQuick(i) && gty.maternalChrom.getQuick(j)) {
                            twoDProb.setQuick(i, j, twoDProb.getQuick(i, j) + 1.0);
                        }
                    }
                }
            }
        }

        for (int i = 0; i < snpNum; i++) {
            for (int j = i; j < snpNum; j++) {
                twoDProb.setQuick(i, j, twoDProb.getQuick(i, j) / chromSize);
            }
        }
        //LocalString.print2DRealMatrix(twoDProb);
        return twoDProb;
    }

    public DoubleMatrix2D calculateHillRobertsonLD(List<Individual> allIndiv, int snpNum, int phenotypeCode) throws Exception {
        int indivSize = allIndiv.size();

        DoubleMatrix2D corMat = new DenseDoubleMatrix2D(snpNum, snpNum);
        DoubleMatrix2D jointProb00 = new DenseDoubleMatrix2D(snpNum, snpNum);
        DoubleMatrix2D jointProb11 = new DenseDoubleMatrix2D(snpNum, snpNum);
        DoubleMatrix2D jointProb01 = new DenseDoubleMatrix2D(snpNum, snpNum);
        DoubleMatrix2D jointProb10 = new DenseDoubleMatrix2D(snpNum, snpNum);
        DoubleMatrix1D margProb0 = new DenseDoubleMatrix1D(snpNum);

        for (int i = 0; i < snpNum; i++) {
            margProb0.setQuick(i, 0);
            for (int j = i; j < snpNum; j++) {
                corMat.setQuick(i, j, 0.0);
            }
        }
        int chromSize = 0;
        for (int k = 0; k < indivSize; k++) {
            StatusGtySet gty = allIndiv.get(k).markerGtySet;
            if (allIndiv.get(k).getAffectedStatus() != phenotypeCode) {
                continue;
            }
            chromSize++;
        }
        chromSize = chromSize * 2;
        //note the chromsome number is alwalys eaqual to  indivSize * 2 at every locus due to the missing genotypes
        for (int k = 0; k < indivSize; k++) {
            StatusGtySet gty = allIndiv.get(k).markerGtySet;
            if (allIndiv.get(k).getAffectedStatus() != phenotypeCode) {
                continue;
            }
            for (int i = 0; i < snpNum; i++) {
                if (gty.existence.getQuick(i)) {
                    if (!gty.paternalChrom.getQuick(i)) {
                        margProb0.setQuick(i, margProb0.getQuick(i) + 1);
                    }

                    if (!gty.maternalChrom.getQuick(i)) {
                        margProb0.setQuick(i, margProb0.getQuick(i) + 1);
                    }
                }

                for (int j = i; j < snpNum; j++) {
                    if (gty.existence.getQuick(i) && gty.existence.getQuick(j)) {
                        if (!gty.paternalChrom.getQuick(i) && !gty.paternalChrom.getQuick(j)) {
                            jointProb00.setQuick(i, j, jointProb00.getQuick(i, j) + 1.0);
                        } else if (!gty.paternalChrom.getQuick(i) && gty.paternalChrom.getQuick(j)) {
                            jointProb01.setQuick(i, j, jointProb01.getQuick(i, j) + 1.0);
                        } else if (gty.paternalChrom.getQuick(i) && !gty.paternalChrom.getQuick(j)) {
                            jointProb10.setQuick(i, j, jointProb10.getQuick(i, j) + 1.0);
                        } else {
                            jointProb11.setQuick(i, j, jointProb11.getQuick(i, j) + 1.0);
                        }

                        if (!gty.maternalChrom.getQuick(i) && !gty.maternalChrom.getQuick(j)) {
                            jointProb00.setQuick(i, j, jointProb00.getQuick(i, j) + 1.0);
                        } else if (!gty.maternalChrom.getQuick(i) && gty.maternalChrom.getQuick(j)) {
                            jointProb01.setQuick(i, j, jointProb01.getQuick(i, j) + 1.0);
                        } else if (gty.maternalChrom.getQuick(i) && !gty.maternalChrom.getQuick(j)) {
                            jointProb10.setQuick(i, j, jointProb10.getQuick(i, j) + 1.0);
                        } else {
                            jointProb11.setQuick(i, j, jointProb11.getQuick(i, j) + 1.0);
                        }
                    }
                }
            }
        }

        for (int i = 0; i < snpNum; i++) {
            for (int j = i; j < snpNum; j++) {
                jointProb00.setQuick(i, j, jointProb00.getQuick(i, j) / chromSize);
                jointProb01.setQuick(i, j, jointProb01.getQuick(i, j) / chromSize);
                jointProb10.setQuick(i, j, jointProb10.getQuick(i, j) / chromSize);
                jointProb11.setQuick(i, j, jointProb11.getQuick(i, j) / chromSize);
            }
            margProb0.setQuick(i, margProb0.getQuick(i) / chromSize);
        }
        double tmpSq = indivSize * 2;
        for (int i = 0; i < snpNum; i++) {
            for (int j = i; j < snpNum; j++) {
                // double maxRS = Math.min(margProb0.get(i) * (1 - margProb0.get(j)), margProb0.get(j) * (1 - margProb0.get(i)));
                // double minRS = Math.max(-margProb0.get(i) * margProb0.get(j), -(1 - margProb0.get(i)) * (1 - margProb0.get(j)));
                double sq = (jointProb00.getQuick(i, j) * jointProb11.getQuick(i, j) - jointProb01.getQuick(i, j) * jointProb10.getQuick(i, j));
                // sq = sq / Math.sqrt(margProb0.getQuick(i) * (1 - margProb0.getQuick(i)) * margProb0.getQuick(j) * (1 - margProb0.getQuick(j)));
                sq = sq * sq / (margProb0.getQuick(i) * (1 - margProb0.getQuick(i)) * margProb0.getQuick(j) * (1 - margProb0.getQuick(j)));
                //correction
                //sq = (tmpSq * sq - 1) / (tmpSq - 3);
                //johny's correction                
                sq = 1 - (tmpSq - 3) / (tmpSq - 2) * (1 - sq) * (1 + 2 * (1 - sq) / (tmpSq - 3.3));
                corMat.setQuick(i, j, sq);
                corMat.setQuick(j, i, sq);
            }
        }
        // System.out.println(corMat);

        return corMat;
    }

    public DoubleMatrix2D calculateHillRobertsonLD(List<Individual> allIndiv, int snpNum) throws Exception {
        int indivSize = allIndiv.size();
        int chromSize = indivSize * 2;
        DoubleMatrix2D corMat = new DenseDoubleMatrix2D(snpNum, snpNum);
        DoubleMatrix2D jointProb00 = new DenseDoubleMatrix2D(snpNum, snpNum);
        DoubleMatrix2D jointProb11 = new DenseDoubleMatrix2D(snpNum, snpNum);
        DoubleMatrix2D jointProb01 = new DenseDoubleMatrix2D(snpNum, snpNum);
        DoubleMatrix2D jointProb10 = new DenseDoubleMatrix2D(snpNum, snpNum);
        DoubleMatrix1D margProb0 = new DenseDoubleMatrix1D(snpNum);

        for (int i = 0; i < snpNum; i++) {
            margProb0.setQuick(i, 0);
            for (int j = i; j < snpNum; j++) {
                corMat.setQuick(i, j, 0.0);
            }
        }

        //note the chromsome number is alwalys eaqual to  indivSize * 2 at every locus due to the missing genotypes
        for (int k = 0; k < indivSize; k++) {
            StatusGtySet gty = allIndiv.get(k).markerGtySet;
            for (int i = 0; i < snpNum; i++) {
                if (gty.existence.getQuick(i)) {
                    if (!gty.paternalChrom.getQuick(i)) {
                        margProb0.setQuick(i, margProb0.getQuick(i) + 1);
                    }

                    if (!gty.maternalChrom.getQuick(i)) {
                        margProb0.setQuick(i, margProb0.getQuick(i) + 1);
                    }
                }

                for (int j = i; j < snpNum; j++) {
                    if (gty.existence.getQuick(i) && gty.existence.getQuick(j)) {
                        if (!gty.paternalChrom.getQuick(i) && !gty.paternalChrom.getQuick(j)) {
                            jointProb00.setQuick(i, j, jointProb00.getQuick(i, j) + 1.0);
                        } else if (!gty.paternalChrom.getQuick(i) && gty.paternalChrom.getQuick(j)) {
                            jointProb01.setQuick(i, j, jointProb01.getQuick(i, j) + 1.0);
                        } else if (gty.paternalChrom.getQuick(i) && !gty.paternalChrom.getQuick(j)) {
                            jointProb10.setQuick(i, j, jointProb10.getQuick(i, j) + 1.0);
                        } else {
                            jointProb11.setQuick(i, j, jointProb11.getQuick(i, j) + 1.0);
                        }

                        if (!gty.maternalChrom.getQuick(i) && !gty.maternalChrom.getQuick(j)) {
                            jointProb00.setQuick(i, j, jointProb00.getQuick(i, j) + 1.0);
                        } else if (!gty.maternalChrom.getQuick(i) && gty.maternalChrom.getQuick(j)) {
                            jointProb01.setQuick(i, j, jointProb01.getQuick(i, j) + 1.0);
                        } else if (gty.maternalChrom.getQuick(i) && !gty.maternalChrom.getQuick(j)) {
                            jointProb10.setQuick(i, j, jointProb10.getQuick(i, j) + 1.0);
                        } else {
                            jointProb11.setQuick(i, j, jointProb11.getQuick(i, j) + 1.0);
                        }
                    }
                }
            }
        }

        for (int i = 0; i < snpNum; i++) {
            for (int j = i; j < snpNum; j++) {
                jointProb00.setQuick(i, j, jointProb00.getQuick(i, j) / chromSize);
                jointProb01.setQuick(i, j, jointProb01.getQuick(i, j) / chromSize);
                jointProb10.setQuick(i, j, jointProb10.getQuick(i, j) / chromSize);
                jointProb11.setQuick(i, j, jointProb11.getQuick(i, j) / chromSize);
            }
            margProb0.setQuick(i, margProb0.getQuick(i) / chromSize);
        }
        double tmpSq = indivSize;
        for (int i = 0; i < snpNum; i++) {
            for (int j = i; j < snpNum; j++) {
                // double maxRS = Math.min(margProb0.get(i) * (1 - margProb0.get(j)), margProb0.get(j) * (1 - margProb0.get(i)));
                // double minRS = Math.max(-margProb0.get(i) * margProb0.get(j), -(1 - margProb0.get(i)) * (1 - margProb0.get(j)));
                double sq = (jointProb00.getQuick(i, j) * jointProb11.getQuick(i, j) - jointProb01.getQuick(i, j) * jointProb10.getQuick(i, j));
                // sq = sq / Math.sqrt(margProb0.getQuick(i) * (1 - margProb0.getQuick(i)) * margProb0.getQuick(j) * (1 - margProb0.getQuick(j)));
                sq = sq * sq / (margProb0.getQuick(i) * (1 - margProb0.getQuick(i)) * margProb0.getQuick(j) * (1 - margProb0.getQuick(j)));
                //correction
                // sq = (tmpSq * sq - 1) / (tmpSq - 3);
                //johny's correction                
               // sq = 1 - (tmpSq - 3) / (tmpSq - 2) * (1 - sq) * (1 + 2 * (1 - sq) / (tmpSq - 3.3));

                corMat.setQuick(i, j, sq);
                corMat.setQuick(j, i, sq);
            }
        }
        // System.out.println(corMat);

        return corMat;
    }

    public DoubleMatrix2D calculateHillRobertsonLD(List<Individual> allIndiv, List<AnnotSNP> snpList) throws Exception {
        int indivSize = allIndiv.size();
        int chromSize = indivSize * 2;
        int snpNum = snpList.size();
        DoubleMatrix2D corMat = new DenseDoubleMatrix2D(snpNum, snpNum);
        DoubleMatrix2D jointProb00 = new DenseDoubleMatrix2D(snpNum, snpNum);
        DoubleMatrix2D jointProb11 = new DenseDoubleMatrix2D(snpNum, snpNum);
        DoubleMatrix2D jointProb01 = new DenseDoubleMatrix2D(snpNum, snpNum);
        DoubleMatrix2D jointProb10 = new DenseDoubleMatrix2D(snpNum, snpNum);
        DoubleMatrix1D margProb0 = new DenseDoubleMatrix1D(snpNum);

        for (int i = 0; i < snpNum; i++) {
            margProb0.setQuick(i, 0);
            for (int j = i; j < snpNum; j++) {
                corMat.setQuick(i, j, 0.0);
            }
        }

        int pos1, pos2;
        //note the chromsome number is alwalys eaqual to  indivSize * 2 at every locus due to the missing genotypes
        for (int k = 0; k < indivSize; k++) {
            StatusGtySet gty = allIndiv.get(k).markerGtySet;
            for (int i = 0; i < snpNum; i++) {
                AnnotSNP snp1 = snpList.get(i);
                pos1 = snp1.order;
                if (gty.existence.getQuick(pos1)) {
                    if (!gty.paternalChrom.getQuick(pos1)) {
                        margProb0.setQuick(i, margProb0.getQuick(i) + 1);
                    }

                    if (!gty.maternalChrom.getQuick(pos1)) {
                        margProb0.setQuick(i, margProb0.getQuick(i) + 1);
                    }
                }

                for (int j = i; j < snpNum; j++) {
                    AnnotSNP snp2 = snpList.get(j);
                    pos2 = snp2.order;

                    if (gty.existence.getQuick(pos1) && gty.existence.getQuick(pos2)) {
                        if (!gty.paternalChrom.getQuick(pos1) && !gty.paternalChrom.getQuick(pos2)) {
                            jointProb00.setQuick(i, j, jointProb00.getQuick(i, j) + 1.0);
                        } else if (!gty.paternalChrom.getQuick(pos1) && gty.paternalChrom.getQuick(pos2)) {
                            jointProb01.setQuick(i, j, jointProb01.getQuick(i, j) + 1.0);
                        } else if (gty.paternalChrom.getQuick(pos1) && !gty.paternalChrom.getQuick(pos2)) {
                            jointProb10.setQuick(i, j, jointProb10.getQuick(i, j) + 1.0);
                        } else {
                            jointProb11.setQuick(i, j, jointProb11.getQuick(i, j) + 1.0);
                        }

                        if (!gty.maternalChrom.getQuick(pos1) && !gty.maternalChrom.getQuick(pos2)) {
                            jointProb00.setQuick(i, j, jointProb00.getQuick(i, j) + 1.0);
                        } else if (!gty.maternalChrom.getQuick(pos1) && gty.maternalChrom.getQuick(pos2)) {
                            jointProb01.setQuick(i, j, jointProb01.getQuick(i, j) + 1.0);
                        } else if (gty.maternalChrom.getQuick(pos1) && !gty.maternalChrom.getQuick(pos2)) {
                            jointProb10.setQuick(i, j, jointProb10.getQuick(i, j) + 1.0);
                        } else {
                            jointProb11.setQuick(i, j, jointProb11.getQuick(i, j) + 1.0);
                        }
                    }
                }
            }
        }

        for (int i = 0; i < snpNum; i++) {
            for (int j = i; j < snpNum; j++) {
                jointProb00.setQuick(i, j, jointProb00.getQuick(i, j) / chromSize);
                jointProb01.setQuick(i, j, jointProb01.getQuick(i, j) / chromSize);
                jointProb10.setQuick(i, j, jointProb10.getQuick(i, j) / chromSize);
                jointProb11.setQuick(i, j, jointProb11.getQuick(i, j) / chromSize);
            }
            margProb0.setQuick(i, margProb0.getQuick(i) / chromSize);
        }
        double tmpSq = 0;
        for (int i = 0; i < snpNum; i++) {
            for (int j = i; j < snpNum; j++) {
                // double maxRS = Math.min(margProb0.get(i) * (1 - margProb0.get(j)), margProb0.get(j) * (1 - margProb0.get(i)));
                // double minRS = Math.max(-margProb0.get(i) * margProb0.get(j), -(1 - margProb0.get(i)) * (1 - margProb0.get(j)));
                double sq = (jointProb00.getQuick(i, j) * jointProb11.getQuick(i, j) - jointProb01.getQuick(i, j) * jointProb10.getQuick(i, j));
                sq = sq / Math.sqrt(margProb0.getQuick(i) * (1 - margProb0.getQuick(i)) * margProb0.getQuick(j) * (1 - margProb0.getQuick(j)));
                corMat.setQuick(i, j, sq);
                corMat.setQuick(j, i, sq);
            }
        }
        // System.out.println(corMat);

        return corMat;
    }

    public DoubleMatrix2D calculateGenotypeCorrelation(List<Individual> allIndiv, int snpNum) throws Exception {
        int indivSize = allIndiv.size();
        DoubleMatrix2D corMat = new DenseDoubleMatrix2D(snpNum, snpNum);
        corMat.assign(0.0);
        DoubleArrayList genotypeCodeList1 = new DoubleArrayList();
        DoubleArrayList genotypeCodeList2 = new DoubleArrayList();

        //note the chromsome number is alwalys eaqual to  indivSize * 2 at every locus due to the missing genotypes
        for (int j = 0; j < snpNum; j++) {
            corMat.setQuick(j, j, 1);
            for (int t = j + 1; t < snpNum; t++) {
                genotypeCodeList1.clear();
                genotypeCodeList2.clear();
                for (int k = 0; k < indivSize; k++) {
                    StatusGtySet gty = allIndiv.get(k).markerGtySet;
                    if (gty.existence.getQuick(j) && gty.existence.getQuick(t)) {
                        if (!gty.paternalChrom.getQuick(j) && !gty.maternalChrom.getQuick(j)) {
                            genotypeCodeList1.add(0);
                        } else if (gty.paternalChrom.getQuick(j) && gty.maternalChrom.getQuick(j)) {
                            genotypeCodeList1.add(2);
                        } else {
                            genotypeCodeList1.add(1);
                        }

                        if (!gty.paternalChrom.getQuick(t) && !gty.maternalChrom.getQuick(t)) {
                            genotypeCodeList2.add(0);
                        } else if (gty.paternalChrom.getQuick(t) && gty.maternalChrom.getQuick(t)) {
                            genotypeCodeList2.add(2);
                        } else {
                            genotypeCodeList2.add(1);
                        }
                    }
                }
                double mean1 = Descriptive.mean(genotypeCodeList1);
                double mean2 = Descriptive.mean(genotypeCodeList2);
                double sd1 = Descriptive.sampleVariance(genotypeCodeList1, mean1);
                double sd2 = Descriptive.sampleVariance(genotypeCodeList2, mean2);
                double r = Descriptive.correlation(genotypeCodeList1, Math.sqrt(sd1), genotypeCodeList2, Math.sqrt(sd2));
                // r = r * r;
                // r = (indivSize * r - 1) / (indivSize - 3);
                corMat.setQuick(j, t, r);
                corMat.setQuick(t, j, r);
            }
        }
        return corMat;
    }

    public DoubleMatrix2D calculateGenotypeCorrelation(List<Individual> allIndiv, int snpNum, int diseaseStatus) throws Exception {
        int indivSize = allIndiv.size();
        DoubleMatrix2D corMat = new DenseDoubleMatrix2D(snpNum, snpNum);
        corMat.assign(0.0);
        DoubleArrayList genotypeCodeList1 = new DoubleArrayList();
        DoubleArrayList genotypeCodeList2 = new DoubleArrayList();
        double tmpSq = indivSize * 2;
        //note the chromsome number is alwalys eaqual to  indivSize * 2 at every locus due to the missing genotypes
        for (int j = 0; j < snpNum; j++) {
            corMat.setQuick(j, j, 1);
            for (int t = j + 1; t < snpNum; t++) {
                genotypeCodeList1.clear();
                genotypeCodeList2.clear();
                for (int k = 0; k < indivSize; k++) {
                    if (allIndiv.get(k).getAffectedStatus() != diseaseStatus) {
                        continue;
                    }
                    StatusGtySet gty = allIndiv.get(k).markerGtySet;
                    if (gty.existence.getQuick(j) && gty.existence.getQuick(t)) {
                        if (!gty.paternalChrom.getQuick(j) && !gty.maternalChrom.getQuick(j)) {
                            genotypeCodeList1.add(2);
                        } else if (gty.paternalChrom.getQuick(j) && gty.maternalChrom.getQuick(j)) {
                            genotypeCodeList1.add(0);
                        } else {
                            genotypeCodeList1.add(1);
                        }

                        if (!gty.paternalChrom.getQuick(t) && !gty.maternalChrom.getQuick(t)) {
                            genotypeCodeList2.add(2);
                        } else if (gty.paternalChrom.getQuick(t) && gty.maternalChrom.getQuick(t)) {
                            genotypeCodeList2.add(0);
                        } else {
                            genotypeCodeList2.add(1);
                        }
                    }
                }
                double mean1 = Descriptive.mean(genotypeCodeList1);
                double mean2 = Descriptive.mean(genotypeCodeList2);
                double sd1 = Descriptive.sampleVariance(genotypeCodeList1, mean1);
                double sd2 = Descriptive.sampleVariance(genotypeCodeList2, mean2);
                double r = Descriptive.correlation(genotypeCodeList1, Math.sqrt(sd1), genotypeCodeList2, Math.sqrt(sd2));
                /*
                sd1 = 0;
                mean1 = 0;
                for (int a = 0; a < genotypeCodeList1.size(); a++) {
                sd1 += (genotypeCodeList1.getQuick(a) * genotypeCodeList1.getQuick(a));
                mean1 += genotypeCodeList1.getQuick(a);
                }
                mean1 = mean1 / genotypeCodeList1.size();
                sd1 = ((sd1 - genotypeCodeList1.size() * mean1 * mean1) / (genotypeCodeList1.size()));
                // sd1 = Math.sqrt(sd1);
                sd2 = 0;
                mean2 = 0;
                for (int a = 0; a < genotypeCodeList2.size(); a++) {
                sd2 += (genotypeCodeList2.getQuick(a) * genotypeCodeList2.getQuick(a));
                mean2 += genotypeCodeList2.getQuick(a);
                }
                mean2 = mean2 / genotypeCodeList2.size();
                sd2 = ((sd2 - genotypeCodeList2.size() * mean2 * mean2) / (genotypeCodeList2.size()));
                // sd2 = Math.sqrt(sd2);
                r = Descriptive.correlation(genotypeCodeList1, Math.sqrt(sd1), genotypeCodeList2, Math.sqrt(sd2));
                double cov = 0;
                for (int a = 0; a < genotypeCodeList2.size(); a++) {
                cov += (genotypeCodeList1.getQuick(a) * genotypeCodeList2.getQuick(a));
                }
                r = (cov - genotypeCodeList2.size() * mean1 * mean2) / genotypeCodeList2.size();                
                r = r / (Math.sqrt(sd1) * Math.sqrt(sd2));
                 */
                r = r * r;
                //correction
                //sq = (tmpSq * sq - 1) / (tmpSq - 3);
                //johny's correction                
                r = 1 - (tmpSq - 3) / (tmpSq - 2) * (1 - r) * (1 + 2 * (1 - r) / (tmpSq - 3.3));
                corMat.setQuick(j, t, r);
                corMat.setQuick(t, j, r);
            }
        }
        return corMat;
    }

    public List<Individual> extractGenotypes(List<Individual> allIndiv, int diseaseStatus) throws Exception {
        int indivSize = allIndiv.size();
        List<Individual> newIndvis = new ArrayList<Individual>();

        for (int k = 0; k < indivSize; k++) {
            if (allIndiv.get(k).getAffectedStatus() != diseaseStatus) {
                continue;
            }
            newIndvis.add(allIndiv.get(k));
        }
        return newIndvis;
    }

    public DoubleMatrix2D calculateGenotypeCorrelation(List<Individual> allIndiv, List<AnnotSNP> snpList) throws Exception {
        int indivSize = allIndiv.size();
        int snpNum = snpList.size();
        DoubleMatrix2D corMat = new DenseDoubleMatrix2D(snpNum, snpNum);
        corMat.assign(0.0);
        DoubleArrayList genotypeCodeList1 = new DoubleArrayList();
        DoubleArrayList genotypeCodeList2 = new DoubleArrayList();


        //note the chromsome number is alwalys eaqual to  indivSize * 2 at every locus due to the missing genotypes
        for (int j = 0; j < snpNum; j++) {
            int pos1 = snpList.get(j).order;
            corMat.setQuick(j, j, 1);
            for (int t = j + 1; t < snpNum; t++) {
                int pos2 = snpList.get(t).order;
                genotypeCodeList1.clear();
                genotypeCodeList2.clear();
                for (int k = 0; k < indivSize; k++) {
                    StatusGtySet gty = allIndiv.get(k).markerGtySet;
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
                double mean1 = Descriptive.mean(genotypeCodeList1);
                double mean2 = Descriptive.mean(genotypeCodeList2);
                double sd1 = Descriptive.sampleVariance(genotypeCodeList1, mean1);
                double sd2 = Descriptive.sampleVariance(genotypeCodeList2, mean2);
                double r = Descriptive.correlation(genotypeCodeList1, Math.sqrt(sd1), genotypeCodeList2, Math.sqrt(sd2));
                corMat.setQuick(j, t, r);
                corMat.setQuick(t, j, r);
            }
        }
        return corMat;
    }

    public DoubleMatrix2D calculate2DProbability(List<Individual> allIndiv, int startIndex, int endIndex) throws Exception {
        int indivSize = allIndiv.size();
        int chromSize = indivSize * 2;
        int snpNum = endIndex - startIndex + 1;
        int index1, index2;
        DoubleMatrix2D twoDProb = new DenseDoubleMatrix2D(snpNum, snpNum);
        for (int i = 0; i < snpNum; i++) {
            for (int j = i; j < snpNum; j++) {
                twoDProb.setQuick(i, j, 0.0);
            }
        }

        for (int k = 0; k < indivSize; k++) {
            StatusGtySet gty = allIndiv.get(k).markerGtySet;
            for (int i = 0; i < snpNum; i++) {
                for (int j = i; j < snpNum; j++) {
                    index1 = i + startIndex;
                    index2 = j + startIndex;
                    if (gty.existence.getQuick(index1) && gty.existence.getQuick(index2)) {
                        if (gty.paternalChrom.getQuick(index1) && gty.paternalChrom.getQuick(index2)) {
                            twoDProb.setQuick(i, j, twoDProb.getQuick(i, j) + 1.0);
                        }
                        if (gty.maternalChrom.getQuick(index1) && gty.maternalChrom.getQuick(index2)) {
                            twoDProb.setQuick(i, j, twoDProb.getQuick(i, j) + 1.0);
                        }
                    }
                }
            }
        }

        for (int i = 0; i < snpNum; i++) {
            for (int j = i; j < snpNum; j++) {
                twoDProb.setQuick(i, j, twoDProb.getQuick(i, j) / chromSize);
            }
        }
        // LocalString.print2DRealMatrix(twoDProb);
        return twoDProb;
    }

    public void summarizePhenotypeProperty(List<Individual> allIndiv, List<LiabilityDiseaseSNP> snpList) throws Exception {
        int suscepNum = snpList.size();

        int[] suscepLoci = new int[suscepNum];
        boolean[] riskAlleles = new boolean[suscepNum];
        int indivSize = allIndiv.size();

        for (int j = 0; j < suscepNum; j++) {
            suscepLoci[j] = snpList.get(j).getLDMarkerPosition();
            riskAlleles[j] = snpList.get(j).isRiskAlleleLable();
        }


        double prevalence = 0.0;
        double[][] penerances = new double[suscepNum][3];
        int[][] accounts = new int[suscepNum][3];
        for (int i = 0; i < suscepNum; i++) {
            Arrays.fill(penerances[i], 0);
            Arrays.fill(accounts[i], 0);
        }

        int affected = 0;
        for (int k = 0; k < indivSize; k++) {
            Individual indiv = allIndiv.get(k);
            affected = indiv.getAffectedStatus();
            if (affected == 2) {
                prevalence += 1.0;
            }

            StatusGtySet gty = indiv.markerGtySet;
            for (int i = 0; i < suscepNum; i++) {
                if (gty.existence.getQuick(suscepLoci[i])) {
                    //ensure 2 denote 2 risk alleles and 1 1 risk allele and 0 no risk allele
                    if (gty.paternalChrom.getQuick(suscepLoci[i]) == riskAlleles[i] && gty.maternalChrom.getQuick(suscepLoci[i]) == riskAlleles[i]) {
                        accounts[i][2] += 1;
                        if (affected == 2) {
                            penerances[i][2] += 1.0;
                        }
                    } else if (gty.paternalChrom.getQuick(suscepLoci[i]) != riskAlleles[i] && gty.maternalChrom.getQuick(suscepLoci[i]) != riskAlleles[i]) {
                        accounts[i][0] += 1;
                        if (affected == 2) {
                            penerances[i][0] += 1.0;
                        }
                    } else {
                        accounts[i][1] += 1;
                        if (affected == 2) {
                            penerances[i][1] += 1.0;
                        }
                    }
                }
            }
        }
        /*
        int count = 0;
        StatusGtySet gty1 = indiv.traitGtySet;
        for (int i = 0; i < suscepNum; i++) {
        if (gty1.paternalChrom.getQuick(i) != gty.paternalChrom.getQuick(i)) {
        count++;
        }
        if (gty1.maternalChrom.getQuick(i) != gty.maternalChrom.getQuick(i)) {
        count++;
        }
        }
        System.out.println(count);
         */
        DecimalFormat df = new DecimalFormat("0.0000000", new DecimalFormatSymbols(Locale.US));
        StringBuffer infor = new StringBuffer();
        infor.append("Prevalence: ");
        infor.append(df.format(prevalence / indivSize));
        infor.append("\n");
        infor.append("riskRatio1\triskRatio2\n");
        BufferedWriter genePBw = new BufferedWriter(new FileWriter("geneticRisk.txt", true));
        for (int i = 0; i < suscepNum; i++) {
            for (int j = 0; j < 3; j++) {
                penerances[i][j] /= accounts[i][j];
            }
            infor.append(df.format(penerances[i][1] / penerances[i][0]));
            infor.append("\t");
            infor.append(df.format(penerances[i][2] / penerances[i][0]));
            infor.append("\n");
            //save in  hard disk
            genePBw.write(df.format(penerances[i][1] / penerances[i][0]));
            genePBw.write("\t");
            genePBw.write(df.format(penerances[i][1] / penerances[i][0]));
            genePBw.write("\t");
        }
        genePBw.newLine();
        genePBw.close();

        System.out.println(infor);
    }

    public void summarizePhenotypePropertyContinuousTrait(List<Individual> allIndiv, int[] suscepLoci, boolean[] riskAlleles) throws Exception {
        int suscepNum = suscepLoci.length;
        int indivSize = allIndiv.size();

        double prevalence = 0.0;
        double[][] penerances = new double[suscepNum][3];
        int[][] accounts = new int[suscepNum][3];
        for (int i = 0; i < suscepNum; i++) {
            Arrays.fill(penerances[i], 0);
            Arrays.fill(accounts[i], 0);
        }
        int traitNum = allIndiv.get(0).getMainTrait().length;

        DoubleArrayList[] traitValues = new DoubleArrayList[traitNum];
        for (int i = 0; i < traitNum; i++) {
            traitValues[i] = new DoubleArrayList();
        }
        for (int k = 0; k < indivSize; k++) {
            Individual indiv = allIndiv.get(k);
            for (int i = 0; i < traitNum; i++) {
                traitValues[i].add(indiv.getMainTrait()[i]);
            }

            StatusGtySet gty = indiv.markerGtySet;
            for (int i = 0; i < suscepNum; i++) {
                if (gty.existence.getQuick(suscepLoci[i])) {
                    //ensure 2 denote 2 risk alleles and 1 1 risk allele and 0 no risk allele
                    if (gty.paternalChrom.getQuick(suscepLoci[i]) == riskAlleles[i] && gty.maternalChrom.getQuick(suscepLoci[i]) == riskAlleles[i]) {
                        accounts[i][2] += 1;

                    } else if (gty.paternalChrom.getQuick(suscepLoci[i]) != riskAlleles[i] && gty.maternalChrom.getQuick(suscepLoci[i]) != riskAlleles[i]) {
                        accounts[i][0] += 1;

                    } else {
                        accounts[i][1] += 1;

                    }
                }
            }
        }

        for (int i = 0; i < traitNum; i++) {
            double meanT = Descriptive.mean(traitValues[i]);
            double varT = Descriptive.sampleVariance(traitValues[i], meanT);
           System.out.println("Mean: " + meanT + "   Variance: " + varT);
            for (int k = 0; k < indivSize; k++) {
                Individual indiv = allIndiv.get(k);
                indiv.getMainTrait()[i] -= meanT;
            }
        }
        if (traitNum > 1) {
            for (int i = 0; i < traitNum; i++) {
                for (int j = i + 1; j < traitNum; j++) {
                    double cov = Descriptive.covariance(traitValues[i], traitValues[j]);
                    double mean1 = Descriptive.mean(traitValues[i]);
                    double mean2 = Descriptive.mean(traitValues[j]);
                    double sd1 = Descriptive.sampleVariance(traitValues[i], mean1);
                    double sd2 = Descriptive.sampleVariance(traitValues[j], mean2);
                    double corr = Descriptive.correlation(traitValues[i], sd1, traitValues[j], sd2);
                    System.out.println("Covariance: " + cov + "   Correlation: " + corr);
                }
            }
        }

        /*
        DecimalFormat df = new DecimalFormat("0.0000000", new DecimalFormatSymbols(Locale.US));
        StringBuffer infor = new StringBuffer();
        infor.append("Prevalence: ");
        infor.append(df.format(prevalence / indivSize));
        infor.append("\n");
        infor.append("riskRatio1\triskRatio2\n");
        BufferedWriter genePBw = new BufferedWriter(new FileWriter("geneticRisk.txt", true));
        for (int i = 0; i < suscepNum; i++) {
        }
        
        genePBw.close();
        System.out.println(infor);
         */
    }

    public void convertQuantiativeTrait2BinaryTrait(List<Individual> allIndiv, double cutoff, boolean isOver) throws Exception {
        int indivSize = allIndiv.size();
        double prevalence = 0.0;

        DoubleArrayList traitValues = new DoubleArrayList();
        for (int k = 0; k < indivSize; k++) {
            Individual indiv = allIndiv.get(k);
            traitValues.add(indiv.getMainTrait()[0]);
            if (isOver) {
                if (indiv.getMainTrait()[0] >= cutoff) {
                    indiv.setAffectedStatus(2);
                    prevalence += 1;
                } else {
                    indiv.setAffectedStatus(1);
                }
            } else {
                if (indiv.getMainTrait()[0] <= cutoff) {
                    indiv.setAffectedStatus(2);
                    prevalence += 1;
                } else {
                    indiv.setAffectedStatus(1);
                }
            }
        }

        prevalence = prevalence / indivSize;
      //  System.out.println("Prevelance: " + prevalence);

        double meanT = Descriptive.mean(traitValues);
        double varT = Descriptive.sampleVariance(traitValues, meanT);
        //System.out.println(meanT + " " + varT + " Var " + varT);
        BufferedWriter testBw = new BufferedWriter(new FileWriter("test.txt", true));
        for (int i = 0; i < indivSize; i++) {
            testBw.write(String.valueOf(traitValues.getQuick(i)));
            testBw.write("\n");
        }
        testBw.close();
    }

    public void convertQuantiativeTrait2TwoBinaryTrait(List<Individual> allIndiv, double cutoff, boolean isOver, int traitNum) throws Exception {
        int indivSize = allIndiv.size();
        double prevalence = 0.0;

        DoubleArrayList traitValues = new DoubleArrayList();
        for (int t = 0; t < traitNum; t++) {
            prevalence = 0.0;
            for (int k = 0; k < indivSize; k++) {
                Individual indiv = allIndiv.get(k);
                traitValues.add(indiv.getMainTrait()[t]);
                if (isOver) {
                    if (indiv.getMainTrait()[t] >= cutoff) {
                        indiv.getMainTrait()[t] = 2;
                        prevalence += 1;
                    } else {
                        indiv.getMainTrait()[t] = 1;
                    }
                } else {
                    if (indiv.getMainTrait()[t] <= cutoff) {
                        indiv.getMainTrait()[t] = 2;
                        prevalence += 1;
                    } else {
                        indiv.getMainTrait()[t] = 1;
                    }
                }
            }
            prevalence = prevalence / indivSize;
            System.out.println("Prevelance: " + prevalence);

            double meanT = Descriptive.mean(traitValues);
            double varT = Descriptive.sampleVariance(traitValues, meanT);
            System.out.println(meanT + " " + varT + " Var " + varT);
        }

    }

    public void summarizePhenotypeProperty1(List<Individual> allIndiv, List<DiseaseSNP> snpList) throws Exception {
        int suscepNum = snpList.size();

        int[] suscepLoci = new int[suscepNum];
        boolean[] riskAlleles = new boolean[suscepNum];
        int indivSize = allIndiv.size();

        for (int j = 0; j < suscepNum; j++) {
            suscepLoci[j] = snpList.get(j).getLDMarkerPosition();
            riskAlleles[j] = snpList.get(j).isRiskAlleleLable();
        }


        double prevalence = 0.0;
        double[][] penerances = new double[suscepNum][3];
        int[][] accounts = new int[suscepNum][3];
        for (int i = 0; i < suscepNum; i++) {
            Arrays.fill(penerances[i], 0);
            Arrays.fill(accounts[i], 0);
        }

        int affected = 0;
        for (int k = 0; k < indivSize; k++) {
            Individual indiv = allIndiv.get(k);
            affected = indiv.getAffectedStatus();
            if (affected == 2) {
                prevalence += 1.0;
            }

            StatusGtySet gty = indiv.markerGtySet;
            for (int i = 0; i < suscepNum; i++) {
                if (gty.existence.getQuick(suscepLoci[i])) {
                    //ensure 2 denote 2 risk alleles and 1 1 risk allele and 0 no risk allele
                    if (gty.paternalChrom.getQuick(suscepLoci[i]) == riskAlleles[i] && gty.maternalChrom.getQuick(suscepLoci[i]) == riskAlleles[i]) {
                        accounts[i][2] += 1;
                        if (affected == 2) {
                            penerances[i][2] += 1.0;
                        }
                    } else if (gty.paternalChrom.getQuick(suscepLoci[i]) == riskAlleles[i] || gty.maternalChrom.getQuick(suscepLoci[i]) == riskAlleles[i]) {
                        accounts[i][1] += 1;
                        if (affected == 2) {
                            penerances[i][1] += 1.0;
                        }
                    } else {
                        accounts[i][0] += 1;
                        if (affected == 2) {
                            penerances[i][0] += 1.0;
                        }
                    }
                }
            }
        }
        /*
        int count = 0;
        StatusGtySet gty1 = indiv.traitGtySet;
        for (int i = 0; i < suscepNum; i++) {
        if (gty1.paternalChrom.getQuick(i) != gty.paternalChrom.getQuick(i)) {
        count++;
        }
        if (gty1.maternalChrom.getQuick(i) != gty.maternalChrom.getQuick(i)) {
        count++;
        }
        }
        System.out.println(count);
         */
        DecimalFormat df = new DecimalFormat("0.0000000", new DecimalFormatSymbols(Locale.US));
        StringBuffer infor = new StringBuffer();
        infor.append("Prevalence: ");
        infor.append(df.format(prevalence / indivSize));
        infor.append("\n");
        infor.append("riskRatio1\triskRatio2\n");
        BufferedWriter genePBw = new BufferedWriter(new FileWriter("geneticRisk.txt", true));
        for (int i = 0; i < suscepNum; i++) {
            for (int j = 0; j < 3; j++) {
                penerances[i][j] /= accounts[i][j];
            }
            infor.append(df.format(penerances[i][1] / penerances[i][0]));
            infor.append("\t");
            infor.append(df.format(penerances[i][2] / penerances[i][0]));
            infor.append("\n");
            //save in  hard disk
            genePBw.write(df.format(penerances[i][1] / penerances[i][0]));
            genePBw.write("\t");
            genePBw.write(df.format(penerances[i][1] / penerances[i][0]));
            genePBw.write("\t");
        }
        genePBw.newLine();
        genePBw.close();

        System.out.println(infor);
    }

    public void summarizePhenotypeProperty1(List<Individual> allIndiv, List<DiseaseSNP> snpList, IntArrayList caseIndexes, IntArrayList controlIndexes) throws Exception {
        int suscepNum = 0;
        if (snpList != null) {
            suscepNum = snpList.size();
        }

        int[] suscepLoci = new int[suscepNum];
        boolean[] riskAlleles = new boolean[suscepNum];
        int indivSize = allIndiv.size();

        for (int j = 0; j < suscepNum; j++) {
            suscepLoci[j] = snpList.get(j).getLDMarkerPosition();
            riskAlleles[j] = snpList.get(j).isRiskAlleleLable();
        }


        double prevalence = 0.0;
        double[][] penerances = new double[suscepNum][3];
        int[][] accounts = new int[suscepNum][3];
        for (int i = 0; i < suscepNum; i++) {
            Arrays.fill(penerances[i], 0);
            Arrays.fill(accounts[i], 0);
        }

        int affected = 0;
        for (int k = 0; k < indivSize; k++) {
            Individual indiv = allIndiv.get(k);
            affected = indiv.getAffectedStatus();
            if (affected == 2) {
                prevalence += 1.0;
                caseIndexes.add(k);
            } else {
                controlIndexes.add(k);
            }



            StatusGtySet gty = indiv.markerGtySet;
            for (int i = 0; i < suscepNum; i++) {
                if (gty.existence.getQuick(suscepLoci[i])) {
                    //ensure 2 denote 2 risk alleles and 1 1 risk allele and 0 no risk allele
                    if (gty.paternalChrom.getQuick(suscepLoci[i]) == riskAlleles[i] && gty.maternalChrom.getQuick(suscepLoci[i]) == riskAlleles[i]) {
                        accounts[i][2] += 1;
                        if (affected == 2) {
                            penerances[i][2] += 1.0;
                        }
                    } else if (gty.paternalChrom.getQuick(suscepLoci[i]) == riskAlleles[i] || gty.maternalChrom.getQuick(suscepLoci[i]) == riskAlleles[i]) {
                        accounts[i][1] += 1;
                        if (affected == 2) {
                            penerances[i][1] += 1.0;
                        }
                    } else {
                        accounts[i][0] += 1;
                        if (affected == 2) {
                            penerances[i][0] += 1.0;
                        }
                    }
                }
            }
        }
        /*
        int count = 0;
        StatusGtySet gty1 = indiv.traitGtySet;
        for (int i = 0; i < suscepNum; i++) {
        if (gty1.paternalChrom.getQuick(i) != gty.paternalChrom.getQuick(i)) {
        count++;
        }
        if (gty1.maternalChrom.getQuick(i) != gty.maternalChrom.getQuick(i)) {
        count++;
        }
        }
        System.out.println(count);
         */
        DecimalFormat df = new DecimalFormat("0.0000000", new DecimalFormatSymbols(Locale.US));
        StringBuffer infor = new StringBuffer();
        infor.append("Prevalence: ");
        infor.append(df.format(prevalence / indivSize));
        infor.append("\n");
        infor.append("riskRatio1\triskRatio2\n");
        BufferedWriter genePBw = new BufferedWriter(new FileWriter("geneticRisk.txt", true));
        for (int i = 0; i < suscepNum; i++) {
            for (int j = 0; j < 3; j++) {
                penerances[i][j] /= accounts[i][j];
            }
            infor.append(df.format(penerances[i][1] / penerances[i][0]));
            infor.append("\t");
            infor.append(df.format(penerances[i][2] / penerances[i][0]));
            infor.append("\n");
            //save in  hard disk
            genePBw.write(df.format(penerances[i][1] / penerances[i][0]));
            genePBw.write("\t");
            genePBw.write(df.format(penerances[i][1] / penerances[i][0]));
            genePBw.write("\t");
        }
        genePBw.newLine();
        genePBw.close();

         System.out.println(infor);
    }

    public void summarizePhenotypeProperty2(List<Individual> allIndiv, List<DiseaseSNP> snpList, IntArrayList caseIndexes, IntArrayList controlIndexes) throws Exception {
        int suscepNum = 0;
        if (snpList != null) {
            suscepNum = snpList.size();
        }

        int[] suscepLoci = new int[suscepNum];
        boolean[] riskAlleles = new boolean[suscepNum];
        int indivSize = allIndiv.size();

        for (int j = 0; j < suscepNum; j++) {
            suscepLoci[j] = snpList.get(j).getLDMarkerPosition();
            riskAlleles[j] = snpList.get(j).isRiskAlleleLable();
        }


        double prevalence = 0.0;
        double[][] penerances = new double[suscepNum][3];
        int[][] accounts = new int[suscepNum][3];
        for (int i = 0; i < suscepNum; i++) {
            Arrays.fill(penerances[i], 0);
            Arrays.fill(accounts[i], 0);
        }

        int affected = 0;
        for (int k = 0; k < indivSize; k++) {
            Individual indiv = allIndiv.get(k);
            affected = (int) indiv.getMainTrait()[0];
            if (affected == 2) {
                prevalence += 1.0;
                caseIndexes.add(k);
            } else {
                controlIndexes.add(k);
            }



            StatusGtySet gty = indiv.markerGtySet;
            for (int i = 0; i < suscepNum; i++) {
                if (gty.existence.getQuick(suscepLoci[i])) {
                    //ensure 2 denote 2 risk alleles and 1 1 risk allele and 0 no risk allele
                    if (gty.paternalChrom.getQuick(suscepLoci[i]) == riskAlleles[i] && gty.maternalChrom.getQuick(suscepLoci[i]) == riskAlleles[i]) {
                        accounts[i][2] += 1;
                        if (affected == 2) {
                            penerances[i][2] += 1.0;
                        }
                    } else if (gty.paternalChrom.getQuick(suscepLoci[i]) == riskAlleles[i] || gty.maternalChrom.getQuick(suscepLoci[i]) == riskAlleles[i]) {
                        accounts[i][1] += 1;
                        if (affected == 2) {
                            penerances[i][1] += 1.0;
                        }
                    } else {
                        accounts[i][0] += 1;
                        if (affected == 2) {
                            penerances[i][0] += 1.0;
                        }
                    }
                }
            }
        }
        /*
        int count = 0;
        StatusGtySet gty1 = indiv.traitGtySet;
        for (int i = 0; i < suscepNum; i++) {
        if (gty1.paternalChrom.getQuick(i) != gty.paternalChrom.getQuick(i)) {
        count++;
        }
        if (gty1.maternalChrom.getQuick(i) != gty.maternalChrom.getQuick(i)) {
        count++;
        }
        }
        System.out.println(count);
         */
        DecimalFormat df = new DecimalFormat("0.0000000", new DecimalFormatSymbols(Locale.US));
        StringBuffer infor = new StringBuffer();
        infor.append("Prevalence: ");
        infor.append(df.format(prevalence / indivSize));
        infor.append("\n");
        infor.append("riskRatio1\triskRatio2\n");
        BufferedWriter genePBw = new BufferedWriter(new FileWriter("geneticRisk.txt", true));
        for (int i = 0; i < suscepNum; i++) {
            for (int j = 0; j < 3; j++) {
                penerances[i][j] /= accounts[i][j];
            }
            infor.append(df.format(penerances[i][1] / penerances[i][0]));
            infor.append("\t");
            infor.append(df.format(penerances[i][2] / penerances[i][0]));
            infor.append("\n");
            //save in  hard disk
            genePBw.write(df.format(penerances[i][1] / penerances[i][0]));
            genePBw.write("\t");
            genePBw.write(df.format(penerances[i][1] / penerances[i][0]));
            genePBw.write("\t");
        }
        genePBw.newLine();
        genePBw.close();

        System.out.println(infor);
    }

    public boolean checkNullHypothesis(List<Individual> allIndiv, List<DiseaseSNP> snpList, double tolerance) throws Exception {
        int suscepNum = snpList.size();

        int[] suscepLoci = new int[suscepNum];
        boolean[] riskAlleles = new boolean[suscepNum];
        int indivSize = allIndiv.size();

        for (int j = 0; j < suscepNum; j++) {
            suscepLoci[j] = snpList.get(j).getLDMarkerPosition();
            riskAlleles[j] = snpList.get(j).isRiskAlleleLable();
        }


        double prevalence = 0.0;
        double[][] penerances = new double[suscepNum][3];
        int[][] accounts = new int[suscepNum][3];
        for (int i = 0; i < suscepNum; i++) {
            Arrays.fill(penerances[i], 0);
            Arrays.fill(accounts[i], 0);
        }

        int affected = 0;
        for (int k = 0; k < indivSize; k++) {
            Individual indiv = allIndiv.get(k);
            affected = indiv.getAffectedStatus();
            if (affected == 2) {
                prevalence += 1.0;
            }

            StatusGtySet gty = indiv.markerGtySet;
            for (int i = 0; i < suscepNum; i++) {
                if (gty.existence.getQuick(suscepLoci[i])) {
                    //ensure 2 denote 2 risk alleles and 1 1 risk allele and 0 no risk allele
                    if (gty.paternalChrom.getQuick(suscepLoci[i]) == riskAlleles[i] && gty.maternalChrom.getQuick(suscepLoci[i]) == riskAlleles[i]) {
                        accounts[i][2] += 1;
                        if (affected == 2) {
                            penerances[i][2] += 1.0;
                        }
                    } else if (gty.paternalChrom.getQuick(suscepLoci[i]) == riskAlleles[i] || gty.maternalChrom.getQuick(suscepLoci[i]) == riskAlleles[i]) {
                        accounts[i][1] += 1;
                        if (affected == 2) {
                            penerances[i][1] += 1.0;
                        }
                    } else {
                        accounts[i][0] += 1;
                        if (affected == 2) {
                            penerances[i][0] += 1.0;
                        }
                    }
                }
            }
        }
        /*
        int count = 0;
        StatusGtySet gty1 = indiv.traitGtySet;
        for (int i = 0; i < suscepNum; i++) {
        if (gty1.paternalChrom.getQuick(i) != gty.paternalChrom.getQuick(i)) {
        count++;
        }
        if (gty1.maternalChrom.getQuick(i) != gty.maternalChrom.getQuick(i)) {
        count++;
        }
        }
        System.out.println(count);
         */
        boolean isNull = true;
        DecimalFormat df = new DecimalFormat("0.0000000", new DecimalFormatSymbols(Locale.US));
        StringBuffer infor = new StringBuffer();
        infor.append("Prevalence: ");
        infor.append(df.format(prevalence / indivSize));
        infor.append("\n");
        infor.append("riskRatio1\triskRatio2\n");
        BufferedWriter genePBw = new BufferedWriter(new FileWriter("geneticRisk.txt", true));
        for (int i = 0; i < suscepNum; i++) {
            for (int j = 0; j < 3; j++) {
                penerances[i][j] /= accounts[i][j];
            }
            double ratio0 = penerances[i][1] / penerances[i][0];
            infor.append(df.format(ratio0));
            infor.append("\t");
            genePBw.write(df.format(ratio0));
            genePBw.write("\t");
            if (Math.abs(ratio0 - 1) > tolerance) {
                isNull = false;
            }
            double ratio1 = penerances[i][2] / penerances[i][0];
            infor.append(df.format(ratio1));
            infor.append("\n");
            genePBw.write(df.format(ratio1));
            genePBw.write("\n");
            if (Math.abs(ratio1 - 1) > tolerance) {
                isNull = false;
            }
        }
        genePBw.close();

        System.out.println(infor);
        return isNull;
    }

    public void summarizePhenotypePropertyforLogitDiseaseSNP(List<Individual> allIndiv, List<LogitDiseaseSNP> snpList) throws Exception {
        int suscepNum = snpList.size();
        int[] suscepLoci = new int[suscepNum];
        boolean[] riskAlleles = new boolean[suscepNum];
        int indivSize = allIndiv.size();
        for (int j = 0; j < suscepNum; j++) {
            suscepLoci[j] = snpList.get(j).getLDMarkerPosition();
            riskAlleles[j] = snpList.get(j).isRiskAlleleLable();
        }
        double prevalence = 0.0;
        double[][] penerances = new double[suscepNum][3];
        int[][] accounts = new int[suscepNum][3];
        for (int i = 0; i < suscepNum; i++) {
            Arrays.fill(penerances[i], 0);
            Arrays.fill(accounts[i], 0);
        }

        for (int k = 0; k < indivSize; k++) {
            Individual indiv = allIndiv.get(k);
            if (indiv.getAffectedStatus() == 2) {
                prevalence += 1.0;
            }

            StatusGtySet gty = indiv.markerGtySet;
            for (int i = 0; i < suscepNum; i++) {
                if (gty.existence.getQuick(suscepLoci[i])) {
                    //ensure 2 denote 2 risk alleles and 1 1 risk allele and 0 no risk allele
                    if (gty.paternalChrom.getQuick(suscepLoci[i]) == riskAlleles[i] && gty.maternalChrom.getQuick(suscepLoci[i]) == riskAlleles[i]) {
                        accounts[i][2] += 1;
                        if (indiv.getAffectedStatus() == 2) {
                            penerances[i][2] += 1.0;
                        }
                    } else if (gty.paternalChrom.getQuick(suscepLoci[i]) != riskAlleles[i] && gty.maternalChrom.getQuick(suscepLoci[i]) != riskAlleles[i]) {
                        accounts[i][0] += 1;
                        if (indiv.getAffectedStatus() == 2) {
                            penerances[i][0] += 1.0;
                        }
                    } else {
                        accounts[i][1] += 1;
                        if (indiv.getAffectedStatus() == 2) {
                            penerances[i][1] += 1.0;
                        }
                    }
                }
            }
        }
        /*
        int count = 0;
        StatusGtySet gty1 = indiv.traitGtySet;
        for (int i = 0; i < suscepNum; i++) {
        if (gty1.paternalChrom.getQuick(i) != gty.paternalChrom.getQuick(i)) {
        count++;
        }
        if (gty1.maternalChrom.getQuick(i) != gty.maternalChrom.getQuick(i)) {
        count++;
        }
        }
        System.out.println(count);
         */
        DecimalFormat df = new DecimalFormat("0.0000000", new DecimalFormatSymbols(Locale.US));
        StringBuffer infor = new StringBuffer();
        infor.append("Prevalence: ");
        infor.append(df.format(prevalence / indivSize));
        infor.append("\n");
        infor.append("riskRatio1\triskRatio2\n");
        for (int i = 0; i < suscepNum; i++) {
            for (int j = 0; j < 3; j++) {
                penerances[i][j] /= accounts[i][j];
            }
            infor.append(df.format(penerances[i][1] / penerances[i][0]));
            infor.append("\t");
            infor.append(df.format(penerances[i][2] / penerances[i][0]));
            infor.append("\n");
        }
        System.out.println(infor);
    }

    /*
    public void summarizePhenotypeProperty(List<Individual> allIndiv, List<LogitDiseaseSNP> snpList) throws Exception {
    int suscepNum = snpList.size();
    int[] suscepLoci = new int[suscepNum];
    boolean[] riskAlleles = new boolean[suscepNum];
    int indivSize = allIndiv.size();
    for (int j = 0; j < suscepNum; j++) {
    suscepLoci[j] = snpList.get(j).getLDMarkerPosition();
    riskAlleles[j] = snpList.get(j).isRiskAlleleLable();
    }
    double prevalence = 0.0;
    double[][] penerances = new double[suscepNum][3];
    int[][] accounts = new int[suscepNum][3];
    for (int i = 0; i < suscepNum; i++) {
    Arrays.fill(penerances[i], 0);
    Arrays.fill(accounts[i], 0);
    }
    
    for (int k = 0; k < indivSize; k++) {
    Individual indiv = allIndiv.get(k);
    if (indiv.getAffectedStatus() == 2) {
    prevalence += 1.0;
    }
    StatusGtySet gty = indiv.markerGtySet;
    for (int i = 0; i < suscepNum; i++) {
    if (gty.existence.getQuick(suscepLoci[i])) {
    //ensure 2 denote 2 risk alleles and 1 1 risk allele and 0 no risk allele
    if (gty.paternalChrom.getQuick(suscepLoci[i]) == riskAlleles[i] && gty.maternalChrom.getQuick(suscepLoci[i]) == riskAlleles[i]) {
    accounts[i][2] += 1;
    if (indiv.getAffectedStatus() == 2) {
    penerances[i][2] += 1.0;
    }
    
    } else if (!gty.paternalChrom.getQuick(suscepLoci[i]) == riskAlleles[i] && gty.maternalChrom.getQuick(suscepLoci[i]) != riskAlleles[i]) {
    accounts[i][0] += 1;
    if (indiv.getAffectedStatus() == 2) {
    penerances[i][0] += 1.0;
    }
    } else {
    accounts[i][1] += 1;
    if (indiv.getAffectedStatus() == 2) {
    penerances[i][1] += 1.0;
    }
    }
    }
    }
    }
    DecimalFormat df = new DecimalFormat("0.0000000", new DecimalFormatSymbols(Locale.US));
    StringBuffer infor = new StringBuffer();
    infor.append("Prevalence: ");
    infor.append(df.format(prevalence / indivSize));
    infor.append("\n");
    infor.append("riskRatio1\triskRatio2\n");
    for (int i = 0; i < suscepNum; i++) {
    for (int j = 0; j < 3; j++) {
    penerances[i][j] /= accounts[i][j];
    }
    infor.append(df.format(penerances[i][1] / penerances[i][0]));
    infor.append("\t");
    infor.append(df.format(penerances[i][2] / penerances[i][0]));
    infor.append("\n");
    }
    System.out.println(infor);
    }
     */
}
