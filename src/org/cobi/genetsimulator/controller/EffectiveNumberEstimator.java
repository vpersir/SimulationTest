/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.genetsimulator.controller;

import cern.colt.list.IntArrayList;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import java.io.BufferedWriter;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

/**
 *
 * @author mxli
 */
public class EffectiveNumberEstimator {

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

    public static DoubleMatrix2D removeRedundantItems(DoubleMatrix2D corrMat, double maxCorr, IntArrayList indexes) {
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
            IntArrayList tmpList = new IntArrayList();
            tmpList.addAllOf(indexes);
            indexes.clear();
            DoubleMatrix2D poweredCorrMat = new DenseDoubleMatrix2D(newSampleSize, newSampleSize);
            int incRow = 0;
            int incCol = 0;
            for (int i = 0; i < originalSampleSize; i++) {
                if (highlyCorrIndexes.contains(i)) {
                    continue;
                }
                indexes.add(tmpList.getQuick(i));
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

    public static double calculateEffectSampleSizeColtMatrixMyMethodByPCov(DoubleMatrix2D corrMat, Set<Integer> selectedSampleIndex) throws Exception {
        DoubleMatrix2D poweredCorrMat = corrMat.copy();
        int originalSampleSize = poweredCorrMat.columns();
        // DoubleMatrix2D corrMat = ColtMatrixBasic.readMatrixFromFile("test.txt", originalSampleSize, originalSampleSize);
        int newSampleSize = originalSampleSize;

        if (selectedSampleIndex != null && selectedSampleIndex.size() > 0) {
            // System.out.println("Removed columns and rows: " + highlyCorrIndexes.toString());
            newSampleSize = selectedSampleIndex.size();

            DoubleMatrix2D tmpCorMat = new DenseDoubleMatrix2D(newSampleSize, newSampleSize);
            int incRow = 0;
            int incCol = 0;
            for (int i = 0; i < originalSampleSize; i++) {
                if (!selectedSampleIndex.contains(i)) {
                    continue;
                }
                incCol = 0;
                for (int j = 0; j < originalSampleSize; j++) {
                    if (!selectedSampleIndex.contains(j)) {
                        continue;
                    }
                    tmpCorMat.setQuick(incRow, incCol, poweredCorrMat.getQuick(i, j));
                    incCol++;
                }
                incRow++;
            }
            poweredCorrMat = tmpCorMat;
            // System.out.println(corrMat.toString());
        } else if (selectedSampleIndex != null && selectedSampleIndex.isEmpty()) {
            return 0;
        }

        if (newSampleSize == 1) {
            return 1;
        }

        //I found this function is less error-prone  than the  EigenDecompositionImpl 2.0 and slightly faster
        // System.out.println(poweredCorrMat.toString());
        EigenvalueDecomposition ed = new EigenvalueDecomposition(poweredCorrMat);

        DoubleMatrix1D eVR = ed.getRealEigenvalues();

        //DoubleMatrix1D eVI = ed.getImagEigenvalues();
        // System.out.println(eVR.toString());
        // System.out.println(eVI.toString());
        //double effectSampleSize = newSampleSize;
        double effectSampleSize = newSampleSize;

        for (int i = 0; i < newSampleSize; i++) {
            if (Double.isNaN(eVR.get(i))) {
                System.err.println("NaN error for eigen values!");
            }
            if (eVR.getQuick(i) > 1) {
                effectSampleSize -= (eVR.getQuick(i) - 1);
            }
        }
        return (effectSampleSize);
    }

    public static double calculateEigenvalueSum(DoubleMatrix2D corrMat) throws Exception {

        //DoubleMatrix1D eVI = ed.getImagEigenvalues();
        // System.out.println(eVR.toString());
        // System.out.println(eVI.toString());
        //double effectSampleSize = newSampleSize;
        double effectSampleSize = 0;

        return (effectSampleSize);
    }

    public static double calculateEffectSampleSizeLargeMatrixMyMethodByRCov(DoubleMatrix2D corrMat) throws Exception {
        int originalSampleSize = corrMat.columns();
        double effectNumSum = 0;
        DoubleMatrix2D poweredCorrMat = new DenseDoubleMatrix2D(originalSampleSize, originalSampleSize);
        for (int i = 0; i < originalSampleSize; i++) {
            poweredCorrMat.setQuick(i, i, 1);
            for (int j = i + 1; j < originalSampleSize; j++) {
                double x = corrMat.getQuick(i, j);
                x = x * x;
                //when r2
                //y = 0.7723x6 - 1.5659x5 + 1.201x4 - 0.2355x3 + 0.2184x2 + 0.6086x
                x = (((((0.7723 * x - 1.5659) * x + 1.201) * x - 0.2355) * x + 0.2184) * x + 0.6086) * x;

                if (x > 1E-6) {
                    poweredCorrMat.setQuick(i, j, x);
                    poweredCorrMat.setQuick(j, i, x);
                }
            }
        }
        // System.out.println(poweredCorrMat.toString());
        int maxExpandLen = 50;
        int forceStopLen = 250;
        int initialSite = 0;
        int startIndex = 0;
        int endIndex = startIndex + 1;
        double minThreshold = 0.05;
        int expandIndex = endIndex + 1;
        int expandLen = 1;
        boolean isSucceed = false;
        initialSite = startIndex;
        while (startIndex < originalSampleSize) {
            //find a site whose LD beteen startIndex start to be less than weakCorrelationThreshold
            while (poweredCorrMat.getQuick(startIndex, endIndex) >= minThreshold) {
                endIndex++;
                if (endIndex >= originalSampleSize) {
                    break;
                }
            }
            //start from next  stopIndex
            expandIndex = endIndex + 1;
            if (expandIndex >= originalSampleSize) {
                DoubleMatrix2D par = poweredCorrMat.viewPart(initialSite, initialSite, originalSampleSize - initialSite, originalSampleSize - initialSite);
                double partEffectiveNum = calculateEffectSampleSizeColtMatrixMyMethodByPCov(par, null);
                //System.out.println(par.toString());
                effectNumSum += partEffectiveNum;
                System.out.println("Real:" + (endIndex - initialSite) + " Estimated: " + partEffectiveNum);
                break;
            }
            expandLen = 1;
            isSucceed = false;
            //check more sites whose LD beteen startIndex is still less than weakCorrelationThreshold
            while (poweredCorrMat.getQuick(startIndex, expandIndex) < minThreshold) {
                expandLen++;
                expandIndex++;
                if (expandLen >= maxExpandLen || expandIndex >= originalSampleSize) {
                    isSucceed = true;
                    break;
                }
            }
            if (isSucceed) {
                //if it gots 
                startIndex++;
                if (startIndex == endIndex) {
                    DoubleMatrix2D par = poweredCorrMat.viewPart(initialSite, initialSite, endIndex - initialSite, endIndex - initialSite);
                    double partEffectiveNum = calculateEffectSampleSizeColtMatrixMyMethodByPCov(par, null);
                    //System.out.println(par.toString());
                    effectNumSum += partEffectiveNum;
                    System.out.println("Real:" + (endIndex - initialSite) + " Estimated: " + partEffectiveNum);
                    initialSite = startIndex;
                    //search at a new position
                    endIndex = startIndex + 2;
                }
            } else if (endIndex - initialSite >= forceStopLen) {
                DoubleMatrix2D par = poweredCorrMat.viewPart(initialSite, initialSite, endIndex - initialSite, endIndex - initialSite);
                double partEffectiveNum = calculateEffectSampleSizeColtMatrixMyMethodByPCov(par, null);
                System.out.println("Real:" + (endIndex - initialSite) + " Estimated: " + partEffectiveNum);
                effectNumSum += partEffectiveNum;
                initialSite = endIndex;
                startIndex = endIndex;
                endIndex = startIndex + 1;
            } else {
                //re-check
                startIndex = initialSite;
                endIndex = expandIndex;
            }

        }

        return effectNumSum;
    }

    public static double fastCalculateEffectSampleSizeLDSparseMatrixMyMethodByRCov(LDSparseMatrix ldRsMatrix,
            int maxBlockLen, double weakCorrelationThreshold,
            int maxCheckingNum, int maxCheckingDistance, BufferedWriter bwLog, String chromName) throws Exception {
        Set<Integer> snpPositionSet = ldRsMatrix.getAllUniqueIndexes();
        int totalObservedSize = snpPositionSet.size();
        if (totalObservedSize == 0) {
            return 0;
        } else if (totalObservedSize == 1) {
            return 1;
        }

        int[] snpPositionArray = new int[totalObservedSize];
        int count = 0;
        for (Iterator<Integer> iter = snpPositionSet.iterator(); iter.hasNext();) {
            snpPositionArray[count] = iter.next();
            count++;
        }

        Arrays.sort(snpPositionArray);
        IntArrayList indexInBlock = new IntArrayList();

        int originalSampleSize = snpPositionArray.length;
        double effectNumSum = 0;

        // System.out.println(poweredCorrMat.toString());
        int inBlockFirst = 0;
        int inBlockLast = 0;
        int outBlockFirst = 1;

        int checkingIndex = outBlockFirst;
        int checkingLen = 1;
        boolean allAreLess = false;

        boolean removeRedundant = true;
        double redundanctLD = 0.9988;
        StringBuilder info = new StringBuilder();
        boolean debug = false;
        double partEffectiveNum = 0;
        int movedRow;
        while (inBlockLast <= originalSampleSize) {
            //find a site whose LD beteen inBlockFirst start to be less than weakCorrelationThreshold
            while ((outBlockFirst < originalSampleSize) && (outBlockFirst - inBlockFirst < maxBlockLen) && (ldRsMatrix.getLDAt(snpPositionArray[inBlockFirst], snpPositionArray[outBlockFirst]) >= weakCorrelationThreshold)) {
                outBlockFirst++;
            }
            inBlockLast = outBlockFirst - 1;

            if (outBlockFirst >= originalSampleSize) {
                if (debug) {
                    System.out.print(chromName + "\t" + snpPositionArray[inBlockFirst] + "\t" + snpPositionArray[originalSampleSize - 1] + "\t" + (originalSampleSize - inBlockFirst) + "\t");
                }
                for (int i = inBlockFirst; i < originalSampleSize; i++) {
                    indexInBlock.add(snpPositionArray[i]);
                }
                DoubleMatrix2D par = ldRsMatrix.subDenseLDMatrix(indexInBlock);
                ldRsMatrix.releaseLDData();
                indexInBlock.clear();
                if (removeRedundant) {
                    par = EffectiveNumberEstimator.removeRedundantItems(par, redundanctLD);
                }

                //partEffectiveNum = calculateEffectSampleSizeColtMatrixMyMethodByPCov(par, null);
                partEffectiveNum = calculateEffectSampleSizeLargeMatrixMyMethodByRCov(par);
                // System.out.println(par.toString());
                effectNumSum += partEffectiveNum;
                info.append(chromName).append("\t").append(snpPositionArray[inBlockFirst]).append("\t").append(snpPositionArray[originalSampleSize - 1]).append("\t").append(originalSampleSize - inBlockFirst).append("\t").append(partEffectiveNum).append("\n");
                if (debug) {
                    System.out.println(partEffectiveNum);
                }
                break;
            }

            if (outBlockFirst - inBlockFirst >= maxBlockLen) {
                if (debug) {
                    System.out.print(chromName + "\t" + snpPositionArray[inBlockFirst] + "\t" + snpPositionArray[inBlockLast] + "\t" + (outBlockFirst - inBlockFirst) + "\t");
                }

                for (int i = inBlockFirst; i < outBlockFirst; i++) {
                    indexInBlock.add(snpPositionArray[i]);
                }
                DoubleMatrix2D par = ldRsMatrix.subDenseLDMatrix(indexInBlock);
                ldRsMatrix.releaseLDData();
                indexInBlock.clear();
                if (removeRedundant) {
                    par = EffectiveNumberEstimator.removeRedundantItems(par, redundanctLD);
                }
                //partEffectiveNum = calculateEffectSampleSizeColtMatrixMyMethodByPCov(par, null);
                partEffectiveNum = calculateEffectSampleSizeLargeMatrixMyMethodByRCov(par);
                info.append(chromName).append("\t").append(snpPositionArray[inBlockFirst]).append("\t").append(snpPositionArray[inBlockLast]).append("\t").append(outBlockFirst - inBlockFirst).append("\t").append(partEffectiveNum).append("\n");
                if (debug) {
                    System.out.println(partEffectiveNum);
                }
                effectNumSum += partEffectiveNum;

                inBlockFirst = outBlockFirst;
                continue;
            }

            movedRow = inBlockLast;
            //check LD  beteen inBlockLast and checkingIndex
            while (movedRow >= inBlockFirst) {
                allAreLess = false;
                checkingIndex = outBlockFirst;
                checkingLen = 1;

                while ((checkingIndex < originalSampleSize)
                        && (ldRsMatrix.getLDAt(snpPositionArray[movedRow], snpPositionArray[checkingIndex]) < weakCorrelationThreshold)) {
                    //stop the search early to avoid unncessary compare between 2 snps far away.
                    if ((checkingLen >= maxCheckingNum) && ((snpPositionArray[checkingIndex] - snpPositionArray[inBlockLast]) >= maxCheckingDistance)) {
                        allAreLess = true;
                        break;
                    }
                    checkingLen++;
                    checkingIndex++;
                }
                //it must reach the end of a relatively short region and all are passing the check
                if (checkingIndex >= originalSampleSize) {
                    allAreLess = true;
                }

                if (!allAreLess) {
                    break;
                }
                movedRow--;
            }

            if (allAreLess) {
                if (debug) {
                    System.out.print(chromName + "\t" + snpPositionArray[inBlockFirst] + "\t" + snpPositionArray[inBlockLast] + "\t" + (outBlockFirst - inBlockFirst) + "\t");
                }

                for (int i = inBlockFirst; i < outBlockFirst; i++) {
                    indexInBlock.add(snpPositionArray[i]);
                }
                DoubleMatrix2D par = ldRsMatrix.subDenseLDMatrix(indexInBlock);
                ldRsMatrix.releaseLDData();
                indexInBlock.clear();
                if (removeRedundant) {
                    par = EffectiveNumberEstimator.removeRedundantItems(par, redundanctLD);
                }
                partEffectiveNum = calculateEffectSampleSizeColtMatrixMyMethodByPCov(par, null);
                //System.out.println(par.toString());
                effectNumSum += partEffectiveNum;

                info.append(chromName).append("\t").append(snpPositionArray[inBlockFirst]).append("\t").append(snpPositionArray[inBlockLast]).append("\t").append(outBlockFirst - inBlockFirst).append("\t").append(partEffectiveNum).append("\n");
                if (debug) {
                    System.out.println(partEffectiveNum);
                }
                inBlockFirst = outBlockFirst;

            } else {
                // System.out.println(ldRsMatrix.getLDAt(snpPositionArray[startIndex], snpPositionArray[checkingIndex]) + "\t" + snpPositionArray[inBlockFirst] + "\t" + snpPositionArray[stopIndex - 1] + "\t" + (stopIndex - inBlockFirst));
                //go to the minExtendLen and re-check
                inBlockLast = checkingIndex;
                outBlockFirst = inBlockLast + 1;
                //force to cut into a block with maxBlockLen
                if (outBlockFirst - inBlockFirst >= maxBlockLen) {
                    outBlockFirst = inBlockFirst + maxBlockLen;
                }
            }
        }

        // GlobalManager.logger.info(info);
        if (bwLog != null) {
            bwLog.write(info.toString());
        }
        return effectNumSum;
    }

    public static double calculateEffectSampleSizeColtMatrixMyMethodByRCov(DoubleMatrix2D corrMat, Set<Integer> selectedSampleIndex) throws Exception {
        int originalSampleSize = corrMat.columns();
        DoubleMatrix2D poweredCorrMat = new DenseDoubleMatrix2D(originalSampleSize, originalSampleSize);
        for (int i = 0; i < originalSampleSize; i++) {
            poweredCorrMat.setQuick(i, i, 1);
            for (int j = i + 1; j < originalSampleSize; j++) {
                double x = corrMat.getQuick(i, j);
                x = x * x;
                //when r2
                //y = 0.7723x6 - 1.5659x5 + 1.201x4 - 0.2355x3 + 0.2184x2 + 0.6086x
                x = (((((0.7723 * x - 1.5659) * x + 1.201) * x - 0.2355) * x + 0.2184) * x + 0.6086) * x;

                if (x > 1E-6) {
                    poweredCorrMat.setQuick(i, j, x);
                    poweredCorrMat.setQuick(j, i, x);
                }
            }
        }
        // System.out.println(poweredCorrMat.toString());

        // DoubleMatrix2D corrMat = ColtMatrixBasic.readMatrixFromFile("test.txt", originalSampleSize, originalSampleSize);
        int newSampleSize = originalSampleSize;
        if (selectedSampleIndex != null && selectedSampleIndex.size() > 0) {
            // System.out.println("Removed columns and rows: " + highlyCorrIndexes.toString());
            newSampleSize = selectedSampleIndex.size();

            DoubleMatrix2D tmpCorMat = new DenseDoubleMatrix2D(newSampleSize, newSampleSize);
            int incRow = 0;
            int incCol = 0;
            for (int i = 0; i < originalSampleSize; i++) {
                if (!selectedSampleIndex.contains(i)) {
                    continue;
                }
                incCol = 0;
                for (int j = 0; j < originalSampleSize; j++) {
                    if (!selectedSampleIndex.contains(j)) {
                        continue;
                    }
                    tmpCorMat.setQuick(incRow, incCol, poweredCorrMat.getQuick(i, j));
                    incCol++;
                }
                incRow++;
            }
            poweredCorrMat = tmpCorMat;
            // System.out.println(corrMat.toString());
        } else if (selectedSampleIndex != null && selectedSampleIndex.isEmpty()) {
            return 0;
        }

        if (newSampleSize == 1) {
            return 1;
        }

        //I found this function is less error-prone  than the  EigenDecompositionImpl 2.0 and slightly faster
        //System.out.println(poweredCorrMat.toString());
        EigenvalueDecomposition ed = new EigenvalueDecomposition(poweredCorrMat);

        DoubleMatrix1D eVR = ed.getRealEigenvalues();

        //DoubleMatrix1D eVI = ed.getImagEigenvalues();
        // System.out.println(eVR.toString());
        // System.out.println(eVI.toString());
        //double effectSampleSize = newSampleSize;
        double effectSampleSize = newSampleSize;

        for (int i = 0; i < newSampleSize; i++) {
            if (Double.isNaN(eVR.get(i))) {
                System.err.println("NaN error for eigen values!");
            }
            if (eVR.getQuick(i) > 1) {
                effectSampleSize -= (eVR.getQuick(i) - 1);
            }
        }
        return (effectSampleSize);
    }
}
