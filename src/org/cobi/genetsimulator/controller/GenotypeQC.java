/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.genetsimulator.controller;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.cobi.genetsimulator.entity.AnnotSNP;
import org.cobi.genetsimulator.entity.GenotypeBasedLDSparseMatrix;
import org.cobi.genetsimulator.entity.Individual;
import org.cobi.genetsimulator.entity.StatusGtySet;

/**
 *
 * @author MX Li
 */
public class GenotypeQC {

    //not thread theft
    public String removeMissingCaseControl(List<AnnotSNP> mainSnpMap, List<Individual> indList, double missingRate, int[] counts, int chrID) throws Exception {
        StringBuilder msg = new StringBuilder();
        int snpPos;
        boolean allMiss;
        List<AnnotSNP> tmpSNPMap = new ArrayList<AnnotSNP>();

        tmpSNPMap.addAll(mainSnpMap);
        mainSnpMap.clear();
        double nonmissingRate = 1 - missingRate;
        int totalCounterCase, totalCountercontrol;
        int counterCase, countercontrol;
        int counter = 0;
        int listSize = tmpSNPMap.size();
        int indivSize = indList.size();
        for (int i = 0; i < listSize; i++) {
            AnnotSNP snp = tmpSNPMap.get(i);
            snpPos = snp.order;
            allMiss = true;
            totalCounterCase = 0;
            totalCountercontrol = 0;
            counterCase = 0;
            countercontrol = 0;
            for (int j = 0; j < indivSize; j++) {
                Individual mIndivi = indList.get(j);
                totalCountercontrol++;
                StatusGtySet gtySet = mIndivi.markerGtySets[chrID];
                if (gtySet.existence.getQuick(snpPos)) {
                    countercontrol++;
                }
            }

            if ((((double) countercontrol) / totalCountercontrol < nonmissingRate) || (((double) counterCase) / totalCounterCase < nonmissingRate)) {

                msg.append(" ");
                counter++;
            } else {
                mainSnpMap.add(snp);
            }

        }

        counts[0] = mainSnpMap.size();
        counts[1] = listSize;
        String info = (mainSnpMap.size() + " SNPs (out of " + listSize + ") pass the fliter with missing call rate over " + missingRate);

        if (msg.length() > 0) {
            //return info+"\n"+counter + " SNP are excluded due to no genotypes:\n" + msg.toString();
            return info;
        } else {
            return info;
        }
    }

    //not thread theft
    public String removeMissingCaseControl(List<AnnotSNP> mainSnpMap, List<Individual> indList, double missingRate, int[] counts) throws Exception {
        StringBuilder msg = new StringBuilder();
        int snpPos;
        boolean allMiss;
        List<AnnotSNP> tmpSNPMap = new ArrayList<AnnotSNP>();

        tmpSNPMap.addAll(mainSnpMap);
        mainSnpMap.clear();
        double nonmissingRate = 1 - missingRate;
        int totalCounterCase, totalCountercontrol;
        int counterCase, countercontrol;
        int counter = 0;
        int listSize = tmpSNPMap.size();
        int indivSize = indList.size();
        for (int i = 0; i < listSize; i++) {
            AnnotSNP snp = tmpSNPMap.get(i);
            snpPos = snp.order;
            allMiss = true;
            totalCounterCase = 0;
            totalCountercontrol = 0;
            counterCase = 0;
            countercontrol = 0;
            for (int j = 0; j < indivSize; j++) {
                Individual mIndivi = indList.get(j);
                totalCountercontrol++;
                StatusGtySet gtySet = mIndivi.markerGtySet;
                if (gtySet.existence.getQuick(snpPos)) {
                    countercontrol++;
                }
            }

            if ((((double) countercontrol) / totalCountercontrol < nonmissingRate) || (((double) counterCase) / totalCounterCase < nonmissingRate)) {

                msg.append(" ");
                counter++;
            } else {
                mainSnpMap.add(snp);
            }

        }

        counts[0] = mainSnpMap.size();
        counts[1] = listSize;
        String info = (mainSnpMap.size() + " SNPs (out of " + listSize + ") pass the fliter with missing call rate over " + missingRate);

        if (msg.length() > 0) {
            //return info+"\n"+counter + " SNP are excluded due to no genotypes:\n" + msg.toString();
            return info;
        } else {
            return info;
        }
    }

    public String removeByMAF(List<AnnotSNP> mainSnpMap, List<Individual> indList, double minAF, int chrID) throws Exception {
        StringBuilder msg = new StringBuilder();
        int snpPos;
        boolean allMiss;
        List<AnnotSNP> tmpSNPMap = new ArrayList<AnnotSNP>();

        tmpSNPMap.addAll(mainSnpMap);
        mainSnpMap.clear();
        double alleleAFreq = 0.0;
        int counter0, counter1;

        int counter = 0;
        int listSize = tmpSNPMap.size();
        int indivSize = indList.size();
        for (int i = 0; i < listSize; i++) {
            AnnotSNP snp = tmpSNPMap.get(i);
            snpPos = snp.order;
            allMiss = true;
            counter0 = 0;
            counter1 = 0;
            for (int j = 0; j < indivSize; j++) {
                Individual mIndivi = indList.get(j);
                StatusGtySet gtySet = mIndivi.markerGtySets[chrID];

                if (gtySet.existence.getQuick(snpPos)) {
                    if (gtySet.paternalChrom.getQuick(snpPos)) {
                        counter1++;
                    } else {
                        counter0++;
                    }
                    if (gtySet.maternalChrom.getQuick(snpPos)) {
                        counter1++;
                    } else {
                        counter0++;
                    }
                }
            }
            alleleAFreq = counter0 * 1.0 / (counter1 + counter0);
            if (alleleAFreq > 0.5) {
                alleleAFreq = 1 - alleleAFreq;
            }
            if (alleleAFreq < minAF) {
                msg.append(snp.getRSID());
                msg.append(" ");
                counter++;
            } else {
                mainSnpMap.add(snp);
            }
        }
        System.out.println(mainSnpMap.size() + " SNPs (out of " + listSize + ") are lefted after the filtering!");
        tmpSNPMap = null;
        if (msg.length() > 0) {
            return counter + " SNP are excluded due to no genotypes:\n" + msg.toString();
        } else {
            return "";
        }
    }

    public String ldPruning(List<AnnotSNP> mainSnpMap, GenotypeBasedLDSparseMatrix ldCorr, double maxCorr, boolean ingoreNOGty) throws Exception {
        List<AnnotSNP> tmpSNPMap = new ArrayList<AnnotSNP>();
        tmpSNPMap.addAll(mainSnpMap);
        mainSnpMap.clear();
        int listSize = tmpSNPMap.size();
        Set<Integer> highlyCorrIndexes = new HashSet<Integer>();
        int windowSize = 50;
        int stepLen = 5;
        double r, c;
        int[] counts = new int[2];
        for (int s = 0; s < listSize; s += stepLen) {
            for (int i = s; (i - s <= windowSize) && (i < listSize); i++) {
                AnnotSNP snp1 = tmpSNPMap.get(i);

                if (highlyCorrIndexes.contains(i)) {
                    continue;
                }
                for (int j = i + 1; (j - i <= windowSize) && (j < listSize); j++) {
                    if (highlyCorrIndexes.contains(j)) {
                        continue;
                    }
                    AnnotSNP snp2 = tmpSNPMap.get(j);

                    r = ldCorr.getLDAt(snp1.getPhysicalPosition(), snp2.getPhysicalPosition());
                    /*
                     r = Math.sqrt(ldRMatrix.getQuick(i, j));
                     //for R
                     c = (0.6065 * r - 1.033) * r + 1.7351;
                     if (c > 2) {
                     c = 2;
                     }
                     */

                    //R2
                    //R2
                    //y = -35.741x6 + 111.16x5 - 128.42x4 + 66.906x3 - 14.641x2 + 0.6075x + 0.8596
                    //c = (((((-35.741 * r + 111.16) * r - 128.42) * r + 66.906) * r - 14.641) * r + 0.6075) * r + 0.8596;
                    //y = 0.2725x2 - 0.3759x + 0.8508
                    //c = (0.2725 * r - 0.3759) * r + 0.8508;
                    // y = 0.2814x2 - 0.4308x + 0.86
                    //c = (0.2814 * r - 0.4308) * r + 0.86;
                    //y = -0.155x + 0.8172
                    //c = -0.155 * r + 0.8172;
                    r = r * r;
                    // r = Math.pow(r, c);
                    if (r >= maxCorr) {
                        highlyCorrIndexes.add(j);
                    }
                }
            }
        }

        counts[0] = listSize - highlyCorrIndexes.size();
        counts[1] = listSize;
        String info = (listSize - (highlyCorrIndexes.size()) + " SNPs (out of " + listSize + ") passed LD pruning (r2>=" + maxCorr + ").");

        for (int s = 0; s < listSize; s++) {
            AnnotSNP snp1 = tmpSNPMap.get(s);
            if (!highlyCorrIndexes.contains(s)) {
                mainSnpMap.add(snp1);
            }
        }
        return info;
    }

    public void outCorrelation(List<AnnotSNP> mainSnpMap, GenotypeBasedLDSparseMatrix ldCorr, int startIndex, int endIndex, String outPath) throws Exception {
        List<AnnotSNP> tmpSNPMap = new ArrayList<AnnotSNP>();
        tmpSNPMap.addAll(mainSnpMap);
        mainSnpMap.clear();
        int listSize = tmpSNPMap.size();

        double r, c;
        BufferedWriter snpPBw = new BufferedWriter(new FileWriter(outPath));
        for (int i = startIndex; i < endIndex; i++) {
            AnnotSNP snp1 = tmpSNPMap.get(i);
            snpPBw.write(String.valueOf(snp1.getaAlleleFreq()));
            for (int j = startIndex; j < endIndex; j++) {
                AnnotSNP snp2 = tmpSNPMap.get(j);
                r = ldCorr.getLDAt(snp1.getPhysicalPosition(), snp2.getPhysicalPosition());
                snpPBw.write("\t" + r);
            }
            snpPBw.write("\n");
        }
        snpPBw.close();
    }

    public String removeByMAF(List<AnnotSNP> mainSnpMap, List<Individual> indList, double minAF) throws Exception {
        StringBuilder msg = new StringBuilder();
        int snpPos;
        boolean hasMissing;
        List<AnnotSNP> tmpSNPMap = new ArrayList<AnnotSNP>();

        tmpSNPMap.addAll(mainSnpMap);
        mainSnpMap.clear();
        double alleleAFreq = 0.0;
        int counter0, counter1;

        int counter = 0;
        int listSize = tmpSNPMap.size();
        int indivSize = indList.size();
        double numMissing = 0;
        for (int i = 0; i < listSize; i++) {
            AnnotSNP snp = tmpSNPMap.get(i);

            snpPos = snp.order;
            hasMissing = false;
            counter0 = 0;
            counter1 = 0;
            numMissing = 0;
            for (int j = 0; j < indivSize; j++) {
                Individual mIndivi = indList.get(j);
                StatusGtySet gtySet = mIndivi.markerGtySet;

                if (gtySet.existence.getQuick(snpPos)) {
                    if (gtySet.paternalChrom.getQuick(snpPos)) {
                        counter1++;
                    } else {
                        counter0++;
                    }
                    if (gtySet.maternalChrom.getQuick(snpPos)) {
                        counter1++;
                    } else {
                        counter0++;
                    }
                } else {
                    //    hasMissing = true;
                    //    break;
                    numMissing++;
                }
            }
            if (hasMissing) {
                // continue;
            }
            if (numMissing / indivSize > 0.5) {
                continue;
            }
            alleleAFreq = counter0 * 1.0 / (counter1 + counter0);
            if (alleleAFreq > 0.5) {
                alleleAFreq = 1 - alleleAFreq;
            }
            if (alleleAFreq < minAF) {
                msg.append(snp.getRSID());
                msg.append(" ");
                counter++;
            } else {
                mainSnpMap.add(snp);
            }
        }
        System.out.println(mainSnpMap.size() + " SNPs (out of " + listSize + ") are lefted after the MAF filtering (" + minAF + ")");
        tmpSNPMap = null;
        if (msg.length() > 0) {
            return counter + " SNP are excluded due to no genotypes:\n" + msg.toString();
        } else {
            return "";
        }
    }

    public void calculateMAF(List<AnnotSNP> mainSnpMap, List<Individual> indList, boolean[] maLabels, double[] mafs) throws Exception {
        int snpPos;
        double alleleAFreq = 0.0;
        int counter0, counter1;

        int counter = 0;
        int listSize = mainSnpMap.size();
        int indivSize = indList.size();
        for (int i = 0; i < listSize; i++) {
            AnnotSNP snp = mainSnpMap.get(i);
            snpPos = snp.order;

            counter0 = 0;
            counter1 = 0;
            for (int j = 0; j < indivSize; j++) {
                Individual mIndivi = indList.get(j);
                if (mIndivi.markerGtySet.existence.getQuick(snpPos)) {
                    if (mIndivi.markerGtySet.paternalChrom.getQuick(snpPos)) {
                        counter1++;
                    } else {
                        counter0++;
                    }
                    if (mIndivi.markerGtySet.maternalChrom.getQuick(snpPos)) {
                        counter1++;
                    } else {
                        counter0++;
                    }
                }
            }

            alleleAFreq = counter0 * 1.0 / (counter1 + counter0);
            if (alleleAFreq > 0.5) {
                maLabels[i] = false;
                mafs[i] = 1 - alleleAFreq;
            } else {
                maLabels[i] = true;
                mafs[i] = alleleAFreq;
            }
            //  System.out.println(mafs[i] + " " + maLabels[i]);
        }
    }

    public void calculateMAF(List<AnnotSNP> mainSnpMap, List<Individual> indList, int chrID) throws Exception {
        int snpPos;
        double alleleAFreq = 0.0;
        int counter0, counter1;

        int counter = 0;
        int listSize = mainSnpMap.size();
        int indivSize = indList.size();
        for (int i = 0; i < listSize; i++) {
            AnnotSNP snp = mainSnpMap.get(i);
            snpPos = snp.order;

            counter0 = 0;
            counter1 = 0;
            for (int j = 0; j < indivSize; j++) {
                Individual mIndivi = indList.get(j);
                StatusGtySet gtySet = mIndivi.markerGtySets[chrID];
                if (gtySet.existence.getQuick(snpPos)) {
                    if (gtySet.paternalChrom.getQuick(snpPos)) {
                        counter1++;
                    } else {
                        counter0++;
                    }
                    if (gtySet.maternalChrom.getQuick(snpPos)) {
                        counter1++;
                    } else {
                        counter0++;
                    }
                }
            }

            alleleAFreq = counter0 * 1.0 / (counter1 + counter0);
            //  System.out.println(mafs[i] + " " + maLabels[i]);
            snp.setAAlleleFreq(alleleAFreq);
        }
    }

    public void calculateMAF(List<AnnotSNP> mainSnpMap, List<Individual> indList) throws Exception {
        int snpPos;
        double alleleAFreq = 0.0;
        int counter0, counter1;

        int counter = 0;
        int listSize = mainSnpMap.size();
        int indivSize = indList.size();
        for (int i = 0; i < listSize; i++) {
            AnnotSNP snp = mainSnpMap.get(i);
            snpPos = snp.order;

            counter0 = 0;
            counter1 = 0;
            for (int j = 0; j < indivSize; j++) {
                Individual mIndivi = indList.get(j);
                StatusGtySet gtySet = mIndivi.markerGtySet;
                if (gtySet.existence.getQuick(snpPos)) {
                    if (gtySet.paternalChrom.getQuick(snpPos)) {
                        counter1++;
                    } else {
                        counter0++;
                    }
                    if (gtySet.maternalChrom.getQuick(snpPos)) {
                        counter1++;
                    } else {
                        counter0++;
                    }
                }
            }

            alleleAFreq = counter0 * 1.0 / (counter1 + counter0);
            //  System.out.println(mafs[i] + " " + maLabels[i]);
            snp.setAAlleleFreq(alleleAFreq);
        }
    }
}
