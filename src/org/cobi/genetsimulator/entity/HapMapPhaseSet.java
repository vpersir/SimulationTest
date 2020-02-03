/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.genetsimulator.entity;

import cern.colt.bitvector.BitVector;
import cern.colt.matrix.DoubleMatrix2D;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.TreeSet;
import org.apache.commons.math.linear.RealMatrix;
import org.cobi.eqtlsimulator.Constants;
import static org.cobi.eqtlsimulator.Constants.FILE_READER_BUFFER_SIZE;
import static org.cobi.eqtlsimulator.Constants.HAP_GENOTYPE_LOAD_FACTOR;
import static org.cobi.eqtlsimulator.Constants.MISSING_ALLELE_NAME;
 
 
import org.cobi.genetsimulator.controller.MultilocusLogitModel;
import org.cobi.genetsimulator.controller.MultilocusLogitModel.LogitDiseaseSNP;
import org.cobi.genetsimulator.controller.PopuStatSummarizer;

/**
 *
 * @author mxli
 */
public class HapMapPhaseSet implements Constants {

    List<AnnotSNP> snpList;
    List<Individual> allIndiv;
    PopuStatSummarizer stater;

    public HapMapPhaseSet() {
        snpList = new ArrayList<AnnotSNP>();
        allIndiv = new ArrayList<Individual>();
        stater = new PopuStatSummarizer();
    }

    public DoubleMatrix2D get2DProbability() throws Exception {
        return stater.calculate2DProbability(allIndiv, snpList.size());
    }

    public void assignRiskLociFreq(MultilocusLogitModel logModel) throws Exception {
        List<LogitDiseaseSNP> susbSnpList = logModel.getSnpList();
        int len = susbSnpList.size();
        int[] index = new int[len];
        for (int i = 0; i < len; i++) {
            index[i] = susbSnpList.get(i).getLDMarkerPosition();
        }
        RealMatrix twoDProb = stater.calculate2DProbability(allIndiv, index);
        for (int i = 0; i < len; i++) {
            LogitDiseaseSNP ssnp = susbSnpList.get(i);
            if (ssnp.isRiskAlleleLable()) {
                ssnp.setRiskAlleleFrequency(twoDProb.getEntry(i, i));
            } else {
                ssnp.setRiskAlleleFrequency(1.0 - twoDProb.getEntry(i, i));
            }
        }
      //  logModel.setJointGtyProb(stater.calculateJointSusceptibGtyProbability(allIndiv, susbSnpList));
    }

    public int getLociNum() {
        return snpList.size();
    }

    public void readHapMapPhase(String filePath, String chroStr, long[] physicalRegions) throws Exception {
        snpList.clear();
        allIndiv.clear();
        Map<Long, AnnotSNP> tmpSnpMap = new HashMap<Long, AnnotSNP>();
        Map<String, Individual> tmpIndivGty = new HashMap<String, Individual>();
        int chroNum;
        boolean isForward = true;
        StringBuffer strRegion = new StringBuffer();
        String line = null;
        StringTokenizer tokenizer = null;
        String delmilit = "\t ,";
        Map<String, Integer> enteredHapIndiv = new HashMap<String, Integer>();

        StringBuffer runningInfo = new StringBuffer();
        TreeSet<String> unexpectedSNPSet = new TreeSet<String>();
        int totalColNum = -1;
        int startGtyCol = -1;
        int rsIndex = -1;
        int posIndex = -1;

        int availabeSNPSpace = 100;
        List<String> indiviChipLabelList = new ArrayList<String>();
        String tmpStr;
        char strandLabel = '+';
        if (!isForward) {
            strandLabel = '-';
        }

        if (chroStr.equalsIgnoreCase("X")) {
            chroNum = 23;
        } else if (chroStr.equalsIgnoreCase("Y")) {
            chroNum = 24;
        } else if (chroStr.equalsIgnoreCase("MT")) {
            chroNum = 25;
        } else if (chroStr.equalsIgnoreCase("XY")) {
            chroNum = 26;
        } else {
            chroNum = Integer.parseInt(chroStr);
        }

        int regionNum = physicalRegions.length / 2;
        boolean checkRegion = false;
        int index = 0;
        long pos = 0;

        boolean invalid = true;
        int curSNPNum = 0;

        char aAllele, bAllele;
        boolean noABAllele = false;
        for (int j = 0; j < regionNum; j++) {
            if (physicalRegions[j * 2] != -9 || physicalRegions[j * 2 + 1] != -9) {
                checkRegion = true;
                if (physicalRegions[j * 2 + 1] == -9) {
                    physicalRegions[1] = Long.MAX_VALUE;
                }
                strRegion.append(" [");
                strRegion.append(physicalRegions[0]);
                strRegion.append(",");
                strRegion.append(physicalRegions[1]);
                strRegion.append("]bp");
            }
        }

        File file = new File(filePath);
        String info = "Encoding HapMap File " + file.getName() + " ...";
        System.out.println(info);

        BufferedReader br = new BufferedReader(new FileReader(file), FILE_READER_BUFFER_SIZE);
        //read individual label
        if ((line = br.readLine()) != null) {
            tokenizer = new StringTokenizer(line, delmilit);
            totalColNum = 0;
            startGtyCol = -1;
            rsIndex = -1;
            posIndex = -1;
            indiviChipLabelList.clear();
            while (tokenizer.hasMoreTokens()) {
                tmpStr = tokenizer.nextToken().trim();
                if (tmpStr.startsWith("NA")) {
                    if (startGtyCol < 0) {
                        startGtyCol = totalColNum;
                    }
                    //some HapMap file has duplicated sample
                    int awfulIndex = tmpStr.indexOf(".");
                    if (awfulIndex >= 0) {
                        tmpStr = tmpStr.substring(0, awfulIndex);
                    }
                    tmpStr = tmpStr.substring(0, tmpStr.length() - 2);
                    if (tmpIndivGty.get(tmpStr) == null) {
                        StatusGtySet genotypes = new StatusGtySet();
                        genotypes.existence = new BitVector(availabeSNPSpace);
                        genotypes.paternalChrom = new BitVector(availabeSNPSpace);
                        genotypes.maternalChrom = new BitVector(availabeSNPSpace);
                        Individual indiv = new Individual(tmpStr);
                        indiv.markerGtySet = genotypes;
                        tmpIndivGty.put(tmpStr, indiv);
                    }
                    indiviChipLabelList.add(tmpStr);
                    if (enteredHapIndiv.get(tmpStr) != null) {
                        int times = enteredHapIndiv.get(tmpStr);
                        if (times >= 2) {
                            System.out.println("The individual " + tmpStr + " is duplicated!");
                        }
                        enteredHapIndiv.put(tmpStr, (int) 2);
                    } else {
                        enteredHapIndiv.put(tmpStr, 1);
                    }
                } else if (tmpStr.equals("rsID")) {
                    rsIndex = totalColNum;
                } else if (tmpStr.indexOf("position") >= 0) {
                    posIndex = totalColNum;
                }
                totalColNum++;
            }
        }

        String[] cells = new String[totalColNum];
        int counter = 0;
        AnnotSNP snp = null;
        //start to read genotypes
        while ((line = br.readLine()) != null) {
            line = line.toUpperCase();
            index = 0;
            tokenizer = new StringTokenizer(line, delmilit);
            invalid = false;
            aAllele = MISSING_ALLELE_NAME;
            bAllele = MISSING_ALLELE_NAME;
            noABAllele = true;

            //decide allele names
            while (tokenizer.hasMoreTokens()) {
                cells[index] = tokenizer.nextToken().trim();
                if (index == posIndex) {
                    //filter physical region
                    pos = Long.parseLong(cells[posIndex]);
                    if (checkRegion) {
                        for (int k = 0; k < regionNum; k++) {
                            if (pos < physicalRegions[k * 2] || pos > physicalRegions[k * 2 + 1]) {
                                invalid = true;
                                break;
                            }
                        }
                    }
                    snp = tmpSnpMap.get(pos);
                    if (snp != null) {
                        aAllele = snp.getAAllele();
                        bAllele = snp.getBAllele();
                        if (aAllele != MISSING_ALLELE_NAME && bAllele != MISSING_ALLELE_NAME) {
                            noABAllele = false;
                        }
                    }
                }
                if (noABAllele) {
                    if (index >= startGtyCol) {
                        if (aAllele == MISSING_ALLELE_NAME) {
                            aAllele = cells[index].charAt(0);
                            if (snp != null) {
                                snp.setAAllele(aAllele);
                            }
                        } else if (bAllele == MISSING_ALLELE_NAME && aAllele != cells[index].charAt(0)) {
                            bAllele = cells[index].charAt(0);
                            if (snp != null) {
                                snp.setBAllele(bAllele);
                            }
                            noABAllele = false;
                        }
                    }
                }
                index++;
            }
            if (invalid) {
                continue;
            }

            counter++;
            if (snp == null) {
                snp = new AnnotSNP("rs" + cells[rsIndex].substring(2), chroNum, (int)pos, 0.0, aAllele, bAllele);
                snp.setStrand(strandLabel);
                snp.order = curSNPNum;

                curSNPNum++;
                if (curSNPNum > availabeSNPSpace) {
                    //adjust the size of bitvectors
                    availabeSNPSpace += HAP_GENOTYPE_LOAD_FACTOR;
                    for (Map.Entry<String, Individual> m1 : tmpIndivGty.entrySet()) {
                        m1.getValue().markerGtySet.paternalChrom.setSize(availabeSNPSpace);
                        m1.getValue().markerGtySet.maternalChrom.setSize(availabeSNPSpace);
                        m1.getValue().markerGtySet.existence.setSize(availabeSNPSpace);
                    }
                }
                tmpSnpMap.put(pos, snp);
            }

            aAllele = snp.getAAllele();
            bAllele = snp.getBAllele();

            for (int k = startGtyCol; k < totalColNum; k += 2) {
                // make sure the genotype is made up of two alleles. otherwise regarded it as a missing data
                String indivLabel = indiviChipLabelList.get((k - startGtyCol));
                Individual indiv = tmpIndivGty.get(indivLabel);
                StatusGtySet gty = indiv.markerGtySet;

                //re-code: A->0; B->1;
                if (cells[k].charAt(0) == bAllele) {
                    gty.paternalChrom.putQuick(snp.order, true);
                    gty.existence.putQuick(snp.order, true);
                } else if (cells[k].charAt(0) == aAllele) {
                    gty.existence.putQuick(snp.order, true);
                } else if (cells[k].charAt(0) != 'N' && cells[k].charAt(0) != '0') {
                    // throw new Exception("Unexpected genotype " + cells[k] + " at " + snp.getRSID() + " for " + indivLabel);
                    String unexpectInfo = "Unexpected genotype " + cells[k] + " at " + snp.getRSID() + " for " + indivLabel;
                    //((IGG3View) GlobalVariables.currentApplication.getMainView()).setBriefRunningInfor(unexpectInfo);
                    //GlobalVariables.addInforLog(unexpectInfo);
                    unexpectedSNPSet.add(snp.getRSID());
                }
                if (cells[k + 1].charAt(0) == bAllele) {
                    gty.maternalChrom.putQuick(snp.order, true);
                } else if (cells[k + 1].charAt(0) != aAllele && cells[k + 1].charAt(0) != 'N' && cells[k + 1].charAt(0) != '0') {
                    //throw new Exception("Unexpected genotype " + cells[k + 1] + " at " + snp.getRSID() + " for " + indivLabel);
                    String unexpectInfo = "Unexpected genotype " + cells[k + 1] + " at " + snp.getRSID() + " for " + indivLabel;
                    //((IGG3View) GlobalVariables.currentApplication.getMainView()).setBriefRunningInfor(unexpectInfo);
                    //GlobalVariables.addInforLog(unexpectInfo);
                    unexpectedSNPSet.add(snp.getRSID());
                }
            }
        }
        br.close();
        if (unexpectedSNPSet.size() > 0) {
            System.out.println(unexpectedSNPSet.size() + " SNPs have unexpected genotypes in " + file);
            //GlobalVariables.addInforLog(unexpectedSNPSet.size() + " SNPs have unexpected genotypes in " + file + ":\n" + unexpectedSNPSet.toString());
            System.out.println(unexpectedSNPSet.toString());
            unexpectedSNPSet.clear();
        }

        runningInfo.delete(0, runningInfo.length());
        runningInfo.append(file + " encoded!\n");
        runningInfo.append("The number of SNPs on chromosome ");
        runningInfo.append(chroStr);
        runningInfo.append(" ");
        runningInfo.append(strRegion);
        runningInfo.append(" is ");
        runningInfo.append(counter);
        runningInfo.append(".");

        System.out.println(runningInfo.toString());
       // GlobalVariables.addInforLog(runningInfo.toString());

        int indivSize = indiviChipLabelList.size();
        for (int i = 0; i < indivSize; i += 2) {
            allIndiv.add(tmpIndivGty.get(indiviChipLabelList.get(i)));
        }
        tmpIndivGty = null;
        snpList.addAll(tmpSnpMap.values());
        Collections.sort(snpList, new SNPOrderComparator());
        tmpSnpMap = null;
    }

    public void removeMonotonyGenotype() {
        int indivSize = allIndiv.size();
        int snpSize = snpList.size();
        int polyNum = 0;
        List<Integer> indexList = new ArrayList<Integer>();

        for (int i = snpSize - 1; i >= 0; i--) {
            AnnotSNP snp = snpList.get(i);
            if (snp.getAAllele() != MISSING_ALLELE_NAME && snp.getBAllele() != MISSING_ALLELE_NAME) {
                polyNum++;
                indexList.add(snp.order);
            } else {
                snpList.remove(i);
            }
        }
        Collections.sort(indexList);

        if (polyNum < snpSize) {
            for (int j = 0; j < indivSize; j++) {
                Individual indiv = allIndiv.get(j);
                BitVector exist = new BitVector(polyNum);

                BitVector chrom1 = new BitVector(polyNum);
                BitVector chrom2 = new BitVector(polyNum);
                for (int i = 0; i < polyNum; i++) {
                    exist.putQuick(i, true);
                    chrom1.putQuick(i, indiv.markerGtySet.paternalChrom.getQuick(indexList.get(i)));
                    chrom2.putQuick(i, indiv.markerGtySet.maternalChrom.getQuick(indexList.get(i)));
                }
                indiv.markerGtySet.existence = exist;
                indiv.markerGtySet.paternalChrom = chrom1;
                indiv.markerGtySet.maternalChrom = chrom2;
            }
        }
    }
}
