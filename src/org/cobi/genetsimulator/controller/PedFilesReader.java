/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.genetsimulator.controller;

import cern.colt.bitvector.BitVector;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;

import java.util.TreeSet;
import org.cobi.eqtlsimulator.Constants;
import static org.cobi.eqtlsimulator.Constants.CHROM_NAME;
import static org.cobi.eqtlsimulator.Constants.FILE_READER_BUFFER_SIZE;
import static org.cobi.eqtlsimulator.Constants.STANDARD_CHROM_NUM;
 
 
import org.cobi.genetsimulator.entity.AnnotSNP;
import org.cobi.genetsimulator.entity.Individual;
import org.cobi.genetsimulator.entity.PedFileSet;
import static org.cobi.genetsimulator.entity.PlinkDataset.DEFAULT_MISSING_GTY_NAME;
import static org.cobi.genetsimulator.entity.PlinkDataset.MISSING_ALLELE_NAME;
import static org.cobi.genetsimulator.entity.PlinkDataset.MISSING_STRAND_NAME;
import org.cobi.genetsimulator.entity.SNPFileIndexComparator;
import org.cobi.genetsimulator.entity.SNPPosiComparator;
import org.cobi.genetsimulator.entity.StatusGtySet;

/**
 *
 * @author MX Li
 */
public class PedFilesReader implements Constants {

    public void loadPedGenoypes(PedFileSet fileSet, List<AnnotSNP> snpList, List<Individual> indList, int startChom, int endChrom) throws Exception {

        Map<Integer, Integer> snpNumPerChrom = new HashMap<Integer, Integer>();
        for (int i = startChom; i < endChrom; i++) {
            List<AnnotSNP> tmpSnpList = new ArrayList<AnnotSNP>();
            readSNPsinMapFileByPysic(fileSet, tmpSnpList, CHROM_NAME[i], null);
            if (tmpSnpList.size() == 0) {
                continue;
            }
            Collections.sort(tmpSnpList, new SNPPosiComparator());
            snpNumPerChrom.put(STANDARD_CHROM_NUM[i], tmpSnpList.size());
            snpList.addAll(tmpSnpList);
            tmpSnpList = null;
        }
        if (snpList.size() == 0) {
            return;
        }
        Collections.sort(snpList, new SNPFileIndexComparator());
        readGenotypeinPedigreeFile(fileSet, indList, snpList);
    }

    protected boolean readGenotypeinPedigreeFile(PedFileSet fileSet, List<Individual> indList, List<AnnotSNP> snpList) throws Exception {
        String pedFileName = fileSet.getPedigreeFileName();
        BufferedReader br = new BufferedReader(new FileReader(pedFileName), FILE_READER_BUFFER_SIZE);
        TreeSet<String> unexpectedSNPSet = new TreeSet<String>();
        TreeSet<String> duplicatedIndiv = new TreeSet<String>();
        String line;
        String delimiter = "\t\" \",/";
        int totalSNPSize = snpList.size();
        int gtyStaringCol = -1;

        String sigGty1;
        String sigGty2;
        int unexpectedGtyNum = 0;

        String info = "Encoding Pedigree File " + pedFileName + " ...";
        System.out.println(info);
      //  GlobalVariables.addInforLog(info);

        List<String> traitNames = new ArrayList<String>();
        traitNames.add("AffectedStatus");
        //guess genotype starting column
        gtyStaringCol = 5 + traitNames.size();
        int snpOrder = 0;

        char aAllele, bAllele;
        boolean[] needAlleleNames = new boolean[totalSNPSize];
        for (int i = 0; i < totalSNPSize; i++) {
            AnnotSNP snp = snpList.get(i);
            if (snp.getAAllele() == MISSING_ALLELE_NAME || snp.getBAllele() == MISSING_ALLELE_NAME) {
                needAlleleNames[i] = true;
            } else {
                needAlleleNames[i] = false;
            }
        }
        String indivLabel = null;
        StringBuffer tmpBuffer = new StringBuffer();
        //long t3 = System.currentTimeMillis();
        while ((line = br.readLine()) != null) {
            line = line.toUpperCase();
            StringTokenizer tokenizer = new StringTokenizer(line, delimiter);

            Individual indiv = new Individual();
            tmpBuffer.delete(0, tmpBuffer.length());
            tmpBuffer.append(tokenizer.nextToken().trim());
            indiv.setFamilyID(tmpBuffer.toString());
            tmpBuffer.delete(0, tmpBuffer.length());
            tmpBuffer.append(tokenizer.nextToken().trim());
            indiv.setIndividualID(tmpBuffer.toString());
            tmpBuffer.delete(0, tmpBuffer.length());
            tmpBuffer.append(tokenizer.nextToken().trim());
            indiv.setDadID(tmpBuffer.toString());
            tmpBuffer.delete(0, tmpBuffer.length());
            tmpBuffer.append(tokenizer.nextToken().trim());
            indiv.setMomID(tmpBuffer.toString());

            indiv.setGender(Integer.valueOf(tokenizer.nextToken().trim()));
            indiv.setAffectedStatus(Integer.valueOf(tokenizer.nextToken().trim()));
            indiv.setLabelInChip(indiv.getFamilyID() + "@*@" + indiv.getIndividualID());

            /*
             * assume there is only one tait in the pedigree file
            for (int i = 5; i < gtyStaringCol; i++) {
            tmpBuffer.delete(0, tmpBuffer.length());
            tmpBuffer.append(tokenizer.nextToken().trim());
            indiv.addTrait(tmpBuffer.toString());
            }
             */
            indivLabel = indiv.getLabelInChip();
            if (duplicatedIndiv.contains(indivLabel)) {
                String unexpectInfo = "Duplicated Individuals  in " +
                        pedFileName + " for PedID " + indiv.getFamilyID() + " IndivID " + indiv.getIndividualID();
              //  GlobalVariables.addInforLog(unexpectInfo);
                throw new Exception(unexpectInfo);
            }
            StatusGtySet sGty = new StatusGtySet();
            sGty.paternalChrom = new BitVector(totalSNPSize);
            sGty.maternalChrom = new BitVector(totalSNPSize);
            sGty.existence = new BitVector(totalSNPSize);
            duplicatedIndiv.add(indivLabel);

            int j = -1;
            int halfJ = 0;
            sigGty1 = null;
            snpOrder = 0;

            for (int i = 0; i < totalSNPSize; i++) {
                AnnotSNP snp = snpList.get(i);
                do {
                    j++;
                    halfJ = j / 2;
                    if (snp.probeSetID == halfJ) {
                        tmpBuffer.delete(0, tmpBuffer.length());
                        tmpBuffer.append(tokenizer.nextToken().trim());
                        sigGty1 = tmpBuffer.toString();
                    } else {
                        tokenizer.nextToken();
                    }
                } while (snp.probeSetID > halfJ);

                snp.order = snpOrder;
                aAllele = snp.getAAllele();
                bAllele = snp.getBAllele();
                tmpBuffer.delete(0, tmpBuffer.length());
                tmpBuffer.append(tokenizer.nextToken().trim());
                sigGty2 = tmpBuffer.toString();

                if (needAlleleNames[i]) {
                    //decide the allele names
                    if (aAllele != MISSING_ALLELE_NAME && bAllele == MISSING_ALLELE_NAME) {
                        //decide the allele names
                        if (sigGty1.charAt(0) != aAllele && sigGty1.charAt(0) != DEFAULT_MISSING_GTY_NAME) {
                            bAllele = sigGty1.charAt(0);
                            snp.setBAllele(bAllele);
                            needAlleleNames[i] = false;
                        } else if (sigGty2.charAt(0) != aAllele && sigGty2.charAt(0) != DEFAULT_MISSING_GTY_NAME) {
                            bAllele = sigGty2.charAt(0);
                            snp.setBAllele(bAllele);
                            needAlleleNames[i] = false;
                        }
                    } else if (aAllele == MISSING_ALLELE_NAME && bAllele == MISSING_ALLELE_NAME) {
                        //decide the allele names
                        if (sigGty1.charAt(0) != DEFAULT_MISSING_GTY_NAME) {
                            aAllele = sigGty1.charAt(0);
                            snp.setAAllele(aAllele);
                            if (sigGty2.charAt(0) != DEFAULT_MISSING_GTY_NAME && sigGty2.charAt(0) != aAllele) {
                                bAllele = sigGty2.charAt(0);
                                snp.setBAllele(bAllele);
                                needAlleleNames[i] = false;
                            }
                        } else if (sigGty2.charAt(0) != DEFAULT_MISSING_GTY_NAME) {
                            aAllele = sigGty2.charAt(0);
                            snp.setAAllele(aAllele);
                        }
                    } else if (aAllele == MISSING_ALLELE_NAME && bAllele != MISSING_ALLELE_NAME) {
                        //decide the allele names
                        if (sigGty1.charAt(0) != bAllele && sigGty1.charAt(0) != DEFAULT_MISSING_GTY_NAME) {
                            aAllele = sigGty1.charAt(0);
                            snp.setAAllele(aAllele);
                            needAlleleNames[i] = false;
                        } else if (sigGty2.charAt(0) != bAllele && sigGty2.charAt(0) != DEFAULT_MISSING_GTY_NAME) {
                            aAllele = sigGty2.charAt(0);
                            snp.setAAllele(aAllele);
                            needAlleleNames[i] = false;
                        }
                    }
                }

                // anyhow the  aAllele and bAllele must be known
                //re-code: A->0; B->1;
                if (sigGty1.charAt(0) == bAllele) {
                    sGty.existence.putQuick(snpOrder, true);
                    sGty.paternalChrom.putQuick(snpOrder, true);
                } else if (sigGty1.charAt(0) == aAllele) {
                    sGty.existence.putQuick(snpOrder, true);
                } else if (sigGty1.charAt(0) != DEFAULT_MISSING_GTY_NAME) {
                    //throw new Exception("Unexpected genotype in " + filePath + " for " + probID);
                    String unexpectInfo = "Unexpected genotype " + sigGty1 + " in " + pedFileName + " for " + snp.getRSID() + " of PedID " + indiv.getFamilyID() + " IndivID " + indiv.getIndividualID();
                    unexpectedSNPSet.add(snp.getRSID());
                    //((IGG3View) GlobalVariables.currentApplication.getMainView()).setBriefRunningInfor(unexpectInfo);
                   // GlobalVariables.addInforLog(unexpectInfo);
                    unexpectedGtyNum++;
                }

                j++;

                if (sigGty2.charAt(0) == bAllele) {
                    sGty.maternalChrom.putQuick(snpOrder, true);
                } else if (sigGty2.charAt(0) != aAllele && sigGty2.charAt(0) != DEFAULT_MISSING_GTY_NAME) {
                    //throw new Exception("Unexpected genotype in " + filePath + " for " + probID);
                    String unexpectInfo = "Unexpected genotype " + sigGty2 + " in " + pedFileName + " for " + snp.getRSID() + " of PedID " + indiv.getFamilyID() + " IndivID " + indiv.getIndividualID();
                    unexpectedSNPSet.add(snp.getRSID());
                    //((IGG3View) GlobalVariables.currentApplication.getMainView()).setBriefRunningInfor(unexpectInfo);
                   // GlobalVariables.addInforLog(unexpectInfo);
                    unexpectedGtyNum++;
                }
                snpOrder++;
            }
            //System.out.println(indiv.getLabelInChip());
            indiv.setMarkerGtySet(sGty);
            indList.add(indiv);
            tokenizer = null;
            line = null;
        }
        br.close();
        if (unexpectedGtyNum > 0) {
           // GlobalVariables.addInforLog(unexpectedSNPSet.size() + " SNPs have unexpected genotypes in " + pedFileName + ":\n" + unexpectedSNPSet.toString());
            unexpectedSNPSet.clear();
        }
        return false;
    }

    /**
     * retrieve SNP information from the Map files
     * @param FileName
     * @param indices
     * @param arry
     * @param genetRegions pair of physical regions
     * @throws java.lang.Exception
     * @return
     */
    protected int readSNPsinMapFileByPysic(PedFileSet fileSet, List<AnnotSNP> snpList,
            String chromosome, long[] physicalRegions) throws Exception {
        File mapFile = new File(fileSet.getMapFileName());
        BufferedReader br = new BufferedReader(new FileReader(mapFile), FILE_READER_BUFFER_SIZE);
        String line = null;

        String delmilit = ", \t";
        int chromIndex = 0;
        int rsIndex = 1;
        int posIndex = 3;
        int maxIndex = Math.max(chromIndex, posIndex);
        maxIndex = Math.max(maxIndex, rsIndex);
        boolean checkRegion = false;

        int j = 0;
        StringBuffer strRegion = new StringBuffer();
        int regionNum = -1;
        if (physicalRegions != null) {
            regionNum = physicalRegions.length / 2;
            for (j = 0; j < regionNum; j++) {
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
        }

        int counter = 0;
        int chroNum;
        if (chromosome.equalsIgnoreCase("X")) {
            chroNum = 23;
        } else if (chromosome.equalsIgnoreCase("Y")) {
            chroNum = 24;
        } else if (chromosome.equalsIgnoreCase("MT")) {
            chroNum = 25;
        } else if (chromosome.equalsIgnoreCase("XY")) {
            chroNum = 26;
        } else {
            chroNum = Integer.parseInt(chromosome);
        }

        if (checkRegion) {
            regionNum = physicalRegions.length / 2;
        }
        boolean invalid = false;
        long pos = 0;

        //skip the head line
        // br.readLine();

        boolean diffChrom = false;
        String tmpStr;
        String rsID = null;

        int index;
        Map<Long, Integer> physicalPos = new HashMap<Long, Integer>(10000, 0.9f);
        double decodeGenet = -9;
        char strand = MISSING_STRAND_NAME;
        int lineCounter = -1;
        StringBuffer tmpBuffer = new StringBuffer();
        try {
            while ((line = br.readLine()) != null) {
                if (line.trim().length() == 0) {
                    continue;
                }
                lineCounter++;
                StringTokenizer tokenizer = new StringTokenizer(line, delmilit);
                index = 0;
                diffChrom = false;
                invalid = false;
                while (tokenizer.hasMoreTokens()) {
                    tmpBuffer.delete(0, tmpBuffer.length());
                    tmpBuffer.append(tokenizer.nextToken().trim());
                    tmpStr = tmpBuffer.toString();
                    //filter other chromosomes, which is ordered in the annotation file
                    if (index == chromIndex) {
                        if (tmpStr.compareTo(chromosome) != 0) {
                            diffChrom = true;
                            break;
                        }
                    } else if (index == rsIndex) {
                        rsID = tmpStr;
                    } else if (index == posIndex) {
                        //filter items without physical postion
                        if (tmpStr.startsWith("-")) {
                            invalid = true;
                            break;
                        }

                        pos = Long.parseLong(tmpStr);
                        if (physicalPos.get(pos) != null) {
                            invalid = true;
                            break;
                        }

                        //filter physical region
                        if (checkRegion) {
                            for (int i = 0; i < regionNum; i++) {
                                if (pos < physicalRegions[i * 2] || pos > physicalRegions[i * 2 + 1]) {
                                    invalid = true;
                                    break;
                                }
                            }
                        }
                    }


                    if (index == maxIndex) {
                        break;
                    }
                    index++;
                }

                if (diffChrom) {
                    continue;
                }
                if (invalid) {
                    continue;
                }

                //split is much slower than Tokenizer
                // linecells = line.split(delmilit);
                /*
                if (!LocalString.isNumeric(linecells[3])) {
                //System.markerOut.println(line);
                //continue;
                linecells[3] = "-1";
                }*/

                if (!rsID.startsWith("rs")) {
                    rsID = "Chr" + chromosome + "Pos" + pos;
                    //continue;
                }
                physicalPos.put(pos, 0);

                // System.out.println(linecells[1]); 
                AnnotSNP aSNP = new AnnotSNP(rsID, chroNum, (int)pos, decodeGenet);
                aSNP.probeSetID = lineCounter;

                aSNP.setStrand(strand);
                //MISSING_STRAND_NAME indicates unknown strand information
                strand = MISSING_STRAND_NAME;
                snpList.add(aSNP);
            }


            StringBuffer runningInfo = new StringBuffer();
            runningInfo.append("The number of SNPs on chromosome ");
            runningInfo.append(chromosome);
            runningInfo.append(" ");
            runningInfo.append(strRegion);
            runningInfo.append(" in map file ");
            runningInfo.append(mapFile.getName());
            runningInfo.append(" is ");
            runningInfo.append(snpList.size());
            runningInfo.append(".");

           // GlobalVariables.addInforLog(runningInfo.toString());
            return (lineCounter + 1);
        } finally {
            br.close();
        }
    }

    /**
     * retrieve SNP information from the Map files
     * @param FileName
     * @param indices
     * @param arry
     * @param genetRegions pair of physical regions
     * @throws java.lang.Exception
     * @return
     */
    protected int readSNPsinMapFile(PedFileSet fileSet, List<AnnotSNP> snpList) throws Exception {
        File mapFile = new File(fileSet.getMapFileName());
        BufferedReader br = new BufferedReader(new FileReader(mapFile), FILE_READER_BUFFER_SIZE);
        String line = null;

        String delmilit = ", \t";
        int chromIndex = 0;
        int rsIndex = 1;
        int posIndex = 3;
        int maxIndex = Math.max(chromIndex, posIndex);
        maxIndex = Math.max(maxIndex, rsIndex);

        int j = 0;
        StringBuffer strRegion = new StringBuffer();

        boolean invalid = false;
        long pos = 0;

        //skip the head line
        // br.readLine();

        String tmpStr;
        String rsID = null;
        int chromosome = -1;

        int index;
        Map<Long, Integer> physicalPos = new HashMap<Long, Integer>(10000, 0.9f);
        double decodeGenet = -9;
        char strand = MISSING_STRAND_NAME;
        int lineCounter = -1;
        StringBuilder tmpBuffer = new StringBuilder();
        try {
            while ((line = br.readLine()) != null) {
                if (line.trim().length() == 0) {
                    continue;
                }
                lineCounter++;
                StringTokenizer tokenizer = new StringTokenizer(line, delmilit);
                index = 0;

                invalid = false;
                while (tokenizer.hasMoreTokens()) {
                    tmpBuffer.delete(0, tmpBuffer.length());
                    tmpBuffer.append(tokenizer.nextToken().trim());
                    tmpStr = tmpBuffer.toString();
                    //filter other chromosomes, which is ordered in the annotation file
                    if (index == chromIndex) {
                        if (tmpStr.equalsIgnoreCase("X")) {
                            chromosome = 23;
                        } else if (tmpStr.equalsIgnoreCase("Y")) {
                            chromosome = 24;
                        } else if (tmpStr.equalsIgnoreCase("---")) {
                            chromosome = 0;
                        } else {
                            chromosome = Integer.parseInt(tmpStr);
                        }
                    } else if (index == rsIndex) {
                        rsID = tmpStr;
                    } else if (index == posIndex) {
                        //filter items without physical postion
                        if (tmpStr.startsWith("-")) {
                            invalid = true;
                            break;
                        }
                        pos = Long.parseLong(tmpStr);
                        if (physicalPos.get(pos) != null) {
                            invalid = true;
                            break;
                        }
                    }

                    if (index == maxIndex) {
                        break;
                    }
                    index++;
                }

                if (invalid) {
                    continue;
                }


                if (!rsID.startsWith("rs")) {
                    rsID = "Chr" + chromosome + "Pos" + pos;
                    //continue;
                }
                physicalPos.put(pos, 0);

                // System.out.println(linecells[1]); 
                AnnotSNP aSNP = new AnnotSNP(rsID, chromosome, (int) pos, decodeGenet);
                aSNP.probeSetID = lineCounter;

                aSNP.setStrand(strand);
                //MISSING_STRAND_NAME indicates unknown strand information
                strand = MISSING_STRAND_NAME;
                snpList.add(aSNP);
            }
 
            StringBuilder runningInfo = new StringBuilder();
            runningInfo.append("The number of SNPs on chromosome ");
            runningInfo.append(chromosome);
            runningInfo.append(" ");
            runningInfo.append(strRegion);
            runningInfo.append(" in map file ");
            runningInfo.append(mapFile.getName());
            runningInfo.append(" is ");
            runningInfo.append(snpList.size());
            runningInfo.append(".");

           // GlobalVariables.addInforLog(runningInfo.toString());
            return (lineCounter + 1);
        } finally {
            br.close();
        }
    }

    public void readObjects(PedFileSet fileSet, List<Individual> indList, List<AnnotSNP> snpList) throws Exception {
        //default names of objects
        File pedFile = new File(fileSet.getPedigreeFileName() + ".obj");
        File mapFile = new File(fileSet.getMapFileName() + ".obj");
        if (!pedFile.exists()) {
            readSNPsinMapFile(fileSet, snpList);
            readGenotypeinPedigreeFile(fileSet, indList, snpList);
            FileOutputStream objFOut = new FileOutputStream(pedFile.getCanonicalPath());
            BufferedOutputStream objOBfs = new BufferedOutputStream(objFOut);
            ObjectOutputStream localObjOut = new ObjectOutputStream(objOBfs);
            localObjOut.writeObject(indList);
            localObjOut.flush();
            localObjOut.close();
            objOBfs.flush();
            objOBfs.close();
            objFOut.close();

            objFOut = new FileOutputStream(mapFile.getCanonicalPath());
            objOBfs = new BufferedOutputStream(objFOut);
            localObjOut = new ObjectOutputStream(objOBfs);
            localObjOut.writeObject(snpList);
            localObjOut.flush();
            localObjOut.close();
            objOBfs.flush();
            objOBfs.close();
            objFOut.close();
        } else {
            FileInputStream objFIn = new FileInputStream(pedFile);
            BufferedInputStream objIBfs = new BufferedInputStream(objFIn);
            ObjectInputStream localObjIn = new ObjectInputStream(objIBfs);
            indList.addAll((List<Individual>) localObjIn.readObject());
            localObjIn.close();
            objIBfs.close();
            objFIn.close();
            System.out.println(indList.size() + " individuals loaded");
            objFIn = new FileInputStream(mapFile);
            objIBfs = new BufferedInputStream(objFIn);
            localObjIn = new ObjectInputStream(objIBfs);
            snpList.addAll((List<AnnotSNP>) localObjIn.readObject());
            localObjIn.close();
            objIBfs.close();
            objFIn.close();
            System.out.println(snpList.size() + " markers loaded");
        }
    }
}
