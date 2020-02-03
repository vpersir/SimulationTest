/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.genetsimulator.entity;

import cern.colt.bitvector.BitVector;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.TreeSet;

/**
 *
 * @author mxli
 */
public class PlinkDataset {

    public static char MISSING_ALLELE_NAME = 'X';
    public static char MISSING_STRAND_NAME = '0';
    public static char DEFAULT_MISSING_GTY_NAME = '0';
    protected String pedigreeFileName;
    protected String mapFileName;
    //if plinkBinaryFileName is not null, it must be a plink binary file
    protected String plinkBinaryFileName;

    public PlinkDataset() {
    }

    public PlinkDataset(String pedigreeFileName, String mapFileName, String plinkBinaryFileName) {
        this.pedigreeFileName = pedigreeFileName;
        this.mapFileName = mapFileName;
        this.plinkBinaryFileName = plinkBinaryFileName;
    }

    public boolean avaibleFiles() {
        File file = new File(pedigreeFileName);
        if (!file.exists()) {
            return false;
        }
        file = new File(mapFileName);
        if (!file.exists()) {
            return false;
        }
        file = new File(plinkBinaryFileName);
        if (!file.exists()) {
            return false;
        }
        return true;
    }

    /**
     * This function can generate genetic map for merlin
     *
     */
    public boolean writeMapFile(String outFileString, List<AnnotSNP> snpInfo) throws Exception {
        DecimalFormat df = new DecimalFormat();
        df.setMaximumFractionDigits(6);
        int snpInfoSize = 0;
        snpInfoSize = snpInfo.size();
        AnnotSNP snp = null;
        AnnotSNP tmpSnpInfoRow = null;
        double lastPosotion = -1.0;
        double currentPostion = -1.0;
        double incTime = 1.0;
        BufferedWriter bw = new BufferedWriter(new FileWriter(outFileString));
        for (int k = 0; k < snpInfoSize; k++) {
            snp = snpInfo.get(k);
            bw.write(String.valueOf(snp.getChromosomeNum()));
            bw.write(" ");
            bw.write(snp.getRSID());
            bw.write(" ");
            if (snp.getDecodeGeneticMap() < 0) {
                bw.write("0.0");
            } else {
                currentPostion = snp.getDecodeGeneticMap();
                //make sure the genetic map of the markers in ascending order
                if (lastPosotion > currentPostion) {
                    int s = k + 2;
                    while (true) {
                        if (s >= snpInfoSize) {
                            currentPostion = lastPosotion + 0.00001 * incTime;
                            break;
                        }

                        tmpSnpInfoRow = snpInfo.get(s);
                        if (tmpSnpInfoRow.getDecodeGeneticMap() < 0) {
                            s++;
                            continue;
                        }
                        currentPostion = tmpSnpInfoRow.getDecodeGeneticMap();
                        if (lastPosotion < currentPostion) {
                            currentPostion = lastPosotion + (currentPostion - lastPosotion) / (s - k);
                            break;
                        } else {
                            s++;
                        }
                    }
                }
                //make sure the genetic map of the markers are unique
                if (lastPosotion == currentPostion) {
                    bw.write(df.format(lastPosotion + 0.00001 * incTime));
                    //System.bedOut.println(snp[1]+" : "+incTime);
                    incTime++;
                } else {
                    bw.write(df.format(currentPostion));
                }
            }
            bw.write(" ");
            bw.write(String.valueOf(snp.getPhysicalPosition()));
            bw.write("\n");
            lastPosotion = currentPostion;
        }
        bw.close();
        return true;
    }
//, List<AnnotSNP> snpInfo

    public void generalLinkPedGenotypeNumberAllele(String outFileString, List<Individual> indivList, List<AnnotSNP> snpInfo) throws Exception {
        String missingGty = "0 0";
        int snpNum = snpInfo.size();
        BufferedWriter bw = new BufferedWriter(new FileWriter(outFileString));
        int[] orders = new int[snpNum];

        for (int i = 0; i < snpNum; i++) {
            orders[i] = snpInfo.get(i).order;
            // orders[i] = i;
        }

        List<String> avaiTrais = new ArrayList<String>();
        // avaiTrais.add("Disease");
        int traitSize = avaiTrais.size();
        traitSize--;
        int indiviSize = indivList.size();
        for (int k = 0; k < indiviSize; k++) {
            Individual indiv = indivList.get(k);
            bw.write(indiv.getFamilyID());

            bw.write(" ");
            bw.write(indiv.getIndividualID());
            bw.write(" ");
            bw.write(indiv.getDadID());
            bw.write(" ");
            bw.write(indiv.getMomID());
            bw.write(" ");
            bw.write(String.valueOf(indiv.getGender()));
            bw.write(" ");
            bw.write(String.valueOf(indiv.getAffectedStatus()));
            if (traitSize >= 0) {
                bw.write(" ");
                //write trait value
                for (int j = 0; j < traitSize; j++) {
                    bw.write(indiv.getTraits().get(j));
                    bw.write(" ");
                }
                bw.write(indiv.getTraits().get(traitSize));
            }
            StatusGtySet indivGty = indivList.get(k).markerGtySet;
            if (indivGty != null) {
                for (int j = 0; j < snpNum; j++) {
                    bw.write(" ");
                    if (indivGty.existence.getQuick(orders[j])) {
                        if (indivGty.paternalChrom.getQuick(orders[j])) {
                            bw.write("2");
                        } else {
                            bw.write("1");
                        }
                        bw.write(" ");
                        if (indivGty.maternalChrom.getQuick(orders[j])) {
                            bw.write("2");
                        } else {
                            bw.write("1");
                        }
                    } else {
                        bw.write(missingGty);
                    }
                }
            } else {
                for (int j = 0; j < snpNum; j++) {
                    bw.write(" ");
                    bw.write(missingGty);
                }
            }

            bw.write("\n");
        }
        bw.close();
    }

    public void generalLinkPhenotype(String outFileString, List<Individual> indivList) throws Exception {
        BufferedWriter bw = new BufferedWriter(new FileWriter(outFileString));
        int traitSize = indivList.get(0).getMainTrait().length;
        traitSize--;
        int indiviSize = indivList.size();
        for (int k = 0; k < indiviSize; k++) {
            Individual indiv = indivList.get(k);
            bw.write(indiv.getFamilyID());

            bw.write(" ");
            bw.write(indiv.getIndividualID());
            if (traitSize >= 0) {
                bw.write(" ");
                //write trait value
                for (int j = 0; j < traitSize; j++) {
                    bw.write(String.valueOf(indiv.getMainTrait()[j]));
                    bw.write(" ");
                }
                bw.write(String.valueOf(indiv.getMainTrait()[traitSize]));
            }
            bw.write("\n");
        }
        bw.close();
    }

    public void generalLinkDiseaseTrait(String outFileString, List<Individual> indivList) throws Exception {
        BufferedWriter bw = new BufferedWriter(new FileWriter(outFileString));

        int indiviSize = indivList.size();
        for (int k = 0; k < indiviSize; k++) {
            Individual indiv = indivList.get(k);
            bw.write(indiv.getFamilyID());

            bw.write(" ");
            bw.write(indiv.getIndividualID());

            bw.write(" ");

            bw.write(String.valueOf(indiv.getMainTrait()[0]));
            // bw.write(String.valueOf(indiv.getAffectedStatus()));

            bw.write("\n");
        }
        bw.close();
    }

    public void generalLinkPedGenotypeNumberAllele(String outFileString, List<Individual> indivList) throws Exception {
        String missingGty = "0 0";
        int snpNum = indivList.get(0).markerGtySet.paternalChrom.size();
        BufferedWriter bw = new BufferedWriter(new FileWriter(outFileString));
        int[] orders = new int[snpNum];

        for (int i = 0; i < snpNum; i++) {
            //orders[i] = snpInfo.get(i).order;
            orders[i] = i;
        }

        List<String> avaiTrais = new ArrayList<String>();
        // avaiTrais.add("Disease");
        int traitSize = avaiTrais.size();
        traitSize--;
        int indiviSize = indivList.size();
        for (int k = 0; k < indiviSize; k++) {
            Individual indiv = indivList.get(k);
            bw.write(indiv.getFamilyID());

            bw.write(" ");
            bw.write(indiv.getIndividualID());
            bw.write(" ");
            bw.write(indiv.getDadID());
            bw.write(" ");
            bw.write(indiv.getMomID());
            bw.write(" ");
            bw.write(String.valueOf(indiv.getGender()));
            bw.write(" ");
            bw.write(String.valueOf(indiv.getAffectedStatus()));
            if (traitSize >= 0) {
                bw.write(" ");
                //write trait value
                for (int j = 0; j < traitSize; j++) {
                    bw.write(indiv.getTraits().get(j));
                    bw.write(" ");
                }
                bw.write(indiv.getTraits().get(traitSize));
            }
            StatusGtySet indivGty = indivList.get(k).markerGtySet;
            if (indivGty != null) {
                for (int j = 0; j < snpNum; j++) {
                    bw.write(" ");
                    if (indivGty.existence.getQuick(orders[j])) {
                        if (indivGty.paternalChrom.getQuick(orders[j])) {
                            bw.write("2");
                        } else {
                            bw.write("1");
                        }
                        bw.write(" ");
                        if (indivGty.maternalChrom.getQuick(orders[j])) {
                            bw.write("2");
                        } else {
                            bw.write("1");
                        }
                    } else {
                        bw.write(missingGty);
                    }
                }
            } else {
                for (int j = 0; j < snpNum; j++) {
                    bw.write(" ");
                    bw.write(missingGty);
                }
            }

            bw.write("\n");
        }
        bw.close();
    }

    // public void calcualte
    public void readPlinkBinaryFormatPedigreeGenotype(List<Individual> indList, List<AnnotSNP> snpList) throws Exception {
        Map<String, StatusGtySet> indivGtyMap = new HashMap<String, StatusGtySet>();

        String pedFileName = pedigreeFileName;
        BufferedReader br = new BufferedReader(new FileReader(pedFileName));
        TreeSet<String> unexpectedSNPSet = new TreeSet<String>();
        String line;
        String delimiter = "\t\" \",/";
        int totalSNPSize = snpList.size();
        int gtyStaringCol = -1;

        String indivLabel;

        int unexpectedGtyNum = 0;
        List<String> indiviIDInPed = new ArrayList<String>();

        //assume only one trait int the file
        int traitNum = 1;
        //guess genotype starting column
        gtyStaringCol = 5 + traitNum;
        int snpOrder = 0;

        StringBuilder tmpBuffer = new StringBuilder();
        //long t3 = System.currentTimeMillis();
        while ((line = br.readLine()) != null) {

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
            tmpBuffer.delete(0, tmpBuffer.length());
            tmpBuffer.append(tokenizer.nextToken().trim());
            indiv.setGender(Integer.valueOf(tmpBuffer.toString()));
            indiv.setLabelInChip(indiv.getFamilyID() + "@*@" + indiv.getIndividualID());
            for (int i = 5; i < gtyStaringCol; i++) {
                tmpBuffer.delete(0, tmpBuffer.length());
                tmpBuffer.append(tokenizer.nextToken().trim());
                indiv.addTrait(tmpBuffer.toString());
            }
            indiv.mainTrait = new double[2];
            //System.out.println(indiv.getLabelInChip());
            indList.add(indiv);

            indivLabel = indiv.getLabelInChip();
            StatusGtySet sGty = indivGtyMap.get(indivLabel);
            if (sGty == null) {
                sGty = new StatusGtySet();
                //////////////////////////////
                // Allocate space for SNPs
                sGty.paternalChrom = new BitVector(totalSNPSize);
                sGty.maternalChrom = new BitVector(totalSNPSize);
                sGty.existence = new BitVector(totalSNPSize);
                indivGtyMap.put(indivLabel, sGty);
                indiviIDInPed.add(indivLabel);
            } else {
                String unexpectInfo = "Duplicated Individuals  in " + pedFileName + " for PedID " + indiv.getFamilyID() + " IndivID " + indiv.getIndividualID();
                System.out.println(unexpectInfo);
            }
            tokenizer = null;
            line = null;
        }
        br.close();

        //start to read binary genotypes
//openPlinkPedFormat
        DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(plinkBinaryFileName)));
        boolean bfile_SNP_major = false;

        StringBuilder logInfor = new StringBuilder();
        byte bt = in.readByte();

        boolean v1_bfile = true;
        boolean[] b = new boolean[8];
        for (int i = 7; i >= 0; i--) {
            if (((1 << i) & bt) != 0) {
                b[i] = true;
            } else {
                b[i] = false;
            }
        }

        // printBitSet(b);
        // If v1.00 file format
        // Magic numbers for .bed file: 00110110 11011000 = v1.00 bed file
        if ((b[2] && b[3] && b[5] && b[6]) && !(b[0] || b[1] || b[4] || b[7])) {
            // Next number
            bt = in.readByte();
            for (int i = 7; i >= 0; i--) {
                if (((1 << i) & bt) != 0) {
                    b[i] = true;
                } else {
                    b[i] = false;
                }
            }
            // printBitSet(b);
            if ((b[0] && b[1] && b[3] && b[4]) && !(b[2] || b[5] || b[6] || b[7])) {
                // Read SNP/Ind major coding
                bt = in.readByte();
                for (int i = 7; i >= 0; i--) {
                    if (((1 << i) & bt) != 0) {
                        b[i] = true;
                    } else {
                        b[i] = false;
                    }
                }
                if (b[0]) {
                    bfile_SNP_major = true;
                } else {
                    bfile_SNP_major = false;
                }

                if (bfile_SNP_major) {
                    logInfor.append("Detected that binary PED file is v1.00 SNP-major mode\n");
                } else {
                    logInfor.append("Detected that binary PED file is v1.00 individual-major mode\n");
                }

            } else {
                v1_bfile = false;
            }

        } else {
            v1_bfile = false;
        }

        // Reset file if < v1
        if (!bfile_SNP_major) {
            logInfor.append("Warning, old BED file <v1.00 : will try to recover...\n");
            logInfor.append("  bs.get(t you should --make-bs.get(d from PED )\n");
            in.close();
            in = new DataInputStream(new BufferedInputStream(new FileInputStream(plinkBinaryFileName)));
            bt = in.readByte();
            for (int i = 7; i >= 0; i--) {
                if (((1 << i) & bt) != 0) {
                    b[i] = true;
                } else {
                    b[i] = false;
                }
            }
        }

        // If 0.99 file format
        if ((!v1_bfile) && (b[1] || b[2] || b[3] || b[4] || b[5] || b[6] || b[7])) {
            logInfor.append("\n *** Possibs.get(e probs.get(em: guessing that BED is < v0.99      *** \n");
            logInfor.append(" *** High chance of data corruption, spurious results    *** \n");
            logInfor.append(" *** Unles you are _sure_ this really is an old BED file *** \n");
            logInfor.append(" *** you should recreate PED -> BED                      *** \n\n");

            bfile_SNP_major = false;
            in.close();
            in = new DataInputStream(new BufferedInputStream(new FileInputStream(plinkBinaryFileName)));
        } else if (!v1_bfile) {
            if (b[0]) {
                bfile_SNP_major = true;
            } else {
                bfile_SNP_major = false;
            }

            logInfor.append("Binary PED file is v0.99\n");
            if (bfile_SNP_major) {
                logInfor.append("Detected that binary PED file is in SNP-major mode\n");
            } else {
                logInfor.append("Detected that binary PED file is in individual-major mode\n");
            }
        }
        System.out.println(logInfor);

//http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml
        String info = ("Reading genotype bit-file from [ " + plinkBinaryFileName + " ] \n");

        System.out.println(info);
        int indiviNum = indiviIDInPed.size();

        ///////////////////////////
        // SNP-major mode
        if (bfile_SNP_major) {
            snpOrder = 0;
            int byteNumPerSNP = indiviNum / 4;
            if (indiviNum % 4 > 0) {
                byteNumPerSNP++;
            }
            int currentBytePoint = 0;
            int shouldBytePoint = 0;
            for (int i = 0; i < totalSNPSize; i++) {
                AnnotSNP snp = snpList.get(i);
                shouldBytePoint = snp.probeSetID * byteNumPerSNP;
                if (currentBytePoint < shouldBytePoint) {
                    in.skipBytes(shouldBytePoint - currentBytePoint);
                    currentBytePoint = shouldBytePoint;
                }
                snp.order = snpOrder;
                int k = 0;
                while (k < indiviNum) {
                    bt = in.readByte();
                    currentBytePoint++;
                    for (int s = 7; s >= 0; s--) {
                        if (((1 << s) & bt) != 0) {
                            b[s] = true;
                        } else {
                            b[s] = false;
                        }
                    }

                    // printBitSet(b);
                    for (int s = 0; s < 8; s += 2) {
                        StatusGtySet sGty = indivGtyMap.get(indiviIDInPed.get(k));
                        if (!b[s] && b[s + 1]) {
                            sGty.existence.putQuick(i, true);
                            sGty.paternalChrom.putQuick(i, true);
                        } else if (!b[s] && !b[s + 1]) {
                            sGty.existence.putQuick(i, true);
                        } else if (b[s] && b[s + 1]) {
                            sGty.existence.putQuick(i, true);
                            sGty.maternalChrom.putQuick(i, true);
                            sGty.paternalChrom.putQuick(i, true);
                        }
                        k++;
                        if (k == indiviNum) {
                            break;
                        }
                    }
                }
                snpOrder++;
            }
        } else {
            int j = -1;
            int currentSNPNum = 0;

            snpOrder = 0;
            int s = 0;
            //This code has not been tested yet because the latest version generate SNP-major model
            ///////////////////////////
            // individual-major mode
            for (int k = 0; k < indiviNum; k++) {
                currentSNPNum = 0;
                StatusGtySet sGty = indivGtyMap.get(indiviIDInPed.get(k));
                for (int i = 0; i < totalSNPSize; i++) {
                    AnnotSNP snp = snpList.get(i);
                    while (currentSNPNum < snp.probeSetID) {
                        bt = in.readByte();
                        for (int t = 7; t >= 0; t--) {
                            if (((1 << t) & bt) != 0) {
                                b[t] = true;
                            } else {
                                b[t] = false;
                            }
                        }
                        currentSNPNum += 4;
                    }
                    s = (snp.probeSetID % 4) * 2;
                    if (!b[s] && b[s + 1]) {
                        sGty.existence.putQuick(i, true);
                        sGty.paternalChrom.putQuick(i, true);
                    } else if (!b[s] && !b[s + 1]) {
                        sGty.existence.putQuick(i, true);
                    } else if (b[s] && b[s + 1]) {
                        sGty.existence.putQuick(i, true);
                        sGty.maternalChrom.putQuick(i, true);
                        sGty.paternalChrom.putQuick(i, true);
                    }
                    snp.order = snpOrder;
                    snpOrder++;
                }
            }
        }

        in.close();
        if (unexpectedGtyNum > 0) {
            System.out.println((unexpectedGtyNum / 2) + " Unexpected genotypes in " + pedFileName + " for " + unexpectedSNPSet.size() + " SNPs (detailed in the log).");
            System.out.println(unexpectedSNPSet.size() + " SNPs have unexpected genotypes in " + pedFileName + ":\n" + unexpectedSNPSet.toString());

            System.out.println(unexpectedSNPSet);
            unexpectedSNPSet.clear();
        }

        int indivNum = indList.size();
        for (int i = 0; i < indivNum; i++) {
            Individual indiv = indList.get(i);
            indiv.markerGtySet = indivGtyMap.get(indiv.getLabelInChip());
        }

        indivGtyMap.clear();
    }

    // public void calcualte
    public void readPlinkBinaryFormatFamily(List<Individual> indList) throws Exception {
        Set<String> indivLabelSet = new HashSet<String>();
        String pedFileName = pedigreeFileName;
        BufferedReader br = new BufferedReader(new FileReader(pedFileName));
        String line;
        String delimiter = "\t\" \",/";

        int gtyStaringCol = -1;
        String indivLabel;
        //assume only one trait int the file
        int traitNum = 1;
        //guess genotype starting column
        gtyStaringCol = 5 + traitNum;

        StringBuilder tmpBuffer = new StringBuilder();
        //long t3 = System.currentTimeMillis();
        while ((line = br.readLine()) != null) {

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
            tmpBuffer.delete(0, tmpBuffer.length());
            tmpBuffer.append(tokenizer.nextToken().trim());
            indiv.setGender(Integer.valueOf(tmpBuffer.toString()));
            indiv.setLabelInChip(indiv.getFamilyID() + "@*@" + indiv.getIndividualID());
            for (int i = 5; i < gtyStaringCol; i++) {
                tmpBuffer.delete(0, tmpBuffer.length());
                tmpBuffer.append(tokenizer.nextToken().trim());
                indiv.addTrait(tmpBuffer.toString());
            }
            //System.out.println(indiv.getLabelInChip());
            indList.add(indiv);

            indivLabel = indiv.getLabelInChip();
            if (indivLabelSet.contains(indivLabel)) {
                String unexpectInfo = "Duplicated Individuals  in " + pedFileName + " for PedID " + indiv.getFamilyID() + " IndivID " + indiv.getIndividualID();
                System.out.println(unexpectInfo);
            } else {
                indivLabelSet.add(indivLabel);
                StatusGtySet sGty = new StatusGtySet();
                //////////////////////////////0
                sGty.paternalChrom = new BitVector(0);
                sGty.maternalChrom = new BitVector(0);
                sGty.existence = new BitVector(0);
                indiv.markerGtySet = sGty;
            }

            tokenizer = null;
            line = null;
        }
        br.close();
    }

    // public void calcualte
    public void readPlinkBinaryFormatGenotype(List<Individual> indList, List<AnnotSNP> snpList) throws Exception {
        TreeSet<String> unexpectedSNPSet = new TreeSet<String>();
        int totalSNPSize = snpList.size();
        int unexpectedGtyNum = 0;
        int snpOrder = 0;

        //start to read binary genotypes
//openPlinkPedFormat
        DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(plinkBinaryFileName)));
        boolean bfile_SNP_major = false;

        StringBuilder logInfor = new StringBuilder();
        byte bt = in.readByte();

        boolean v1_bfile = true;
        boolean[] b = new boolean[8];
        for (int i = 7; i >= 0; i--) {
            if (((1 << i) & bt) != 0) {
                b[i] = true;
            } else {
                b[i] = false;
            }
        }

        // printBitSet(b);
        // If v1.00 file format
        // Magic numbers for .bed file: 00110110 11011000 = v1.00 bed file
        if ((b[2] && b[3] && b[5] && b[6]) && !(b[0] || b[1] || b[4] || b[7])) {
            // Next number
            bt = in.readByte();
            for (int i = 7; i >= 0; i--) {
                if (((1 << i) & bt) != 0) {
                    b[i] = true;
                } else {
                    b[i] = false;
                }
            }
            // printBitSet(b);
            if ((b[0] && b[1] && b[3] && b[4]) && !(b[2] || b[5] || b[6] || b[7])) {
                // Read SNP/Ind major coding
                bt = in.readByte();
                for (int i = 7; i >= 0; i--) {
                    if (((1 << i) & bt) != 0) {
                        b[i] = true;
                    } else {
                        b[i] = false;
                    }
                }
                if (b[0]) {
                    bfile_SNP_major = true;
                } else {
                    bfile_SNP_major = false;
                }

                if (bfile_SNP_major) {
                    logInfor.append("Detected that binary PED file is v1.00 SNP-major mode\n");
                } else {
                    logInfor.append("Detected that binary PED file is v1.00 individual-major mode\n");
                }

            } else {
                v1_bfile = false;
            }

        } else {
            v1_bfile = false;
        }

        // Reset file if < v1
        if (!bfile_SNP_major) {
            logInfor.append("Warning, old BED file <v1.00 : will try to recover...\n");
            logInfor.append("  bs.get(t you should --make-bs.get(d from PED )\n");
            in.close();
            in = new DataInputStream(new BufferedInputStream(new FileInputStream(plinkBinaryFileName)));
            bt = in.readByte();
            for (int i = 7; i >= 0; i--) {
                if (((1 << i) & bt) != 0) {
                    b[i] = true;
                } else {
                    b[i] = false;
                }
            }
        }

        // If 0.99 file format
        if ((!v1_bfile) && (b[1] || b[2] || b[3] || b[4] || b[5] || b[6] || b[7])) {
            logInfor.append("\n *** Possibs.get(e probs.get(em: guessing that BED is < v0.99      *** \n");
            logInfor.append(" *** High chance of data corruption, spurious results    *** \n");
            logInfor.append(" *** Unles you are _sure_ this really is an old BED file *** \n");
            logInfor.append(" *** you should recreate PED -> BED                      *** \n\n");

            bfile_SNP_major = false;
            in.close();
            in = new DataInputStream(new BufferedInputStream(new FileInputStream(plinkBinaryFileName)));
        } else if (!v1_bfile) {
            if (b[0]) {
                bfile_SNP_major = true;
            } else {
                bfile_SNP_major = false;
            }

            logInfor.append("Binary PED file is v0.99\n");
            if (bfile_SNP_major) {
                logInfor.append("Detected that binary PED file is in SNP-major mode\n");
            } else {
                logInfor.append("Detected that binary PED file is in individual-major mode\n");
            }
        }
        System.out.println(logInfor);

        int indivNum = indList.size();
        for (int i = 0; i < indivNum; i++) {
            Individual indiv = indList.get(i);
            indiv.markerGtySet.existence.setSize(totalSNPSize);
            indiv.markerGtySet.existence.replaceFromToWith(0, totalSNPSize - 1, false);
            indiv.markerGtySet.maternalChrom.setSize(totalSNPSize);
            indiv.markerGtySet.maternalChrom.replaceFromToWith(0, totalSNPSize - 1, false);
            indiv.markerGtySet.paternalChrom.setSize(totalSNPSize);
            indiv.markerGtySet.paternalChrom.replaceFromToWith(0, totalSNPSize - 1, false);
        }
//http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml
        String info = ("Reading genotype bit-file from [ " + plinkBinaryFileName + " ] \n");
        System.out.println(info);
        ///////////////////////////
        // SNP-major mode
        if (bfile_SNP_major) {
            snpOrder = 0;
            int byteNumPerSNP = indivNum / 4;
            if (indivNum % 4 > 0) {
                byteNumPerSNP++;
            }
            int currentBytePoint = 0;
            int shouldBytePoint = 0;
            for (int i = 0; i < totalSNPSize; i++) {
                AnnotSNP snp = snpList.get(i);
                shouldBytePoint = snp.probeSetID * byteNumPerSNP;
                if (currentBytePoint < shouldBytePoint) {
                    in.skipBytes(shouldBytePoint - currentBytePoint);
                    currentBytePoint = shouldBytePoint;
                }
                snp.order = snpOrder;
                int k = 0;
                while (k < indivNum) {
                    bt = in.readByte();
                    currentBytePoint++;
                    for (int s = 7; s >= 0; s--) {
                        if (((1 << s) & bt) != 0) {
                            b[s] = true;
                        } else {
                            b[s] = false;
                        }
                    }

                    // printBitSet(b);
                    for (int s = 0; s < 8; s += 2) {
                        StatusGtySet sGty = indList.get(k).markerGtySet;
                        if (!b[s] && b[s + 1]) {
                            sGty.existence.putQuick(i, true);
                            sGty.paternalChrom.putQuick(i, true);
                        } else if (!b[s] && !b[s + 1]) {
                            sGty.existence.putQuick(i, true);
                        } else if (b[s] && b[s + 1]) {
                            sGty.existence.putQuick(i, true);
                            sGty.maternalChrom.putQuick(i, true);
                            sGty.paternalChrom.putQuick(i, true);
                        }
                        k++;
                        if (k == indivNum) {
                            break;
                        }
                    }
                }
                snpOrder++;
            }
        } else {
            throw new Exception("Cannot handle individual-major mode");
        }

        in.close();
        if (unexpectedGtyNum > 0) {
            System.out.println((unexpectedGtyNum / 2) + " Unexpected genotypes in " + plinkBinaryFileName + " for " + unexpectedSNPSet.size() + " SNPs (detailed in the log).");
            System.out.println(unexpectedSNPSet.size() + " SNPs have unexpected genotypes in " + plinkBinaryFileName + ":\n" + unexpectedSNPSet.toString());
            System.out.println(unexpectedSNPSet);
            unexpectedSNPSet.clear();
        }
    }

    // public void calcualte
    public void readPlinkBinaryFormatGenotype(List<Individual> indList, List<AnnotSNP> snpList, int chromID) throws Exception {
        TreeSet<String> unexpectedSNPSet = new TreeSet<String>();
        int totalSNPSize = snpList.size();
        int unexpectedGtyNum = 0;
        int snpOrder = 0;

        //start to read binary genotypes
//openPlinkPedFormat
        DataInputStream in = new DataInputStream(new BufferedInputStream(new FileInputStream(plinkBinaryFileName)));
        boolean bfile_SNP_major = false;

        StringBuilder logInfor = new StringBuilder();
        byte bt = in.readByte();

        boolean v1_bfile = true;
        boolean[] b = new boolean[8];
        for (int i = 7; i >= 0; i--) {
            if (((1 << i) & bt) != 0) {
                b[i] = true;
            } else {
                b[i] = false;
            }
        }

        // printBitSet(b);
        // If v1.00 file format
        // Magic numbers for .bed file: 00110110 11011000 = v1.00 bed file
        if ((b[2] && b[3] && b[5] && b[6]) && !(b[0] || b[1] || b[4] || b[7])) {
            // Next number
            bt = in.readByte();
            for (int i = 7; i >= 0; i--) {
                if (((1 << i) & bt) != 0) {
                    b[i] = true;
                } else {
                    b[i] = false;
                }
            }
            // printBitSet(b);
            if ((b[0] && b[1] && b[3] && b[4]) && !(b[2] || b[5] || b[6] || b[7])) {
                // Read SNP/Ind major coding
                bt = in.readByte();
                for (int i = 7; i >= 0; i--) {
                    if (((1 << i) & bt) != 0) {
                        b[i] = true;
                    } else {
                        b[i] = false;
                    }
                }
                if (b[0]) {
                    bfile_SNP_major = true;
                } else {
                    bfile_SNP_major = false;
                }

                if (bfile_SNP_major) {
                    logInfor.append("Detected that binary PED file is v1.00 SNP-major mode\n");
                } else {
                    logInfor.append("Detected that binary PED file is v1.00 individual-major mode\n");
                }

            } else {
                v1_bfile = false;
            }

        } else {
            v1_bfile = false;
        }

        // Reset file if < v1
        if (!bfile_SNP_major) {
            logInfor.append("Warning, old BED file <v1.00 : will try to recover...\n");
            logInfor.append("  bs.get(t you should --make-bs.get(d from PED )\n");
            in.close();
            in = new DataInputStream(new BufferedInputStream(new FileInputStream(plinkBinaryFileName)));
            bt = in.readByte();
            for (int i = 7; i >= 0; i--) {
                if (((1 << i) & bt) != 0) {
                    b[i] = true;
                } else {
                    b[i] = false;
                }
            }
        }

        // If 0.99 file format
        if ((!v1_bfile) && (b[1] || b[2] || b[3] || b[4] || b[5] || b[6] || b[7])) {
            logInfor.append("\n *** Possibs.get(e probs.get(em: guessing that BED is < v0.99      *** \n");
            logInfor.append(" *** High chance of data corruption, spurious results    *** \n");
            logInfor.append(" *** Unles you are _sure_ this really is an old BED file *** \n");
            logInfor.append(" *** you should recreate PED -> BED                      *** \n\n");

            bfile_SNP_major = false;
            in.close();
            in = new DataInputStream(new BufferedInputStream(new FileInputStream(plinkBinaryFileName)));
        } else if (!v1_bfile) {
            if (b[0]) {
                bfile_SNP_major = true;
            } else {
                bfile_SNP_major = false;
            }

            logInfor.append("Binary PED file is v0.99\n");
            if (bfile_SNP_major) {
                logInfor.append("Detected that binary PED file is in SNP-major mode\n");
            } else {
                logInfor.append("Detected that binary PED file is in individual-major mode\n");
            }
        }
        System.out.println(logInfor);

        int indivNum = indList.size();
        for (int i = 0; i < indivNum; i++) {
            Individual indiv = indList.get(i);
            if (indiv.markerGtySets[chromID] == null) {
                indiv.markerGtySets[chromID] = new StatusGtySet(totalSNPSize);
            } else {
                indiv.markerGtySets[chromID].existence.setSize(totalSNPSize);
                indiv.markerGtySets[chromID].existence.replaceFromToWith(0, totalSNPSize - 1, false);
                indiv.markerGtySets[chromID].maternalChrom.setSize(totalSNPSize);
                indiv.markerGtySets[chromID].maternalChrom.replaceFromToWith(0, totalSNPSize - 1, false);
                indiv.markerGtySets[chromID].paternalChrom.setSize(totalSNPSize);
                indiv.markerGtySets[chromID].paternalChrom.replaceFromToWith(0, totalSNPSize - 1, false);
            }
        }
//http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml
        String info = ("Reading genotype bit-file from [ " + plinkBinaryFileName + " ] \n");
        System.out.println(info);
        ///////////////////////////
        // SNP-major mode
        if (bfile_SNP_major) {
            snpOrder = 0;
            int byteNumPerSNP = indivNum / 4;
            if (indivNum % 4 > 0) {
                byteNumPerSNP++;
            }
            int currentBytePoint = 0;
            int shouldBytePoint = 0;
            for (int i = 0; i < totalSNPSize; i++) {
                AnnotSNP snp = snpList.get(i);
                shouldBytePoint = snp.probeSetID * byteNumPerSNP;
                if (currentBytePoint < shouldBytePoint) {
                    in.skipBytes(shouldBytePoint - currentBytePoint);
                    currentBytePoint = shouldBytePoint;
                }
                snp.order = snpOrder;
                int k = 0;
                while (k < indivNum) {
                    bt = in.readByte();
                    currentBytePoint++;
                    for (int s = 7; s >= 0; s--) {
                        if (((1 << s) & bt) != 0) {
                            b[s] = true;
                        } else {
                            b[s] = false;
                        }
                    }

                    // printBitSet(b);
                    for (int s = 0; s < 8; s += 2) {
                        StatusGtySet sGty = indList.get(k).markerGtySets[chromID];
                        if (!b[s] && b[s + 1]) {
                            sGty.existence.putQuick(i, true);
                            sGty.paternalChrom.putQuick(i, true);
                        } else if (!b[s] && !b[s + 1]) {
                            sGty.existence.putQuick(i, true);
                        } else if (b[s] && b[s + 1]) {
                            sGty.existence.putQuick(i, true);
                            sGty.maternalChrom.putQuick(i, true);
                            sGty.paternalChrom.putQuick(i, true);
                        }
                        k++;
                        if (k == indivNum) {
                            break;
                        }
                    }
                }
                snpOrder++;
            }
        } else {
            throw new Exception("Cannot handle individual-major mode");
        }

        in.close();
        if (unexpectedGtyNum > 0) {
            System.out.println((unexpectedGtyNum / 2) + " Unexpected genotypes in " + plinkBinaryFileName + " for " + unexpectedSNPSet.size() + " SNPs (detailed in the log).");
            System.out.println(unexpectedSNPSet.size() + " SNPs have unexpected genotypes in " + plinkBinaryFileName + ":\n" + unexpectedSNPSet.toString());
            System.out.println(unexpectedSNPSet);
            unexpectedSNPSet.clear();
        }

    }

    /**
     * retrieve SNP information from the Map files
     *
     * @param FileName
     * @param indices
     * @param arry
     * @param genetRegions pair of physical regions
     * @throws java.lang.Exception
     * @return
     */
    public int readPlinkBinaryFormatMapByPysic(List<AnnotSNP> snpList, String chromosome, int[] physicalRegions) throws Exception {
        File mapFile = new File(mapFileName);
        if (!mapFile.exists()) {
            return -1;
        }
        BufferedReader br = new BufferedReader(new FileReader(mapFile));
        String line = null;

        String delmilit = ", \t";
        int chromIndex = 0;
        int posIndex = 3;
        int rsIndex = 1;
        int strandIndex = -1;
        int maxIndex = Math.max(chromIndex, posIndex);

        maxIndex = Math.max(maxIndex, rsIndex);
        maxIndex = Math.max(maxIndex, strandIndex);
        boolean checkRegion = false;
        int regionNum = physicalRegions.length / 2;

        int j = 0;
        StringBuilder strRegion = new StringBuilder();
        for (j = 0; j < regionNum; j++) {
            if (physicalRegions[j * 2] != -9 || physicalRegions[j * 2 + 1] != -9) {
                checkRegion = true;
                if (physicalRegions[j * 2 + 1] == -9) {
                    physicalRegions[1] = Integer.MAX_VALUE;
                }
                strRegion.append(" [");
                strRegion.append(physicalRegions[j * 2]);
                strRegion.append(",");
                strRegion.append(physicalRegions[j * 2 + 1]);
                strRegion.append("]bp");
            }
        }

        int lineCounter = -1;

        int chroNum;
        if (chromosome.equalsIgnoreCase("X")) {
            chroNum = 23;
        } else if (chromosome.equalsIgnoreCase("Y")) {
            chroNum = 24;
        } else if (chromosome.equalsIgnoreCase("MT")) {
            chroNum = 26;
        } else if (chromosome.equalsIgnoreCase("XY")) {
            chroNum = 25;
        } else {
            chroNum = Integer.parseInt(chromosome);
        }

        if (checkRegion) {
            regionNum = physicalRegions.length / 2;
        }
        boolean invalid = false;
        int pos = 0;

        boolean diffChrom = false;
        String tmpStr;
        String rsID = null;
        StringBuilder tmpBuffer = new StringBuilder();
        int index;
        Map<Integer, Integer> physicalPos = new HashMap<Integer, Integer>(10000, 0.9f);
        double decodeGenet = -9;
        char strand = MISSING_STRAND_NAME;
        int alleleAIndex = 4;
        int alleleBIndex = 5;
        maxIndex = Math.max(maxIndex, alleleAIndex);
        maxIndex = Math.max(maxIndex, alleleBIndex);
        int currentChromNum = 0;
        char alleleA = MISSING_ALLELE_NAME, alleleB = MISSING_ALLELE_NAME;
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
                alleleA = MISSING_ALLELE_NAME;
                alleleB = MISSING_ALLELE_NAME;
                while (tokenizer.hasMoreTokens()) {
                    tmpBuffer.delete(0, tmpBuffer.length());
                    tmpBuffer.append(tokenizer.nextToken().trim());
                    tmpStr = tmpBuffer.toString();
                    if (index == chromIndex) {
                        //plink binary format only have numbers for all chromsomes
                        if (tmpStr.equals("X")) {
                            currentChromNum = 23;
                        } else {
                            currentChromNum = Integer.parseInt(tmpStr);
                        }
                        if (currentChromNum != chroNum) {
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

                        pos = Integer.parseInt(tmpStr);
                        if (physicalPos.get(pos) != null) {
                            invalid = true;
                            break;
                        }

                        //filter physical region
                        if (checkRegion) {
                            invalid = true;
                            for (int i = 0; i < regionNum; i++) {
                                if (pos >= physicalRegions[i * 2] && pos <= physicalRegions[i * 2 + 1]) {
                                    invalid = false;
                                    break;
                                }
                            }
                        }
                    } else if (index == strandIndex) {
                        strand = tmpStr.charAt(0);
                    } else if (index == alleleAIndex) {
                        alleleA = tmpStr.charAt(0);
                    } else if (index == alleleBIndex) {
                        alleleB = tmpStr.charAt(0);
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
                AnnotSNP aSNP = new AnnotSNP(rsID, chroNum, pos, -1);
                aSNP.setAAllele(alleleA);
                aSNP.setBAllele(alleleB);

                //a global index in the map file
                aSNP.probeSetID = lineCounter;
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

            System.out.println(runningInfo.toString());
            return (snpList.size());
        } finally {
            br.close();
        }
    }

    //___________________________________________TEST______________________________________________________________
    public int readPlinkBinaryFormatMapByPysicTEST(List<AnnotSNP> snpList, String chromosome, int[] physicalRegions) throws Exception {
        File mapFile = new File(mapFileName);
        if (!mapFile.exists()) {
            return -1;
        }
        BufferedReader br = new BufferedReader(new FileReader(mapFile));
        String line = null;

        String delmilit = ", \t";
        int chromIndex = 0;
        int posIndex = 3;
        int rsIndex = 1;
        int strandIndex = -1;
        int maxIndex = Math.max(chromIndex, posIndex);

        maxIndex = Math.max(maxIndex, rsIndex);
        maxIndex = Math.max(maxIndex, strandIndex);
        boolean checkRegion = false;
        int regionNum = physicalRegions.length / 2;

        int j = 0;
        StringBuilder strRegion = new StringBuilder();
        for (j = 0; j < regionNum; j++) {
            if (physicalRegions[j * 2] != -9 || physicalRegions[j * 2 + 1] != -9) {
                checkRegion = true;
                if (physicalRegions[j * 2 + 1] == -9) {
                    physicalRegions[1] = Integer.MAX_VALUE;
                }
                strRegion.append(" [");
                strRegion.append(physicalRegions[j * 2]);
                strRegion.append(",");
                strRegion.append(physicalRegions[j * 2 + 1]);
                strRegion.append("]bp");
            }
        }

        int lineCounter = -1;

        int chroNum;
        if (chromosome.equalsIgnoreCase("X")) {
            chroNum = 23;
        } else if (chromosome.equalsIgnoreCase("Y")) {
            chroNum = 24;
        } else if (chromosome.equalsIgnoreCase("MT")) {
            chroNum = 26;
        } else if (chromosome.equalsIgnoreCase("XY")) {
            chroNum = 25;
        } else {
            chroNum = Integer.parseInt(chromosome);
        }

        if (checkRegion) {
            regionNum = physicalRegions.length / 2;
        }
        boolean invalid = false;
        int pos = 0;

        boolean diffChrom = false;
        String tmpStr;
        String rsID = null;
        StringBuilder tmpBuffer = new StringBuilder();
        int index;
        Map<Integer, Integer> physicalPos = new HashMap<Integer, Integer>(10000, 0.9f);
        double decodeGenet = -9;
        char strand = MISSING_STRAND_NAME;
        int alleleAIndex = 4;
        int alleleBIndex = 5;
        maxIndex = Math.max(maxIndex, alleleAIndex);
        maxIndex = Math.max(maxIndex, alleleBIndex);
        int currentChromNum = 0;
        char alleleA = MISSING_ALLELE_NAME, alleleB = MISSING_ALLELE_NAME;
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
                alleleA = MISSING_ALLELE_NAME;
                alleleB = MISSING_ALLELE_NAME;
                while (tokenizer.hasMoreTokens()) {
                    tmpBuffer.delete(0, tmpBuffer.length());
                    tmpBuffer.append(tokenizer.nextToken().trim());
                    tmpStr = tmpBuffer.toString();
                    if (index == chromIndex) {
                        //plink binary format only have numbers for all chromsomes
                        if (tmpStr.equals("X")) {
                            currentChromNum = 23;
                        } else {
                            currentChromNum = Integer.parseInt(tmpStr);
                        }
                        if (currentChromNum != chroNum) {
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

                        pos = Integer.parseInt(tmpStr);
                        if (physicalPos.get(pos) != null) {
                            invalid = true;
                            break;
                        }

                        //filter physical region
                        if (checkRegion) {
                            invalid = true;
                            for (int i = 0; i < regionNum; i++) {
                                if (pos >= physicalRegions[i * 2] && pos <= physicalRegions[i * 2 + 1]) {
                                    invalid = false;
                                    break;
                                }
                            }
                        }
                    } else if (index == strandIndex) {
                        strand = tmpStr.charAt(0);
                    } else if (index == alleleAIndex) {
                        alleleA = tmpStr.charAt(0);
                    } else if (index == alleleBIndex) {
                        alleleB = tmpStr.charAt(0);
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
                AnnotSNP aSNP = new AnnotSNP(rsID, chroNum, pos, -1);
                aSNP.setAAllele(alleleA);
                aSNP.setBAllele(alleleB);

                //a global index in the map file
                aSNP.probeSetID = lineCounter;
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

            //System.out.println(runningInfo.toString());
            return (snpList.size());
        } finally {
            br.close();
        }
    }
    //___________________________________________________TEST___________________________________________________________________________________________


    public int readPlinkBinaryFormatMapByPysic(List<AnnotSNP> snpList) throws Exception {
        File mapFile = new File(mapFileName);
        if (!mapFile.exists()) {
            return -1;
        }
        BufferedReader br = new BufferedReader(new FileReader(mapFile));
        String line = null;

        String delmilit = ", \t";
        int chromIndex = 0;
        int posIndex = 3;
        int rsIndex = 1;
        int strandIndex = -1;
        int maxIndex = Math.max(chromIndex, posIndex);

        maxIndex = Math.max(maxIndex, rsIndex);
        maxIndex = Math.max(maxIndex, strandIndex);

        int j = 0;

        int lineCounter = -1;

        int chroNum = 0;

        boolean invalid = false;
        int pos = 0;

        boolean diffChrom = false;
        String tmpStr;
        String rsID = null;
        StringBuilder tmpBuffer = new StringBuilder();
        int index;
        Map<Integer, Integer> physicalPos = new HashMap<Integer, Integer>(10000, 0.9f);
        double decodeGenet = -9;
        char strand = MISSING_STRAND_NAME;
        int alleleAIndex = 4;
        int alleleBIndex = 5;
        maxIndex = Math.max(maxIndex, alleleAIndex);
        maxIndex = Math.max(maxIndex, alleleBIndex);
        int currentChromNum = 0;
        char alleleA = MISSING_ALLELE_NAME, alleleB = MISSING_ALLELE_NAME;
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
                alleleA = MISSING_ALLELE_NAME;
                alleleB = MISSING_ALLELE_NAME;
                while (tokenizer.hasMoreTokens()) {
                    tmpBuffer.delete(0, tmpBuffer.length());
                    tmpBuffer.append(tokenizer.nextToken().trim());
                    tmpStr = tmpBuffer.toString();
                    if (index == chromIndex) {
                        //plink binary format only have numbers for all chromsomes
                        if (tmpStr.equals("X")) {
                            currentChromNum = 23;
                        } else {
                            currentChromNum = Integer.parseInt(tmpStr);
                        }

                    } else if (index == rsIndex) {
                        rsID = tmpStr;
                    } else if (index == posIndex) {
                        //filter items without physical postion
                        if (tmpStr.startsWith("-")) {
                            invalid = true;
                            break;
                        }

                        pos = Integer.parseInt(tmpStr);
                        if (physicalPos.get(pos) != null) {
                            invalid = true;
                            break;
                        }

                    } else if (index == strandIndex) {
                        strand = tmpStr.charAt(0);
                    } else if (index == alleleAIndex) {
                        alleleA = tmpStr.charAt(0);
                    } else if (index == alleleBIndex) {
                        alleleB = tmpStr.charAt(0);
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
                physicalPos.put(pos, 0);

                // System.out.println(linecells[1]);
                AnnotSNP aSNP = new AnnotSNP(rsID, chroNum, pos, -1);
                aSNP.setAAllele(alleleA);
                aSNP.setBAllele(alleleB);

                //a global index in the map file
                aSNP.probeSetID = lineCounter;
                //MISSING_STRAND_NAME indicates unknown strand information
                strand = MISSING_STRAND_NAME;
                snpList.add(aSNP);

            }

            StringBuilder runningInfo = new StringBuilder();
            runningInfo.append("The number of SNPs  map file ");
            runningInfo.append(mapFile.getName());
            runningInfo.append(" is ");
            runningInfo.append(snpList.size());
            runningInfo.append(".");

            System.out.println(runningInfo.toString());
            return (snpList.size());
        } finally {
            br.close();
        }
    }

    public int countSNPNuminChrom(String chromosome) throws Exception {
        File mapFile = new File(mapFileName);
        if (!mapFile.exists()) {
            return -1;
        }
        BufferedReader br = new BufferedReader(new FileReader(mapFile));
        String line = null;

        String delmilit = ", \t";

        int chroNum;
        if (chromosome.equalsIgnoreCase("X")) {
            chroNum = 23;
        } else if (chromosome.equalsIgnoreCase("Y")) {
            chroNum = 24;
        } else if (chromosome.equalsIgnoreCase("MT")) {
            chroNum = 26;
        } else if (chromosome.equalsIgnoreCase("XY")) {
            chroNum = 25;
        } else {
            chroNum = Integer.parseInt(chromosome);
        }
        int chromSNPCount = 0;

        int currentChromNum;
        while ((line = br.readLine()) != null) {
            if (line.trim().length() == 0) {
                continue;
            }
            StringTokenizer tokenizer = new StringTokenizer(line, delmilit);

            if (tokenizer.hasMoreTokens()) {
                currentChromNum = Integer.parseInt(tokenizer.nextToken().trim());
                if (currentChromNum == chroNum) {
                    chromSNPCount++;
                }
            }
        }
        return chromSNPCount;
    }

    /**
     * retrieve SNP information from the Map files
     *
     * @param FileName
     * @param indices
     * @param arry
     * @param genetRegions pair of physical regions
     * @throws java.lang.Exception
     * @return
     */
    public int readPlinkBinaryFormatMapByBlock(List<AnnotSNP> snpList, String chromosome, int blockStart, int blockLen) throws Exception {
        File mapFile = new File(mapFileName);
        if (!mapFile.exists()) {
            return -1;
        }
        BufferedReader br = new BufferedReader(new FileReader(mapFile));
        String line = null;

        String delmilit = ", \t";
        int chromIndex = 0;
        int posIndex = 3;
        int rsIndex = 1;
        int strandIndex = -1;
        int maxIndex = Math.max(chromIndex, posIndex);

        maxIndex = Math.max(maxIndex, rsIndex);
        maxIndex = Math.max(maxIndex, strandIndex);

        int j = 0;
        int lineCounter = -1;

        int chroNum;
        if (chromosome.equalsIgnoreCase("X")) {
            chroNum = 23;
        } else if (chromosome.equalsIgnoreCase("Y")) {
            chroNum = 24;
        } else if (chromosome.equalsIgnoreCase("MT")) {
            chroNum = 26;
        } else if (chromosome.equalsIgnoreCase("XY")) {
            chroNum = 25;
        } else {
            chroNum = Integer.parseInt(chromosome);
        }

        boolean invalid = false;
        int pos = 0;

        boolean diffChrom = false;
        String tmpStr;
        String rsID = null;
        StringBuilder tmpBuffer = new StringBuilder();
        int index;
        Map<Integer, Integer> physicalPos = new HashMap<Integer, Integer>(10000, 0.9f);
        double decodeGenet = -9;
        char strand = MISSING_STRAND_NAME;
        int alleleAIndex = 4;
        int alleleBIndex = 5;
        maxIndex = Math.max(maxIndex, alleleAIndex);
        maxIndex = Math.max(maxIndex, alleleBIndex);
        int currentChromNum = 0;
        char alleleA = MISSING_ALLELE_NAME, alleleB = MISSING_ALLELE_NAME;
        int chromSNPCount = 0;
        int reservedSNPCount = 0;
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
                alleleA = MISSING_ALLELE_NAME;
                alleleB = MISSING_ALLELE_NAME;
                while (tokenizer.hasMoreTokens()) {
                    tmpBuffer.delete(0, tmpBuffer.length());
                    tmpBuffer.append(tokenizer.nextToken().trim());
                    tmpStr = tmpBuffer.toString();
                    if (index == chromIndex) {
                        //plink binary format only have numbers for all chromsomes
                        currentChromNum = Integer.parseInt(tmpStr);
                        if (currentChromNum != chroNum) {
                            diffChrom = true;
                            break;
                        } else {
                            chromSNPCount++;
                        }
                    } else if (index == rsIndex) {
                        rsID = tmpStr;
                    } else if (index == posIndex) {
                        //filter items without physical postion
                        if (tmpStr.startsWith("-")) {
                            invalid = true;
                            break;
                        }

                        pos = Integer.parseInt(tmpStr);
                        if (physicalPos.get(pos) != null) {
                            invalid = true;
                            break;
                        }
                    } else if (index == strandIndex) {
                        strand = tmpStr.charAt(0);
                    } else if (index == alleleAIndex) {
                        alleleA = tmpStr.charAt(0);
                    } else if (index == alleleBIndex) {
                        alleleB = tmpStr.charAt(0);
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

                if (chromSNPCount < blockStart) {
                    continue;
                }

                reservedSNPCount++;
                if (reservedSNPCount > blockLen) {
                    break;
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
                AnnotSNP aSNP = new AnnotSNP(rsID, chroNum, pos, -1);
                aSNP.setAAllele(alleleA);
                aSNP.setBAllele(alleleB);
                aSNP.probeSetID = lineCounter;
                //MISSING_STRAND_NAME indicates unknown strand information
                strand = MISSING_STRAND_NAME;
                snpList.add(aSNP);
            }

            return (snpList.size());
        } finally {
            br.close();
        }
    }

    /**
     * retrieve SNP information from the Map files
     *
     * @param FileName
     * @param indices
     * @param arry
     * @param genetRegions pair of physical regions
     * @throws java.lang.Exception
     * @return
     */
    public int readPlinkBinaryFormatMapBySNPList(List<AnnotSNP> snpList, String chromosome, HashMap<Integer, AnnotSNP> chromSNPPValue) throws Exception {
        File mapFile = new File(mapFileName);
        if (!mapFile.exists()) {
            return -1;
        }
        BufferedReader br = new BufferedReader(new FileReader(mapFile));
        String line = null;

        String delmilit = ", \t";
        int chromIndex = 0;
        int posIndex = 3;
        int rsIndex = 1;
        int strandIndex = -1;
        int maxIndex = Math.max(chromIndex, posIndex);

        maxIndex = Math.max(maxIndex, rsIndex);
        maxIndex = Math.max(maxIndex, strandIndex);
        int j = 0;
        StringBuilder strRegion = new StringBuilder();

        int lineCounter = -1;

        int chroNum;
        if (chromosome.equalsIgnoreCase("X")) {
            chroNum = 23;
        } else if (chromosome.equalsIgnoreCase("Y")) {
            chroNum = 24;
        } else if (chromosome.equalsIgnoreCase("MT")) {
            chroNum = 26;
        } else if (chromosome.equalsIgnoreCase("XY")) {
            chroNum = 25;
        } else {
            chroNum = Integer.parseInt(chromosome);
        }

        boolean invalid = false;
        int pos = 0;

        boolean diffChrom = false;
        String tmpStr;
        String rsID = null;
        StringBuilder tmpBuffer = new StringBuilder();
        int index;
        Map<Integer, Integer> physicalPos = new HashMap<Integer, Integer>(10000, 0.9f);
        double decodeGenet = -9;
        char strand = MISSING_STRAND_NAME;
        int alleleAIndex = 4;
        int alleleBIndex = 5;
        maxIndex = Math.max(maxIndex, alleleAIndex);
        maxIndex = Math.max(maxIndex, alleleBIndex);
        int currentChromNum = 0;
        char alleleA = MISSING_ALLELE_NAME, alleleB = MISSING_ALLELE_NAME;
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
                alleleA = MISSING_ALLELE_NAME;
                alleleB = MISSING_ALLELE_NAME;
                while (tokenizer.hasMoreTokens()) {
                    tmpBuffer.delete(0, tmpBuffer.length());
                    tmpBuffer.append(tokenizer.nextToken().trim());
                    tmpStr = tmpBuffer.toString();
                    if (index == chromIndex) {
                        //plink binary format only have numbers for all chromsomes
                        currentChromNum = Integer.parseInt(tmpStr);
                        if (currentChromNum != chroNum) {
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

                        pos = Integer.parseInt(tmpStr);
                        if (physicalPos.get(pos) != null) {
                            invalid = true;
                            break;
                        }
                        //filter by physical  
                        if (!chromSNPPValue.containsKey(pos)) {
                            invalid = true;
                            break;
                        }

                    } else if (index == strandIndex) {
                        strand = tmpStr.charAt(0);
                    } else if (index == alleleAIndex) {
                        alleleA = tmpStr.charAt(0);
                    } else if (index == alleleBIndex) {
                        alleleB = tmpStr.charAt(0);
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
                AnnotSNP aSNP = new AnnotSNP(rsID, chroNum, pos, -1);
                aSNP.setAAllele(alleleA);
                aSNP.setBAllele(alleleB);
                aSNP.probeSetID = lineCounter;
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

            System.out.println(runningInfo.toString());
            return (snpList.size());
        } finally {
            br.close();
        }
    }
}
