/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.eqtlsimulator;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashMap;
import java.util.Map;
import java.util.StringTokenizer;

/**
 *
 * @author MX Li
 */
public class Options {

    public String ldRFilePath = "geneSymbols.txt";
    public int qtlNum = 3;
    public int[] diseaseLoci;
    public int[] ldBlockBoundary;
    public boolean isIndependent = false;
    public double rrHeterozygous = 0;
    public double rrHomozygous = 0;
    public double allelicRR = 1.0;
    public int permutationTime = 0;
    public int markerBlockNum = 1;
    public String diseaseModel = "ADD";
    Map<String, String> optionMap;
    public String pedFileName;
    public String mapFileName;
    public String plinFileName;
    public String chromosome;
    public long[] regions;
    public String inputFolder = "./";
    public String outputFolder = "./";
    public double genetVarPercentage = 0.05;

    public String[] geneList;
    public double[][] determineCoeffs;

    public boolean parseOptions() throws Exception {
        pedFileName = optionMap.get("--ped-file");
        if (pedFileName == null) {
            String infor = "No --map-file option for gene symbols";
            //  throw new Exception(infor);
            //return false;
        } else {
            pedFileName = pedFileName.trim();
            System.out.println("Ped File: " + pedFileName);
        }

        mapFileName = optionMap.get("--map-file");
        if (mapFileName == null) {
            String infor = "No --map-file option for gene symbols";
            //  throw new Exception(infor);
            //return false;
        } else {
            mapFileName = mapFileName.trim();
            System.out.println("Map File: " + mapFileName);
        }

        plinFileName = optionMap.get("--plink-binary");
        if (plinFileName == null) {
            String infor = "No --plink-binary option for Plink binnary file";
            //  throw new Exception(infor);
            //return false;
        } else {
            plinFileName = plinFileName.trim();
            System.out.println("Plink binnary file: " + plinFileName);
        }

        inputFolder = optionMap.get("--input-folder");
        if (inputFolder == null) {
            String infor = "No ---input-folder option for Plink binnary file";
            //  throw new Exception(infor);
            //return false;
            inputFolder = "./";
        } else {
            inputFolder = inputFolder.trim();
            if (inputFolder.endsWith("\\")) {
                inputFolder = inputFolder + "\\";
            } else if (!inputFolder.endsWith("/")) {
                inputFolder = inputFolder + "/";
            }
            System.out.println("Input folder: " + inputFolder);
        }
        outputFolder = optionMap.get("--output-folder");
        if (outputFolder == null) {
            String infor = "No --output-folder option for Plink binnary file";
            //  throw new Exception(infor);
            //return false;
            outputFolder = "./";
        } else {
            outputFolder = outputFolder.trim();
            if (!outputFolder.endsWith("/")) {
                outputFolder = outputFolder;
            }
            System.out.println("Output folder: " + outputFolder);
        }

        chromosome = optionMap.get("--chromsome");
        if (chromosome == null) {
            String infor = "No --chromosome option for Plink binnary file";
            //  throw new Exception(infor);
            //return false;
        } else {
            chromosome = chromosome.trim();
            System.out.println("Chromosome: " + chromosome);
        }

        String tmp = optionMap.get("--regions");
        if (tmp == null) {
            String infor = "No --chromosome option for Plink binnary file";
            //  throw new Exception(infor);
            //return false;
        } else {
            System.out.println("Regions: " + tmp);
            String[] r1 = tmp.split(";");
            regions = new long[r1.length * 2];
            for (int i = 0; i < r1.length; i++) {
                String[] b1 = r1[i].split("-");
                regions[i * 2] = Long.parseLong(b1[0].trim());
                regions[i * 2 + 1] = Long.parseLong(b1[1].trim());
            }
        }
        tmp = optionMap.get("--gene-list");
        if (tmp == null) {
            String infor = "No --chromosome option for Plink binnary file";
            //  throw new Exception(infor);
            //return false;
        } else {
            System.out.println("Genes: " + tmp);
            geneList = tmp.trim().split(",");
        }

        tmp = optionMap.get("--gene-r2");
        if (tmp == null) {
            String infor = "No --gene-r2 option for Plink binnary file";
            //  throw new Exception(infor);
            //return false;
        } else {
            System.out.println("Gene's detrmination coefficients: " + tmp);
            String[] r1 = tmp.trim().split(";");
            determineCoeffs = new double[r1.length][r1.length];
            for (int i = 0; i < r1.length; i++) {
                String[] b1 = r1[i].split(",");
                for (int j = 0; j < b1.length; j++) {
                    determineCoeffs[i][j] = Double.parseDouble(b1[j]);
                }
            }
        }

        String para = optionMap.get("--qtl-num");
        if (para == null) {
            String infor = "No --qtl-num option for the marker number";
            throw new Exception(infor);
            //return false;
        } else {
            para = para.trim();
            System.out.println("eQTL Num: " + para);
            qtlNum = Integer.parseInt(para);
        }

        para = optionMap.get("--disease-loci");
        if (para == null) {
            String infor = "No --disease-loci option for the disease loci";
            //throw new Exception(infor);
            //return false;
        } else {
            para = para.trim();
            System.out.println("Disease Loci: " + para);
            StringTokenizer tokenizer = new StringTokenizer(para);
            diseaseLoci = new int[tokenizer.countTokens()];
            int i = 0;
            while (tokenizer.hasMoreTokens()) {
                diseaseLoci[i] = Integer.parseInt(tokenizer.nextToken().trim());
                i++;
            }
        }
        para = optionMap.get("--ld-block-boundary");
        if (para == null) {
            String infor = "No --ld-block-boundary option for the disease loci";
            //throw new Exception(infor);
            //System.out.println(infor);
            //return false;
        } else {
            para = para.trim();
            System.out.println("LD block boundary: " + para);
            StringTokenizer tokenizer = new StringTokenizer(para);
            ldBlockBoundary = new int[tokenizer.countTokens()];
            int i = 0;
            while (tokenizer.hasMoreTokens()) {
                ldBlockBoundary[i] = Integer.parseInt(tokenizer.nextToken().trim());
                i++;
            }
        }
        para = optionMap.get("--marker-block-num");
        if (para == null) {
            String infor = "No --marker-block-num option for  marker block number and use the default " + markerBlockNum;
            //throw new Exception(infor);
            //System.out.println(infor);
            //return false;
        } else {
            para = para.trim();
            System.out.println("Marker block number: " + para);
            markerBlockNum = Integer.parseInt(para);
        }

        para = optionMap.get("--independent");
        if (para == null) {
            String infor = "No --independent option for the independence of markers";
           // throw new Exception(infor);
            //return false;
        } else {
            para = para.trim();
            System.out.println("Independent: " + para);
            if (para.equals("N")) {
                isIndependent = false;
            } else {
                isIndependent = true;
            }
        }

        para = optionMap.get("--genet-var");
        if (para == null) {
            String infor = "No --marker-block-num option for  marker block number and use the default " + markerBlockNum;
            //throw new Exception(infor);
            //  System.out.println(infor);
            //return false;
        } else {
            para = para.trim();
            System.out.println("Genetic variance: " + para);
            genetVarPercentage = Double.parseDouble(para);
        }


        /*
        para = optionMap.get("--rr-heterozygous").trim();
        if (para == null) {
        String infor = "No --rr-heterozygous option for the relative risk of heterozygous";
        throw new Exception(infor);
        //return false;
        } else {
        System.out.println("Relative risk of heterozygous: " + para);
        rrHeterozygous = Double.parseDouble(para);
        }
        
        para = optionMap.get("--rr-homozygous").trim();
        if (para == null) {
        String infor = "No --rr-homozygous option for the relative risk of homozygous";
        throw new Exception(infor);
        //return false;
        } else {
        System.out.println("Relative risk of homozygous: " + para);
        rrHomozygous = Double.parseDouble(para);
        }
         */
        para = optionMap.get("--allelic-relative-risk");
        if (para == null) {
            String infor = "No --allelic-relative-risk option for the allelic relative risk";
            //throw new Exception(infor);
            //return false;
        } else {
            para = para.trim();
            System.out.println("Allelic relative risk: " + para);
            allelicRR = Double.parseDouble(para);
        }

        para = optionMap.get("--disease-model");
        if (para == null) {
            String infor = "No --disease-model option for the disease model";
            //throw new Exception(infor);
            //return false;
        } else {
            para = para.trim();
            System.out.println("Disease model: " + para);
            diseaseModel = para;
        }
        para = optionMap.get("--permutation-time").trim();
        if (para == null) {
            String infor = "No --permutation-time option for the permutation time";
            throw new Exception(infor);
            //return false;
        } else {
            para = para.trim();
            System.out.println("Permuation time: " + para);
            permutationTime = Integer.parseInt(para);
        }
        return true;
    }

    public void readOptions(String fileName) throws Exception {
        BufferedReader br = new BufferedReader(new FileReader(fileName));
        String line = "";
        int lineNumber = 0;
        StringBuilder tmpStr = new StringBuilder();
        optionMap = new HashMap<String, String>();
        String name;

        //assume every parameter has a line
        while ((line = br.readLine()) != null) {
            if (line.trim().length() == 0) {
                continue;
            }
            int endNum = line.indexOf('#');
            if (endNum > 0) {
                line = line.substring(0, endNum);
            }

            StringTokenizer tokenizer = new StringTokenizer(line);
            //sometimes tokenizer.nextToken() can not release memory
            //parameter Name
            name = tmpStr.append(tokenizer.nextToken().trim()).toString();
            tmpStr.delete(0, tmpStr.length());
            //parameter value
            while (tokenizer.hasMoreTokens()) {
                tmpStr.append(tokenizer.nextToken().trim());
                tmpStr.append(" ");
            }
            optionMap.put(name, tmpStr.toString());
            tmpStr.delete(0, tmpStr.length());
            lineNumber++;
        }
        br.close();

    }
}
