/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.genetsimulator.controller;

import cern.colt.list.DoubleArrayList;
import cern.jet.stat.Descriptive;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;
import org.cobi.util.text.LocalString;
 
import umontreal.ssj.randvar.NormalGen;
import umontreal.ssj.randvarmulti.MultinormalCholeskyGen;
import umontreal.ssj.randvarmulti.MultinormalGen;
import umontreal.ssj.rng.MT19937;
import umontreal.ssj.rng.WELL607;
 

/**
 *
 * @author mxli
 */
public class MultiplePhenotypeSimulator {

    public void mulitQuantitiveTraits(String pedFileName, double[][] corr, List<String> titles) throws Exception {
        BufferedReader br = new BufferedReader(new FileReader(pedFileName));

        String line;
        String delimiter = "\t\" \",/";


        int traitNum = corr.length;
        double[] mean = new double[traitNum];
        Arrays.fill(mean, 0);
        double[] values = new double[traitNum];
        List<String[]> indiviIDInPed = new ArrayList<String[]>();

        StringBuilder tmpBuffer = new StringBuilder();
        //long t3 = System.currentTimeMillis();
        while ((line = br.readLine()) != null) {
            line = line.toUpperCase();
            String[] cells = new String[2];
            StringTokenizer tokenizer = new StringTokenizer(line, delimiter);
            tmpBuffer.append(tokenizer.nextToken());
            cells[0] = tmpBuffer.toString();
            tmpBuffer.delete(0, tmpBuffer.length());

            tmpBuffer.append(tokenizer.nextToken());
            cells[1] = tmpBuffer.toString();
            tmpBuffer.delete(0, tmpBuffer.length());
            indiviIDInPed.add(cells);
        }
        br.close();
        int indivNum = indiviIDInPed.size();
        BufferedWriter bw = new BufferedWriter(new FileWriter(pedFileName + ".sim"));
        int[] seeds = new int[19];
        for (int i = 0; i < seeds.length; i++) {
            seeds[i] = (int) (Math.random() * 1000000);
            System.out.println(seeds[i]);
        }
        WELL607 we = new WELL607();
        we.setSeed(seeds);

        NormalGen ng = new NormalGen(new MT19937(we));
        MultinormalGen multiG = new MultinormalCholeskyGen(ng, mean, corr);

        DoubleArrayList[] valueLists = new DoubleArrayList[traitNum];
        for (int j = 0; j < traitNum; j++) {
            valueLists[j] = new DoubleArrayList();
        }
        bw.write("FID\tIID");
        for (String str : titles) {
            bw.write("\t");
            bw.write(str);
        }
        bw.write("\n");
        for (int i = 0; i < indivNum; i++) {
            multiG.nextPoint(values);
            bw.write(indiviIDInPed.get(i)[0]);
            bw.write("\t");
            bw.write(indiviIDInPed.get(i)[1]);

            for (int j = 0; j < traitNum; j++) {
                // indiviIDInPed.get(i)[j + 2] = String.valueOf(values[j]);
                bw.write("\t");
                bw.write(String.valueOf(values[j]));
                valueLists[j].add(values[j]);
            }
            bw.write("\n");
        }
        bw.close();

        for (int j = 0; j < traitNum; j++) {
            for (int t = 0; t < traitNum; t++) {
                double mean1 = Descriptive.mean(valueLists[j]);
                double mean2 = Descriptive.mean(valueLists[t]);
                double sd1 = Descriptive.sampleVariance(valueLists[j], mean1);
                double sd2 = Descriptive.sampleVariance(valueLists[t], mean2);
                double averagePCorre = (Descriptive.correlation(valueLists[j], Math.sqrt(sd1), valueLists[t], Math.sqrt(sd2)));
                System.out.print(averagePCorre);
                System.out.print(" ");
            }
            System.out.println();
        }
    }

    public void mergePhenoResults(String folder, List<String> names, int id) throws Exception {
        BufferedWriter bw = new BufferedWriter(new FileWriter(folder + "plink.all.qassoc." + id));
        Map<String, String[]> allMapedPVales = new HashMap<String, String[]>();
        String line;
        String delimiter = "\t\" \",/";
        bw.write("CHR\tSNP\tBP");
        for (int i = 0; i < names.size(); i++) {
            bw.write("\t");
            bw.write(names.get(i));

            BufferedReader br = new BufferedReader(new FileReader(folder + id+"." + names.get(i) + ".qassoc"));
            br.readLine();
            while ((line = br.readLine()) != null) {
                StringTokenizer tokenizer = new StringTokenizer(line, delimiter);
                String[] cells = new String[9];
                int t = 0;
                while (tokenizer.hasMoreTokens()) {
                    cells[t] = tokenizer.nextToken();
                    t++;
                }
                if (allMapedPVales.containsKey(cells[1])) {
                    String[] reuls = allMapedPVales.get(cells[1]);
                    reuls[3 + i] = cells[8];
                } else {
                    String[] reuls = new String[names.size() + 3];
                    reuls[0] = cells[0];
                    reuls[1] = cells[1];
                    reuls[2] = cells[2];
                    reuls[3 + i] = cells[8];
                    allMapedPVales.put(cells[1], reuls);
                }

            }
            br.close();
        }
        bw.write("\n");

        for (Map.Entry<String, String[]> m : allMapedPVales.entrySet()) {
            String[] reuls = m.getValue();
            bw.write(reuls[0]);

            for (int t = 1; t < reuls.length; t++) {
                bw.write("\t");
                bw.write(reuls[t]);
            }
            bw.write("\n");
        }
        bw.close();
    }

    public double[][] readCorrelationMatrix(String filePath, List<String> traitIDs) throws Exception {
        BufferedReader br = new BufferedReader(new FileReader(filePath));
        String line = null;
        int colNum = -1;

        line = br.readLine();
        line = line.trim();
        StringTokenizer tokenizer = new StringTokenizer(line);


        colNum = 0;
        while (tokenizer.hasMoreTokens()) {
            traitIDs.add(tokenizer.nextToken().trim());
        }
        colNum = traitIDs.size();
        double[][] corrMatrix = new double[colNum][colNum];
        int rowNum = 0;
        String tmpStr = null;
        while ((line = br.readLine()) != null) {
            line = line.trim();
            if (line.trim().length() == 0) {
                continue;
            }
            tokenizer = new StringTokenizer(line);

            colNum = 0;

            while (tokenizer.hasMoreTokens()) {
                tmpStr = tokenizer.nextToken().trim();
                if (LocalString.isNumeric(tmpStr)) {
                    // corrMatrix[rowNum][colNum] = Math.abs(Double.parseDouble(tmpStr));
                    corrMatrix[rowNum][colNum] = (Double.parseDouble(tmpStr));
                    colNum++;
                }
            }
            rowNum++;
        }
        br.close();
        return corrMatrix;
    }

    public static void main(String[] args) {
        try {
            MultiplePhenotypeSimulator mps = new MultiplePhenotypeSimulator();
            List<String> traitIDs = new ArrayList<String>();
            //double[][] corrMatrix = mps.readCorrelationMatrix("E:\\home\\mxli\\MyPapers\\Gene- and pathway-based genetic associaton analysis\\Documents\\multi phenotype\\data\\20131005\\PHENO_CORR_DEF.txt", traitIDs);
            double[][] corrMatrix = mps.readCorrelationMatrix("E:\\home\\mxli\\MyJava\\WAT\\PHENO_CORR_DEF.txt", traitIDs);
            // mps.mulitQuantitiveTraits("E:/home/mxli/MyJava/PI/58c_nbs.fam", corrMatrix, traitIDs); 
            //E:\home\mxli\MyJava\PI>plink --bfile 58c_nbs --assoc --pheno 58c_nbs.fam.sim --all-pheno

            //mps.mergePhenoResults("E:/home/mxli/MyJava/PI/", traitIDs);

            Runtime rt = Runtime.getRuntime();
            String line;
            for (int i = 2; i < 10; i++) {
                mps.mulitQuantitiveTraits("E:/home/mxli/MyJava/WAT/flip_cleanedTW_SC_hwe.fam", corrMatrix, traitIDs);
                Process p = rt.exec("E:/home/mxli/MyJava/WAT/plink --bfile E:/home/mxli/MyJava/WAT/flip_cleanedTW_SC_hwe --assoc --pheno E:/home/mxli/MyJava/WAT/flip_cleanedTW_SC_hwe.fam.sim --all-pheno --allow-no-sex --out E:/home/mxli/MyJava/WAT/"+i);

                BufferedReader input = new BufferedReader(new InputStreamReader(p.getInputStream()));
                while ((line = input.readLine()) != null) {
                    System.out.println(line);
                }
                input.close();

                System.out.println("Done.");
                p.waitFor();
                p.destroy();
                mps.mergePhenoResults("E:/home/mxli/MyJava/WAT/", traitIDs, i);
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
}
