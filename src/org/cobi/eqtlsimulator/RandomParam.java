package org.cobi.eqtlsimulator;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class RandomParam {
    /*
        --input-folder ./
        --output-folder ./out
        --plink-binary flip_cleanedTW_SC_hwe
        --chromsome 1
        --qtl-num 4
        --gene-list PRDM16,PRKCZ,GNB1
        --gene-r2 1,0,0;0,1,0.1;0,0,1
        --permutation-time 1
        --genet-var 0.2
     */
    static String[] option = {"--input-folder",
                            "--output-folder",
                            "--plink-binary",
                            "--chromsome",
                            "--qtl-num",
                            "--gene-list",
                            "--gene-r2",
                            "--permutation-time",
                            "--genet-var"};


    public RandomParam(String chr, int gene_num, int qtl_num, int permutationTime, double genet) {
        try {
            String[] info = {
                    " ./",
                    " ./out",
                    " flip_cleanedTW_SC_hwe",
                    " ",
                    " ",
                    " ",
                    " ",
                    " ",
                    " "};

            info[3] += chr;
            info[4] += Integer.toString(qtl_num);

            BufferedReader br = null;
            br = new BufferedReader(new FileReader(new File("F:\\java\\SimulationTestCode\\SeqGeneB36.txt")));
            String line;
            int gene_count = 0;
            boolean region_flag = false;
            List<String> gene_list = new ArrayList<String>();
            while ((line = br.readLine()) != null) {
                if (gene_count++ == 0) {
                    continue;
                }
                String[] str = line.split("\t");
                if(str[2].equals(chr)){
                    region_flag = true;
                    gene_list.add(str[1]);
                } else if(!region_flag){
                    continue;
                } else if(region_flag){
                    break;
                }
            }
            br.close();
            Collections.shuffle(gene_list);
            for(int i = 0; i < gene_num; i++){
                info[5] += gene_list.get(i);
                if(i != gene_num - 1){
                    info[5] += ",";
                }
            }

            int r2_num = gene_num * (gene_num - 1) / 2;
            String[] r2 = new String[r2_num];
            int r2_count = 0;
            for (int i = 0; i < gene_num; i++) {
                for (int j = 0; j < gene_num; j++) {
                    if (i == j) {
                        info[6] += "1";
                    } else if (i + j - 1 < r2_count) {
                        info[6] += r2[i + j - 1];
                    } else {
                        r2[r2_count] = Double.toString(Math.random());
                        info[6] += r2[r2_count];
                        r2_count++;
                    }

                    if (j != gene_num - 1) {
                        info[6] += ",";
                    }
                }
                if (i != gene_num - 1) {
                    info[6] += ";";
                }
            }

            info[7] += Integer.toString(permutationTime);
            info[8] += Double.toString(genet);

            BufferedWriter bw = new BufferedWriter(new FileWriter(new File("F:\\java\\SimulationTestCode\\param.txt")));
            for (int i = 0; i < info.length; i++) {
                bw.write(option[i] + info[i] + "\r\n");
            }
            bw.flush();
            bw.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        RandomParam param = new RandomParam("1",1,2, 1, 0.005);

    }
}
