package org.cobi.eqtlsimulator;

import cern.colt.matrix.DoubleMatrix2D;
import org.cobi.genetsimulator.entity.AnnotSNP;

import java.io.*;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class CAVIAR_FILE {
    String FilePath = "/home/ds/FineMapping/caviar/test/";

    public String z_files(List<AnnotSNP> FullSnpList, double[] zscores) throws IOException {
        String ZFilePath = FilePath + "data.Z";
        BufferedWriter ZWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(ZFilePath)));
        for(int i = 0; i < FullSnpList.size(); i++) {
            AnnotSNP snp = FullSnpList.get(i);
            ZWriter.write(snp.getRSID() + " " + zscores[i]);
            ZWriter.newLine();
        }
        ZWriter.flush();
        ZWriter.close();
        return ZFilePath;
    }


    public String ld_file(DoubleMatrix2D LDMatrix) throws IOException {
        String LDFilePath = FilePath + "data.LD";
        BufferedWriter LDWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(LDFilePath)));
        for(int i = 0; i < LDMatrix.rows(); i++){
            for(int j =0; j < LDMatrix.columns(); j++){
                double ld = LDMatrix.getQuick(i, j);
                LDWriter.write(Double.toString(ld));
                if(j < LDMatrix.columns() - 1){
                    LDWriter.write(" ");
                }
            }
            LDWriter.newLine();
        }
        LDWriter.flush();
        LDWriter.close();
        return LDFilePath;
    }

    public void CAVIAR(){
        try {
            String[] cmd = new String[]{"/home/ds/FineMapping/caviar/CAVIAR-C++/CAVIAR", "-o", "/home/ds/FineMapping/caviar/test/data.out", "-l", "/home/ds/FineMapping/caviar/test/data.LD", "-z", "/home/ds/FineMapping/caviar/test/data.Z", "-c", "2"};
            Process ps = Runtime.getRuntime().exec(cmd);
            BufferedReader br = new BufferedReader(new InputStreamReader(ps.getInputStream()));
            StringBuffer sb = new StringBuffer();
            String line;
            while ((line = br.readLine()) != null) {
                sb.append(line).append("\n");
            }
            String result = sb.toString();
            System.out.println(result);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public int[] out_file() throws IOException {
        String SNPFilePath = FilePath + "data.out_post";
        BufferedReader SNPReader = new BufferedReader(new InputStreamReader(new FileInputStream(SNPFilePath)));
        String line;

        int count = 0;
        Map<String, Double> snp = new HashMap<>();
        while((line = SNPReader.readLine()) != null){
            if(count ++ == 0){
                continue;
            }
            String[] str = line.split("\\s+");
            snp.put(str[0], Double.valueOf(str[2]));
        }

        String[] snpid = new String[2];
        double max1 = 0, max2 = 0;
        for(Map.Entry<String, Double> entry : snp.entrySet()){
            if(entry.getValue() > max1){
                snpid[1] = snpid[0];
                max2 = max1;
                snpid[0] = entry.getKey();
                max1 = entry.getValue();
            }else if(entry.getValue() > max2){
                snpid[1] = entry.getKey();
                max2 = entry.getValue();
            }
        }

        int[] index = new int[2];
        index[0] = -1;
        index[1] = -1;
        String ZFilePath = FilePath + "data.Z";
        BufferedReader ZFileReader = new BufferedReader(new InputStreamReader(new FileInputStream(ZFilePath)));
        count = -1;
        boolean flag1 = false, flag2 = false;
        while ((line = ZFileReader.readLine()) != null){
            count ++;
            String[] str = line.split("\\s+");
            if(snpid[0].equals(str[0])){
                index[0] = count;
                flag1 = true;
            }else if(snpid[1].equals(str[0])) {
                index[1] = count;
                flag2 = true;
            }

            if(flag1 && flag2){
                break;
            }
        }
        return index;
    }
}
