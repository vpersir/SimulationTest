package org.cobi.eqtlsimulator;

import cern.colt.matrix.DoubleMatrix2D;
import org.cobi.genetsimulator.entity.AnnotSNP;

import java.io.*;
import java.util.List;

public class DAP_FILE {
    String FilePath = "/home/ds/FineMapping/dap/test/";

    public String est_file(List<AnnotSNP> FullSnpList, double[] beta, double[] se) throws IOException {
        String estFilePath = FilePath + "sim.est.dat";
        BufferedWriter estWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(estFilePath)));
        for(int i = 0; i < FullSnpList.size(); i++) {
            AnnotSNP snp = FullSnpList.get(i);
            estWriter.write(snp.getRSID() + " " + beta[i] + " " + se[i]);
            estWriter.newLine();
        }
        estWriter.flush();
        estWriter.close();
        return estFilePath;
    }
    public String z_file(List<AnnotSNP> FullSnpList, double[] Zscores) throws IOException {
        String ZFilePath = FilePath + "data.z";
        BufferedWriter ZWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(ZFilePath)));
        for(int i = 0; i < FullSnpList.size(); i++) {
            AnnotSNP snp = FullSnpList.get(i);
            ZWriter.write(snp.getRSID() + " " + Zscores[i]);
            ZWriter.newLine();
        }
        ZWriter.flush();
        ZWriter.close();
        return ZFilePath;
    }

    public String ld_file(DoubleMatrix2D LDMatrix) throws IOException {
        String LDFilePath = FilePath + "sim.LD.dat";
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

    public void dap_g(double[][] phenos){
        int indivSize = phenos.length;
        double SST = 0;
        for(int i = 0; i < indivSize; i++){
            SST += phenos[i][0];
        }
        try {
            //dap-g -d_est sample_data/sim.est.dat -d_ld sample_data/sim.1.LD.dat -d_n 343 -d_syy 515.6
            String[] cmd = new String[]{"/home/ds/FineMapping/dap/dap_src/dap-g", "-d_est", "/home/ds/FineMapping/dap/test/sim.est.dat", "-d_ld", "/home/ds/FineMapping/dap/test/sim.LD.dat", "-d_n", Integer.toString(indivSize), "-d_syy", Double.toString(SST), "-o", "/home/ds/FineMapping/dap/test/output.ss"};
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
        String SNPFilePath = FilePath + "output.ss";
        BufferedReader SNPReader = new BufferedReader(new InputStreamReader(new FileInputStream(SNPFilePath)));
        String line;
        String[] snpid = new String[2];
        boolean snpflag = false;
        int count = 0;
        while((line = SNPReader.readLine()) != null){
            String[] str = line.split("\\s+");
            if(str[0].equals("((1))")){
                snpflag = true;
            }

            if(snpflag){
                if(count == 0) {
                    snpid[0] = str[1];
                }else if (count == 1) {
                    snpid[1] = str[1];
                }else {
                    break;
                }
                count++;
            }
        }

        int[] index = new int[2];
        index[0] = -1;
        index[1] = -1;
        String estFilePath = FilePath + "sim.est.dat";
        BufferedReader estFileReader = new BufferedReader(new InputStreamReader(new FileInputStream(estFilePath)));
        count = -1;
        boolean flag1 = false, flag2 = false;
        while ((line = estFileReader.readLine()) != null){
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
