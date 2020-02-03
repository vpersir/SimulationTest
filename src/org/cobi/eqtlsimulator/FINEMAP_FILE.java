package org.cobi.eqtlsimulator;

import cern.colt.list.IntArrayList;
import cern.colt.matrix.DoubleMatrix2D;
import org.cobi.genetsimulator.entity.AnnotSNP;
import org.cobi.util.stat.LogisticRegression;

import java.io.*;
import java.util.List;

public class FINEMAP_FILE {
    String FilePath = "/home/ds/FineMapping/finemap_v1.3_x86_64/test/";

    public void master_file(List<AnnotSNP> FullSnpList, DoubleMatrix2D LDMatrix, double[] beta, double[] se) throws IOException {
        String MasterFilePath = FilePath + "master";
        BufferedWriter MasterWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(MasterFilePath)));
        MasterWriter.write("z;ld;snp;config;cred;log;n_samples");
        MasterWriter.newLine();
        String str = z_file(FullSnpList, beta, se) + ";" + ld_file(LDMatrix) + ";" + FilePath +
                "data.snp;" + FilePath + "data.config;" + FilePath + "data.cred;" + FilePath + "data.log;" + FullSnpList.size();
        MasterWriter.write(str);
        MasterWriter.flush();
        MasterWriter.close();
    }
    public String z_file(List<AnnotSNP> FullSnpList, double[] beta, double[] se) throws IOException {
        String ZFilePath = FilePath + "data.z";
        BufferedWriter ZWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(ZFilePath)));
        ZWriter.write("rsid chromosome position allele1 allele2 maf beta se");
        ZWriter.newLine();
        for(int i = 0; i < FullSnpList.size(); i++) {
            AnnotSNP snp = FullSnpList.get(i);
            ZWriter.write(snp.getRSID() + " " + snp.getChromosomeNum() + " " + snp.getPhysicalPosition() + " "
            + snp.getAAllele() + " " + snp.getBAllele() + " " + snp.getAAlleleFreq() + " " + beta[i] + " " + se[i]);
            ZWriter.newLine();
        }
        ZWriter.flush();
        ZWriter.close();
        return ZFilePath;
    }

    public String ld_file(DoubleMatrix2D LDMatrix) throws IOException {
        String LDFilePath = FilePath + "data.ld";
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

    public String k_file(double[] arrayW) throws IOException {
        String KFilePath = FilePath + "data.k";
        BufferedWriter KWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(KFilePath)));

        double sumW = 0.0;
        for(int i = 0; i < arrayW.length; i++){
            sumW += arrayW[i];
        }

        for(int i = 0; i < arrayW.length; i++){
            KWriter.write(Double.toString(arrayW[i] / sumW));
            if(i < arrayW.length - 1){
                KWriter.write(" ");
            }
        }
        KWriter.flush();
        KWriter.close();
        return KFilePath;
    }

    public void fineMap(){
        try {
            // /home/ds/FineMapping/finemap_v1.3_x86_64/finemap_v1.3.1_x86_64 --sss --in-files example/master
            String[] cmd = new String[]{"/home/ds/FineMapping/finemap_v1.3_x86_64/finemap_v1.3_x86_64", "--sss", "--in-files", "/home/ds/FineMapping/finemap_v1.3_x86_64/test/master", "--n-causal-snps", "2"};
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

    public int[] snp_file() throws IOException {
        String SNPFilePath = FilePath + "data.snp";
        BufferedReader SNPReader = new BufferedReader(new InputStreamReader(new FileInputStream(SNPFilePath)));
        String line;
        int count = -1;
        int[] index = new int[2];
        index[0] = -1;
        index[1] = -1;
        while((line = SNPReader.readLine()) != null){
            if(count == -1){
                count++;
                continue;
            }
            if(count < 2) {
                String[] str = line.split(" ");
                index[count] = Integer.parseInt(str[0]);
                count++;
                continue;
            }
            break;
        }
        return index;
    }
}
