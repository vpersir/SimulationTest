/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.genetsimulator.controller;

import java.util.ArrayList;
import java.util.List;
import org.cobi.eqtlsimulator.Constants;
 
import org.cobi.genetsimulator.entity.DiseaseSNP;
import org.cobi.util.text.LocalNumber;
 
 
/**
 *
 * @author mxli
 */
//an implementation of the m-locus model at http://biostat.mc.vanderbilt.edu/twiki/pub/Main/GWAsimulator/GWAsimulator_v2.0.pdf
public class MultilocusLogitModel implements Constants {
    double prevalence = 0.1;
    List<LogitDiseaseSNP> snpList = null;
    double alpha = 0.0;

    public MultilocusLogitModel(double prevalence) {
        this.prevalence = prevalence;
        snpList = new ArrayList<LogitDiseaseSNP>();
    }

    public double getPrevalence() {
        return prevalence;
    }

    public double getAlpha() {
        return alpha;
    }

    public void setAlpha(double alpha) {
        this.alpha = alpha;
    }



    public void setPrevalence(double prevalence) {
        this.prevalence = prevalence;
    }

    public List<LogitDiseaseSNP> getSnpList() {
        return snpList;
    }

    public void setSnpList(List<LogitDiseaseSNP> snpList) {
        this.snpList = snpList;
    }

  
 

    public void addSusceptibSNPLoci(double freq, int index, boolean riskAlleleLabel, double genotypicRelativeRiskAa, double genotypicRelativeRiskAA, boolean isXosome) throws Exception {
        LogitDiseaseSNP snp = new LogitDiseaseSNP(freq, genotypicRelativeRiskAa, genotypicRelativeRiskAA,riskAlleleLabel,index, 0.2,0.8, isXosome);
       
        snpList.add(snp);
    }


    private double calculateAveragePenetranceHWD(double alpha) {
        int snpSize = snpList.size();
        int times = (int) Math.pow(3, snpSize);
        double sum = 0.0;
        double jointProb = 1.0;
        double alphaBeta;
        for (int i = 0; i < times; i++) {
            String gtyLable = LocalNumber.decimal2XSystem(i, 3, snpSize);
            alphaBeta = alpha;
            jointProb = 1.0;
            for (int j = 0; j < snpSize; j++) {
                if (gtyLable.charAt(j) == '1') {
                    alphaBeta += snpList.get(j).getBeta1();
                    jointProb = jointProb * 2 * snpList.get(j).riskAlleleFrequency * (1 - snpList.get(j).riskAlleleFrequency);
                } else if (gtyLable.charAt(j) == '2') {
                    alphaBeta += snpList.get(j).getBeta2();
                    jointProb = jointProb * snpList.get(j).riskAlleleFrequency * snpList.get(j).riskAlleleFrequency;
                } else {
                    jointProb = jointProb * (1 - snpList.get(j).riskAlleleFrequency) * (1 - snpList.get(j).riskAlleleFrequency);
                }
            }
            alphaBeta = Math.exp(alphaBeta);
            sum += (jointProb / (1 + alphaBeta));
        }
        return (1 - sum);
    }
/*
    public void calculateJointGtyPenerance() {
        int snpSize = snpList.size();
        int times = (int) Math.pow(3, snpSize);
        jointGtyPenerance = new double[times];
        for (int i = 0; i < times; i++) {
            String gtyLable = LocalNumber.decimal2XSystem(i, 3, snpSize);
            double alphaBeta = alpha;
            for (int j = 0; j < snpSize; j++) {
                if (gtyLable.charAt(j) == '1') {
                    alphaBeta += snpList.get(j).getBeta1();
                } else if (gtyLable.charAt(j) == '2') {
                    alphaBeta += snpList.get(j).getBeta2();
                }
            }
            alphaBeta = Math.exp(alphaBeta);
            jointGtyPenerance[i] = 1 - 1 / (1 + alphaBeta);
        }
    }

    public double getJointGtyPenerance(StringBuffer gtyLable) {
        int order = LocalNumber.xSystem2decimal(gtyLable, 3);
        return jointGtyPenerance[order];
    }
*/
    public double exploreBetaAndAlphaHWD() {
        int snpSize = snpList.size();
        for (int i = 0; i < snpSize; i++) {
            snpList.get(i).calculateBetaHWD();
        }
        double apha2 = 0.0, apha1 = 0;
        double searchLen = 1.0;
        int grid = 4;
        boolean notYetFound = true;
        double precise = 1e-6;
        double tmp = 0;
        //this search algorithm is designed for increasing fuction
        while (notYetFound && searchLen >= precise) {
            tmp = calculateAveragePenetranceHWD(apha1);
            if (tmp < prevalence) {
                do {
                    apha1 += searchLen;
                    tmp = calculateAveragePenetranceHWD(apha1);
                } while (tmp < prevalence);
                apha2 = apha1 - searchLen;
            } else if (tmp == prevalence) {
                notYetFound = false;
            } else {
                do {
                    apha1 -= searchLen;
                    tmp = calculateAveragePenetranceHWD(apha1);
                } while (tmp > prevalence);
                apha2 = apha1 + searchLen;
            }
            if (tmp == prevalence) {
                notYetFound = false;
            } else {
                searchLen /= grid;
                apha1 = (apha1 + apha2) / 2;
            }
        }
        alpha = apha1;
        return alpha;
    }

    public double exploreBetaAndAlphaGty() {
        int snpSize = snpList.size();
        for (int i = 0; i < snpSize; i++) {
            snpList.get(i).calculateBetaHWD();
        }
        double apha2 = 0.0, apha1 = 0;
        double searchLen = 1.0;
        int grid = 4;
        boolean notYetFound = true;
        double precise = 1e-6;
        double tmp = 0;
        //this search algorithm is designed for increasing fuction
        while (notYetFound && searchLen >= precise) {
            tmp = calculateAveragePenetranceHWD(apha1);
            if (tmp < prevalence) {
                do {
                    apha1 += searchLen;
                    tmp = calculateAveragePenetranceHWD(apha1);
                } while (tmp < prevalence);
                apha2 = apha1 - searchLen;
            } else if (tmp == prevalence) {
                notYetFound = false;
            } else {
                do {
                    apha1 -= searchLen;
                    tmp = calculateAveragePenetranceHWD(apha1);
                } while (tmp > prevalence);
                apha2 = apha1 + searchLen;
            }
            if (tmp == prevalence) {
                notYetFound = false;
            } else {
                searchLen /= grid;
                apha1 = (apha1 + apha2) / 2;
            }
        }
        alpha = apha1;
        return alpha;
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        MultilocusLogitModel mulit = new MultilocusLogitModel(0.1);
        try {
            for (int i = 0; i < 3; i++) {
               // mulit.addSusceptibSNPLoci(i, true, 1.4, 1.4, false);
            }
            mulit.getSnpList().get(0).setRiskAlleleFrequency(0.4167);
            mulit.getSnpList().get(1).setRiskAlleleFrequency(0.0750);
            mulit.getSnpList().get(2).setRiskAlleleFrequency(0.2167);
            for (int i = 0; i < 3; i++) {
                mulit.getSnpList().get(i).calculateBetaHWD();
                System.out.println(mulit.getSnpList().get(i).beta1 + "\t" + mulit.getSnpList().get(i).beta2);
            }
            System.out.println(mulit.exploreBetaAndAlphaHWD());
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    public class LogitDiseaseSNP extends DiseaseSNP{
 
        double beta1 = 0.0;
        double beta2 = 0.0;      
        boolean isOnXChrom;
      
   
        public LogitDiseaseSNP(double riskAlleleFrequency, double genotypicRelativeRiskAa, double genotypicRelativeRiskAA,  boolean riskAlleleLable,int originalOrder, double markerAlleleFreq, double rSquare,boolean isOnXChrom) {
            super(riskAlleleFrequency, genotypicRelativeRiskAa, genotypicRelativeRiskAa, riskAlleleLable, originalOrder);            
                       this.isOnXChrom = isOnXChrom;
        }

        public double getBeta1() {
            return beta1;
        }

        public void setBeta1(double beta1) {
            this.beta1 = beta1;
        }

        public double getBeta2() {
            return beta2;
        }

        public void setBeta2(double beta2) {
            this.beta2 = beta2;
        }

   

        public void calculateBetaHWD() {
            double expAlpha = 0;
            if (isOnXChrom) {
                expAlpha = (0.5 * (1 - riskAlleleFrequency) * (2 - riskAlleleFrequency) + riskAlleleFrequency * (1 - riskAlleleFrequency) * genotypicRelativeRiskAa + 0.5 * riskAlleleFrequency * (1 + riskAlleleFrequency) * genotypicRelativeRiskAA) / prevalence - 1;
            } else {
                expAlpha = ((1 - riskAlleleFrequency) * (1 - riskAlleleFrequency) + 2 * riskAlleleFrequency * (1 - riskAlleleFrequency) * genotypicRelativeRiskAa + riskAlleleFrequency * riskAlleleFrequency * genotypicRelativeRiskAA) / prevalence - 1;
            }
            beta1 = -Math.log(((expAlpha + 1) / genotypicRelativeRiskAa - 1) / expAlpha);
            beta2 = -Math.log(((expAlpha + 1) / genotypicRelativeRiskAA - 1) / expAlpha);
            System.out.println("SNP "+LDMarkerPosition+" beta1 "+beta1+" beta2 "+beta2);
        }

        public void calculateBetaGty() {
            double expAlpha = 0;
            if (isOnXChrom) {
                expAlpha = (0.5 * (1 - riskAlleleFrequency) * (2 - riskAlleleFrequency) + riskAlleleFrequency * (1 - riskAlleleFrequency) * genotypicRelativeRiskAa + 0.5 * riskAlleleFrequency * (1 + riskAlleleFrequency) * genotypicRelativeRiskAA) / prevalence - 1;
            } else {
                expAlpha = ((1 - riskAlleleFrequency) * (1 - riskAlleleFrequency) + 2 * riskAlleleFrequency * (1 - riskAlleleFrequency) * genotypicRelativeRiskAa + riskAlleleFrequency * riskAlleleFrequency * genotypicRelativeRiskAA) / prevalence - 1;
            }
            beta1 = -Math.log(((expAlpha + 1) / genotypicRelativeRiskAa - 1) / expAlpha);
            beta2 = -Math.log(((expAlpha + 1) / genotypicRelativeRiskAA - 1) / expAlpha);
        }
    }
}
