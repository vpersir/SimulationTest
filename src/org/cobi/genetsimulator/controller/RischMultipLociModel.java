/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.genetsimulator.controller;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.cobi.genetsimulator.entity.DiseaseSNP;

import org.cobi.util.math.CombinationGenerator;

/**
 *
 * @author mxli
 * This class only works for a small number of DSL as it will calculate the joint probablity of genotypes
 */
public class RischMultipLociModel {

    public double baselineRisk = 0.0;
    public double prevalance = 0.1;
    int lociNum;
    //assume every allele has the same effect
    double[] jointGenotypePenerances;
    List<DiseaseSNP> diseaseSNPs = null;

    public void addDiseaseSNPLoci(double riskAlleleFreq, boolean riskAlleleLabel, int LDMarkerPosition) throws Exception {
        DiseaseSNP snp = new DiseaseSNP(riskAlleleFreq, riskAlleleLabel, LDMarkerPosition);
        diseaseSNPs.add(snp);
    }

    public void cleanDiseaseSNPLoci() {
        if (diseaseSNPs != null) {
            diseaseSNPs.clear();
        }
    }

    public List<DiseaseSNP> getDiseaseSNPs() {
        return diseaseSNPs;
    }

    public void setDiseaseSNPs(List<DiseaseSNP> diseaseSNPs) {
        this.diseaseSNPs = diseaseSNPs;
    }

    public RischMultipLociModel(double pre) {
        prevalance = pre;
        diseaseSNPs = new ArrayList<DiseaseSNP>();
    }

    public double getBaselineRisk() {
        return baselineRisk;
    }

    public void setBaselineRisk(double baselineRisk) {
        this.baselineRisk = baselineRisk;
    }

    public int getLociNum() {
        return lociNum;
    }

    public void setLociNum(int lociNum) {
        this.lociNum = lociNum;
    }

    public double[] getJointGenotypePenerances() {
        return jointGenotypePenerances;
    }

    public void setJointGenotypePenerances(double[] jointGenotypePenerances) {
        this.jointGenotypePenerances = jointGenotypePenerances;
    }

    public double getPrevalance() {
        return prevalance;
    }

    public void setPrevalance(double prevalance) {
        this.prevalance = prevalance;
    }

    public void calculateAdditiveModelParams(double allelicRisk) {
        lociNum = diseaseSNPs.size();
        if (lociNum == 0) {
            System.err.println("No disease SNPs!");
            return;
        }
        double[] riskAlleleFreq = new double[lociNum];
        for (int i = 0; i < lociNum; i++) {
            riskAlleleFreq[i] = diseaseSNPs.get(i).getRiskAlleleFrequency();
        }

        //assume each locus has two alleles
        int alleleNum = riskAlleleFreq.length * 2;

        int[] indices;
        Set<Integer> selectedIDs = new HashSet<Integer>();
        double[] nonRiskAlleleFreq = new double[lociNum];

        double accumulatedProb = 0;
        double localProb = 0;
        double tmpProb = 1;

        for (int i = 0; i < lociNum; i++) {
            nonRiskAlleleFreq[i] = 1 - riskAlleleFreq[i];
        }

        for (int a = 1; a <= alleleNum; a++) {
            CombinationGenerator x = new CombinationGenerator(alleleNum, a);
            localProb = 0;
            while (x.hasMore()) {
                indices = x.getNext();
                selectedIDs.clear();
                for (int i = 0; i < indices.length; i++) {
                    selectedIDs.add(indices[i]);
                }
                tmpProb = 1;
                for (int i = 0; i < alleleNum; i++) {
                    if (selectedIDs.contains(i)) {
                        tmpProb *= riskAlleleFreq[i / 2];
                    } else {
                        tmpProb *= nonRiskAlleleFreq[i / 2];
                    }
                }
                localProb += tmpProb;
            }
            accumulatedProb += (a * localProb);
        }
        baselineRisk = prevalance / (1 + accumulatedProb * (allelicRisk - 1));
        jointGenotypePenerances = new double[alleleNum + 1];
        for (int i = 0; i <= alleleNum; i++) {
            jointGenotypePenerances[i] = baselineRisk * (1 + i * (allelicRisk - 1));
            //System.out.println(jointGenotypePenerances[i]);
        }
    }

    public void calculateMultiplicativeModelParams(double allelicRisk) {
        lociNum = diseaseSNPs.size();
        if (lociNum == 0) {
            System.err.println("No disease SNPs!");
            return;
        }
        double[] riskAlleleFreq = new double[lociNum];
        for (int i = 0; i < lociNum; i++) {
            riskAlleleFreq[i] = diseaseSNPs.get(i).getRiskAlleleFrequency();
        }
        //assume each locus has two alleles
        int alleleNum = riskAlleleFreq.length * 2;

        int[] indices;
        Set<Integer> selectedIDs = new HashSet<Integer>();
        double[] nonRiskAlleleFreq = new double[lociNum];
        double testProb = 0;

        double accumulatedProb = 0;
        double localProb = 0;
        double tmpProb = 1;

        for (int i = 0; i < lociNum; i++) {
            nonRiskAlleleFreq[i] = 1 - riskAlleleFreq[i];
        }

        for (int i = 0; i < alleleNum; i++) {
            tmpProb *= nonRiskAlleleFreq[i / 2];
        }
        accumulatedProb += tmpProb;
        testProb = tmpProb;
        // System.out.println(tmpProb);
        for (int a = 1; a <= alleleNum; a++) {
            CombinationGenerator x = new CombinationGenerator(alleleNum, a);
            localProb = 0;
            while (x.hasMore()) {
                indices = x.getNext();
                selectedIDs.clear();
                for (int i = 0; i < indices.length; i++) {
                    selectedIDs.add(indices[i]);
                }
                tmpProb = 1;
                for (int i = 0; i < alleleNum; i++) {
                    if (selectedIDs.contains(i)) {
                        tmpProb *= riskAlleleFreq[i / 2];
                    } else {
                        tmpProb *= nonRiskAlleleFreq[i / 2];
                    }
                }
                localProb += tmpProb;
            }
            testProb += localProb;
            // System.out.println(localProb);
            accumulatedProb += (Math.pow(allelicRisk, a) * localProb);
        }

        System.out.println(testProb);
        baselineRisk = prevalance / accumulatedProb;
        jointGenotypePenerances = new double[alleleNum + 1];
        for (int i = 0; i <= alleleNum; i++) {
            jointGenotypePenerances[i] = baselineRisk * (Math.pow(allelicRisk, i));
            // System.out.println(jointGenotypePenerances[i]);
        }
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        double[] riskAlleleFreq = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
        double[] allelicRisk = {1.14, 1.14, 1.14, 1.14, 1.14, 1.14, 1.14};
        int lociNum = allelicRisk.length;

        RischMultipLociModel rsModel = new RischMultipLociModel(0.1);
        try {
            for (int i = 0; i < riskAlleleFreq.length; i++) {
                rsModel.addDiseaseSNPLoci(riskAlleleFreq[i], true, i);
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        rsModel.calculateAdditiveModelParams(1.14);
        double[] accumlateRiks = rsModel.jointGenotypePenerances;
        for (int i = 0; i < accumlateRiks.length; i++) {
            System.out.println(accumlateRiks[i]);
        }
    }
}
