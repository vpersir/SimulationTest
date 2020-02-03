/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.genetsimulator.controller;

import java.util.ArrayList;
import java.util.List;
import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;
 
 
import org.cobi.genetsimulator.entity.DiseaseSNP;

/**
 *
 * @author Miaoxin Li
 */
public class LiabilityThresholdModel {

    double overallDiseasePrevalence;
    double overallThreshold;
    double liabilityThreshold;
    double backgroundStandardDeviation = 13.0;
    List<PredisposingSNP> liabiSNPs = new ArrayList<PredisposingSNP>();
    NormalDistribution normDis = new NormalDistributionImpl(0.0, 1.0);

    public LiabilityThresholdModel(double overallDiseasePrevalence) {
        this.overallDiseasePrevalence = overallDiseasePrevalence;
    }

    public void addLiabilitySNPLoci(double freq, int index, boolean riskAlleleLabel, double genotypicRelativeRisk1, double genotypicRelativeRisk2) throws Exception {
        PredisposingSNP snp = new PredisposingSNP(freq,  genotypicRelativeRisk1, genotypicRelativeRisk2,riskAlleleLabel, index);
        snp.calculateUnderlyingParamters();
        liabiSNPs.add(snp);
    }

    public void addLiabilitySNPLoci(double riskAlleleFrequency, double genotypicRelativeRisk1, double genotypicRelativeRisk2) throws Exception {
         PredisposingSNP snp = new PredisposingSNP(riskAlleleFrequency,  genotypicRelativeRisk1, genotypicRelativeRisk2,true, 0);
        snp.calculateUnderlyingParamters();
        liabiSNPs.add(snp);
    }

    public class PredisposingSNP extends DiseaseSNP {

        double dominanceAdditiveRiskEffectRatio;
        double additiveRiskEffectSize;

        public PredisposingSNP(double riskAlleleFrequency, double genotypicRelativeRisk1, double genotypicRelativeRisk2, boolean riskAlleleLable, int originalOrder) {
            super(riskAlleleFrequency, genotypicRelativeRisk1, genotypicRelativeRisk2, riskAlleleLable, originalOrder);

        }

        public double getDominanceAdditiveRiskEffectRatio() {
            return dominanceAdditiveRiskEffectRatio;
        }

        public void setDominanceAdditiveRiskEffectRatio(double dominanceAdditiveRiskEffectRatio) {
            this.dominanceAdditiveRiskEffectRatio = dominanceAdditiveRiskEffectRatio;
        }

        public void calculateUnderlyingParamters() throws Exception {
            double relativePrevalence = (1 - riskAlleleFrequency) * (1 - riskAlleleFrequency) + 2 * (1 - riskAlleleFrequency) * riskAlleleFrequency * genotypicRelativeRiskAa + riskAlleleFrequency * riskAlleleFrequency * genotypicRelativeRiskAA;
            relativePrevalence = overallDiseasePrevalence / relativePrevalence;
            normDis.setMean(0);
            normDis.setStandardDeviation(1);
            double x1 = normDis.inverseCumulativeProbability(1 - relativePrevalence);
            double x2 = normDis.inverseCumulativeProbability(1 - genotypicRelativeRiskAa * relativePrevalence);
            double x3 = normDis.inverseCumulativeProbability(1 - genotypicRelativeRiskAA * relativePrevalence);
            additiveRiskEffectSize = 0.5 * backgroundStandardDeviation * (x1 - x3);
            dominanceAdditiveRiskEffectRatio = (2 * x2 - x1 - x3) / (x3 - x1);
            liabilityThreshold = 0.5 * backgroundStandardDeviation * (x1 + x3);
            double effectVar = 2 * additiveRiskEffectSize * additiveRiskEffectSize * riskAlleleFrequency * (1 - riskAlleleFrequency) *
                    ((1 + dominanceAdditiveRiskEffectRatio * (1 - 2 * riskAlleleFrequency)) * (1 + dominanceAdditiveRiskEffectRatio * (1 - 2 * riskAlleleFrequency)) + 2 * riskAlleleFrequency * (1 - riskAlleleFrequency) * dominanceAdditiveRiskEffectRatio * riskAlleleFrequency);
            double heritablity = effectVar / (effectVar + backgroundStandardDeviation * backgroundStandardDeviation);
            StringBuffer reults = new StringBuffer();

            normDis.setMean(-additiveRiskEffectSize);
            normDis.setStandardDeviation(backgroundStandardDeviation);
            reults.append("Penetrance xx: ");
            x1 = 1 - normDis.cumulativeProbability(liabilityThreshold);
            reults.append(x1);
            reults.append("\n");

            normDis.setMean(additiveRiskEffectSize * dominanceAdditiveRiskEffectRatio);
            x2 = 1 - normDis.cumulativeProbability(liabilityThreshold);
            reults.append("Penetrance xX: ");
            reults.append(x2 / x1);
            reults.append("\n");
            normDis.setMean(additiveRiskEffectSize);
            reults.append("Penetrance XX: ");
            x3 = 1 - normDis.cumulativeProbability(liabilityThreshold);
            reults.append(x3 / x1);
            reults.append("\n");
            reults.append("Dominanc Effect Size: ");
            reults.append(dominanceAdditiveRiskEffectRatio * additiveRiskEffectSize);
            reults.append("\n");
            reults.append("Liability Threshold: ");
            reults.append(liabilityThreshold);
            reults.append("\n");
            reults.append("Heritablity: ");
            reults.append(heritablity * 100);
            reults.append("%");
            System.out.println(reults);
        }
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        System.out.println("  ************************************************************");
        System.out.println("  *                        Hello!                            *");
        System.out.println("  *          Thank you for using my program!                 *");
        System.out.println("  *       If there are any bugs, please inform me.           *");
        System.out.println("  *              MX Li; Email: limx54@yahoo.com              *");
        System.out.println("  ************************************************************");

        if (args.length < 4) {
            System.out.println("java -jar GenetSimulator.jar <prevalence> <risk allele frequency in general population> <heterozygous genetic risk> <homozygous genetic risk>");
            return;
        }

        double prevalence = Double.parseDouble(args[0]);
        double risk_allele_frequency = Double.parseDouble(args[1]);
        double heterozygous_risk = Double.parseDouble(args[2]);
        double homozygous_risk = Double.parseDouble(args[3]);

        LiabilityThresholdModel liability = new LiabilityThresholdModel(prevalence);
        System.out.println("prevalence: " + prevalence);
        System.out.println("risk allele frequency in general population: " + risk_allele_frequency);
        System.out.println("heterozygous relative genetic risk: " + heterozygous_risk);
        System.out.println("homozygous relative genetic risk: " + homozygous_risk);
        System.out.println();

        try {
            liability.addLiabilitySNPLoci(risk_allele_frequency, heterozygous_risk, homozygous_risk);
        // liability.addLiabilitySNPLoci(0.3, 2000, 5000);
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
}
