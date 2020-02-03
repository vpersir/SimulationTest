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
 * @author mxli
 */
public class MultiLocusLiabilityThresholdModel {

    final double DEFAULT_MISSING_NUMBER = -9.99;
    final double TOLERATED_ERROR = 1e-6;
    double overallDiseasePrevalence;
    double overallThreshold;
    double accumulatedVariance = 0.0;
    List<LiabilityDiseaseSNP> liabiSNPs = new ArrayList<LiabilityDiseaseSNP>();
    NormalDistribution normDis = new NormalDistributionImpl(0.0, 1.0);

    public MultiLocusLiabilityThresholdModel(double overallDiseasePrevalence) throws Exception {
        this.overallDiseasePrevalence = overallDiseasePrevalence;
        normDis.setMean(0);
        normDis.setStandardDeviation(1);
        overallThreshold = normDis.inverseCumulativeProbability(1 - overallDiseasePrevalence);
        accumulatedVariance = 0.0;
    }

    public void addLiabilitySNPLoci(double riskAlleleFreq, boolean riskAlleleLabel, double genotypicRelativeRisk1, double genotypicRelativeRisk2,
            int LDMarkerPosition, double markerAlleleFreq, double rSquare) throws Exception {
        LiabilityDiseaseSNP snp = new LiabilityDiseaseSNP(riskAlleleFreq,
                genotypicRelativeRisk1, genotypicRelativeRisk2, riskAlleleLabel, LDMarkerPosition, markerAlleleFreq, rSquare);
        snp.calculateLiabilities();
        liabiSNPs.add(snp);
        accumulatedVariance += snp.explainedVariance;
        if (accumulatedVariance >= 1) {
            throw new Exception("Error! The overal explianed variance is larger than 1!");
        }
    }

    public double getAccumulatedVariance() {
        return accumulatedVariance;
    }

    public void setAccumulatedVariance(double accumulatedVariance) {
        this.accumulatedVariance = accumulatedVariance;
    }

    public List<LiabilityDiseaseSNP> getLiabiSNPs() {
        return liabiSNPs;
    }

    public void setLiabiSNPs(List<LiabilityDiseaseSNP> liabiSNPs) {
        this.liabiSNPs = liabiSNPs;
    }

    public double getOverallThreshold() {
        return overallThreshold;
    }

    public void setOverallThreshold(double overallThreshold) {
        this.overallThreshold = overallThreshold;
    }

    public class LiabilityDiseaseSNP extends DiseaseSNP {

        public double liabilityAA; //homozygous genotype of risk alleles
        public double liabilityAa; //heterozygous genotype of risk alleles
        public double liabilityaa; //homozygous genotype of wide-type alleles   
        double explainedVariance;
        public double probAcM;
        public double probAcm;

        public LiabilityDiseaseSNP(double riskAlleleFrequency, double genotypicRelativeRisk1, double genotypicRelativeRisk2, boolean riskAlleleLable, 
                int LDMarkerPosition, double markerAlleleFreq, double rSquare) throws Exception {
            super(riskAlleleFrequency, genotypicRelativeRisk1, genotypicRelativeRisk2, riskAlleleLable, LDMarkerPosition);
            this.probAcM = rSquare * Math.sqrt(riskAlleleFrequency * (1 - riskAlleleFrequency) * markerAlleleFreq * (1 - markerAlleleFreq)) + riskAlleleFrequency * markerAlleleFreq;
            this.probAcm = riskAlleleFrequency - this.probAcM;
            this.probAcM = this.probAcM / markerAlleleFreq;
            this.probAcm = this.probAcm / (1 - markerAlleleFreq);
            //System.out.println("P(A|M)= " + probAcM);
            //System.out.println("P(A|m)= " + probAcm);
            if (probAcM < 0 && probAcM > 1) {
                throw new Exception("P(A|M)= " + probAcM+" P(A|M) is not proper!");
            }
            if (probAcm < 0 && probAcm > 1) {
                throw new Exception("P(A|m)= " + probAcm+" P(A|m) is not proper!");
            }
        }

        public void calculateLiabilities() throws Exception {
            double relativePrevalence = (1 - riskAlleleFrequency) * (1 - riskAlleleFrequency) + 2 * (1 - riskAlleleFrequency) * riskAlleleFrequency * genotypicRelativeRiskAa + riskAlleleFrequency * riskAlleleFrequency * genotypicRelativeRiskAA;
            StringBuffer reults = new StringBuffer();
            double x1 = overallDiseasePrevalence / (relativePrevalence);
            solveQuadraticEquation(x1);
            double quantile1 = overallThreshold + (Math.pow((1 - riskAlleleFrequency), 2) * liabilityaa + 2 * (1 - riskAlleleFrequency) * riskAlleleFrequency * liabilityAa) / (riskAlleleFrequency * riskAlleleFrequency);
            quantile1 = quantile1 / (Math.sqrt(1 - explainedVariance));
            double diff1 = 1 - normDis.cumulativeProbability(quantile1) - genotypicRelativeRiskAA * x1;
            if (Math.abs(diff1) < TOLERATED_ERROR) {
                // System.out.println(diff1);
                liabilityAA = -(1 - riskAlleleFrequency) * (liabilityaa * (1 - riskAlleleFrequency) / (riskAlleleFrequency * riskAlleleFrequency) + 2 * liabilityAa / riskAlleleFrequency);
                System.out.println("Explained Variance:" + explainedVariance);
                return;
            }
            double diff2;

            double searchLen = 0.1;
            double startPoint = x1 + searchLen;
            double direction = 1;
            do {
                if (startPoint <= 0 || startPoint >= 1) {
                    break;
                }
                if (solveQuadraticEquation(startPoint)) {
                    quantile1 = overallThreshold + (Math.pow((1 - riskAlleleFrequency), 2) * liabilityaa + 2 * (1 - riskAlleleFrequency) * riskAlleleFrequency * liabilityAa) / (riskAlleleFrequency * riskAlleleFrequency);
                    quantile1 = quantile1 / (Math.sqrt(1 - explainedVariance));
                    diff2 = 1 - normDis.cumulativeProbability(quantile1) - genotypicRelativeRiskAA * startPoint;
                    if (Math.abs(diff2) < TOLERATED_ERROR) {
                        // System.out.println(diff2 + " " + startPoint);
                        liabilityAA = -(1 - riskAlleleFrequency) * (liabilityaa * (1 - riskAlleleFrequency) / (riskAlleleFrequency * riskAlleleFrequency) + 2 * liabilityAa / riskAlleleFrequency);
                        //System.out.println("Explained Variance:" + explainedVariance);
                        System.out.println("Explained Variance:" + explainedVariance);
                        return;
                    } else if (diff1 > 0 && diff2 > 0) {
                        if (diff1 < diff2) {
                            direction = -1;
                        }
                        startPoint += searchLen * direction;
                    } else if (diff1 < 0 && diff2 < 0) {
                        if (diff1 > diff2) {
                            direction = -1;
                        }
                        startPoint += searchLen * direction;
                    } else {
                        // otherwise the symbols of diffs must be opersite
                        //go back to the startPoint
                        startPoint = startPoint - searchLen * direction;
                        searchLen /= 10;
                        startPoint += searchLen * direction;
                    }
                } else {
                    //go back to the startPoint and narrow down the search length
                    startPoint = startPoint - searchLen * direction;
                    searchLen /= 10;
                    startPoint += searchLen * direction;
                }
            } while (true);


        /*
        reults.append(coefficients[2]);
        reults.append("\n");
        reults.append("B: ");
        reults.append(coefficients[1]);
        reults.append("\n");
        reults.append("C: ");
        reults.append(coefficients[0]);
        reults.append("\n");
        
        reults.append("\n");
        reults.append("Heritablity: ");
        //  reults.append(heritablity * 100);
        reults.append("%");
        System.out.println(reults);
        normDis.setMean(0);
        normDis.setStandardDeviation(1);
        double tt = (overallThreshold - liabilityaa) / Math.sqrt(1 - explainedVariance);
        double prob0 = 1 - normDis.cumulativeProbability(tt);
        double tt1 = (overallThreshold - liabilityAa) / Math.sqrt(1 - explainedVariance);
        double prob1 = 1 - normDis.cumulativeProbability(tt1);
        double tt2 = (overallThreshold - liabilityAA) / Math.sqrt(1 - explainedVariance);
        double prob2 = 1 - normDis.cumulativeProbability(tt2);
        System.out.println(prob1 / prob0);
        System.out.println(prob2 / prob0);
        
        normDis.setMean(liabilityaa);
        normDis.setStandardDeviation(Math.sqrt(1 - explainedVariance));
        prob0 = 1 - normDis.cumulativeProbability(overallThreshold);
        normDis.setMean(liabilityAa);
        prob1 = 1 - normDis.cumulativeProbability(overallThreshold);
        normDis.setMean(liabilityAA);
        prob2 = 1 - normDis.cumulativeProbability(overallThreshold);
        System.out.println(prob1 / prob0);
        System.out.println(prob2 / prob0);
         */
        }

        private double[] obtainCoefficientsQuadraticEquation(double x1) throws Exception {
            double[] coefficients = new double[3];
            double r1x1 = normDis.inverseCumulativeProbability(1 - genotypicRelativeRiskAa * x1);
            x1 = normDis.inverseCumulativeProbability(1 - x1);
            double af = (1 - riskAlleleFrequency) * (1 - riskAlleleFrequency) * (1 + ((1 - riskAlleleFrequency) / riskAlleleFrequency) * (1 - riskAlleleFrequency) / riskAlleleFrequency) + 1 / (x1 * x1);
            coefficients[2] = af * Math.pow(x1 / r1x1, 2) + 2 * (1 - riskAlleleFrequency) * (2 - riskAlleleFrequency) + 4 * Math.pow(1 - riskAlleleFrequency, 3) * x1 / (riskAlleleFrequency * r1x1);
            coefficients[1] = 2 * af * overallThreshold * (1 - x1 / r1x1) * x1 / r1x1 + 4 * Math.pow(1 - riskAlleleFrequency, 3) / riskAlleleFrequency * overallThreshold * (1 - x1 / r1x1) - 2 * overallThreshold / (x1 * r1x1);
            coefficients[0] = af * Math.pow(overallThreshold * (1 - x1 / r1x1), 2) - 1 - Math.pow(overallThreshold / x1, 2) + 2 * overallThreshold * overallThreshold / (x1 * r1x1);
            return coefficients;
        }

        private boolean solveQuadraticEquation(double x1) throws Exception {
            double[] coefficients = obtainCoefficientsQuadraticEquation(x1); 
            double A,      B,      C,      root1,      root2,      discriminant;
            A = coefficients[2];
            B = coefficients[1];
            C = coefficients[0];
            double r1x1 = normDis.inverseCumulativeProbability(1 - genotypicRelativeRiskAa * x1);
            x1 = normDis.inverseCumulativeProbability(1 - x1);
            //Info.
            //System.out.println("This snippet computes the quadratic formula." + "\nIt works for quadratic equations in the Standard Form: Ax^2 + Bx + C = 0" + "\nIt even determines whether the roots are imaginary!\n");
            /*Now, we'll check the discriminant to see if the roots are imaginary or not.
             *Discriminant = B^2 - 4(A)(C)
             *If negative, then the roots are imaginary, otherwise the roots are real (if zero, there is one real root): */
            discriminant = (Math.pow(B, 2)) - (4 * A * C);
            if (discriminant == 0.0) {           //One root.

                root1 = (-B + Math.sqrt(discriminant)) / (2 * A);
                liabilityAa = root1;
                liabilityaa = overallThreshold * (1 - x1 / r1x1) + liabilityAa * x1 / r1x1;
                explainedVariance = Math.pow(1 - riskAlleleFrequency, 2) * (1 + Math.pow((1 - riskAlleleFrequency) / riskAlleleFrequency, 2)) * Math.pow(liabilityaa, 2);
                explainedVariance += 2 * (1 - riskAlleleFrequency) * (2 - riskAlleleFrequency) * Math.pow(liabilityAa, 2);
                explainedVariance += 4 * Math.pow(1 - riskAlleleFrequency, 3) * liabilityAa * liabilityaa / riskAlleleFrequency;

                System.out.println(String.format("There is one root at: %.4f\n\n", root1));
                if (explainedVariance > 1) {
                    liabilityAa = liabilityaa = liabilityAA = explainedVariance = DEFAULT_MISSING_NUMBER;
                    //System.out.println(String.format("There is one root at: %.4f\n\n", root1));
                    //System.out.println("No proper solution");
                    return false;
                }
                return true;
            } else {
                if (discriminant > 0.0) {        //Two real roots.

                    root1 = (-B + Math.sqrt(discriminant)) / (2 * A);
                    double liabilityAa1 = root1;
                    double liabilityaa1 = overallThreshold * (1 - x1 / r1x1) + liabilityAa1 * x1 / r1x1;
                    double Vg1 = Math.pow(1 - riskAlleleFrequency, 2) * (1 + Math.pow((1 - riskAlleleFrequency) / riskAlleleFrequency, 2)) * Math.pow(liabilityaa1, 2);
                    Vg1 += 2 * (1 - riskAlleleFrequency) * (2 - riskAlleleFrequency) * Math.pow(liabilityAa1, 2);
                    Vg1 += 4 * Math.pow(1 - riskAlleleFrequency, 3) * liabilityAa1 * liabilityaa1 / riskAlleleFrequency;
                    //double tmp1 = 1 - Math.pow((overallThreshold - liabilityaa1) / x1, 2);
                    //double tt = overallThreshold - liabilityaa1 - x1 * Math.sqrt(1 - Vg1);

                    root2 = (-B - Math.sqrt(discriminant)) / (2 * A);
                    double liabilityAa2 = root2;
                    double liabilityaa2 = overallThreshold * (1 - x1 / r1x1) + liabilityAa2 * x1 / r1x1;
                    double Vg2 = Math.pow(1 - riskAlleleFrequency, 2) * (1 + Math.pow((1 - riskAlleleFrequency) / riskAlleleFrequency, 2)) * Math.pow(liabilityaa2, 2);
                    Vg2 += 2 * (1 - riskAlleleFrequency) * (2 - riskAlleleFrequency) * Math.pow(liabilityAa2, 2);
                    Vg2 += 4 * Math.pow(1 - riskAlleleFrequency, 3) * liabilityAa2 * liabilityaa2 / riskAlleleFrequency;
                    //double tmp2 = 1 - Math.pow((overallThreshold - liabilityaa2) / x1, 2);
                    //double tt2 = overallThreshold - liabilityaa2 - x1 * Math.sqrt(1 - Vg2);
                    // the variance is unlikely to be very large
                    if (Vg1 < 1 && Vg2 < 1) {
                        if (Vg1 < Vg2) {
                            liabilityAa = root1;
                            liabilityaa = liabilityaa1;
                            explainedVariance = Vg1;
                        } else {
                            liabilityAa = root2;
                            liabilityaa = liabilityaa2;
                            explainedVariance = Vg2;
                        }
                        //System.out.println("Two possible liabilites");
                        return true;
                    } else if (Vg1 < 1) {
                        liabilityAa = root1;
                        liabilityaa = liabilityaa1;
                        explainedVariance = Vg1;
                    } else if (Vg2 < 1) {
                        liabilityAa = root2;
                        liabilityaa = liabilityaa2;
                        explainedVariance = Vg2;
                    } else {
                        liabilityAa = liabilityaa = liabilityAA = explainedVariance = DEFAULT_MISSING_NUMBER;
                        //System.out.println("No proper solution");
                        return false;
                    }
                    //System.out.println(String.format("There are two real roots at: %.4f and %.4f\n\n", root1, root2));
                    return true;
                } else {
                    if (discriminant < 0.0) {    //Two imaginary roots.
                        //System.out.println("No proper solution");
                        root1 = (-B + Math.sqrt(-discriminant)) / (2 * A);
                        root2 = (-B - Math.sqrt(-discriminant)) / (2 * A);
                    //System.out.println(String.format("There are two imaginary roots at: %.4fi and %.4fi\n\n", root1, root2));
                    }
                    liabilityAa = liabilityaa = liabilityAA = explainedVariance = DEFAULT_MISSING_NUMBER;
                    return false;
                }
            }
        }
    }
}
