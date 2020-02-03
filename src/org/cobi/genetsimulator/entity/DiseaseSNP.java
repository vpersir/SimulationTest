/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.cobi.genetsimulator.entity;

import org.cobi.eqtlsimulator.Constants;
import static org.cobi.eqtlsimulator.Constants.DEFAULT_MISSING_DOUBLE;

 

/**
 *
 * @author mxli
 */
public class DiseaseSNP implements Constants {

    public double riskAlleleFrequency = DEFAULT_MISSING_DOUBLE;
    public double genotypicRelativeRiskAa = DEFAULT_MISSING_DOUBLE;; //relative risk of heterozygous genotype of risk alleles
    public double genotypicRelativeRiskAA = DEFAULT_MISSING_DOUBLE;;//relative risk of homozygous genotype of risk alleles
    public boolean riskAlleleLable = true;
    public int LDMarkerPosition;

    public DiseaseSNP(double riskAlleleFrequency, double genotypicRelativeRiskAa, double genotypicRelativeRiskAA, boolean riskAlleleLable, int LDMarkerPosition) {
        this.genotypicRelativeRiskAa = genotypicRelativeRiskAa;
        this.genotypicRelativeRiskAA = genotypicRelativeRiskAA;
        this.riskAlleleFrequency = riskAlleleFrequency;
        this.riskAlleleLable = riskAlleleLable;
        this.LDMarkerPosition = LDMarkerPosition;
    }

    public DiseaseSNP(double riskAlleleFrequency,  boolean riskAlleleLable, int LDMarkerPosition) {
        this.riskAlleleFrequency = riskAlleleFrequency;
        this.riskAlleleLable = riskAlleleLable;
        this.LDMarkerPosition = LDMarkerPosition;
    }

    public int getLDMarkerPosition() {
        return LDMarkerPosition;
    }

    public void setLDMarkerPosition(int LDMarkerPosition) {
        this.LDMarkerPosition = LDMarkerPosition;
    }

    public double getGenotypicRelativeRiskAA() {
        return genotypicRelativeRiskAA;
    }

    public void setGenotypicRelativeRiskAA(double genotypicRelativeRiskAA) {
        this.genotypicRelativeRiskAA = genotypicRelativeRiskAA;
    }

    public double getGenotypicRelativeRiskAa() {
        return genotypicRelativeRiskAa;
    }

    public void setGenotypicRelativeRiskAa(double genotypicRelativeRiskAa) {
        this.genotypicRelativeRiskAa = genotypicRelativeRiskAa;
    }

    public double getRiskAlleleFrequency() {
        return riskAlleleFrequency;
    }

    public void setRiskAlleleFrequency(double riskAlleleFrequency) {
        this.riskAlleleFrequency = riskAlleleFrequency;
    }

    public boolean isRiskAlleleLable() {
        return riskAlleleLable;
    }

    public void setRiskAlleleLable(boolean riskAlleleLable) {
        this.riskAlleleLable = riskAlleleLable;
    }
}
